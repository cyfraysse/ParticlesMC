"""
ParticlesMC: Monte Carlo simulation framework for particle systems.

Provides core types, utilities, and the `particlesmc` command implemented with Comonicon.
Exports commonly-used types (e.g., `Particles`, `Model`) and helper functions for simulation control, I/O, and moves.
"""
module ParticlesMC

using Arianna, StaticArrays, Transducers
using Comonicon, TOML
using Comonicon: @main

export Particles
abstract type Particles <: AriannaSystem end


include("utils.jl")
include("neighbours.jl")
include("models.jl")
include("molecules.jl")
include("atoms.jl")
include("moves.jl")
include("rotation.jl")

"""Return the position of particle `i` in `system`.

# Arguments
- `system::Particles`: the particle system
- `i::Int`: particle index

# Returns
- Coordinates of particle `i` (e.g., an `SVector` or array).
"""
get_position(system::Particles, i::Int) = @inbounds system.position[i]

"""Return the species index (type) of particle `i`.

# Arguments
- `system::Particles`: the particle system
- `i::Int`: particle index

# Returns
- `Int`: species identifier of particle `i`.
"""
get_species(system::Particles, i::Int) = @inbounds system.species[i]

"""Return the interaction `Model` between the species of particles `i` and `j`.

# Arguments
- `system::Particles`: the particle system
- `i::Int`, `j::Int`: particle indices

# Returns
- `Model` object or callable describing pair interactions for the two species.
"""
get_model(system::Particles, i::Int, j::Int) = @inbounds system.model_matrix[get_species(system, i), get_species(system, j)]

"""Return the simulation box of `system`.

# Returns
- Box description (usually vector or struct) representing periodic box extents.
"""
get_box(system::Particles) = system.box

"""Return the neighbour list of `system`.

# Returns
- The neighbour list object used for pair evaluations (e.g., `NeighbourList`, `LinkedList`).
"""
get_neighbour_list(system::Particles) = system.neighbour_list

"""Return the number of particles in `system`.

Overloads `Base.length` for `Particles`.
"""
Base.length(system::Particles) = system.N

"""Return a proper index range for `system`.

Overloads `Base.eachindex` to allow fast indexing.
"""
Base.eachindex(system::Particles) = Base.OneTo(length(system))

"""Return the (position, species) tuple for atom `i`.

This overload supports indexing `atoms[i]` to get coordinates and species.
"""
Base.getindex(system::Atoms, i::Int) = system.position[i], system.species[i]

"""Iterate over `Atoms` or `Molecules` returning the position and next state.

Conforms to Julia iterator interface; yields the position of the current index.
"""
function Base.iterate(system::Union{Atoms,Molecules}, state=1)
    state > length(system) && return nothing  # Stop iteration
    return (system.position[state], state + 1)  # Return element & next state
end

"""Compute the energy contribution of particle `i` in `system` using its neighbour list.

If a neighbour list is provided it will be used for the evaluation.
"""
function compute_energy_particle(system::Particles, i::Int)
    return compute_energy_particle(system, i, system.neighbour_list)
end

"""Compute the energy contribution for each particle index in `ids`.

Returns an array with per-particle energies by mapping `compute_energy_particle`.
"""
function compute_energy_particle(system::Particles, ids::AbstractVector)
    return map(i -> compute_energy_particle(system, i), ids)
end

"""
    Helper for building scheduler based on TOML schedul parameters
"""

function parse_schedule(scheduler_params, steps, burn)
    if haskey(scheduler_params, "multi_origins")
        tw   = scheduler_params["multi_origins"]["tw"]
        N    = scheduler_params["multi_origins"]["N"]
        return build_schedule(steps, MultiOrigins(tw, N), burn=burn)
    elseif haskey(scheduler_params, "log_base")
        interval = get(scheduler_params, "linear_interval", 1)
        block    = build_schedule(interval, 0, 2.0)
        return build_schedule(steps, burn, block)
    else
        interval = get(scheduler_params, "linear_interval", 1)
        return build_schedule(steps, burn, interval)
    end
end

export energy
#export nearest_image_distance
export Model, GeneralKG, JBB, BHHP, SoftSpheres, KobAndersen, Trimer
export NeighbourList, LinkedList, CellList, EmptyList, VerletList
export Atoms, Molecules
export Displacement, DiscreteSwap, MoleculeFlip
export fold_back, System
export SimpleGaussian, DoubleUniform, EnergyBias
export sample_action!, log_proposal_density, reward, invert_action!, delta_log_target_density
export perform_action!, revert_action!
include("IO/IO.jl")
using .IO: XYZ, EXYZ, LAMMPS, load_configuration, load_chains
export XYZ, EXYZ, LAMMPS, load_configuration, load_chains
export ComputeRotation, StorePhiTrajectories, StoreLastPhiFrame, load_phi_frame



"""
ParticlesMC implemented in Comonicon.

# Arguments

- `params`: Path to the TOML parameter file.
"""
@main function particlesmc(params::String)
    if !isfile(params)
        error("Parameter file '$params' does not exist in the current path.")
    end
    params = TOML.parsefile(params)

    # Extract system parameters
    system = params["system"]
    temperature = system["temperature"]
    density = system["density"]
    config = system["config"]
    model = get(system, "model", nothing)
    if model === nothing
        model = params["model"]
    end  # optional field
    filename = isdir(config) ? system["filename"] : ""
    list_type = get(system, "list_type", "LinkedList")  # optional field
    list_parameters = get(system, "list_parameters", nothing)  # optional field
    bonds = get(system, "bonds", nothing)

    # Extract simulation parameters
    sim = params["simulation"]
    steps = sim["steps"]
    burn = get(sim, "burn", 0)
    seed = sim["seed"]
    parallel = sim["parallel"]
    output_path = get(sim, "output_path", "./")

    # Checkpoint / restart parameter (absent = single job, no checkpointing)
    t_restart = get(sim, "t_restart", nothing)

    # ── Auto-detect restart ───────────────────────────────────────────────────
    t_start = isnothing(t_restart) ? 0 : Arianna.detect_restart(output_path, EXYZ())

    if t_start >= steps
        println("Simulation already complete (checkpoint at t=$t_start, steps=$steps). Nothing to do.")
        return
    end
    t_start > 0 && println("Checkpoint found at t=$t_start — resuming simulation.")

    # ── Load chains ───────────────────────────────────────────────────────────
    load_args = Dict(
        "temperature"     => temperature,
        "density"         => density,
        "model"           => model,
        "list_type"       => list_type,
        "list_parameters" => list_parameters,
    )
    bonds !== nothing && (load_args["bonds"] = bonds)

    if t_start > 0
        n_chains = length([f for f in readdir(joinpath(output_path, "chains"))
                           if isdir(joinpath(output_path, "chains", f))])
        ckpt_files = [joinpath(output_path, "chains", "$c", "lastframe.exyz")
                      for c in 1:n_chains]
        all(isfile, ckpt_files) || error("lastframe.exyz missing for one or more chains")
        chains = vcat([load_chains(f, args=load_args) for f in ckpt_files]...)
    else
        if !isfile(config) && !isdir(config)
            error("Configuration file '$config' does not exist in the current path.")
        end
        chains = load_chains(config, args=load_args, filename=filename)
    end

    algorithm_list = []

    # ── Setup moves ───────────────────────────────────────────────────────────
    pool = []
    for move in sim["move"]
        prob = move["probability"]
        policy = move["policy"]
        action = move["action"]
        parameters = get(move, "parameters", Dict())
        param_obj = ComponentArray()

        if action == "Displacement"
            action_obj = Displacement(0, zero(chains[1].box), 0.0)
            if "sigma" in keys(parameters)
                param_obj = ComponentArray(σ=parameters["sigma"])
            else
                error("Missing parameter 'sigma' for action: $action")
            end
            if policy == "SimpleGaussian"
                policy_obj = SimpleGaussian()
            else
                error("Unsupported policy: $policy for action: $action")
            end
        elseif action == "MoleculeFlip"
            action_obj = MoleculeFlip(0, 0, 0.0)
            param_obj = Vector{Float64}()
            if policy == "DoubleUniform"
                policy_obj = DoubleUniform()
            else
                error("Unsupported policy: $policy for action: $action")
            end
        elseif action == "DiscreteSwap"
            if "species" in keys(parameters)
                species = parameters["species"]
                if length(species) != 2 || eltype(species) != Int
                    error("'species' for action: $action must be an array of two ints")
                end
            else
                error("Missing parameter 'species' for action: $action")
            end
            action_obj = DiscreteSwap(species, chains[1])
            param_obj = Vector{Float64}()
            if policy == "DoubleUniform"
                policy_obj = DoubleUniform()
            else
                error("Unsupported policy: $policy for action: $action")
            end
        else
            error("Unsupported action: $action")
        end
        move_obj = Move(action_obj, policy_obj, param_obj, prob)
        push!(pool, move_obj)
    end
    push!(algorithm_list, (algorithm=Metropolis, pool=pool, seed=seed, parallel=parallel, sweepstep=length(chains[1])))

    # ── Setup observables ─────────────────────────────────────────────────────
    has_compute_rotation = false
    for observable in get(sim, "observable", [])
        alg = observable["algorithm"]
        scheduler_params = observable["scheduler_params"]
        sched = parse_schedule(scheduler_params, steps, burn)
        if alg == "ComputeRotation"
            parameters = get(observable, "parameters", Dict())
            theta_T    = Float64.(get(parameters, "theta_T", [π/4]))
            algorithm  = (
                algorithm=ComputeRotation,
                scheduler=sched,
                theta_T=theta_T,
            )
            has_compute_rotation = true
        else
            error("Unsupported observable algorithm: $alg")
        end
        push!(algorithm_list, algorithm)
    end

    # ── Setup outputs ─────────────────────────────────────────────────────────
    for output in sim["output"]
        alg = output["algorithm"]
        scheduler_params = output["scheduler_params"]
        dependencies = get(output, "dependencies", nothing)
        callbacks = get(output, "callbacks", [])
        fmt = get(output, "fmt", "XYZ")
        sched = parse_schedule(scheduler_params, steps, burn)
        if alg == "StoreCallbacks"
            callbacks = map(c -> eval(Meta.parse("$c")), callbacks)
            algorithm = (
                algorithm=eval(Meta.parse(alg)),
                callbacks=callbacks,
                scheduler=sched,
            )
        elseif alg == "StoreAcceptance"
            dependencies = map(d -> eval(Meta.parse("$d")), dependencies)
            algorithm = (
                algorithm=eval(Meta.parse(alg)),
                dependencies=dependencies,
                scheduler=sched,
            )
        elseif alg == "StoreTrajectories" || alg == "StoreLastFrames"
            algorithm = (
                algorithm=eval(Meta.parse(alg)),
                scheduler=sched,
                fmt=eval(Meta.parse("$(fmt)()")),
            )
        elseif alg == "StorePhiTrajectories" || alg == "StoreLastPhiFrame"
            algorithm = (
                algorithm=eval(Meta.parse(alg)),
                scheduler=sched,
                path=output_path,
            )
        elseif alg == "PrintTimeSteps"
            algorithm = (
                algorithm=eval(Meta.parse(alg)),
                scheduler=sched,
            )
        else
            error("Unsupported output algorithm: $alg")
        end
        push!(algorithm_list, algorithm)
    end

    # ── Auto-add checkpoint algorithms ───────────────────────────────────────
    # StoreLastFrames is required by Arianna when t_restart is set; add it with
    # EXYZ so detect_restart(output_path, EXYZ()) finds the right file.
    if !isnothing(t_restart)
        has_lf = any(c -> haskey(c, :algorithm) && c.algorithm === StoreLastFrames, algorithm_list)
        if !has_lf
            push!(algorithm_list, (algorithm=StoreLastFrames, fmt=EXYZ()))
        end
    end
    if has_compute_rotation
        push!(algorithm_list, (
            algorithm = StoreLastPhiFrame,
            path      = output_path,
        ))
    end

    path = joinpath(output_path)
    simulation = Simulation(chains, algorithm_list, steps;
                            t_start=t_start, t_restart=t_restart, path=path, verbose=true)

    # Run the simulation; exit 1 if more jobs are needed, 0 if complete
    done = run!(simulation)
    exit(done ? 0 : 1)

end

end
