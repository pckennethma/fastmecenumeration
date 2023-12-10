using BenchmarkTools, Graphs, Statistics

include("new_enumerate.jl")
include("utils.jl")

"""
    run_eval(instance_path::String, timeout=600.0, verbose=true)

Run the experiment for a specific instance. Parameters:
- instance_path: Path to the instance file
- timeout: The time in seconds after which a single benchmark execution should be cancelled
- verbose: `true` for detailed logging, `false` otherwise
"""
function run_eval(instance_path::String, output_folder::String, timeout=60.0, verbose=true)
    G = readgraph(instance_path, false)  # Assuming the graph is directed
    n = nv(G)
    mdir, mundir = ndiredges(G), nundiredges(G)
    nedge = ne(G)
    m = mdir + mundir

    if verbose
        @info "Instance: $instance_path (n=$n, m=$m, mdir=$mdir, mundir=$mundir, nedge=$nedge, type=PDAG)"
    end

    ts = time_ns()
    measure = Measurement(ts, ts, Inf, -Inf, 0.0, 0.0, 0, timeout, "")

    # try
    G_ = copy(G)
    res = cpdag_enumerate(G_, measure, output_folder)  # Run cpdag enumeration

    if verbose
        @info "# of DAGs: $res"
    end

    measure.std = sqrt(measure.std / (measure.n - 1))
    mmin = nanosec2millisec(measure.min)
    mmax = nanosec2millisec(measure.max)
    mmean = nanosec2millisec(measure.mean)
    mstd = nanosec2millisec(measure.std)

    if verbose
        @info "Results for 'pdag_enumerate': Min: $mmin ms, Max: $mmax ms, Mean: $mmean ms, Std: $mstd ms"
    end
    # catch err
    #     if !isa(err, TimeoutException)
    #         @error err
    #         exit()
    #     end
    # end
end

# Example usage: julia script.jl "path/to/instance.gr"
if length(ARGS) != 2
    @error "Please provide a single instance path and an output folder as arguments!"
    exit()
end

instance_path = ARGS[1]
output_folder = ARGS[2]
run_eval(instance_path, output_folder, 60.0, true)
