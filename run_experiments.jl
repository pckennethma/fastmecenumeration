using BenchmarkTools, Graphs, Statistics

include("new_enumerate.jl")
include("pdag_enumerate.jl")
include("dfs_enumerate.jl")
include("chickering_enumerate.jl")
include("meek.jl")
include("utils.jl")
include("counting.jl")
include("plots.jl")


@enum InputType begin
    CC
    CPDAG
    PDAG
end

"""
    run_eval(dir="instances/", outf="out.csv", undir=true, itype=CC, plot=false, timeout=600.0, verbose=true)

Run the experiments. Parameters:
- dir: The directory containing the instances
- outf: The file to write the results to
- undir: `true` if the input graph contains only undirected edges and thus encodes every undirected edge as a single edge, `false` if undirected edges are encoded as two arcs
- itype: The input type of the instances, can be `CC` (undirected chordal graphs), `CPDAG`, or `PDAG`
- plot: `true` if the results should be plotted, `false` otherwise
- timeout: The time in seconds after which a single benchmark execution should be cancelled
- verbose: `true` if every single delay should be saved, `false` if only min, max, mean and stddev should be saved
"""
function run_eval(dir="instances/", outf="out.csv", undir=true, itype=CC, plot=false, timeout=60.0, verbose=true)
    algorithms = [enumerate_meek, dfs_enumerate, chickering_enumerate]
    outfile = string(@__DIR__, "/results/", outf)
    fexists = isfile(outfile)

    # Compute MEC size for CC and CPDAG inputs
    mecsz = itype == CC || itype == CPDAG

    open(outfile, "a") do io
        !fexists && write(io, "file,n,m,dir,undir,sz,count,algo,min,max,mean,std\n")
        for (root, dirs, files) in walkdir(string(@__DIR__, "/$dir"))
            for f in files
                (!occursin(".DS_Store", f) && !occursin("README", f) && !occursin(".git", f)) || continue
                fpath = string(root, f)
                G = readgraph(fpath, undir)
                n = nv(G)
                mdir, mundir = ndiredges(G), nundiredges(G)
                m = mdir + mundir
                sz = mecsz ? MECsize(G) : "?"

                @info "$f (n=$n, m=$m, mdir=$mdir, mundir=$mundir, MECsize=$sz, type=$itype)"

                results = mecsz ? Dict("CliquePicking" => sz) : Dict()
                for algo in algorithms
                    if algo == pdag_enumerate && (itype == CC || itype == CPDAG)
                        algo = cpdag_enumerate
                    end

                    @info "\t Running '$algo'..."
                    skipms = itype == CC || itype == CPDAG
                    params = algo == enumerate_meek ? [skipms] : []

                    for i in 1:2 # Do a first run and take result from second run
                        i == 2 && (for _ in 1:5 GC.gc() end) # Call GC five times
                        ts = time_ns()
                        of = string(algo, "-", replace(f, ".gr" => ".csv"))
                        fp = i == 2 && verbose ? string(@__DIR__, "/results", "/$of") : ""
                        measure = Measurement(ts, ts, Inf, -Inf, 0.0, 0.0, 0, timeout, fp)
                        try
                            G_ = copy(G)
                            res = algo(G_, measure, params...)
                            results[string(algo, "-", i)] = res
                            @info "# of DAGs: $res"
                        catch err
                            if !isa(err, TimeoutException)
                                @error err
                                exit()
                            end
                        finally
                            i == 1 && continue
                            measure.std = sqrt(measure.std / (measure.n - 1))
                            mmin = nanosec2millisec(measure.min)
                            mmax = nanosec2millisec(measure.max)
                            mmean = nanosec2millisec(measure.mean)
                            mstd = nanosec2millisec(measure.std)
                            c = measure.n
                            write(io, "$f,$n,$m,$mdir,$mundir,$sz,$c,$algo,$mmin,$mmax,$mmean,$mstd\n")
                            flush(io)
                        end
                    end
                end

                vals = collect(values(results))
                if !isempty(results) && any(r -> r != vals[1], vals)
                    @error "\t Algorithms found different results: $results"
                    exit()
                end
            end
        end
    end
    plot && plotbenchmarks(outfile)
end


l = length(ARGS)
allowed = ["all", "cc", "cpdag", "pdag"]
if l < 1 || (l == 1 && !(ARGS[1] in allowed)) || (l > 1 && !all(a -> a in setdiff(allowed, ["all"]), ARGS))
    @error "Run this file via 'julia $PROGRAM_FILE <TYPE>' with <TYPE> being 'all' or one or more (separated by spaces) of $(join(setdiff(allowed, ["all"]), ", "))."
    exit()
end

types = ARGS[1] == "all" ? setdiff(allowed, ["all"]) : ARGS
for t in types
    if t == "cc"
        run_eval("instances/cc/",    "cc.csv",    true,  CC,    false, 60.0, true)
    elseif t == "cpdag"
        run_eval("instances/cpdag/", "cpdag.csv", false, CPDAG, false, 60.0, true)
    elseif t == "pdag"
        run_eval("instances/pdag/",  "pdag.csv",  false, PDAG,  false, 60.0, true)
    else
        @error "Unsupported input type: $t"
    end
end
