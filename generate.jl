using Graphs
using Random

@isdefined(gencc) || include("generators.jl")

"""
    writegraph(G, f, undir=false)

Save the graph `G` to the file `f`.
Set `undir` to `true` if the graph is fully undirected and edges should not be encoded as two arcs.
"""
function writegraph(G, f, undir=false)
    n = nv(G)
    m = undir ? convert(Int, ne(G) / 2) : ne(G)
    done = Set()
    open(f, "w") do io
        write(io, "$n $m\n\n")
        for u = 1:n, v = 1:n
            if has_edge(G, u, v) && !((u, v) in done)
                write(io, "$u $v\n")
                undir && push!(done, (v, u))
            end
        end
    end
end


allowed = ["cc", "cpdag", "pdag"]
if length(ARGS) != 1 || !(ARGS[1] in allowed)
    @error "Run this file via 'julia $PROGRAM_FILE <TYPE>' with <TYPE> being one of $(join(allowed, ", "))."
    exit()
end

if ARGS[1] == "cc"
    for rep = 1:10
        for n in [16, 32, 64, 128, 256, 512, 1024]
            for d in [3, round(Int, log2(n))]
                m = d*n
                f = string(@__DIR__, "/instances/cc/", "cc-$n-$m-$rep.gr")
                isfile(f) && continue
                G = gencc(n, m, rep)
                writegraph(G, f, true)
            end
        end
    end
elseif ARGS[1] == "cpdag"
    for rep = 1:2
        for n in [2]
            for d in [1]
                m = d*n
                for sf in [true, false]
                    s = sf ? "ba" : "er"
                    f = string(@__DIR__, "/instances/cpdag/", "cpdag-$n-$m-$s-$rep.gr")
                    isfile(f) && continue
                    G = gencpdag(n, m, sf, rep)
                    writegraph(G, f, false)
                end
            end
        end
    end
elseif ARGS[1] == "pdag"
    for rep = 1:10
        Random.seed!(rep)
        for n in [16, 32, 64, 128, 256, 512, 1024]
            for d in [3, round(Int, log2(n))]
                m = d*n
                for sf in [true, false]
                    s = sf ? "ba" : "er"
                    k = d > 3 ? rand(3:7) : d
                    f = string(@__DIR__, "/instances/pdag/", "pdag-$n-$m-$k-$s-$rep.gr")
                    isfile(f) && continue
                    G = genpdag_from_cpdag(n, m, k, sf, rep)
                    writegraph(G, f, false)
                end
            end
        end
    end
else
    @error "Unsupported type: $(ARGS[1])"
end
