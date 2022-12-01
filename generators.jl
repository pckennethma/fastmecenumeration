using Graphs
using Random
using DataFrames
using Distributions
using CausalInference

@isdefined(mcs) || include("counting.jl")
@isdefined(ischordal) || include("utils.jl")
@isdefined(meek!) || include("meek.jl")

"""
    gencc(n, m, seed=123)

Generate an undirected chordal component with `n` vertices and `m` edges.
The procedure works as follows:
1. Generate a Prüfer sequence (https://en.wikipedia.org/wiki/Prüfer_sequence)
2. Convert the Prüfer sequence into a tree
3. Add any edge, check chordality, keep the edge if the graph is still chordal, and repeat until there are `m` edges in the graph

"""
function gencc(n, m, seed=123)
    Random.seed!(seed)
    G = SimpleDiGraph(n)

    # 1. Prüfer sequence of n-2 numbers, each in [1, n-2]
    ps = rand(1:n-2, n-2)

    # 2. Convert the Prüfer sequence into a tree
    deg = fill(1, n)
    for i in ps
        deg[i] += 1
    end

    for i in ps
        for j = 1:n
            if deg[j] == 1
                add_edge!(G, i, j)
                add_edge!(G, j, i)
                deg[i] -= 1
                deg[j] -= 1
                break
            end
        end
    end

    u, v = 0, 0
    for i in 1:n
        if deg[i] == 1
            if u == 0
                u = i
            else
                v = i
                break
            end
        end
    end
    add_edge!(G, u, v)
    add_edge!(G, v, u)
    deg[u] -= 1
    deg[v] -= 1

    # 3. Add remaining edges
    ecount = convert(Int, ne(G) / 2)
    while ecount < m
        u = rand(1:n)
        v = rand(1:n)
        if u != v && !has_edge(G, u, v) && !has_edge(G, v, u)
            add_edge!(G, u, v)
            add_edge!(G, v, u)
            if ischordal(G)
                ecount += 1
            else
                rem_edge!(G, u, v)
                rem_edge!(G, v, u)
            end
        end
    end

    return G
end

"""
    gencpdag(n, m, sf=false, seed=123)

Generate a CPDAG with `n` vertices and `m` edges.
The procedure works as follows:
1. Generate a random DAG D. Set `sf` to `true` to generate a scale-free DAG.
2. Compute the skeleton of D and orient v-structures as in D
3. Apply Meek's rules on the obtained graph

"""
function gencpdag(n, m, sf=false, seed=123)
    Random.seed!(seed)
    D = sf ? genscalefree(n, m, seed) : gendag(n, m, seed)
    S = skeleton(D)
    orient_vstructs!(S, D)
    meek!(S)
    return S
end

"""
    gendag(n, m, seed=123)

Generate a random DAG with `n` vertices and `m` edges by adding
edges randomly until there are `m` edges in the graph. Whenever
an edge is added, check for cycles and remove the edge again if
a cycle has emerged.
"""
function gendag(n, m, seed=123)
    Random.seed!(seed)
    G = SimpleDiGraph(n)
    ecount = 0

    while ecount < m
        u = rand(1:n)
        v = rand(1:n)
        (u != v && !has_edge(G, u, v)) || continue

        add_edge!(G, u, v)

        if !is_cyclic(G)
            ecount += 1
        else
            rem_edge!(G, u, v)
        end
    end

    return G
end

"""
    skeleton(G)

Compute the skeleton of `G` (i.e., the underlying graph with
all edges being undirected).
"""
function skeleton(G)
    S = copy(G)
    for e in edges(G)
        !has_edge(S, e.dst, e.src) && add_edge!(S, e.dst, e.src)
    end
    return S
end

"""
    orient_vstructs!(S, D)

Given a skeleton `S` of a DAG `D`, orient the v-structures of `D` in `S`.
"""
function orient_vstructs!(S, D)
    for u in vertices(D)
        for v in vertices(D)
            (u != v && has_edge(D, u, v) && !has_edge(D, v, u)) || continue
            for w in vertices(D)
                (u != w && v != w && has_edge(D, w, v) && !has_edge(D, v, w)) || continue
                if !has_edge(D, u, w) && !has_edge(D, w, u)
                    rem_edge!(S, v, u)
                    rem_edge!(S, v, w)
                end
            end
        end
    end
end

"""
    genpdag(n, m, s=1000, alpha=0.01, sf=false, seed=123)

Generate a PDAG with `n` vertices and `m` edges.
The procedure works as follows:
1. Generate a random DAG. Set `sf` to `true` to generate a scale-free DAG.
2. Assign weights and sample data from the DAG (`s` is the number of samples)
3. Run the PC algorithm on the data and return the result (`alpha` is the parameter for the PC algorithm)

"""
function genpdag(n, m, s=1000, alpha=0.01, sf=false, seed=123)
    Random.seed!(seed)
    D = sf ? genscalefree(n, m, seed) : gendag(n, m, seed)
    data = sampledata(assignweights(D, seed), s, topological_sort_by_dfs(D), seed)
    G = pcalg(DataFrame(data, :auto), alpha, gausscitest)
    return G
end

"""
    assignweights(G, seed=123)

Assign weights to each edge in `G`.
"""
function assignweights(G, seed=123)
    Random.seed!(seed)
    weights = [Dict() for _ = 1:nv(G)]
    ud = Uniform(1, 2)
    for a in vertices(G), b in outneighbors(G, a)
        weights[a][b] = rand(ud)
    end
    return weights
end

"""
    sampledata(wg, s, ts, seed=123)

Sample `s` data samples from weights `wg` and topological sorting `ts`.
"""
function sampledata(wg, s, ts, seed=123)
    Random.seed!(seed)
    n = size(wg, 1)
    nd = Normal()
    dt = rand(nd, s, n)
    its = zeros(Int64, n)
    for i in 1:n
        its[ts[i]] = i
    end
    # assumes ts 1..n
    for a in its, b in wg[a]
        dt[1:s, b.first] += dt[1:s, a] * b.second
    end
    return dt
end

"""
    genscalefree(n, m, seed=123)

Generate a scale-free DAG (degree distribution following a power law)
with `n` vertices and `m` edges.
"""
function genscalefree(n, m, seed=123)
    Random.seed!(seed)
    G = SimpleDiGraph(n)
    k = convert(Int, floor(1/2 * (n-sqrt(max(n^2 - 4*m, 0)))))

    # Use SimpleDiGraph instead of SimpleGraph
    for e in edges(barabasi_albert(n, k, seed=seed))
        add_edge!(G, e.src, e.dst)
        add_edge!(G, e.dst, e.src)
    end

    # Add random edges if there are not exactly m edges
    # (might happen as we have to specify k instead of m
    # in the Barabasi-Albert model)
    ecount = convert(Int, ne(G) / 2)
    while ecount < m
        u = rand(1:n)
        v = rand(1:n)
        if u != v && !has_edge(G, u, v) && !has_edge(G, v, u)
            add_edge!(G, u, v)
            add_edge!(G, v, u)
            ecount += 1
        end
    end

    ts = randperm(n) # topological sorting
    for a in 1:n, b in inneighbors(G, a)
        ts[b] < ts[a] && rem_edge!(G, a, b)
    end

    @assert !is_cyclic(G) && isextendable(G)
    return G
end

"""
    genpdag_from_cpdag(n, m, k, sf=false, seed=123)

Generate a PDAG with `n` vertices and `m` edges.
The procedure works as follows:
1. Generate a CPDAG based on a random DAG (or on a scale-free DAG if `sf` is `true`)
2. Orient `k` edges randomly one by one and propagate changes via Meek's rules

If there are less than `k` edges oriented after `3*m` tries, the procedure
is cancelled and the current graph is returned.
"""
function genpdag_from_cpdag(n, m, k, sf=false, seed=123)
    Random.seed!(seed)
    G = gencpdag(n, m, sf, seed)

    c = 0
    e = collect(edges(G))
    while k > 0 && c < 3*m
        idx = rand(1:length(e))
        u, v = e[idx].src, e[idx].dst

        if has_edge(G, u, v) && has_edge(G, v, u)
            rem_edge!(G, v, u)
            meek!(G)
            k -= 1
        end

        c += 1
    end

    return G
end
