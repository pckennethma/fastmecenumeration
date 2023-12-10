using Graphs

@isdefined(isextendable) || include("extend.jl")
@isdefined(meek!) || include("meek.jl")
@isdefined(Measurement) || include("utils.jl")
@isdefined(writegraph) || include("utils.jl")

mutable struct PDAGHelper
    parents::Vector{Set{Int64}}
    A::Vector{Set{Int64}}
    invA::Vector{Int64}
    tau::Vector{Int64}
    indeg::Vector{Int64}
    comp::Vector{Int64}
    maxA::Int64
    i::Int64
    count::BigInt
    m::Measurement
end

function pdag_setv(C, H, v)
    H.invA[v] = -H.invA[v] # simple way to mark vertex as visited
    delete!(H.A[H.maxA], v)
    H.tau[H.i] = v         # append v to tau
    H.i += 1

    for w in inneighbors(C, v)
        if H.invA[w] >= 1  # check that w is unvisited
            delete!(H.A[H.invA[w]], w)
            H.invA[w] += 2
            if v in H.parents[w]
                H.indeg[w] -= 1
                H.indeg[w] == 0 && (H.invA[w] += 1)
            end
            push!(H.A[H.invA[w]], w)
        end
    end

    H.maxA += 2
    while H.maxA > 2 && isempty(H.A[H.maxA])
        H.maxA -= 2
    end
end

function pdag_resetv(C, H, v)
    H.invA[v] = -H.invA[v]
    push!(H.A[H.invA[v]], v)
    H.i -= 1
    
    for w in inneighbors(C, v)
        if H.invA[w] >= 1
            delete!(H.A[H.invA[w]], w)
            H.invA[w] -= 2
            if v in H.parents[w]
                H.indeg[w] == 0 && (H.invA[w] -= 1)
                H.indeg[w] += 1
            end
            push!(H.A[H.invA[w]], w)
        end
    end
    
    H.maxA = H.invA[v]
end

function dfs(C, H, v, reachable)
    for w in inneighbors(C, v)
        if (w in H.A[H.maxA]) && !(w in reachable)
            push!(reachable, w)
            dfs(C, H, w, reachable)
        end
    end
end

function construct_DAG(G, H)
    n = nv(G)
    D = copy(G)
    invtau = zeros(Int64, n)
    for i = 1:n
        invtau[H.tau[i]] = i
    end
    D.ne = 0
    for i = 1:n
        filter!(j->H.comp[i] != H.comp[j] || invtau[j] > invtau[i], D.fadjlist[i])
        filter!(j->H.comp[i] != H.comp[j] || invtau[j] < invtau[i], D.badjlist[i])
        D.ne += length(D.fadjlist[i])
    end

    return D
end

function rec_pdag_enum(G, C, H, output_folder="")
    if H.count >= 1024
        @warn "More than 1024 DAGs found, aborting."
        return
    end

    if H.i > nv(G)
        D = construct_DAG(G, H)
        # do something, e.g. print
        index = H.count
        if output_folder != ""
            writegraph(D, joinpath(output_folder, "dag-$index.gr"), false)
        end
        addmeasurement!(H.m)
        H.count += 1
        return
    end

    S = H.A[H.maxA]
    v = first(S)
    pdag_setv(C, H, v)
    rec_pdag_enum(G, C, H, output_folder)
    pdag_resetv(C, H, v)

    reachable = Set{Int64}()
    dfs(C, H, v, reachable)
    
    for x in reachable
        x == v && continue
        pdag_setv(C, H, x)
        rec_pdag_enum(G, C, H, output_folder)
        pdag_resetv(C, H, x)
    end
end

"""
    pdag_enumerate(G, m, output_folder="")

Enumerate all Markov equivalent DAGs of a PDAG `G` with linear-time delay
by first converting `G` into an MPDAG and then making use of the Maximum
Cardinality Search (MCS) algorithm, with a slight modification such that
only vertices with no univisted parent are chosen.

## References
- Theorem 5 in Section 4 (PDAGs and MPDAGs).
- Algorithm 7 (based on MCS-ENUM for buckets) in Appendix B.2.
"""
function pdag_enumerate(G, m, output_folder="")
    # G is orig graph (PDAG), B is the graph containing only buckets, H is helper
    n = nv(G)

    isextendable(G) || return 0
    G = copy(G)
    meek!(G)

    U = copy(G)
    U.ne = 0
    for i = 1:n
        filter!(j->has_edge(G, j, i), U.fadjlist[i])
        filter!(j->has_edge(G, i, j), U.badjlist[i])
        U.ne += length(U.fadjlist[i])
    end

    comp = Vector{Int64}(undef, n)
    components = connected_components(U)
    for i = 1:length(components)
        for x in components[i]
            comp[x] = i
        end
    end

    C = copy(G)
    C.ne = 0
    for i = 1:n
        filter!(j->comp[i] == comp[j], C.fadjlist[i])
        filter!(j->comp[i] == comp[j], C.badjlist[i])
        C.ne += length(C.fadjlist[i])
    end

    parents = [Set{Int64}() for i = 1:n]
    for i = 1:n
        for j in inneighbors(C, i)
            !has_edge(C, i, j) && push!(parents[i], j)
        end
    end
    indeg = [length(parents[i]) for i = 1:n]

    for i = 1:n
        for j in inneighbors(C, i)
            if !has_edge(C, i, j)
                add_edge!(C, i, j)
            end
        end
    end

    A = [Set{Int64}()  for i = 1:2*(n+1)]
    invA = Vector{Int64}(undef, n)
    tau = Vector{Int64}(undef, n)
    for v in vertices(G)
        invA[v] = 1
        indeg[v] == 0 && (invA[v] += 1)
        push!(A[invA[v]], v)
    end
    # has to be 2 else not extendable
    maxA = 2
    i = 1
    count = BigInt(0)

    m.last = time_ns()
    H = PDAGHelper(parents, A, invA, tau, indeg, comp, maxA, i, count, m)


    rec_pdag_enum(G, C, H, output_folder)

    return H.count
end
