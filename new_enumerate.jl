using Graphs

@isdefined(Measurement) || include("utils.jl")

mutable struct CPDAGHelper
    A::Vector{Set{Int64}}
    invA::Vector{Int64}
    tau::Vector{Int64}
    comp::Vector{Int64}
    maxA::Int64
    i::Int64
    count::BigInt
    m::Measurement
end

function cpdag_setv(C, H, v)
    H.invA[v] = -H.invA[v] # simple way to mark vertex as visited
    delete!(H.A[H.maxA], v)
    H.tau[H.i] = v         # append v to tau
    H.i += 1

    for w in inneighbors(C, v)
        if H.invA[w] >= 1  # check that w is unvisited
            delete!(H.A[H.invA[w]], w)
            H.invA[w] += 1
            push!(H.A[H.invA[w]], w)
        end
    end

    H.maxA += 1
    while H.maxA > 1 && isempty(H.A[H.maxA])
        H.maxA -= 1
    end
end

function cpdag_resetv(C, H, v)
    H.invA[v] = -H.invA[v]
    push!(H.A[H.invA[v]], v)
    H.i -= 1
    
    for w in inneighbors(C, v)
        if H.invA[w] >= 1
            delete!(H.A[H.invA[w]], w)
            H.invA[w] -= 1
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

function rec_cpdag_enum(G, C, H)
    if H.i > nv(C)
        D = construct_DAG(G, H)
        # do something, e.g. print
        addmeasurement!(H.m)
        H.count += 1
        return
    end

    S = H.A[H.maxA]
    v = first(S)
    cpdag_setv(C, H, v)
    rec_cpdag_enum(G, C, H)
    cpdag_resetv(C, H, v)

    reachable = Set{Int64}()
    dfs(C, H, v, reachable)
    
    for x in reachable
        x == v && continue
        cpdag_setv(C, H, x)
        rec_cpdag_enum(G, C, H)
        cpdag_resetv(C, H, x)
    end
end

"""
    cpdag_enumerate(G, m)

Enumerate all Markov equivalent DAGs of a CPDAG `G` with linear-time delay
by making use of the Maximum Cardinality Search (MCS) algorithm.

## References
- Theorem 3 in Section 3 (Enumerating AMOs with Linear Delay).
- Algorithm 5 (EnumMCS) in Appendix B.1.
"""
function cpdag_enumerate(G, m)
    # G is CPDAG, C is graph containing only the chordal components, H is helper
    n = nv(G)

    C = copy(G)
    C.ne = 0
    for i = 1:n
        filter!(j->has_edge(G, j, i), C.fadjlist[i])
        filter!(j->has_edge(G, i, j), C.badjlist[i])
        C.ne += length(C.fadjlist[i])
    end

    comp = Vector{Int64}(undef, n)
    components = connected_components(C)
    for i = 1:length(components)
        for x in components[i]
            comp[x] = i
        end
    end

    A = [Set{Int64}() for i = 1:n+1] # 1:n should suffice
    invA = Vector{Int64}(undef, n)
    tau = Vector{Int64}(undef, n)
    for v in vertices(C)
        invA[v] = 1
        push!(A[1], v)
    end
    maxA = 1
    i = 1
    count = BigInt(0)

    m.last = time_ns()
    H = CPDAGHelper(A, invA, tau, comp, maxA, i, count, m)


    rec_cpdag_enum(G, C, H)

    return H.count
end
