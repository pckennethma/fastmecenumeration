using Graphs

@isdefined(isextendable) || include("extend.jl")
@isdefined(Measurement) || include("utils.jl")

"""
    isadjacent(G, u, v)

Check whether vertices `u` and `v` are adjacent in `G`.
"""
function isadjacent(G, u, v)
    return has_edge(G, u, v) || has_edge(G, v, u)
end

"""
    meek!(G)

Apply the four Meek rules exhaustively to `G`.
"""
function meek!(G)
    while true
        changed = false

        # Meek rule 1
        # a -> b - c => b -> c
        for b in vertices(G), c in outneighbors(G, b)
            has_edge(G, c, b) || continue # edge b - c must be undirected
            for a in inneighbors(G, b)
                !has_edge(G, b, a) || continue # edge a -> b must be directed
                if a != b && a != c && !isadjacent(G, a, c)
                    rem_edge!(G, c, b)
                    changed = true
                    break
                end
            end
        end

        # Meek rule 2
        # a -> b -> c and a - c => a -> c
        for a in vertices(G), c in outneighbors(G, a)
            has_edge(G, c, a) || continue # edge a - c must be undirected
            for b in intersect(outneighbors(G, a), inneighbors(G, c))
                # edges a -> b and b -> c must be directed
                if !has_edge(G, b, a) && !has_edge(G, c, b)
                    rem_edge!(G, c, a)
                    changed = true
                    break
                end
            end
        end

        # Meek rule 3
        # a - d -> c <- b with a - b and a - c => a -> c
        for d in vertices(G), c in outneighbors(G, d), b in inneighbors(G, c)
            # edges d -> c and c <- b must be directed
            (!has_edge(G, c, d) && !has_edge(G, c, b)) || continue
            if b != d && !isadjacent(G, b, d)
                for a in intersect(outneighbors(G, b), outneighbors(G, c), outneighbors(G, d))
                    # edges a - b, a - c, and a - d must be undirected
                    (has_edge(G, a, b) && has_edge(G, a, c) && has_edge(G, a, d)) || continue
                    rem_edge!(G, c, a)
                    changed = true
                end
            end
        end

        # Meek rule 4
        # d -> c -> b with a - b, a - c, and a - d => a -> b
        for c in vertices(G), b in outneighbors(G, c), d in inneighbors(G, c)
            # edges d -> c and c -> b must be directed
            (!has_edge(G, c, d) && !has_edge(G, b, c)) || continue
            if b != d && !isadjacent(G, b, d)
                for a in intersect(outneighbors(G, b), outneighbors(G, c), outneighbors(G, d))
                    # edges a - b, a - c, and a - d must be undirected
                    (has_edge(G, a, b) && has_edge(G, a, c) && has_edge(G, a, d)) || continue
                    rem_edge!(G, b, a)
                    changed = true
                end
            end
        end

        !changed && break
    end
end

mutable struct Counter
    c::BigInt
    m::Measurement
end

"""
    enumerate_meek(G, m, skip=false)

Enumerate Markov equivalent DAGs represented by `G` by orienting an edge
in both directions, applying Meek's rules, and recursing for both directions.
`m` is a helper for taking measurements.
Set `skip` to `true` if the input graph is an undirected chordal graph
or a CPDAG, allowing to skip the check whether the input is extendable
and to skip the application of Meek's rules in the first step.

## References
- Algorithm 2 (MEEK-ENUM) in Appendix A.1 (Enumerating Markov Equivalent
DAGs based on Meek's Rules).
"""
function enumerate_meek(G, m, skip=false)
    skip || isextendable(G) || return 0
    m.last = time_ns()
    h = Counter(0, m)
    recenum_meek(G, 1, skip, h)
    return h.c
end

"""
    recenum_meek(G, h)

Enumerate Markov equivalent DAGs represented by `G` recursively by
going recursively through the following steps:

1. Orient an undirected edge
2. Apply Meek's rules
3. Recurse separately with each direction of the selected undirected edge

`lastidx` is the lastly viewed vertex when searching for an undirected edge,
`skipmr` is a flag to skip the application of Meek's rules in the first call,
and `h` is a counter to count the number of outputs.
"""
function recenum_meek(G, lastidx, skipmr, h)
    if !skipmr
        G_ = copy(G)
        meek!(G_)
    else
        G_ = G
    end

    u, v = -1, -1
    for i = lastidx:nv(G_)
        for j = i+1:nv(G_)
            # Select the first undirected edge found
            if has_edge(G_, i, j) && has_edge(G_, j, i)
                u, v = i, j
                lastidx = i
                break
            end
        end
        u != -1 && break
    end

    if u == -1 # No undirected edges left
        h.c += 1
        addmeasurement!(h.m)
        return
    end

    # Recurse for the first orientation
    rem_edge!(G_, u, v)
    recenum_meek(G_, lastidx, false, h)
    add_edge!(G_, u, v)

    # Recurse for the second orientation
    rem_edge!(G_, v, u)
    recenum_meek(G_, lastidx, false, h)
    add_edge!(G_, v, u)
end