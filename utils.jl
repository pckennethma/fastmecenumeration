using Graphs
using Random

struct TimeoutException <: Exception end

mutable struct Measurement
    start::UInt
    last::UInt
    min::AbstractFloat
    max::AbstractFloat
    mean::AbstractFloat
    std::AbstractFloat
    n::Int
    timeout::Int
    f::String
end

"""
    readgraph(file = stdin, undirected = false)

Read a graph from the standard input or a given file and return a directed simple graph.
In case undirected graphs are read (as in the toy examples we give
here), the argument undirected=true may be passed. Then, each edge can be given only once.

In practice, for partially directed graphs like CPDAGs, an undirected
edge is expected to be encoded as two directed edges.

# Examples
```julia-repl
julia> readgraph(stdin, true)
6 11

1 2
1 3
2 3
2 4
2 5
2 6
3 4
3 5
3 6
4 5
5 6
{6, 22} directed simple Int64 graph
```
"""
function readgraph(file = stdin, undirected = false)
    if file != stdin
        infile = open(file, "r")
    else
        infile = stdin
    end
    (n, m) = parse.(Int, split(readline(infile)))
    readline(infile)
    G = SimpleDiGraph(n)
    for i = 1:m
        (a, b) = parse.(Int, split(readline(infile)))
        add_edge!(G, a, b)
        undirected && add_edge!(G, b, a)
    end
    nedge = ne(G)
    @info "Read graph with $n vertices and $nedge edges"
    if file != stdin
        close(infile)
    end
    return G
end

"""
    ischordal(g)

Return true if the given graph is chordal
"""
function ischordal(G)
    mcsorder, invmcsorder, _ = mcs(G, Set())
    
    n = length(mcsorder)
    
    f = zeros(Int, n)
    index = zeros(Int, n)
    for i=n:-1:1
        w = mcsorder[i]
        f[w] = w
        index[w] = i
        for v in neighbors(G, w)
            if invmcsorder[v] > i
                index[v] = i
                if f[v] == v
                    f[v] = w
                end
            end
        end
        for v in neighbors(G, w)
            if invmcsorder[v] > i
                if index[f[v]] > i
                    return false
                end
            end
        end
    end
    return true
end

"""
    nanosec2millisec(t)

Convert `t` from nanoseconds to milliseconds.
"""
function nanosec2millisec(t)
    # Nano /1000 -> Micro /1000 -> Milli /1000 -> Second
    return t / 1000 / 1000
end

"""
    nanosec2sec(t)

Convert `t` from nanoseconds to seconds.
"""
function nanosec2sec(t)
    return t / 1000 / 1000 / 1000
end

"""
    ndiredges(G)

Return the number of directed edges in `G`.
"""
function ndiredges(G)
    return length(filter(e -> !has_edge(G, e.dst, e.src), collect(edges(G))))
end

"""
    nundiredges(G)

Return the number of undirected edges in `G`.
"""
function nundiredges(G)
    return convert(Int, length(filter(e -> has_edge(G, e.dst, e.src), collect(edges(G)))) / 2)
end

"""
    addmeasurement!(m)

Compute the minimum delay, maximum delay, mean delay, and standard deviation
on the fly.
`ts` is the timestamp of the current sample, `m.n` the number of the current
sample, `m.last` the timestamp of the last sample, `m.min` the minimum delay,
`m.max` the maximum delay, `m.mean` the mean delay, and `m.std` the squared
distance from the mean (to obtain the actual standard deviation, divide by
`n-1` and take the square root).
`m.start` contains the start time and `m.timeout` the number of seconds after
which the execution should be ended (regardless of its current state).
"""
function addmeasurement!(m)
    m.n += 1
    ts = time_ns()
    elapsed = ts - m.last
    m.min = m.min < elapsed ? m.min : elapsed
    m.max = m.max > elapsed ? m.max : elapsed
    delta = elapsed - m.mean
    m.mean = m.mean + (elapsed - m.mean) / m.n
    m.std = m.n < 2 ? 0 : m.std + delta * (elapsed - m.mean)

    if !isempty(m.f)
        fexists = isfile(m.f)
        open(m.f, "a") do io
            !fexists && write(io, "number,delay\n")
            write(io, "$(m.n),$(nanosec2millisec(elapsed))\n")
        end
    end

    nanosec2sec(ts - m.start) >= m.timeout && throw(TimeoutException())
    m.last = time_ns()
end


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
