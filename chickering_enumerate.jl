using Graphs

@isdefined(extend) || include("extend.jl")
@isdefined(Measurement) || include("utils.jl")

function edgelist(D)
    el = Vector{Edge}()
    for e in edges(D)
        push!(el, e)
    end
    return el
end

function covered_edges(D, undirected)
    n = nv(D)
    ce = Vector{Edge}()
    for i = 1:n
        for j = 1:n
            (!has_edge(D, i, j) || !(Edge(i,j) in undirected))  && continue
            inneighbors(D, i) == filter(x -> x != i, inneighbors(D, j)) && push!(ce, Edge(i, j))
        end
    end
    return ce
end

function output(D, vis, m)
    push!(vis, edgelist(D))
    # do something, e.g. print
    addmeasurement!(m)
end

function rec_chickering_enumerate(D, vis, turnededges, undirected, i, m)
    output(D, vis, m)
    C = covered_edges(D, undirected)
    for e in C
        (e in turnededges) && continue
        x = src(e)
        y = dst(e)
        rem_edge!(D, x, y)
        add_edge!(D, y, x)
        push!(turnededges, Edge(x,y))
        push!(turnededges, Edge(y,x))
        if !(edgelist(D) in vis)
            rec_chickering_enumerate(D, vis, turnededges, undirected, i+1, m)
        end
        rem_edge!(D, y, x)
        add_edge!(D, x, y)
        delete!(turnededges, Edge(x,y))
        delete!(turnededges, Edge(y,x))
    end
end

"""
    chickering_enumerate(G, m)

Enumerate all Markov equivalent DAGs of `G` based on the transformational
characterization by Chickering (1995), that is, by successive single-edge
reversals starting from an arbitrary DAG in the MEC of `G`.

## References
- Algorithm 3 (CHICKERING-ENUM) in Appendix A.2 (Enumerating Markov
Equivalent DAGs based on Chickering's Transformational Characterization).
"""
function chickering_enumerate(G, m)
    undirected = Set{Edge}()
    for e in edges(G)
        x = src(e)
        y = dst(e)
        has_edge(G, y, x) && push!(undirected, e)
    end
    D = extend(G)
    D == SimpleDiGraph(0) && return 0 # G is not extendable
    vis = Set{Vector{Edge}}()
    turnededges = Set{Edge}()
    m.last = time_ns()
    rec_chickering_enumerate(D, vis, turnededges, undirected, 0, m)
    return length(vis)
end
