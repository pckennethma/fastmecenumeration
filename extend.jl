using Graphs

mutable struct ExtendedGraph
    g::SimpleDiGraph
    alpha::Vector{Int}
    beta::Vector{Int}
    deltaplus_dir::Vector{Int}
    deltaplus_undir::Vector{Int}
    deltaminus_dir::Vector{Int}
    deltaminus_undir::Vector{Int}
end

"""
    setup(G)

Set up the internal data structure for the extension algorithm.
"""
function setup(G)
    n = nv(G)
    return ExtendedGraph(SimpleDiGraph(n), fill(0, n), fill(0, n), fill(0, n), fill(0, n), fill(0, n), fill(0, n))
end

"""
    initialise!(EG, G)

Initialise the internal data structure for the extension algorithm.
"""
function initialise!(EG, G)
    done = Set()
    for v in 1:nv(G)
        for adj in all_neighbors(G, v)
            adj < v || continue
            is_ingoing  = has_edge(G, adj, v)
            is_outgoing = has_edge(G, v, adj)
            if is_ingoing && is_outgoing # edge is undirected
                if !((v, adj) in done)
                    EG.deltaplus_undir[v]    += 1
                    EG.deltaminus_undir[adj] += 1
                    EG.deltaplus_undir[adj]  += 1
                    EG.deltaminus_undir[v]   += 1
                    add_edge!(EG.g, v, adj)
                    add_edge!(EG.g, adj, v)
                    update_alphabeta!(EG, v, adj, +1, false)
                    push!(done, (adj, v))
                end
            else # edge is directed
                if is_ingoing
                    EG.deltaplus_dir[adj] += 1
                    EG.deltaminus_dir[v]  += 1
                    add_edge!(EG.g, adj, v)
                    update_alphabeta!(EG, adj, v, +1, true)
                elseif is_outgoing
                    EG.deltaplus_dir[v]    += 1
                    EG.deltaminus_dir[adj] += 1
                    add_edge!(EG.g, v, adj)
                    update_alphabeta!(EG, v, adj, +1, true)
                end
            end
        end
    end
end

"""
    isps(G, s)

Check whether `s` is a potential sink in `G`.
"""
function isps(G, s)
    return G.deltaplus_dir[s] == 0 &&
        G.beta[s] == G.deltaplus_undir[s] * G.deltaminus_dir[s] &&
        G.alpha[s] == binomial(G.deltaplus_undir[s], 2)
end

"""
    listps(G)

List all potential sinks in `G`.
"""
function listps(G)
    return [v for v in vertices(G.g) if isps(G, v)]
end

"""
    update_alphabeta!(G, u, v, val, is_uv_dir)

Update the alpha- and beta-values of the data structure.
"""
function update_alphabeta!(G, u, v, val, is_uv_dir)
    for x in all_neighbors(G.g, u)
        has_edge(G.g, x, v) || has_edge(G.g, v, x) || continue

        is_ux_undir = has_edge(G.g, u, x) && has_edge(G.g, x, u)
        is_vx_undir = has_edge(G.g, v, x) && has_edge(G.g, x, v)

        !is_uv_dir && is_ux_undir && (G.alpha[u] += val)
        !is_uv_dir && !has_edge(G.g, u, x) && has_edge(G.g, x, u) && (G.beta[u] += val)

        !is_uv_dir && is_vx_undir && (G.alpha[v] += val)
        is_uv_dir  && is_vx_undir && (G.beta[v] += val)
        !is_uv_dir && has_edge(G.g, x, v) && !has_edge(G.g, v, x) && (G.beta[v] += val)

        is_ux_undir && is_vx_undir && (G.alpha[x] += val)
        is_vx_undir && has_edge(G.g, u, x) && !has_edge(G.g, x, u) && (G.beta[x] += val)
        is_ux_undir && !has_edge(G.g, x, v) && has_edge(G.g, v, x) && (G.beta[x] += val)
    end
end

"""
    removeps!(G, s)

Remove the potential sink `s` from `G`, update the data structure and
return a list of previous neighbors of `s` that became a potential sink
after the removal of `s`. 
"""
function removeps!(G, s)
    oldnbrs = copy(all_neighbors(G.g, s))

    for i in copy(inneighbors(G.g, s))
        !has_edge(G.g, s, i) || continue

        # As s is a sink, all outneighbors are undirected
        for u in outneighbors(G.g, s)
            if has_edge(G.g, i, u)
                has_edge(G.g, u, i) ? (G.alpha[u] += -1) : (G.beta[u] += -1)
            end
        end

        rem_edge!(G.g, i, s)
        G.deltaplus_dir[i]  -= 1
        G.deltaminus_dir[s] -= 1
    end

    for u in copy(outneighbors(G.g, s))
        G.deltaplus_undir[s]  -= 1
        G.deltaminus_undir[u] -= 1
        G.deltaplus_undir[u]  -= 1
        G.deltaminus_undir[s] -= 1
        update_alphabeta!(G, s, u, -1, false)
        rem_edge!(G.g, s, u)
        rem_edge!(G.g, u, s)
    end

    return [v for v in oldnbrs if isps(G, v)]
end

"""
    isextendable(G)

Check whether `G` admits a consistent DAG extension.
"""
function isextendable(G)
    G_ = setup(G)
    initialise!(G_, G)

    ps = listps(G_)
    while !isempty(ps)
        s = pop!(ps)
        push!(ps, removeps!(G_, s)...)
    end

    return isempty(edges(G_.g))
end

"""
    extend(G)

Compute a consistent DAG extension of `G`, if it exists. If no consistent
extension exists, an empty graph (`SimpleDiGraph(0)`) is returned.
"""
function extend(G)
    G_ = setup(G)
    initialise!(G_, G)
    D = copy(G)

    ps = listps(G_)
    while !isempty(ps)
        s = pop!(ps)
        for v in outneighbors(G_.g, s)
            rem_edge!(D, s, v)
        end
        push!(ps, removeps!(G_, s)...)
    end

    return isempty(edges(G_.g)) ? D : SimpleDiGraph(0)
end