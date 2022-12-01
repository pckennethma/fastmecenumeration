using Plots
using Statistics
using KernelDensity

"""
    plotbenchmarks(benchfile)

Create a simple plot for the measurements stored in the file `benchfile`.
"""
function plotbenchmarks(benchfile)
    times = Dict()
    open(benchfile, "r") do io
        readline(io) # Remove header
        for l in readlines(io)
            f, n, m, mdir, mundir, sz, c, algo, tmin, tmax, tmean, tstd = split(l, ",")
            tmean == "skipped" && continue
            f = replace(f, r"-\d{1,2}.gr" => "")
            !haskey(times, algo) && (times[algo] = Dict())
            !haskey(times[algo], f) && (times[algo][f] = [])
            push!(times[algo][f], tmean)
        end
    end
    isempty(times) && return
    x = [i for i in sort(collect(keys(collect(values(times))[1])))]
    ys = []
    label = [k for k in sort(collect(keys(times)))]
    for (k, v) in sort(collect(times), by = x -> x[1])
        y = []
        for t in sort(collect(v), by = x -> x[1])
            push!(y, mean(parse.(Float64, t.second)))
        end
        push!(ys, y)
    end
    y = Matrix{AbstractFloat}(undef, length(ys[1]), length(ys))
    for i in eachindex(ys)
        y[:, i] = ys[i]
    end
    p = plot(x, y, label = reshape(label, 1, length(label)), lw = 2)
    savefig(p, replace(benchfile, ".csv" => ".pdf"))
end

"""
    combine_instances(create_plot=false)

Combine measurements from multiple runs into a single measurement.
Set `create_plot` to `true` to create a kernel density based plot afterwards.
Set `count` to `true` to count how many delays are below `k*mean` for different `k`.
"""
function combine_instances(create_plot=false, count=false)
    for s in ["cc", "cpdag", "pdag"]
        outf = string(@__DIR__, "/results/", "$s-combined.csv")
        fexists = isfile(outf)
        open(outf, "a") do io
            !fexists && write(io, "instances,algo,min,max,mean,median,std\n")
            a1 = (s == "pdag" ? "pdag_enumerate" : "cpdag_enumerate")
            for a in [a1, "enumerate_meek", "dfs_enumerate", "chickering_enumerate"]
                for n in [16, 32, 64, 128, 256, 512, 1024]
                    for d in [3, round(Int, log2(n))]
                        m = d*n
                        for sf in (s == "cc" ? [""] : ["ba", "er"])
                            pattern = "$a-$s-$n-$m-$sf"
                            regex = Regex("$a-$s-$n-$m-(\\d-)?$sf")
                            combined_data = Float64[]
                            for f in readdir(string(@__DIR__, "/results/"))
                                (startswith(f, regex) && endswith(f, ".csv")) || continue
                                open(string(@__DIR__, "/results/", f), "r") do iof
                                    readline(iof) # Remove header
                                    data = parse.(Float64, map(l -> split(l, ",")[2], readlines(iof)))
                                    push!(combined_data, data...)
                                end
                            end
                            isempty(combined_data) && continue
                            create_plot && plot_combined(combined_data, pattern)
                            count && count_delays(combined_data, pattern)
                            mind, maxd = minimum(combined_data), maximum(combined_data)
                            meand, mediand = mean(combined_data), median(combined_data)
                            stddev = std(combined_data)
                            write(io, "$pattern*,$a,$mind,$maxd,$meand,$mediand,$stddev\n")
                        end
                    end
                end
            end
        end
    end
end

"""
    plot_combined(combined_data, pattern)

Create a kernel density based plot for measurements in `combined_data`.
The plots are saved in a file named according to `pattern`.
"""
function plot_combined(combined_data, pattern)
    kd = kde(combined_data, npoints=1048576) #2^20
    p = plot(kd.x, kd.density, legend=false)
    savefig(p, string(@__DIR__, "/results/", "$pattern*.pdf"))
end

"""
    combine_delays()

Combine delays from multiple runs into a single file.
"""
function combine_delays()
    resdir = string(@__DIR__, "/results/")
    for s in ["cc", "cpdag", "pdag"]
        a1 = (s == "pdag" ? "pdag_enumerate" : "cpdag_enumerate")
        for a in [a1, "enumerate_meek", "dfs_enumerate", "chickering_enumerate"]
            for n in [16, 32, 64, 128, 256, 512, 1024]
                for d in [3, round(Int, log2(n))]
                    m = d*n
                    for sf in (s == "cc" ? [""] : ["ba", "er"])
                        pattern = "$a-$s-$n-$m-$sf"
                        outf = string(@__DIR__, "/results/", "$pattern*.csv")
                        for f in readdir(resdir)
                            fexists = isfile(outf)
                            (startswith(f, pattern) && endswith(f, ".csv")) || continue
                            open(string(resdir, f), "r") do iof
                                readline(iof) # Remove header
                                open(outf, "a") do io
                                    !fexists && write(io, "number,delay\n")
                                    for line in readlines(iof)
                                        write(io, string(line, "\n"))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

"""
    count_delays(combined_data, pattern)

Create an overview of how many delays in `combined_data` are less or equal
than `k*mean` for different `k`.
The results are saved in a file named according to `pattern`.
"""
function count_delays(combined_data, pattern)
    n = length(combined_data)
    meanval = mean(combined_data)
    counts = Dict()
    open(string(@__DIR__, "/results/", "counts-$pattern*.csv"), "a") do io
        print(io, "mean,$meanval,\n")
        for k in [1, 2, 3, 5, 7, 10]
            counts[k] = 0
            for delay in combined_data
                delay <= k*meanval && (counts[k] += 1)
            end
            percentage = round(counts[k]/n*100, digits=2)
            print(io, "<=$(lpad(k, 2))*mean,$(counts[k])/$n,$percentage%\n")
        end
    end
end

"""
    combine_counts()

Combine counts of multiple parameters into a single file.
"""
function combine_counts()
    resdir = string(@__DIR__, "/results/")
    combined_counts = Dict()

    for s in ["cc", "cpdag", "pdag"]
        a1 = (s == "pdag" ? "pdag_enumerate" : "cpdag_enumerate")
        for a in [a1, "enumerate_meek", "dfs_enumerate", "chickering_enumerate"]
            a_ = a == "cpdag_enumerate" || a == "pdag_enumerate" ? "new_enumerate" : a
            !haskey(combined_counts, a_) && (combined_counts[a_] = Dict())
            for n in [16, 32, 64, 128, 256, 512, 1024]
                for d in [3, round(Int, log2(n))]
                    m = d*n
                    for sf in (s == "cc" ? [""] : ["ba", "er"])
                        pattern = "counts-$a-$s-$n-$m-$sf"
                        for f in readdir(resdir)
                            (startswith(f, pattern) && endswith(f, ".csv")) || continue
                            !haskey(combined_counts[a_], s) && (combined_counts[a_][s] = Dict())
                            open(string(resdir, f), "r") do io
                                readline(io) # remove first line
                                for line in readlines(io)
                                    key, count, _ = split(line, ",")
                                    k = parse(Int, match(r"\d(\d)?", key).match)
                                    c, n = parse.(Int, split(count, "/"))
                                    !haskey(combined_counts[a_][s], k) && (combined_counts[a_][s][k] = [0, 0])
                                    combined_counts[a_][s][k] += [c, n]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    for (a, sd) in combined_counts
        allcounts = Dict()
        for (s, kd) in sd
            for (k, counts) in kd
                !haskey(allcounts, k) && (allcounts[k] = [0, 0])
                allcounts[k] += counts
            end
        end
        combined_counts[a]["all"] = allcounts
    end

    outf = string(@__DIR__, "/results/", "combined-counts.csv")
    open(outf, "a") do io
        write(io, "algorithm,scenario,k,perc\n")
        for (a, sd) in combined_counts
            for (s, kd) in sd
                for (k, counts) in kd
                    c, n = counts
                    percentage = round(c/n*100, digits=2)
                    write(io, "$a,$s,$k,$percentage\n")
                end
            end
        end
    end
end

### How to use
# Start Julia via `julia -i plots.jl`
# 1. Run `combine_instances(false, true)` to combine multiple runs for
# identical parameters into a single measurement and at the same time
# count the number of delays that are less or equal than `k*mean` for
# different `k`. Set the second parameter to `false` if the counting
# part should be skipped.
# 2. Run `combine_counts()` (only if the counting part was not skipped
# before) to combine the number of delays that are less or equal than
# `k*mean` from different inputs into a single average for the table
# in the appendix.
# Note: `combine_delays()` combines delays from multiple runs for
# identical parameters into a single file and is not needed anymore,
# as the delays are now directly counted on the combined data in
# `combine_instances()`.