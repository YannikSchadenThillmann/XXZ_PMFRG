
logrange(min,max,n,a) = (exp(x)+a for x in LinRange(log(min-a),log(max-a),n))

T_max = 1e+5

dense_range = collect(LinRange(0.0,2.,101))
medium_range = collect(logrange(2.,10.,50,1.7))
sparse_range = collect(logrange(10.,T_max,30,8.8))
test = unique!(append!(dense_range,medium_range,sparse_range))

fig=Figure()
ax = Axis(fig[1,1])
scatter!(ax,test,zeros(length(test)))
xlims!(ax,(1.5,2.5))
xlims!(ax,(1.5,100.5))
display("image/png",fig)