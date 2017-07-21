
#list of carbon cycle parameters
p = list()

#time series attributes
p$t0 = 0
p$tf = 6e6
p$dt = 5e4
p$tt = seq(from=p$t0, to=p$tf, by=p$dt )
p$nt = length(p$tt)

#reservoirs
