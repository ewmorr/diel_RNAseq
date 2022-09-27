nls_coefs = read.table("biochem_vars_nls_coefs.txt", header = T)

nls_coefs$phase.dec = nls_coefs$phase - as.integer(nls_coefs$phase)
nls_coefs$time.diff.period = nls_coefs$phase.dec/(2*pi)
nls_coefs$time.diff.hours = nls_coefs$time.diff.period*24
nls_coefs$sin.peak = 10.06 - nls_coefs$time.diff.hours + 6
