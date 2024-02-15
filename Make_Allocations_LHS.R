
set.seed(3011991)

#
# LHS for drawing voccine coverage within each population
#

pops = c(0.151, 0.144, 0.148, 0.164, 0.145, 0.110, 0.085, 0.054)

library(lhs)
library(tictoc)

lims = c(rep(10000, 6), 7500, 5000, 2000)


v = round(randomLHS(1000000, 8), 2); v.cov = v
for (i in 1:ncol(v)) { v[ , i] = v[ , i] * pops[i] }
tot.vax = apply(v, 1, sum)
#v.cov = v.cov[which(floor(tot.vax*10)/10 == n.vax), ]
#v = v[which(floor(tot.vax*10)/10 == n.vax), ]
v.cov = v.cov[which(tot.vax > n.vax & tot.vax < (n.vax+ifelse(n.vax == 0.1, 0.01, 0.05))), ]
v = v[which(tot.vax > n.vax & tot.vax < (n.vax+ifelse(n.vax == 0.1, 0.01, 0.05))), ]
for (i in 1:nrow(v)) { v.cov[i, ] = v.cov[i, ] * n.vax / sum(v[i, ]); v[i, ] = v[i, ] * n.vax / sum(v[i, ]) }

table(floor(tot.vax*10)/10)

if (nrow(v) > 0.2*lims[n.vax*10]) { 
	v = v[1:min(nrow(v), 0.2*lims[n.vax*10]), ]
	v.cov = v.cov[1:min(nrow(v.cov), 0.2*lims[n.vax*10]), ]
}






tic()
# assuring that every age group gets at least 25 observations with 95% coverage (or at max in the case of 10% supply)
if (n.vax > 0.1) { for (a in 1:ncol(v)) {
	cat('*')
	while(sum(v.cov[ , a] > 0.95) < 25) {

		v2 = round(randomLHS(50000, 8), 2); v2[ , a] = round(runif(nrow(v2), 0.97, 1), 2); v.cov2 = v2
#		v2 = round(randomLHS(50000, 8), 2); v.cov2 = v2
		if (n.vax < 0.9) {
			for (i in 1:nrow(v2)) { for (j in c(1:8)[-a]) { v2[i, j] = v2[i, j] * (1 - pops[a]) } } 
		}
		for (i in 1:ncol(v2)) { v2[ , i] = v2[ , i] * pops[i] }
		tot.vax2 = apply(v2, 1, sum)
		
		while(min(tot.vax2) > n.vax) {
			for (i in 1:nrow(v.cov2)) { for (j in c(1:8)[-a]) { v.cov2[i, j] = v.cov2[i, j] * 0.95 } } 
			v2 = v.cov2
			for (i in 1:ncol(v2)) { v2[ , i] = v2[ , i] * pops[i] }
			tot.vax2 = apply(v2, 1, sum)
		}
		
		v.cov2 = v.cov2[which(tot.vax2 > n.vax & tot.vax2 < (n.vax+0.02)), ]
		v2 = v2[which(tot.vax2 > n.vax & tot.vax2 < (n.vax+0.02)), ]
		
		if (!is.null(nrow(v2)) && nrow(v2) > 0) {
			for (i in 1:nrow(v2)) { 
				v.cov2[i, ] = v.cov2[i, ] * n.vax / sum(v2[i, ]) 	
				v2[i, ] = v2[i, ] * n.vax / sum(v2[i, ]) 
			}

			v = rbind(v, v2[1:min(nrow(v2), 25), ]); v.cov = rbind(v.cov, v.cov2[1:min(nrow(v2), 25), ])	
		
		}
	
	}
	
} }

if (n.vax == 0.1) {
for (a in 1:ncol(v)) {
	cat('*')
	while(sum(v.cov[ , a] > min(0.1*0.95/pops[a], 0.95)) < 25) {

		v2 = round(randomLHS(50000, 8), 2); v2[ , a] = round(runif(nrow(v2), min(0.1*0.985/pops[a], 0.985), min(0.1/pops[a], 1)), 2); v.cov2 = v2
#		v2 = round(randomLHS(50000, 8), 2); v.cov2 = v2
		for (i in 1:nrow(v2)) { for (j in c(1:8)[-a]) { v2[i, j] = v2[i, j] * (1 - pops[a]) } } 
		for (i in 1:ncol(v2)) { v2[ , i] = v2[ , i] * pops[i] }
		tot.vax2 = apply(v2, 1, sum)
		
		while(min(tot.vax2) > n.vax) {
			for (i in 1:nrow(v.cov2)) { for (j in c(1:8)[-a]) { v.cov2[i, j] = v.cov2[i, j] * 0.95 } } 
			v2 = v.cov2
			for (i in 1:ncol(v2)) { v2[ , i] = v2[ , i] * pops[i] }
			tot.vax2 = apply(v2, 1, sum)
		}
		
		v.cov2 = v.cov2[which(tot.vax2 > n.vax & tot.vax2 < (n.vax+0.005)), ]
		v2 = v2[which(tot.vax2 > n.vax & tot.vax2 < (n.vax+0.005)), ]
		
		if (!is.null(nrow(v2)) && nrow(v2) > 0) {
			for (i in 1:nrow(v2)) { 
				v.cov2[i, ] = v.cov2[i, ] * n.vax / sum(v2[i, ]) 	
				v2[i, ] = v2[i, ] * n.vax / sum(v2[i, ]) 
			}

			v = rbind(v, v2[1:min(nrow(v2), 25), ]); v.cov = rbind(v.cov, v.cov2[1:min(nrow(v2), 25), ])	
		
		}
	
	}
	
}	
	
	
}
toc()




tic()
while (nrow(v) < 10*lims[10*n.vax]) {
#	cat('*')
	v2 = round(randomLHS(ifelse(n.vax %in% c(0.1, 0.9), 500000, 5000), 8), 2); v.cov2 = v2
	for (i in 1:ncol(v2)) { v2[ , i] = v2[ , i] * pops[i] }
	tot.vax2 = apply(v2, 1, sum)
	v.cov2 = v.cov2[which(floor(tot.vax2*10)/10 == n.vax), ]
	v2 = v2[which(floor(tot.vax2*10)/10 == n.vax), ]
	if (!is.null(nrow(v2)) && nrow(v2) > 0) {
		for (i in 1:nrow(v2)) { 
			v.cov2[i, ] = v.cov2[i, ] * n.vax / sum(v2[i, ]) 
			v2[i, ] = v2[i, ] * n.vax / sum(v2[i, ]) 
		}
	
		v = rbind(v, v2); v.cov = rbind(v.cov, v.cov2)	
	}
	
}
toc()

#v = v[1:lims[10*n.vax], ]
#v.cov = v.cov[1:lims[10*n.vax], ]

allocs = v.cov
for (i in 1:ncol(allocs)) { allocs[ , i] = allocs[ , i] * pops[i] / n.vax }
#allocs = round(allocs, 2)
# for (i in 1:nrow(allocs)) { allocs[i, ] = floor(100*allocs[i, ]) / 100 }

allocs = allocs[which(apply(allocs, 1, sum) == 1), ]
v = v[which(apply(allocs, 1, sum) == 1), ]
v.cov = v.cov[which(apply(allocs, 1, sum) == 1), ]


allocs = allocs[c(1:lims[10*n.vax]), ]
v = v[1:lims[10*n.vax], ]
v.cov = v.cov[1:lims[10*n.vax], ]


rm(v2); rm(v.cov2)

