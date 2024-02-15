


setwd('')

library(Rcpp)
library(tictoc)
library(parallel)


#
# Pull in all needed data/code functions
# 

source('Code/condenseContact.R')
load('Data/VN_Demographics.Rdata')
load('Data/kieshaprem-synthetic-contact-matrices-6e0eebc/output/syntheticcontactmatrices2020/overall/contact_all.rdata')
contact = as.data.frame(contact_all$VNM)
rm(contact_all)


# Parameters for model
load('Params_for_start_full.Rdata') 
best.0 = which(apply(pars[ , -c(1:5)], 1, sum) == max(apply(pars[ , -c(1:5)], 1, sum)))
length(best.0)
best.0 = c(1:50)


#
# Set up model configurations (population(s), R classes)
#

group_labels = c('0s', '10s', '20s', '30s', '40s', '50s', '60s', '70up')
populations = set.pops(1000000, demo$Year, demo$Pop_Total, 
	list(c(0:9), c(10:19), c(20:29), c(30:39), c(40:49), c(50:59), c(60:69), c(70:100)))
load('Data/VN_US_pops.Rdata')
populations = vn.us.pops[[1]]

contact.matrix = condense.contact(contact, 
	list(c(1:2), c(3:4), c(5:6), c(7:8), c(9:10), c(11:12), c(13:14), c(15:16)))
	
num.r = 4

age.bands = c(rep(10, 7), 30)
death.rate.per.1000 = 1
death.share = condense.deathrate(demo$VN_DeathRate_est, demo$Pop_Total/sum(demo$Pop_Total), 
	list(c(0:9), c(10:19), c(20:29), c(30:39), c(40:49), c(50:59), c(60:69), c(70:100))) / 100

starts = c('S' = 0.6, 'I' = 10/1000000, 'II' = 10/1000000, 'R' = round(0.25 / num.r / 3, 3))
starts[1] = 1 - sum(starts[2:3]*3) - (starts[4]*num.r*3)


file.name = 'model_code_multiplegroup.txt'



#
# Flu characteristics
#


flu.hfr = condense.fatality(demo$HFR, demo$Year, demo$Pop_Total, 
	list(c(0:9), c(10:19), c(20:29), c(30:39), c(40:49), c(50:59), c(60:69), c(70:100))) 
flu.ph = condense.fatality(demo$pH, demo$Year, demo$Pop_Total, 
	list(c(0:9), c(10:19), c(20:29), c(30:39), c(40:49), c(50:59), c(60:69), c(70:100)))


tic()
# write model file

source('Code/write_model.R')
toc()



tic()
source('Code/model_code_multiplegroup.txt')
toc()





#
# Making longitidinal beta
#	

load('output/Res_4R8A/betas_for_start_Full.Rdata') 
	

vax_res = list()

for (n.vax in c(seq(0.1, 0.1, by = 0.1))) {


	all.deaths = NULL
	all.cases = NULL
	all.wastes = NULL
	
	cat('\n\n', n.vax, '\n')

	months.vax = c(3); n.days = 31
	while(n.days < 215 * n.vax) {
		months.vax = c(months.vax, max(months.vax)+1)
		n.days = n.days + 30
	}
	months.vax = months.vax %% 12; months.vax[which(months.vax == 0)] = 12; months.vax = sort(months.vax)

	n.days = 31 * sum(c(1, 3, 5, 7, 8, 10, 12) %in% months.vax) + 
		 30 * sum(c(4, 6, 9, 11) %in% months.vax) + 
		 28 * (2 %in% months.vax)


	tic()

	source('~/work/Flu_Vax/Code/Make_Allocations_LHS.R', local = TRUE)
	all_allocs = rbind(rep(0, ncol(all_allocs)), allocs, round(populations/sum(populations), 3))



	for (bb in 1:min(length(best.0), 50)) {

	best = best.0[bb]

	betas = as.numeric(pars[best, 2:4])

	duration = pars[best, 5]

	beta.vec = beta.list[[best]]
	if (length(beta.vec) < 7) { beta.vec[[7]] = 182 }




#
# Making vaccine allocations
#



	tic()
	deaths = rep(NA, nrow(all_allocs))
	wastes = rep(NA, nrow(all_allocs))
	ARs = matrix(NA, nrow = nrow(all_allocs), ncol = 10)

	ve.adj = 0.0 # For attempting different VE values

	res = list()

	res = mclapply(c(1:nrow(all_allocs)), function(x) {
	
		out = sirout_groups(beta_H1 = beta.vec[[1]], beta_B = beta.vec[[2]], beta_H3 = beta.vec[[3]], 
			h_H1 = c(-1.0, flu.ph), h_B = c(-1.0, flu.ph), h_H3 = c(-1.0, flu.ph), 
			d_H1 = c(-1.0, flu.hfr), d_B = c(-1.0, flu.hfr), d_H3 = c(-1.0, flu.hfr), 
			cov = 0, r_i = 1 - (0.5+ve.adj), r_h = 1 - (0.5+ve.adj), r_d = 1 - (0.5+ve.adj), 
			nu = 1/(5/2), rho = 1/(duration/num.r), rho_v = 1/270, 
			sigma12 = 1 - 0.27, sigma13 = 1 - 0.3, sigma23 = 1 - 0.17, 
			N0 = sum(populations), tf = 365*30, 
			vaxmonth = months.vax, vaxgroup = c(as.numeric(all_allocs[x, ])), 
			n_vax = n.vax, days_imp = beta.vec[[7]])
		
		names(out) = all.cols
		mod.out = as.data.frame(out)
	
		inc = get.inc(mod.out, 365*20)
		yrs = rep(1, 365); while(length(yrs) < nrow(inc)) { yrs = c(yrs, rep(max(yrs)+1, 365)) }; yrs = yrs[1:nrow(inc)]
		ars = aggregate(inc, by = list('year' = yrs), FUN = sum)[ , -1] / 1000000
	
		qq = list('deaths' = 10*(mod.out$D[365*30] - mod.out$D[365*20]), 
			'wastes' = sum(mod.out$W[(365*20+1):(365*30)]), 
			'cases' = 10*sum(inc))
	
		return(qq)
	
	}, mc.cores = 24)


	zz = unlist(res)
	deaths = zz[seq(1, length(zz), by = 3)]
	wastes = zz[seq(2, length(zz), by = 3)]
	cases = zz[seq(3, length(zz), by = 3)]


	all.deaths = cbind(all.deaths, deaths)
	all.wastes = cbind(all.wastes, wastes)
	all.cases = cbind(all.cases, cases)


	toc()

	}

	colnames(all.deaths) = paste('deaths', c(1:ncol(all.deaths)), sep = '')
	colnames(all.cases) = paste('cases', c(1:ncol(all.cases)), sep = '')
	colnames(all.wastes) = paste('wastes', c(1:ncol(all.wastes)), sep = '')

	
	vax_res[[length(vax_res)+1]] = cbind(all_allocs, all.cases, all.deaths, all.wastes)
	names(vax_res)[length(vax_res)] = paste('Res', n.vax, sep = '_')
	

}



save(vax_res, file = 'Model_Results/Vax_Ages_Res_allscen_50ve_lhs10pct.Rdata')
save(vax_res, file = 'Model_Results/Vax_Ages_Res_allscen_270dayvd_lhs10pct.Rdata')
save(vax_res, file = 'Model_Results/Vax_Ages_Res_allscen_pop0_lhs10pct.Rdata')


