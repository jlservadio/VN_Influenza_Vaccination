

condense.contact = function(contact, groups) {
	
	new.contact = matrix(NA, nrow = length(groups), ncol = length(groups))
	
	for (i in 1:length(groups)) {
		for (j in 1:length(groups)) {
			
			temp = contact[groups[[i]], groups[[j]]]
			new.contact[i, j] = mean(as.matrix(temp))
			
		}
	}
	
	return(new.contact)
	
}

condense.deathrate = function(death.rates, pop.props, groups){
	
	new.death.rate = rep(NA, length(groups))
	
	for (i in 1:length(groups)) {
		# new.death.rate[i] = sum(death.rates[groups[[i]]] * pop.props[groups[[i]]])
		# new.death.rate[i] = mean(death.rates[1+groups[[i]]])

		new.death.rate[i] = sum(death.rates[1+groups[[i]]] * pop.props[1+groups[[i]]] / sum(pop.props[1+groups[[i]]]))

	}
	
	return(new.death.rate)
	
}


condense.fatality = function(fatality, ages, pops, groups){
	
	new.fatality = rep(NA, length(groups))
	
	for (i in 1:length(groups)) {
		new.fatality[i] = sum(fatality[which(ages %in% groups[[i]])] * 
			pops[which(ages %in% groups[[i]])] / sum(pops[which(ages %in% groups[[i]])]))
	}
	
	return(new.fatality)
	
}


set.pops = function(total.pop, data.year, data.pops, groups) {
	
	pops = NULL
	for (i in 1:length(groups)) {
		pops = c(pops, sum(data.pops[which(data.year %in% groups[[i]])] / sum(data.pops)))
	}
	pops = round(pops * total.pop)
	
	if (sum(pops) > total.pop) {
		while(sum(pops) > total.pop) { z = sample(c(1:length(pops)), 1); pops[z] = pops[z] - 1 }
	} else if (sum(pops) < total.pop) {
		while(sum(pops) < total.pop) { z = sample(c(1:length(pops)), 1); pops[z] = pops[z] + 1 }
	}
	
	return(pops)
	
}



get.inc = function(output, ts) {
	
	cols.h1 = grep('J1', names(output))
	cols.b = grep('J2', names(output))
	cols.h3 = grep('J3', names(output))
	
	inc1.0 = apply(output[which(output$t >= ts-1), cols.h1], 1, sum)
	inc2.0 = apply(output[which(output$t >= ts-1), cols.b], 1, sum)
	inc3.0 = apply(output[which(output$t >= ts-1), cols.h3], 1, sum)
	
	inc = data.frame(
		'inc1' = inc1.0[-1] - inc1.0[-length(inc1.0)], 
		'inc2' = inc2.0[-1] - inc2.0[-length(inc2.0)], 
		'inc3' = inc3.0[-1] - inc3.0[-length(inc3.0)]
		)
		
	return(inc)
	
}


long_beta = function(b0_1, b0_2, b0_3, phi_1, phi_2, phi_3, amp_1, amp_2, amp_3, del_1, del_2, del_3, lgth) {
	
	l.beta.1 = rep(b0_1, 2*lgth)
	l.beta.2 = rep(b0_2, 2*lgth)
	l.beta.3 = rep(b0_3, 2*lgth)
	
	sinu.1 = 1 + (amp_1 * cos((2*pi*((del_1/2) - c(0:(2*del_1 - 1))))/(2 * del_1))); sinu.1[which(sinu.1 < 1)] = 1
	sinu.2 = 1 + (amp_2 * cos((2*pi*((del_2/2) - c(0:(2*del_2 - 1))))/(2 * del_2))); sinu.2[which(sinu.2 < 1)] = 1
	sinu.3 = 1 + (amp_3 * cos((2*pi*((del_3/2) - c(0:(2*del_3 - 1))))/(2 * del_3))); sinu.3[which(sinu.3 < 1)] = 1
	
	for (i in 1:length(phi_1)) {
		l.beta.1[floor(phi_1[i]-(del_1/2)):(floor(phi_1[i]-(del_1/2))+length(sinu.1)-1)] = 
			l.beta.1[floor(phi_1[i]-(del_1/2)):(floor(phi_1[i]-(del_1/2))+length(sinu.1)-1)] * sinu.1
	}
	for (i in 1:length(phi_2)) {
		l.beta.2[floor(phi_2[i]-(del_2/2)):(floor(phi_2[i]-(del_2/2))+length(sinu.2)-1)] = 
			l.beta.2[floor(phi_2[i]-(del_2/2)):(floor(phi_2[i]-(del_2/2))+length(sinu.2)-1)] * sinu.2
	}
	for (i in 1:length(phi_3)) {
		l.beta.3[floor(phi_3[i]-(del_3/2)):(floor(phi_3[i]-(del_3/2))+length(sinu.3)-1)] = 
			l.beta.3[floor(phi_3[i]-(del_3/2)):(floor(phi_3[i]-(del_3/2))+length(sinu.3)-1)] * sinu.3
	}
	
	
	out = list(l.beta.1[1:lgth], l.beta.2[1:lgth], l.beta.3[1:lgth])
	return(out)
	
}


viz.alloc = function(alloc.list, labels, x1 = 0, x2 = 1, cols = rep('white', 8)) {
	
	c.alloc = c(0, cumsum(as.numeric(alloc.list)))
	
	for (i in length(c.alloc):2) {
		polygon(c(x1, x2, x2, x1), c(c.alloc[i], c.alloc[i], 0, 0), col = cols[i-1])
		if (alloc.list[i-1] > 0) {
			text(mean(c(x1, x2)), mean(c.alloc[(i-1):i]), labels[i-1])
		}
	}
	
}

