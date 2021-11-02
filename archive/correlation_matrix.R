library(stats)
library(Hmisc)

# create objects containing a list of continuous, categorical and all variable names
cont.vars <- c("titerincrease","season","prevactiter","age","days_before_vac","bmi")
cat.vars  <- c("subtype","sex","prior_year_vac")
names.all <- c(cont.vars, cat.vars)

# diagonal matrix for bivariate analysis results	
M = diag(nrow=length(names.all))
colnames(M) <- row.names(M) <- names.all

# iterate through all predictors to test the association between the two
for (i in 1:(nrow(M)-1)){
	for ( j in (i+1) : nrow(M)){
		bivars = c(rownames(M)[i],colnames(M)[j])
		ncat = sum(bivars %in% cat.vars)
		
		# both variables are categorical
		if(ncat ==2){ 
			
			# use Chi-square to assess the association between two variables 
			asscn=chisq.test(table(dat_elderly_long[,bivars]))
			#print(asscn)
			# extract p-value 
			M[i,j] <- asscn$p.value
			
			# 1 variable is continuous while another one is categorical (point-biserial)
		} else if(ncat ==1){ 
			
			y = which(bivars %in% cont.vars)
			asscn <- cor.test(dat_elderly_long[,bivars[y]],as.numeric(unlist(dat_elderly_long[,bivars[-y]])))
			#print(summary(asscn))
			# extract p-value 
			M[i,j]  <-asscn$p.value
			
			# both variables are continuous	
		}else{ 
			
			asscn <- cor.test(as.numeric(unlist(dat_elderly_long[,bivars[1]])), as.numeric(unlist(dat_elderly_long[,bivars[2]])),
											 method = "spearman",
											 exact = F)
			#print(summary(asscn))
			# extract p-value 
			M[i,j]<- asscn$p.value
		}
	}
}


M[lower.tri(M, diag = FALSE)] <- t(M)[lower.tri(M, diag = FALSE)]
diag(M) <-0
M.long <- c()
M.long <- data.frame(col = rep(colnames(M),each=nrow(M)),
										 row =  rep(colnames(M),nrow(M)),
										 assocn = c(M))
M.long %>%
	filter(assocn<0.1) %>%
	ggplot(aes(x= factor(col, level = names.bivariate),y=factor(row, level = new.names.all2), fill = assocn)) +
	geom_tile() +
	scale_fill_gradient(low = "steelblue",
											high = "white", space = "Lab")+
	labs(x = NULL, y = NULL) +
	theme(axis.text=element_text(size=15))+   theme(axis.text.x = element_text(angle = 60, hjust = 1))

