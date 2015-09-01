load('RLS_fish.Rdata')

# expanding into long format (row for each individual fish)
.counts.to.rows<-function(fishdat){
	rows.expected<- sum(fishdat$Num)
	fishdat.expanded <- as.data.frame(lapply(fishdat, function(x)rep(x,fishdat$Num)))
	rows.returned<- nrow(fishdat.expanded)
	cat("Expected number of rows in expanded DF: ", rows.expected)
	cat("\nActual number of rows in expanded DF: ", rows.returned)
	if(rows.expected!=rows.returned) stop("error: rows expected not equal to rows returned. check function.\n")
	fishdat.expanded
}

l2lims<-c(5,16)
dat<- fish[fish$l2bm>=l2lims[1]&fish$l2bm<=l2lims[2],]
dat$l2.mid.C<-dat$bin.mid.l2-10.5

# trim un-necessary columns first as it gets very big
	fish.l<- dat[c("SiteCode","Lat_Zone","Sizeclass","Num", "pcbm", "bm.bin","bin.mid.l2","Ecoregion")]
	fish.l<-.counts.to.rows(fish.l)

n.survs <- plyr::ddply(fish, plyr::.(SiteCode), function(x) data.frame(n.surv=length(unique(x$SurveyID))))

