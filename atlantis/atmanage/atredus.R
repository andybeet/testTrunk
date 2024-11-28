# Load configuration
source("redus_conf.R")

library(FLCore)
library(FLa4a)

scaAssess <- function(minYear, maxYear, biomass, params) {

	print("REDUS: SCA assessment")
	print(minYear)	
	print(maxYear)
	print(biomass)
	print(params)

	q <- as.numeric(params[["q"]])

	biomassArr <- biomass[, as.character(minYear:maxYear)]

	dnms <- list(age="all", year=range(biomassArr)["minyear"]:range(biomassArr)["maxyear"])

	bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=dnms))

	index(bioidx) <- stock(biomassArr) * q

	range(bioidx)[c("startf","endf")] <- c(0,0)

	fit <- sca(ple4, FLIndices(bioidx), qmodel=list(~1))

	# TODO: Need to change the return
	return(fit)
}


tpnAssess <- function(minYear, maxYear, biomass, params) {

	print("REDUS: Truth plus noise assessment")

	print(minYear)	
	print(maxYear)
	print(biomass)
	print(params)

	q <- as.numeric(params[["q"]])

	biomassArr <- biomass[, as.character(minYear:maxYear)]
	idx <- stock(biomassArr)
	idx <- idx * exp(rnorm(idx * q, sd = 0.1))[1]
	stock(biomass)[, as.character(minYear:maxYear)] <- idx

	return(biomass)
}


doAssess <- function(startyear, curryear, groupCode) {
	
	# Get variables from global
	assessParameters <- get("assessParameters", envir = .GlobalEnv)
	biomassIndex <- get("biomassIndex", envir = .GlobalEnv)
	assessmentOutput <- get("assessmentOutput", envir = .GlobalEnv)

	assessType <- assessParameters[[groupCode]][["type"]]
	assessParams <- assessParameters[[groupCode]][["params"]]


	if( assessType == "SCA" ) {
		aResult <- scaAssess(startyear, curryear, biomassIndex[[groupCode]], assessParams)
	} else if( assessType == "NON" ) {
		aResult <- tpnAssess(startyear, curryear, biomassIndex[[groupCode]], assessParams)
	}

	print("REDUS: Put back assessment result in global")
	assessmentOutput[[groupCode]] <<- aResult
}

icesHCR <- function(stkSSB, curr_yr, params) {
	print("REDUS: ICES HCR")
	print(stkSSB)	
	print(params)
		
	fmin <- as.numeric(params[["fmin"]])
	ftrg <- as.numeric(params[["ftrg"]])
	blim <- as.numeric(params[["blim"]])
	bsafe <- as.numeric(params[["bsafe"]])
	ssb_lag <- as.numeric(params[["ssb_lag"]])

	ssb <- stkSSB[, ac(curr_yr-ssb_lag)]

	fout <- FLQuant(fmin, dimnames=list(iter=dimnames(stkSSB)$iter))
	
	fout[ssb >= bsafe] <- ftrg

	inbetween <- (ssb < bsafe) & (ssb > blim)
	gradient <- (ftrg - fmin) / (bsafe - blim)
	fout[inbetween] <- (ssb[inbetween] - blim) * gradient + fmin

	return(fout)
}

doHCR <- function(curryear, groupCode) {
	
	# Get variables from global
	hcrParameters <- get("hcrParameters", envir = .GlobalEnv)
        hcrOutput <- get("hcrOutput", envir = .GlobalEnv)
        assessmentOutput <- get("assessmentOutput", envir = .GlobalEnv)

	#print(stock(assessmentOutput[[groupCode]]))

	print(groupCode)
	print(hcrParameters[[groupCode]])

	hcrType <- hcrParameters[[groupCode]][["type"]]
	hcrParams <- hcrParameters[[groupCode]][["params"]]

	# rule
	if( hcrType == "ICS" ) {
		fout <- icesHCR(stock(assessmentOutput[[groupCode]]), curryear, hcrParams)
	}

	# Put it inside the HCR output
	print("REDUS: Put back HCR result in global")
	hcrOutput[[groupCode]][, ac(curryear)] <<- fout[,1]

	print(hcrOutput[[groupCode]])

	return(fout)
}

