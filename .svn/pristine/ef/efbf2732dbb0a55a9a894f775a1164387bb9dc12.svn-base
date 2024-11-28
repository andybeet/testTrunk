# --- DO NOT EDIT THE LINES BELOW
assessParameters = list()
hcrParameters = list()
# --- END

# Global configuration #

# When to start collecting biomass data
collectStart = 0

# When to start the assessments
assessStart = 0

# Species shortnames that are going to be assessed
assessSpecies = c("NCO")

# Assessment parameters 
## SAITHE
assessParameters[["SAI"]] = list(
	type   = "NON",
	params = list(
		q = 0
	)
)

## COD
assessParameters[["NCO"]] = list(
        type   = "NON",
        params = list(
                q = 0
        )
)

# HCR parameters
## SAITHE
hcrParameters[["SAI"]] = list(
	type   = "ICS",
	params = list(
        	ssb_lag = 0,
        	fmin = 0.210,
        	ftrg = 0.536,
        	blim = 107297,
        	bsafe = 149098
	)
)

## COD
hcrParameters[["NCO"]] = list(
        type   = "ICS",
        params = list(
                ssb_lag = 0,
                fmin = 0.4,
                ftrg = 0.6,
                blim = 2500000,
                bsafe = 5000000
        )
)

