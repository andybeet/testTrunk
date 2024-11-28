
include $(top_srcdir)/PreRules.am
lib_LIBRARIES=libatharvest.a

libatharvest_adir=$(includedir)/atharvest

libatharvest_a_SOURCES = atHarvestAnnual.c atHarvest.c atHarvestIndex.c atHarvestIO.c atHarvestTS.c \
	atHarvest.h atHarvestLib.h atHarvestSetup.c atHarvestParamIO.c atHarvestDiscards.c atHarvestDiscards.h \
	atHarvestIO.c atHarvestIO.h atHarvestCatch.c atHarvestImposedCatch.c atImposedCatch.h
    
#libatharvest_a_HEADERS= atHarvest.h atHarvestLib.h 

# overide the -O2 flag
CFLAGS=
WARN += -Werror
