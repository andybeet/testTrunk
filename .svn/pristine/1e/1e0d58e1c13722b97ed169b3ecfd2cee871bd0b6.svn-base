
include $(top_srcdir)/PreRules.am
lib_LIBRARIES=libatimplementation.a

libatimplementation_adir=$(includedir)/atimplementation

c_sources = atImplementationAnnual.c  atImplementationSetup.c atImplementationParamIO.c 
	
h_sources = atImplementationLib.h atImplementation.h 

libatimplementation_a_SOURCES = $(c_sources) $(h_sources)
    
#libimplementation_a_HEADERS= atImplementationLib.h atImplementation.h

# overide the -O2 flag
CFLAGS=
#WARN += -Werror
