
include $(top_srcdir)/PreRules.am
lib_LIBRARIES=libateconomic.a

libateconomic_adir=$(includedir)/ateconomic

libateconomic_a_SOURCES = ateconhelp.c  ateconindicator.c ateconomic.c ateconomicsetup.c \
    ateconresponse.c ateconts.c atquota.c ateconeffort.c \
    ateconio.c ateconParamIO.c  \
    ateconomic.h ateconomicsetup.h atEconomicLib.h\
    ateffortquota.c atEconomicAnnual.c atEconomicUtil.c
    
#libeconomic_a_HEADERS= ateconomic.h atEconomicModule.h ateconomicsetup.h ateconParamIO.h

# overide the -O2 flag
CFLAGS=

WARN += -Wno-unused-parameter
WARN += -Werror
