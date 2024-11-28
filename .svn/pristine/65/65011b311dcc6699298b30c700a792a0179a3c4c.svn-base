
include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libatassess.a

libatassess_adir=$(includedir)/atassess

libatassess_a_SOURCES = atclassical.c atindices.c atsample.c atsurvey.c \
    atdatareading.c atsamplesetup.c \
    atwriting.c atattribute.c atdiet.c atnetwork.c atsolve.c  \
    atassessParamIO.c atAssessTierSetup.c \
    atsample.h atsamplesetup.h doolittle.c unit_lower_triangular.c upper_triangular.c 
 
#atmemory.c   

# overide the -O2 flag
#CFLAGS=

WARN += -Wno-unused-parameter
WARN += -Werror



