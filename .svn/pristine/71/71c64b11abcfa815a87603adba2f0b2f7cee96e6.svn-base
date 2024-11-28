
include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libatecology.a

libatecology_adir=$(includedir)/atecology

libatecology_a_SOURCES = atbiolUtil.c atannualbiology.c atbiolmessage.c atbiology.c atbiolsetup.c \
	atbiophysics.c atcalibtools.c atecology.c atecologyts.c \
	atfluxbreakdown.c atmovement.c atprocess.c atq10.c attime.c \
	atvertprocesses.c atmacrophytes.c \
	atBuildTracer.c atGroupProcesses.c \
	atBiologyXMLParamIO.c atNutrient.c atPhysChemIO.c \
	atdemography.c atdiversity.c atcoral.c atImposeRecruit.c\
    atbiology.h atbiolParamIO.h atbiolsetup.h atbiolUtil.h atecology.h atecologylib.h \
    atecologyModule.h atlog.h atnamelist.h atNutrient.h attime.h additionalTracer.c \
    atExternalScalar.c atExternalScalar.h atdemography.c atIceProcesses.c \
    atLandProcess.c atContaminants.c atmigration.c atexternalpop.c
    
#libatecology_a_HEADERS= atbiology.h atbiolParamIO.h atbiolsetup.h atbiolUtil.h atecology.h atecologylib.h \
#atecologyModule.h atlog.h atnamelist.h atNutrient.h attime.h

# overide the -O2 flag
CFLAGS=

WARN += -Wno-unused-parameter
WARN += -Werror
#WARN += -Wno-unused-but-set-variable
