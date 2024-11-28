include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libatphysics.a

libatphysics_adir=$(includedir)/atphysics

libatphysics_a_SOURCES = atbioirrig.c atbioturb.c atboundary.c atdecay.c atdeposition.c \
	atdiagIO.c atepiIO.c atfishstatIO.c atgas.c atgeomIO.c \
	athydromod.c atindex.c atparamIO.c atphysIO.c atphysics.c \
	atprofile.c  atsaturation.c atsedprops.c atsettle.c \
	atsourcesink.c atsummaryIO.c atsuspension.c attempsalt.c \
	attracerIO.c attransport.c atvdiff.c atvertgeom.c atvmix.c \
	athdiff.c atPhysicsModule.c atPhysics.h atPhysicsModule.h \
    aticeIO.c atLandIO.c atswr.c atagetracerIO.c atdietIO.c
    
#libatphysics_a_HEADERS = atPhysics.h atPhysicsModule.h


# overide the -O2 flag
CFLAGS=
#WARN += -Werror


