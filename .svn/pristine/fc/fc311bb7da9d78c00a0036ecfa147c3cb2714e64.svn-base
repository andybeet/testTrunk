include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libatlink.a

libatlink_adir=$(includedir)/atlink

libatlink_a_SOURCES = atComms.c atlinkconversion.o atlinkexport.c atlinkimport.c atComms.h atlink.h atLinkModule.h
    
#libatlink_a_HEADERS= atComms.h atlink.h atLinkModule.h

# overide the -O2 flag
CFLAGS=
