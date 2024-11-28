include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libatSS3Link.a

c_sources = atSS3Link.c atSS3LinkIO.c atSS3DataGen.c atSS3Util.c atCloseKin.c 

# atSS3Test.c 

h_sources  = $(top_srcdir)/atSS3Link/include/atSS3LinkLib.h \
$(top_srcdir)/atSS3Link/include/atSS3Util.h

libatSS3Link_adir=$(includedir)/atSS3Link

libatSS3Link_a_SOURCES = $(c_sources) $(h_sources)


library_includedir=$(includedir)/$(PACKAGE_NAME)-$(PACKAGE_VERSION)/$(PACKAGE_NAME)
library_include_HEADERS = $(h_sources)

# overide the -O2 flag
CFLAGS=

#WARN += -Wno-unused-parameter
#WARN += -Werror
