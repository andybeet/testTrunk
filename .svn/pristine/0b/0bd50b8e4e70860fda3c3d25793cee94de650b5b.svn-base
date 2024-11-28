include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libatlantisutil.a

libatlantisutil_adir=$(includedir)/atlantisUtil

libatlantisutil_a_SOURCES = atUtilhelp.c atUtil.c atUtilArray.c atUtilUnix.c atUtilIO.c atUtilGroupIO.c atUtilXML.c atUtilFisheryIO.c \
atUtilFisheryXML.c 

h_sources = $(top_srcdir)/atlantisUtil/include/atUtilLib.h $(top_srcdir)/atlantisUtil/include/atTracer.h \
$(top_srcdir)/atlantisUtil/include/atXMLUtil.h $(top_srcdir)/atlantisUtil/include/atFunctGroup.h \
$(top_srcdir)/atlantisUtil/include/atGroupTypes.h

library_includedir=$(includedir)/$(PACKAGE_NAME)-$(PACKAGE_VERSION)/$(PACKAGE_NAME)
library_include_HEADERS = $(h_sources)

# overide the -O2 flag
CFLAGS=

WARN += -Wno-unused-parameter
WARN += -Werror
