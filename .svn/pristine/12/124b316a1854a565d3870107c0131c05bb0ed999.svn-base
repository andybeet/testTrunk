include $(top_srcdir)/PreRules.am

h_sources = $(top_srcdir)/sjwlib/include/datafile.h \
$(top_srcdir)/sjwlib/include/gridlib.h \
$(top_srcdir)/sjwlib/include/mapproj.h \
$(top_srcdir)/sjwlib/include/pointsourcesink.h \
$(top_srcdir)/sjwlib/include/polyline.h \
$(top_srcdir)/sjwlib/include/ptrack.h \
$(top_srcdir)/sjwlib/include/sjwalloc.h \
$(top_srcdir)/sjwlib/include/sjwconfig.h \
$(top_srcdir)/sjwlib/include/sjwlib.h \
$(top_srcdir)/sjwlib/include/timeseries.h \
$(top_srcdir)/sjwlib/include/txt_param.h \
$(top_srcdir)/sjwlib/include/xytoij.h \
$(top_srcdir)/sjwlib/include/colourtable.h \
$(top_srcdir)/sjwlib/include/nrcdflib.h

c_sources = cfft.c colourtable.c contour.c datafile.c decay.c dfcoords.c \
	dfeval.c diffusion.c drandom.c ellipt_coord.c geodetic.c \
	gridangle.c gridmetric.c gridmisc.c keyfile.c mapproj.c \
	memory.c netcdf.c pointsourcesink.c polar_coord.c polyline.c \
	ppbfetch.c ptrack.c quit.c rect_coord.c spline.c text_input.c \
	time.c timeseries.c txt_param.c warn.c waterprops.c \
	weight_fn.c xytoij.c convertFileFormats.c nrcdflib.c erfc.c

lib_LIBRARIES=libsjwlib.a

libsjwlib_adir=$(includedir)/sjwlib

libsjwlib_a_SOURCES = $(h_sources) $(c_sources)
    
library_includedir=$(includedir)/$(PACKAGE_NAME)-$(PACKAGE_VERSION)/$(PACKAGE_NAME)
#library_includedir=$(includedir)/
library_include_HEADERS = $(h_sources)

# overide the -O2 flag
CFLAGS=

#WARN += -Wno-error
#WARN += -Werror


pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libsjwlib.pc
