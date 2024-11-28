include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libatCLAMLink.a

c_sources = atCLAMManage.c atCLAMLinkIO.c atCLAMComm.c atCLAMLinkSetup.c atCLAMImplement.c
h_sources  = $(top_srcdir)/atCLAMLink/include/atCLAMLink.h  

libatCLAMLink_adir=$(includedir)/atlink

libatCLAMLink_a_SOURCES = $(c_sources) $(h_sources)

#atComms.h atlink.h 
    
#libatlink_a_HEADERS= 

library_includedir=$(includedir)/$(PACKAGE_NAME)-$(PACKAGE_VERSION)/$(PACKAGE_NAME)
#library_includedir=$(includedir)/
library_include_HEADERS = $(h_sources)

# overide the -O2 flag
CFLAGS=

WARN += -Wno-unused-parameter
WARN += -Werror