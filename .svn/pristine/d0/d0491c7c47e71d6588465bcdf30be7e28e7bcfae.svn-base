include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libbrokerlink.a

libbrokerlink_adir=$(includedir)/atbrokerlink

libbrokerlink_a_SOURCES = requests.pb-c.h requests.pb-c.c responses.pb-c.h responses.pb-c.c \
network.c network.h NetworkError.c NetworkError.h \
LinkageInterface.c LinkageInterface.h Deserialiser.h Deserialiser.c \
atBrokerLinkInit.c ImportExportData.c ImportExportData.h
   
#atmemory.c   

# overide the -O2 flag
#CFLAGS=

WARN += -Wno-unused-parameter
WARN += -Werror

