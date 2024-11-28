include $(top_srcdir)/PreRules.am

lib_LIBRARIES=libatmanage.a

libatmanage_adir=$(includedir)/atmanage

libatmanage_a_SOURCES = atManageAnnual.c atManage.c \
	atManageSetup.c atManageIndex.c atManageMPATS.c atManageIO.c \
    atManage.h atManageParamIO.h atManagePrivaye.h \
    atManageParamIO.c atManageTier.c atSS3assess.c atCloseKin.c \
    atRlink.h atRlinkRedus.c atRlinkRBC.c atRlinkRAssess.c

# overide the -O2 flag
#WARN += -Wno-unused-parameter
#WARN += -Werror
