FIND_LIBRARY(IntelFCE_Lib NAMES ifcore PATHS /opt/intel/fce/10.1.015/lib)
SET(EXTRALIB ${EXTRALIB} ${IntelFCE_Lib})
