cmake_minimum_required(VERSION 2.6)
SET(INC_DIR .. . /home/michel/Research/Academic/PhysBam /home/michel/Research/Academic/PhysBam/Public_Library /home/icmewang/PhysBam /home/icmewang/PhysBam/Public_Library /home/icmewang/PhysBam /home/icmewang/PhysBam/Public_Library /home/pavery/boost/boost_install/include/boost-1_34_1)
add_definitions(-DTEMPLATE_FIX -D__LINUX -DF_NEEDS_UNDSC -DTYPE_PREC=float)
add_library(PhysBAMIntersect SHARED PhysBAMIntersect.C)
