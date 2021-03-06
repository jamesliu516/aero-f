SET(CMAKE_LIBRARY_PATH /home/pavery/Intel/ARPACK)
## blas and lapack with scalapack and blacs (using intel math kernel library)
SET(LAPACK_LIBRARIES -Wl,--start-group /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_scalapack_lp64.a /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_blacs_openmpi_lp64.a  /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_sequential.a /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_core.a /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_rt.so -Wl,--end-group -lpthread CACHE STRING "Path to a library.")
SET(LAPACK_FOUND true)
SET(SCALAPACK_FOUND true)
SET(BLACS_FOUND true)
SET(BLACSLIB "")
SET(BLACSCLIB "")
SET(SCALAPACKLIB "")
