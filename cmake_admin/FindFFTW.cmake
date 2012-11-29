# copied https://code.google.com/p/gerardus/source/browse/trunk/cpp/src/third-party/cquammen-Clarity/CMake/FindFFTW.cmake?spec=svn671&r=671

# - Find FFTW
# Find the native FFTW includes and library
# This module defines
#  FFTW_INCLUDE_DIR, where to find fftw3.h, etc.
#  FFTW_LIBRARIES, the libraries needed to use FFTW.
#  FFTW_FOUND, If false, do not try to use FFTW.
# also defined, but not for general use are
#  FFTW_LIBRARY, where to find the FFTW library.

FIND_PATH(FFTW_INCLUDE_DIR fftw3.h
/usr/local/include
/usr/include
/opt/local/lib
)

SET(FFTW_NAMES ${FFTW_NAMES} fftw3 fftw3f fftw3-3)
FIND_LIBRARY(FFTW_LIBRARY
  NAMES ${FFTW_NAMES}
  PATHS /usr/lib /usr/local/lib /opt/locala/lib
  )
separate_arguments( FFTW_LIBRARY )
# Find threads part of FFTW

SET(FFTW_THREADS_NAMES ${FFTW_THREADS_NAMES} fftw3f_threads fftw3_threads fftw3-3_threads)
FIND_LIBRARY(FFTW_THREADS_LIBRARY
  NAMES ${FFTW_THREADS_NAMES}
  PATHS /usr/lib /usr/local/lib /opt/local/lib
  )

IF (FFTW_THREADS_LIBRARY AND FFTW_INCLUDE_DIR)
    SET(FFTW_THREADS_LIBRARIES ${FFTW_THREADS_LIBRARY})
    SET(FFTW_THREADS_FOUND "YES")
ELSE (FFTW_THREADS_LIBRARY AND FFTW_INCLUDE_DIR)
  SET(FFTW_THREADS_FOUND "NO")
ENDIF (FFTW_THREADS_LIBRARY AND FFTW_INCLUDE_DIR)


IF (FFTW_THREADS_FOUND)
   IF (NOT FFTW_THREADS_FIND_QUIETLY)
      MESSAGE(STATUS "Found FFTW threads: ${FFTW_THREADS_LIBRARIES}")
   ENDIF (NOT FFTW_THREADS_FIND_QUIETLY)
ELSE (FFTW_THREADS_FOUND)
   IF (FFTW_THREADS_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find FFTW threads library")
   ENDIF (FFTW_THREADS_FIND_REQUIRED)
ENDIF (FFTW_THREADS_FOUND)


IF (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
    SET(FFTW_LIBRARIES ${FFTW_LIBRARY})
    SET(FFTW_FOUND "YES")
ELSE (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
  SET(FFTW_FOUND "NO")
ENDIF (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)


IF (FFTW_FOUND)
   IF (NOT FFTW_FIND_QUIETLY)
      MESSAGE(STATUS "Found FFTW: ${FFTW_LIBRARIES}")
   ENDIF (NOT FFTW_FIND_QUIETLY)
ELSE (FFTW_FOUND)
   IF (FFTW_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find FFTW library")
   ENDIF (FFTW_FIND_REQUIRED)
ENDIF (FFTW_FOUND)

SET (ON_UNIX ${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR
             ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
IF (${ON_UNIX})
   SET (FFTW_EXECUTABLE_LIBRARIES fftw3f fftw3f_threads)
ENDIF (${ON_UNIX})

