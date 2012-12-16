project(FFTWExamples)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Find ITK.
MESSAGE( STATUS "hola"  )
FIND_PACKAGE(FFTW3 REQUIRED)
MESSAGE( STATUS "found: "  ${FFTW3_FOUND})
IF(FFTW3_FOUND)
MESSAGE( STATUS "found it"  )
MESSAGE( STATUS "PATH: "  ${FFTW3_INCLUDE_DIRS})
include_directories(${FFTW3_INCLUDE_DIRS})
set(LIBS ${FFTW3_LIBRARIES})
ENDIF(FFTW3_FOUND)

SET( EXAMPLE_LIST
test_kde test_kde_dct test_extrema
 )


SET( OBJ_LIST_TMP
kde_util
kde
roots
 )
 
SET(OBJ_LIST )

foreach (o ${OBJ_LIST_TMP})
    MESSAGE( STATUS "o: " ${o})
    SET(OBJ_LIST ${OBJ_LIST} ${o}.c )
endforeach()

MESSAGE( STATUS "obj_list: " ${OBJ_LIST})


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99")

FOREACH( var ${EXAMPLE_LIST} )
  ADD_EXECUTABLE( ${var} ${var}.c ${OBJ_LIST})
  TARGET_LINK_LIBRARIES( ${var}  ${LIBS} -lgsl -lgslcblas)
ENDFOREACH( var )



