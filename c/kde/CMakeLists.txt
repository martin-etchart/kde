project(KDE)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

# Find ITK.
MESSAGE( STATUS "hola"  )
FIND_PACKAGE(FFTW3 REQUIRED)
MESSAGE( STATUS "found: "  ${FFTW3_FOUND})
IF(FFTW3_FOUND)
MESSAGE( STATUS "found it"  )
MESSAGE( STATUS "PATH: "  ${FFTW3_INCLUDE_DIRS})
include_directories(${FFTW3_INCLUDE_DIRS})
set(KDE_LIBS ${FFTW3_LIBRARIES})
ENDIF(FFTW3_FOUND)

include_directories("/usr/local/include")
link_directories("/usr/local/lib")

SET( EXAMPLE_LIST test_kde test_kde_dct test_extrema )

SET( OBJ_LIST_TMP roots kde kde_util   )
 
SET(OBJ_LIST )

MESSAGE( STATUS "OBJ LIST")
foreach (var ${OBJ_LIST_TMP})
	MESSAGE( STATUS "o: " ${var})
	set(KDE_HEADERS	${KDE_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/${var}.h PARENT_SCOPE)
	set(KDE_SOURCES	${KDE_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/${var}.c PARENT_SCOPE)

	SET(OBJ_LIST ${OBJ_LIST} ${var}.c )
	if(WIN32)
		SET_SOURCE_FILES_PROPERTIES( ${var}.c PROPERTIES LANGUAGE CXX )
	endif(WIN32)
endforeach()

MESSAGE( STATUS "obj_list: " ${OBJ_LIST})


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99")

set(KDE_LIBS ${KDE_LIBS}  gsl gslcblas)
set(KDE_LIBS ${KDE_LIBS} PARENT_SCOPE)

MESSAGE( STATUS "LIBS: " ${KDE_LIBS})

FOREACH( var ${EXAMPLE_LIST} )
  ADD_EXECUTABLE( ${var} ${var}.c ${OBJ_LIST})
  TARGET_LINK_LIBRARIES( ${var}  ${KDE_LIBS})
ENDFOREACH( var )



