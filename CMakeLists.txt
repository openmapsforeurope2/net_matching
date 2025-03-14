cmake_minimum_required( VERSION 2.8.3 )
cmake_policy(SET CMP0020 NEW)

# Definition du projet
set(PROJECT_NAME net_matching)
project(${PROJECT_NAME})

FIND_PATH(SOCLE_EXTERN_ROOT NAMES SocleExterns.cmake
	PATHS $ENV{SOCLE_EXTERN_ROOT}
	PATHS ../Lib/SOCLE_WORKSPACE/extern-1.0
	PATHS ../Lib/extern-1.0
	PATHS ../Lib/extern
)

IF (NOT SOCLE_EXTERN_ROOT)
	MESSAGE(FATAL_ERROR "Could not locate IGN extern root, define env SOCLE_EXTERN_ROOT")
ENDIF()

include( ${SOCLE_EXTERN_ROOT}/SocleExterns.cmake )


# disable warning substition /W3 par /W4
IF (MSVC)
	SET(CMAKE_CXX_FLAGS_TEMP "")
	IF(${CMAKE_CXX_FLAGS} MATCHES "/W4")
		STRING(REGEX REPLACE "/W4 " "" CMAKE_CXX_FLAGS_TEMP "${CMAKE_CXX_FLAGS}")
		SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_TEMP})
	ENDIF()
	#pour avoir 4G de memoire en 64 bit
	SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE" )
ENDIF (MSVC)

# Maintenant on peut ajouter des lib avec une seule ligne de commande
SET(Boost_USE_STATIC_LIBS ON)
SET(Boost_USE_MULTITHREADED ON)
add_definitions( "-DBOOST_THREAD_USE_LIB" )
add_definitions( "-DBOOST_ALL_NO_LIB" )
find_package(Boost COMPONENTS locale date_time filesystem system serialization thread program_options regex REQUIRED)
find_package(Shapelib REQUIRED)
find_package(Geos REQUIRED)
find_package(CGAL COMPONENTS Core REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package( Pq REQUIRED )
find_package( Qt5 COMPONENTS Widgets Gui Core Xml REQUIRED )
add_definitions(${QT_DEFINITIONS})
find_package( Sqlite REQUIRED )


find_package( IgnSocle COMPONENTS core numeric math tools data filesystem geometry sql feature graph transform shapefile cgal REQUIRED )


# Ajout de LibEPG
FIND_PATH(
    LibEPG
    NAMES LibEPG.cmake
    PATHS
        $ENV{LibEPG_ROOT}
		$ENV{LIB_EPG}
)


IF (NOT LibEPG)
        MESSAGE(FATAL_ERROR "Could not locate LibEPG")
ENDIF()

include (${LibEPG}/LibEPG.cmake)


include_directories(
	${Boost_INCLUDE_DIRS}
	${Qt5Core_INCLUDE_DIRS}
	${Qt5Widgets_INCLUDE_DIRS}
	${Qt5Gui_INCLUDE_DIRS}
	${Qt5Xml_INCLUDE_DIRS}
	${Pq_INCLUDE_DIRS}
	${Geos_INCLUDE_DIRS}
    ${Shapelib_INCLUDE_DIRS}
	${CGAL_INCLUDE_DIRS}
	${GMP_INCLUDE_DIRS}
	${MPFR_INCLUDE_DIRS}
	${Sqlite_INCLUDE_DIRS}
    ${IgnSocle_INCLUDE_DIRS}
    ${LibEPG_INCLUDE_DIRS}
)

link_directories(               
	${Boost_LIBRARY_DIRS}
	${QT_LiBRARY_DIRS}
	${Pq_LIBRARY_DIRS}
	${Geos_LIBRARY_DIRS}
    ${Shapelib_LIBRARY_DIRS}
	${CGAL_LIBRARY_DIRS}
	${GMP_LIBRARY_DIRS}
	${MPFR_LIBRARY_DIRS}
	${Sqlite_LIBRARY_DIRS}
    ${IgnSocle_LIBRARY_DIRS}
    ${LibEPG_LIBRARY_DIRS}
)

#definition de l'output de l'executable	
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

#configuration de l'executable
file(GLOB_RECURSE APP_HEADERS include/*.h)
file(GLOB_RECURSE APP_SOURCES src/*.cpp)

IF(MSVC_IDE)
	SOURCE_GROUP_RECURSE( "${APP_HEADERS}" "${CMAKE_SOURCE_DIR}/include/" "Header Files" )
	SOURCE_GROUP_RECURSE( "${APP_SOURCES}" "${CMAKE_SOURCE_DIR}/src/" "Source Files" )
ENDIF(MSVC_IDE)

#include de mes headers
include_directories(include)

add_executable(
	${PROJECT_NAME} 
	${APP_HEADERS}
	${APP_SOURCES}
)
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

#Configuration de l'edition de liens
target_link_libraries(
	${PROJECT_NAME} 
    ${LibEPG_LIBRARIES}
    ${IgnSocle_LIBRARIES}
	${Pq_LIBRARIES}
	${Sqlite_LIBRARIES}
    ${Shapelib_LIBRARIES}
	${CGAL_LIBRARIES}
	#${CGAL_CORE_LIBRARIES}
	#${CGAL_LIBRARY}
	#${CGAL_CORE_LIBRARY}
	${GMP_LIBRARIES}
	${MPFR_LIBRARIES}
	${Geos_LIBRARIES}
	${Qt5Xml_LIBRARIES}
	${Qt5Widgets_LIBRARIES}
	${Qt5Gui_LIBRARIES}
	${Qt5Core_LIBRARIES}
	${Qt5Core_QTMAIN_LIBRARIES}
	${Boost_LIBRARIES}
)

#-- installation
install(
	TARGETS
	${PROJECT_NAME}
	DESTINATION
	bin
)

if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/version)
    file(READ ${CMAKE_CURRENT_LIST_DIR}/version VERSION)
else()
    message("version file not found")
	set(VERSION 0.0.0)
endif ()

set(CPACK_DEBIAN_PACKAGE_NAME "${PROJECT_NAME}")
set(CPACK_DEBIAN_FILE_NAME "${PROJECT_NAME}-${VERSION}.deb")
set(CPACK_DEBIAN_PACKAGE_VERSION "${VERSION}")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "IGN")

SET(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_INSTALL_PREFIX}/package")

set(CPACK_GENERATOR "DEB")


include(CPack)