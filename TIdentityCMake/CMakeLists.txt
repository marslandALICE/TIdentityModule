
FIND_PACKAGE(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net Matrix)
if ( NOT ROOT_FOUND )
message(FATAL_ERROR "ROOT is required but missing !")
endif ( NOT ROOT_FOUND )

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake/modules)
include(RootInstallDirs)

SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories("inc" ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

FILE(GLOB SRCFILES src/*.cc)
FILE(GLOB HFILES inc/TIdentityFunctions.h inc/TIdentityBase.h inc/TIdentity2D.h)

ROOT_GENERATE_DICTIONARY(TIdentityDict ${HFILES} LINKDEF inc/TIdentityLinkDef.h)

add_library(TIdentity SHARED ${SRCFILES} TIdentityDict.cxx)
target_link_libraries(TIdentity ${ROOT_LIBRARIES})
message(AUTHOR_WARNING "test libraries ${ROOT_LIBRARIES}")

add_subdirectory(test)