
FIND_PACKAGE( ROOT REQUIRED COMPONENTS MathMore )
if ( NOT ROOT_FOUND )
message(FATAL_ERROR "ROOT is required but missing !")
endif ( NOT ROOT_FOUND )
include_directories("./" ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

add_executable(testIdentity_Sign  testIdentity_Sign.C)
add_executable(testIdentityNet_particles testIdentity_NetParticles.C)
target_link_libraries(testIdentity_Sign TIdentity ${ROOT_LIBRARIES})
target_link_libraries(testIdentityNet_particles TIdentity ${ROOT_LIBRARIES})
