#
# Project build configuration.
#



# Set the minimum required CMake version.
cmake_minimum_required( VERSION 3.4 FATAL_ERROR )

# Look for any AnalysisBase-21.2.X release:
find_package( AnalysisBase 21.2 REQUIRED )

option(DO_TRUTHTAGGING "Build with Truth Tagging" OFF)

include( Ext_RestFrames/externals.cmake )

# Set up the usage of CTest:
atlas_ctest_setup()




# Set up the project:
atlas_project( SimpleAnalysisProject 1.0.0
   USE AnalysisBase ${AnalysisBase_VERSION} )

# Generate an environment setup script:
lcg_generate_env( SH_FILE ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh )
install( FILES ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh
   DESTINATION . )

atlas_subdir( SimpleAnalysis )
atlas_add_library( SimpleAnalysisLib
  src/*.cxx
  NO_PUBLIC_HEADERS
  INTERFACE
  LINK_LIBRARIES
  SimpleAnalysisFrameworkLib
)


# Set up the usage of CPack:
atlas_cpack_setup()

