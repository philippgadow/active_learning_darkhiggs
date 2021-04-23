#
# Build configuration for the SimpleAnalysisCodes package
#

# Set the name of the package:
atlas_subdir( SimpleAnalysisCodes )

atlas_add_library( SimpleAnalysisCodesLib
  src/*.cxx
  NO_PUBLIC_HEADERS
  INTERFACE
  PRIVATE_INCLUDE_DIRS
  ${RESTFRAMES_INCLUDE_DIRS}
  LINK_LIBRARIES
  SimpleAnalysisFrameworkLib
)

atlas_install_data(data/*.xml)
atlas_install_data(data/*.root)
atlas_install_data(data/*.json)
