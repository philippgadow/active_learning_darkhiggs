#
# Package building RestFrames
#

# Externals needed by the build of our externals:
find_package( ROOT REQUIRED )

# Temporary directory for the build results:
set( _restFramesBuildDir
   ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/RestFramesBuild )

# Set up the variables that the users can pick up RestFrame with:
set( RESTFRAMES_INCLUDE_DIRS
   $<BUILD_INTERFACE:${_restFramesBuildDir}/include>
   $<INSTALL_INTERFACE:include>
   ${ROOT_INCLUDE_DIRS} )
set( RESTFRAMES_LIBRARIES
   ${_restFramesBuildDir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}RestFrames${CMAKE_SHARED_LIBRARY_SUFFIX} )

# Set up the build of RestFrames for the build area:
ExternalProject_Add( RestFrames
   PREFIX ${CMAKE_BINARY_DIR}
   URL https://github.com/crogan/RestFrames/archive/v1.0.1.tar.gz
   URL_MD5 668e3ca6f301172d7e67b5e85b1ab6d2
   INSTALL_DIR ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}
   CONFIGURE_COMMAND <SOURCE_DIR>/configure
   --prefix=${_restFramesBuildDir} --enable-shared --disable-static
   --with-rootsys=${ROOTSYS}
   BUILD_IN_SOURCE 1
   INSTALL_COMMAND COMMAND make install
   COMMAND ${CMAKE_COMMAND} -E copy_directory
   ${_restFramesBuildDir}/ <INSTALL_DIR>
   BUILD_BYPRODUCTS ${RESTFRAMES_LIBRARIES} )

# Make sure that whoever uses RestFrames, uses ROOT as well:
list( APPEND RESTFRAMES_LIBRARIES ${ROOT_LIBRARIES} )

# Install RestFrames:
install( DIRECTORY ${_restFramesBuildDir}/
   DESTINATION . USE_SOURCE_PERMISSIONS )

# Clean up:
unset( _restFramesBuildDir )
