# Ext_RestFrames

Package is a wrapper class to get the RestFrames package compiled in a RootCore or CMake environment.

## RootCore

### Issues

If you see an issue like

```
compiling Ext_RestFrames
WARNING: package Ext_RestFrames doesn't request pedantic compilation
WARNING: fix this by setting 'PACKAGE_PEDANTIC = 1' in Makefile.RootCore
configuring Ext_RestFrames
  testing location: %%%LOCAL%%%
trying manual install
/share/home/kratsg/RestFramesAnalysis/RootCoreBin/download/Ext_RestFrames/tarball
/cvmfs/atlas.cern.ch/repo/sw/ASG/AnalysisBase/2.3.12/RootCore/scripts/external_download.sh: line 3: cd: crogan-RestFrames-cf9e419: No such file or directory
    failed: installation of Ext_RestFrames failed
RootCore: Error failed to find valid Ext_RestFrames installation in %%%LOCAL%%%
RootCore: Error failed to execute /share/home/kratsg/RestFramesAnalysis/Ext_RestFrames/cmt/precompile.RootCore
```

The way to fix this is to clean out your `download` folder inside `$ROOTCOREBIN` by something like

```
rm -rf RootCoreBin/download/Ext_RestFrames/*
```

## ATLAS CMake

*This section inspired by [this example ATLAS GitLab repo](https://gitlab.cern.ch/akraszna/GitAnalysisTest1)*.

To get this compiled in CMake (for ATLAS), you need to clone this with your other package sources, and then add the dependency correctly in your package's CMakeLists.txt file:

```cmake
# Build a shared library:
atlas_add_library( AnalysisPackage
   AnalysisPackage/*.h Root/*.h Root/*.cxx ${_dictionarySource}
   PUBLIC_HEADERS AnalysisPackage
   INCLUDE_DIRS ${RESTFRAMES_INCLUDE_DIRS}
   LINK_LIBRARIES ${RESTFRAMES_LIBRARIES} EventLoop )
add_dependencies( AnalysisPackage RestFrames )
```

The explicit dependency being added will make sure `RestFrames` gets build first before your package does. Lastly, your top-level project's `CMakeLists.txt` needs to be edited as follows:

```
--- CMakeLists.txt  2017-06-02 08:38:41.538140740 -0500
+++ ../CMakeLists.txt   2017-06-02 08:37:55.548415446 -0500
@@ -28,10 +28,13 @@
 find_package( AtlasCMake QUIET )

 # Find the project that we depend on:
 find_package( AnalysisBase )

+# Include the externals configuration:
+include( Ext_RestFrames/externals.cmake )
+
 # Set up CTest:
 atlas_ctest_setup()

 # Set up a work directory project:
 atlas_project( WorkDir 2.6.3
```

### More details

The logic of the build is the following:

 - `ExternalProject_Add(...)` tells autoconf/automake to install the project under `<build dir>/CMakeFiles/RestFramesBuild/`.
 - The `${RESTFRAMES_INCLUDE_DIRS}` and `${RESTFRAMES_LIBRARIES}` variables pick up the project from that directory, to avoid problems with inclusion orders. Adding `${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/include` to the include path of the build makes things quite complicated, apparently.
 - `ExternalProject_Add(...)` also copies the header and library files under `${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}`, so that they'd be available in the runtime environment when you just use the build directory. As the `setup.sh` file sets up the directories created by RestFrame.
 - The files are also declared to be installed with `install(...)`. This is needed in order to make the headers/libraries available in the RPM/tgz file created by CPack. This is important if you want to be able to run your code on the grid.
 - The package/library/executable that you want to use `RestFrames` in has to explicitly set up a dependency on the `RestFrames` target. And then use the `${RESTFRAMES_INCLUDE_DIRS}` and `${RESTFRAMES_LIBRARIES}` variables to set up its build. The dependency is important to make sure that in my example AnalysisPackage would not start building until the build of RestFrames is finished.

