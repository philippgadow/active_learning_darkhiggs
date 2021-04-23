#
# Build configuration for the SimpleAnalysisFramework package
#

# Set the name of the package:
atlas_subdir( SimpleAnalysisFramework )

# Set up which packages this package depends on:
atlas_depends_on_subdirs(
   PUBLIC
   Control/AthToolSupport/AsgTools
   PhysicsAnalysis/AnalysisCommon/PATCore
   PhysicsAnalysis/D3PDTools/EventLoop
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonFourMomentumCorrection
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonSelectorTools
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonShowerShapeFudgeTool
   PhysicsAnalysis/ElectronPhotonID/IsolationCorrections
   PhysicsAnalysis/ElectronPhotonID/IsolationSelection
   PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonSelectorTools
   PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonMomentumCorrections
   PhysicsAnalysis/Interfaces/AsgAnalysisInterfaces
   PhysicsAnalysis/JetMissingEtID/JetSelectorTools
   PhysicsAnalysis/JetMissingEtID/JetJvtEfficiency
   PhysicsAnalysis/JetMissingEtID/JetMomentTools
   PhysicsAnalysis/JetMissingEtID/JetResolution
   PhysicsAnalysis/JetMissingEtID/JetCalibTools
   PhysicsAnalysis/JetTagging/JetTagPerformanceCalibration/xAODBTaggingEfficiency
   PhysicsAnalysis/TauID/TauAnalysisTools/TauAnalysisTools
   PhysicsAnalysis/SUSYPhys/SUSYTools
   PhysicsAnalysis/UpgradePhys/SmearingFunctions/UpgradePerformanceFunctions
   Reconstruction/MVAUtils
   PRIVATE
   Control/AthContainers
   Control/CxxUtils
   Control/xAODRootAccess
   Tools/PathResolver
   Trigger/TrigAnalysis/TrigDecisionTool
   Trigger/TrigAnalysis/TriggerMatchingTool
   Event/xAOD/xAODBase
   Event/xAOD/xAODCore
   Event/xAOD/xAODJet
   Event/xAOD/xAODMuon
   Event/xAOD/xAODMissingET
   Event/xAOD/xAODEgamma
   Event/xAOD/xAODEventInfo
   Event/xAOD/xAODTracking
   Event/xAOD/xAODTruth
   PhysicsAnalysis/AnalysisCommon/HDF5Utils
)

# External(s) used by the package:
find_package( ROOT COMPONENTS Math Core Hist Physics REQUIRED )
find_package( Boost COMPONENTS filesystem thread system )
find_package(Eigen)
find_package(Lhapdf)
find_package(HDF5 1.10.1 REQUIRED COMPONENTS CXX C)
find_package(lwtnn)
find_package(onnxruntime)


# Generate a dictionary for the library:
atlas_add_root_dictionary( SimpleAnalysisFrameworkLib SimpleAnalysisFrameworkLibDictSrc
   ROOT_HEADERS src/LinkDef.h
   EXTERNAL_PACKAGES ROOT )

add_compile_options(-Werror=int-in-bool-context)

# Build the shared library of the package:
atlas_add_library( SimpleAnalysisFrameworkLib
   SimpleAnalysisFramework/*.h src/*.h src/*.cxx ${SimpleAnalysisFrameworkLibDictSrc}
   PUBLIC_HEADERS SimpleAnalysisFramework
   PRIVATE_INCLUDE_DIRS
   ${ROOT_INCLUDE_DIRS}
   ${EIGEN_INCLUDE_DIRS}
   ${Boost_INCLUDE_DIRS}
   ${LWTNN_INCLUDE_DIRS}
   PUBLIC_INCLUDE_DIRS
   ${RESTFRAMES_INCLUDE_DIRS}
   ${ONNXRUNTIME_INCLUDE_DIRS}
   LINK_LIBRARIES
   ${EIGEN_LIBRARIES}
   ${Boost_LIBRARIES}
   ${LWTNN_LIBRARIES}
   ${ROOT_LIBRARIES}
   ${LHAPDF_LIBRARIES}
   ${RESTFRAMES_LIBRARIES}
   ${ONNXRUNTIME_LIBRARIES}
   AthLinks
   AthContainers
   CxxUtils
   MuonMomentumCorrectionsLib
   IsolationSelectionLib IsolationCorrectionsLib
   ElectronPhotonShowerShapeFudgeToolLib
   ElectronPhotonFourMomentumCorrectionLib
   UpgradePerformanceFunctionsLib
   PathResolver
   JetCalibToolsLib
   JetMomentToolsLib
   JetSelectorToolsLib
   TauAnalysisToolsLib
   MVAUtils

   SUSYToolsLib
   xAODMuon
   xAODBTaggingEfficiencyLib
   xAODRootAccessInterfaces xAODCore xAODEventFormat xAODEventInfo xAODEgamma  xAODTracking xAODTruth
)

if(DO_TRUTHTAGGING)
   target_compile_definitions( SimpleAnalysisFrameworkLib PUBLIC PACKAGE_BTaggingTruthTagging )
endif(DO_TRUTHTAGGING)

add_dependencies(
  SimpleAnalysisFrameworkLib
  RestFrames
)

atlas_add_executable (simpleAnalysis util/simpleAnalysis.cxx
   INCLUDE_DIRS
   ${RESTFRAMES_INCLUDE_DIRS}
   LINK_LIBRARIES
   SimpleAnalysisLib
   SimpleAnalysisCodesLib
   SimpleAnalysisFrameworkLib
)

atlas_add_executable (slimMaker util/slimMaker.cxx
   INCLUDE_DIRS
   ${RESTFRAMES_INCLUDE_DIRS}
   LINK_LIBRARIES
   SimpleAnalysisFrameworkLib
)

atlas_install_data(data/*.json)
atlas_install_data(data/*.xml)
atlas_install_data(data/*.root)
