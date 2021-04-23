#include <cmath>
#include <iostream>

#include "xAODRecoReader.h"
#include <TH2.h>
#include <TStyle.h>
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "SUSYTools/SUSYCrossSection.h"
#include "xAODCore/AuxContainerBase.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h"
#include "ElectronPhotonSelectorTools/AsgPhotonIsEMSelector.h"
#include "ElectronPhotonShowerShapeFudgeTool/ElectronPhotonShowerShapeFudgeTool.h"
#include "ElectronPhotonSelectorTools/egammaPIDdefs.h"
#include "IsolationCorrections/IsolationCorrectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"
#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"
#include "JetJvtEfficiency/JetJvtEfficiency.h"
#include "JetResolution/JERSmearingTool.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "xAODBTaggingEfficiency/BTaggingSelectionTool.h"
#include "TauAnalysisTools/TauOverlappingElectronLLHDecorator.h"
#include "TauAnalysisTools/TauSelectionTool.h"
#include "TauAnalysisTools/TauSmearingTool.h"
#include "TauAnalysisTools/Enums.h"
#include "SimpleAnalysisFramework/TruthEvent.h"
#include "SimpleAnalysisFramework/OutputHandler.h"

using std::vector;

#define CHECK( EXP )                                 \
  do {                                               \
    const StatusCode sc__ = EXP;                     \
    if (!sc__.isSuccess() ) {                        \
      throw std::runtime_error("Failed in call: "#EXP);   \
    }                                                \
  } while(0)

#define WARN_ONCE(warning)          \
  do {                              \
    static bool first=true;         \
    if (first) std::cout<<warning<<std::endl; \
    first=false;                              \
  } while(0)

DefineReader(xAODRecoReader);


static AsgElectronLikelihoodTool* initElecLH(std::string WP) {
  AsgElectronLikelihoodTool* LHTool = new AsgElectronLikelihoodTool(WP);
  CHECK(LHTool->setProperty("WorkingPoint", WP+"Electron"));
  CHECK(LHTool->initialize());
  return LHTool;
}

static CP::IsolationSelectionTool* initIso(std::string WP) {
  CP::IsolationSelectionTool *isoTool = new CP::IsolationSelectionTool(WP);
  if (WP!="FixedCutTight")  //only valid for electrons
    CHECK(isoTool->setProperty("MuonWP", WP));
  CHECK(isoTool->setProperty("ElectronWP", WP));
  if ((WP=="FixedCutTight")||(WP=="FixedCutLoose"))
    CHECK(isoTool->setProperty("PhotonWP", WP));
  CHECK(isoTool->initialize());
  return isoTool;
}

static BTaggingSelectionTool* initBtag(std::string WP) {
  BTaggingSelectionTool* tool = new BTaggingSelectionTool("BTaggingSelection"+WP+"Tool");
  CHECK(tool->setProperty("TaggerName", "MV2c10"));
  CHECK(tool->setProperty("OperatingPoint", "FixedCutBEff_"+WP));
  CHECK(tool->setProperty("JetAuthor", "AntiKt4EMTopoJets"));
  CHECK(tool->setProperty("FlvTagCutDefinitionsFileName", "xAODBTaggingEfficiency/13TeV/2016-20_7-13TeV-MC15-CDI-2017-01-31_v1.root"));
  CHECK(tool->initialize());
  return tool;
}

static TauAnalysisTools::TauSelectionTool* initTauSel(std::string WP) {
  TauAnalysisTools::TauSelectionTool* tool = new TauAnalysisTools::TauSelectionTool("TauSelection"+WP+"Tool" );
  CHECK(tool->setProperty("ConfigPath","SUSYTools/tau_selection_"+WP+".conf"));
  CHECK(tool->initialize());
  return tool;
}

static SG::AuxElement::Accessor<float> acc_filtHT("GenFiltHT");
static SG::AuxElement::Accessor<float> acc_filtMET("GenFiltMET");


xAODRecoReader::xAODRecoReader() : Reader() {}

void xAODRecoReader::AddOptions(po::options_description& desc ) {
  desc.add_options()
    ("readReco,r", "Use reconstructed quantities instead of truth")
    ;   
}

bool xAODRecoReader::isFileType(TFile * fh, po::variables_map&  vm) {
  if (fh->FindKey("CollectionTree") && vm.count("readReco")) {
    std::cout<<"Reading xAOD reco input"<<std::endl;
    return true;
  }
  return false;
}
void xAODRecoReader::Init()
{
  _event=new xAOD::TEvent(xAOD::TEvent::kClassAccess);
  _susytools= new ST::SUSYObjDef_xAOD("mySUSYTools");

  _ElecVeryLooseLH = initElecLH("VeryLooseLH");
  _ElecLooseLH     = initElecLH("LooseLH");
  _ElecLooseBLLH   = initElecLH("LooseBLLH");
  _ElecMediumLH    = initElecLH("MediumLH");
  _ElecTightLH     = initElecLH("TightLH");
  _isoTrackLoose             = initIso("LooseTrackOnly");
  _isoLoose                  = initIso("Loose");
  _isoGradient               = initIso("Gradient");
  _isoGradientLoose          = initIso("GradientLoose");
  _isoFixedCutLoose          = initIso("FixedCutLoose");
  _isoFixedCutTight          = initIso("FixedCutTight");
  _isoFixedCutTightTrackOnly = initIso("FixedCutTightTrackOnly");

  _ElecCalibTool= new CP::EgammaCalibrationAndSmearingTool("ElecCalibTool");
  CHECK(asg::setProperty(_ElecCalibTool, "decorrelationModel", "1NP_v1"));
  // CHECK(asg::setProperty(_ElecCalibTool, "useAFII", 1)); //FIXME: need to be configurable
  CHECK( asg::setProperty(_ElecCalibTool, "ESModel", "es2016data_mc15c"));
  CHECK(asg::setProperty(_ElecCalibTool, "randomRunNumber",EgammaCalibPeriodRunNumbersExample::run_2016));
  CHECK(_ElecCalibTool->initialize());

  _ElecIsoCorrTool = new CP::IsolationCorrectionTool("isoCorrTool");
  CHECK(asg::setProperty( _ElecIsoCorrTool,  "IsMC", true)); //FIXME: to be removed for later versions
  CHECK(_ElecIsoCorrTool->initialize());

  _MuonCalibTool= new CP::MuonCalibrationAndSmearingTool("MuonCalibTool");
  CHECK(_MuonCalibTool->initialize());
  _MuonSelectionTool=new CP::MuonSelectionTool("MuonSelectionTool");
  CHECK(_MuonSelectionTool->setProperty("MaxEta", 2.7));
  CHECK(_MuonSelectionTool->setProperty("MuQuality", 1));
  CHECK(_MuonSelectionTool->initialize());

  _tauSmearingTool=new TauAnalysisTools::TauSmearingTool("TauSmearingTool");
  CHECK(_tauSmearingTool->initialize());
  _tauElORdecorator = new TauAnalysisTools::TauOverlappingElectronLLHDecorator("TauElORDec");
  CHECK(_tauElORdecorator->initialize());
  _TauBDTLooseTool  = initTauSel("loose");
  _TauBDTMediumTool = initTauSel("medium");
  _TauBDTTightTool  = initTauSel("tight");

  _photonTight = new AsgPhotonIsEMSelector ( "PhotonTightIsEMSelector" );
  CHECK(_photonTight->setProperty("isEMMask",egammaPID::PhotonTight));
  CHECK(_photonTight->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMTightSelectorCutDefs.conf"));
  CHECK(_photonTight->initialize());
  _photonLoose = new AsgPhotonIsEMSelector ( "PhotonLooseIsEMSelector" );
  CHECK(_photonLoose->setProperty("isEMMask",egammaPID::PhotonLoose));
  CHECK(_photonLoose->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMTightSelectorCutDefs.conf"));
  CHECK(_photonLoose->initialize());

  _electronPhotonShowerShapeFudgeTool = new ElectronPhotonShowerShapeFudgeTool("FudgeMCTool");
  int FFset = 21; // for MC15 samples, which are based on a geometry derived from GEO-21 from 2015+2016 data
  CHECK(_electronPhotonShowerShapeFudgeTool->setProperty("Preselection",FFset));  //CHECK: is the right setting - is used in SUSYTools, I think
  CHECK(_electronPhotonShowerShapeFudgeTool->initialize());

  _JetCalibrationTool = new JetCalibrationTool("JetCalibTool");
  CHECK(_JetCalibrationTool->setProperty("JetCollection","AntiKt4EMTopo"));
  //FIXME: should maybe use JMS version in next two line
  //FIXME: also no option to use AF-II calibration
  CHECK(_JetCalibrationTool->setProperty("ConfigFile","JES_data2016_data2015_Recommendation_Dec2016.config"));
  CHECK(_JetCalibrationTool->setProperty("CalibSequence","JetArea_Residual_Origin_EtaJES_GSC"));
  CHECK(_JetCalibrationTool->setProperty("IsData",false));
  CHECK(_JetCalibrationTool->initializeTool("JetCalibTool"));
  _JetVertexTaggerTool = new JetVertexTaggerTool("JetVertexTaggerTool");
  CHECK(_JetVertexTaggerTool->setProperty("JVTFileName","JetMomentTools/JVTlikelihood_20140805.root"));
  CHECK(_JetVertexTaggerTool->initialize());
  _JetJvtEfficiencyTool = new CP::JetJvtEfficiency("JetJvtEfficiency");
  CHECK(_JetJvtEfficiencyTool->setProperty("SFFile","JetJvtEfficiency/Moriond2017/JvtSFFile_EM.root"));
  CHECK(_JetJvtEfficiencyTool->initialize());
  //FIXME: do we need to the JER smearing?
  _JetCleaningLooseTool = new JetCleaningTool("JetCleaningLooseTool", JetCleaningTool::LooseBad, false);
  _JetCleaningTightTool = new JetCleaningTool("JetCleaningTightTool", JetCleaningTool::TightBad, false);

  _fatJetCalibrationTool = new JetCalibrationTool("FatJetCalibTool");
  //FIXME: no option to use AF-II calibration
  CHECK(_fatJetCalibrationTool->setProperty("JetCollection","AntiKt10LCTopoTrimmedPtFrac5SmallR20"));
  CHECK(_fatJetCalibrationTool->setProperty("ConfigFile","JES_MC15recommendation_FatJet_Nov2016_QCDCombinationUncorrelatedWeights.config"));
  CHECK(_fatJetCalibrationTool->setProperty("CalibSequence","EtaJES_JMS"));
  CHECK(_fatJetCalibrationTool->setProperty("IsData",false));
  CHECK(_fatJetCalibrationTool->initializeTool("FatJetCalibTool"));

  _Btagging85Tool = initBtag("85");
  _Btagging77Tool = initBtag("77");
  _Btagging70Tool = initBtag("70");
  _Btagging60Tool = initBtag("60");
}


bool xAODRecoReader::processEvent(xAOD::TEvent *xaodEvent,xAOD::TStore */*store*/) {

  const xAOD::EventInfo* eventInfo = 0;
  if ( !xaodEvent->retrieve( eventInfo, "EventInfo").isSuccess() ) {
    throw std::runtime_error("Cannot read EventInfo");
  }
  int eventNumber = eventInfo->eventNumber();
  int mcChannel   = eventInfo->mcChannelNumber();
  int susy_part_id1 = 0;
  int susy_part_id2 = 0;
  int susy_process  = 0;

  const xAOD::TruthParticleContainer* truthparticles = 0;
  if ( xaodEvent->contains<xAOD::TruthParticleContainer>("TruthParticles")) {
    if ( !xaodEvent->retrieve( truthparticles, "TruthParticles").isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key TruthParticles");
    }
  } else {
    if ( !xaodEvent->retrieve( truthparticles, "TruthBSM").isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key TruthBSM");
    }
  }
  _susytools->FindSusyHardProc(truthparticles,susy_part_id1,susy_part_id2);
  if (susy_part_id2==0 && truthparticles->size()>1) susy_part_id2=truthparticles->at(1)->pdgId();
  if (susy_part_id1==0 && truthparticles->size()) susy_part_id1=truthparticles->at(0)->pdgId();
  if ((abs(susy_part_id1)>1000000) && (abs(susy_part_id1)>1000000)) //only consider BSM particles
    susy_process = SUSY::finalState(susy_part_id1,susy_part_id2);

  //FIXME: should recalculate MET
  const xAOD::MissingETContainer* metCont = 0;
  if ( !xaodEvent->retrieve(metCont, "MET_Reference_AntiKt4EMTopo").isSuccess() ){
    throw std::runtime_error("Could not retrieve met with key MET_MET_Reference_AntiKt4EMTopo");
  }
  const xAOD::MissingET* met = (*metCont)["FinalTrk"];
  const xAOD::VertexContainer* vertices(0);
  float primZ=0;
  int vtxIdx=0;
  if ( xaodEvent->retrieve( vertices, "PrimaryVertices" ).isSuccess() ) {
    for ( const auto& vx : *vertices ) {
      if (vx->vertexType() == xAOD::VxType::PriVtx) {
	primZ = vx->z();
	vtxIdx = vx->index();
	break;
      }
    }
  } else throw std::runtime_error("Could not find primary vertex");
  TruthEvent* event=new TruthEvent(met->sumet()/1000.,met->mpx()/1000.,met->mpy()/1000.);
  event->setChannelInfo(mcChannel,susy_process);

  TLorentzVector tlv(0.,0.,0.,0.);

  int idx=0;
  const xAOD::ElectronContainer* electrons(0);
  if ( !xaodEvent->retrieve(electrons, "Electrons") )
     throw std::runtime_error("Could not retrieve electrons");
  for (const auto& electron : *electrons) {
    if (!electron->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON)) continue;
    const xAOD::TrackParticle* track =  electron->trackParticle();
    double el_z0 = track->z0() + track->vz() - primZ;
    double z0sinTheta = el_z0 * TMath::Sin(electron->p4().Theta());
    double d0sig  = 9999;
    try {
      d0sig=xAOD::TrackingHelpers::d0significance( track , eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    } catch(...) {
      d0sig = 9999;
    }
    int id=0;
    if (fabs(d0sig)<5) id |= ED0Sigma5;
    if (fabs(z0sinTheta)<0.5) id |= EZ05mm;
    if (_ElecVeryLooseLH->accept(electron)) id |= EVeryLooseLH;
    if (_ElecLooseLH->accept(electron))     id |= ELooseLH;
    if (_ElecLooseBLLH->accept(electron))   id |= ELooseBLLH;
    if (_ElecMediumLH->accept(electron)) id |= EMediumLH;
    if (_ElecTightLH->accept(electron)) id |= ETightLH;
    xAOD::Electron *corrElec(0);
    if (_ElecCalibTool->correctedCopy(*electron,corrElec) != CP::CorrectionCode::Ok)
      throw std::runtime_error("Failed to apply electron correction");
    if (_ElecIsoCorrTool->applyCorrection(*corrElec)  != CP::CorrectionCode::Ok)
      throw std::runtime_error("Failed to apply electron isolation correction");
    if (_isoTrackLoose->accept(*corrElec)) id |= EIsoLooseTrack|EIsoBoosted; //FIXME: not sure this corresponds to 2015 definition
    if (_isoLoose->accept(*corrElec))      id |= EIsoLoose;
    if (_isoGradient->accept(*corrElec))   id |= EIsoGradient;
    if (_isoGradientLoose->accept(*corrElec))   id |= EIsoGradientLoose;
    if (_isoFixedCutLoose->accept(*corrElec))   id |= EIsoFixedCutLoose;
    if (_isoFixedCutTight->accept(*corrElec))   id |= EIsoFixedCutTight;
    if (_isoFixedCutTightTrackOnly->accept(*corrElec))   id |= EIsoFixedCutTightTrackOnly;
    //FIXME: add isolation for close-by objects?
    tlv.SetPtEtaPhiM(corrElec->pt()/1000.,corrElec->eta(),corrElec->phi(),electron->m()/1000.);
    event->addElectron(tlv,electron->charge(),id,0,idx++);
    delete corrElec;
  }
  idx=0;

  const xAOD::MuonContainer* muons(0);
  if ( !xaodEvent->retrieve(muons, "Muons") )
     throw std::runtime_error("Could not retrieve muons");
  for (const auto& muon : *muons) {
    xAOD::Muon *corrMuon(0);
    if (_MuonCalibTool->correctedCopy(*muon,corrMuon) != CP::CorrectionCode::OutOfValidityRange) {
      //      throw std::runtime_error("Failed to apply muon correction");
    }
    int id = MuNotCosmic|MuQoPSignificance; //FIXME: these flags should be calculated

    const xAOD::TrackParticle* track;
    if (corrMuon->muonType() == xAOD::Muon::SiliconAssociatedForwardMuon) {
      track = corrMuon->trackParticle(xAOD::Muon::CombinedTrackParticle);
      if (!track) continue; // don't treat SAF muons without CB track further (from SUSYTools)
    } else {
      track = corrMuon->primaryTrackParticle();
    }
    double mu_z0 = track->z0() + track->vz() - primZ;
    double z0sinTheta = mu_z0 * TMath::Sin(corrMuon->p4().Theta());
    double d0sig  = 9999;
    try {
      d0sig=xAOD::TrackingHelpers::d0significance( track , eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    } catch(...) {
      d0sig = 9999;
    }
    if (fabs(d0sig)<3) id |= MuD0Sigma3;
    if (fabs(z0sinTheta)<0.5) id |= MuZ05mm;
    xAOD::Muon::Quality muonQuality;
    muonQuality = _MuonSelectionTool->getQuality(*corrMuon);
    bool passesIDRequirements = _MuonSelectionTool->passedIDCuts(*corrMuon);
    bool passMuonCuts = _MuonSelectionTool->passedMuonCuts(*corrMuon);
    if (passMuonCuts&&passesIDRequirements) {
      if(muonQuality <= xAOD::Muon::VeryLoose) id|=MuVeryLoose;
      if(muonQuality <= xAOD::Muon::Loose) id|=MuLoose;
      if(muonQuality <= xAOD::Muon::Medium) id|=MuMedium;
      if(muonQuality <= xAOD::Muon::Tight) id|=MuTight;
      //FIXME: add HighPt muon selection bit
    }
    if (_isoTrackLoose->accept(*corrMuon)) id |= MuIsoLooseTrack|MuIsoBoosted; //FIXME: not sure this corresponds to 2015 definition
    if (_isoLoose->accept(*corrMuon))      id |= MuIsoLoose;
    if (_isoGradient->accept(*corrMuon))   id |= MuIsoGradient;
    if (_isoGradientLoose->accept(*corrMuon))   id |= MuIsoGradientLoose;
    if (_isoFixedCutLoose->accept(*corrMuon))   id |= MuIsoFixedCutLoose;
    if (_isoFixedCutTightTrackOnly->accept(*corrMuon)) id |= MuIsoFixedCutTightTrackOnly;

    tlv.SetPtEtaPhiM(corrMuon->pt()/1000.,corrMuon->eta(),corrMuon->phi(),muon->m()/1000.);
    event->addMuon(tlv,muon->charge(),id,0,idx++);
    delete corrMuon;
  }

  idx=0;
  const xAOD::TauJetContainer* taus = 0;
  if ( !xaodEvent->retrieve( taus, "TauJets").isSuccess() ) {
      throw std::runtime_error("Could not retrieve tau particles with key TauJets");
  }

  for ( const auto& tau : *taus) {
    if (!tau->container()->getConstStore()) continue; //taus were slimmed away

    if (fabs(tau->eta()) <= 2.5 && tau->nTracks() > 0) {
      xAOD::TauJet *corrTau(0);
      // FIXME: from SUSYTools: do truth matching first if required (e.g. for running on primary xAOD)
      // _tauTruthMatch->getTruth(input);

      if (_tauSmearingTool->correctedCopy(*tau,corrTau) != CP::CorrectionCode::Ok)
	throw std::runtime_error("Failed to apply tau correction");
      if (!_tauElORdecorator->decorate(*corrTau).isSuccess())
	throw std::runtime_error("Failed to apply tau decoration");
      int tauId=0;
      if (corrTau->nTracks()==1) tauId |= TauOneProng;
      if (corrTau->nTracks()==3) tauId |= TauThreeProng;
      if (_TauBDTLooseTool->accept(*corrTau)) tauId |= TauLoose;
      if (_TauBDTMediumTool->accept(*corrTau)) tauId |= TauMedium;
      if (_TauBDTTightTool->accept(*corrTau)) tauId |= TauTight;

      tlv.SetPtEtaPhiM(corrTau->pt()/1000.,corrTau->eta(),corrTau->phi(),corrTau->m()/1000.);
      event->addTau(tlv,corrTau->charge(),tauId,0,idx++);
      delete corrTau;
    }
  }

  idx=0;
  const xAOD::PhotonContainer* photons = 0;

  if ( !xaodEvent->retrieve( photons, "Photons").isSuccess() ) {
      throw std::runtime_error("Could not retrieve photon particles with key Photons");
  }

  for ( const auto& photon : *photons) {
    if ( !(photon->author() & (xAOD::EgammaParameters::AuthorPhoton + xAOD::EgammaParameters::AuthorAmbiguous)) ) continue;
    xAOD::Photon *corrPhoton(0);
    if (_ElecCalibTool->correctedCopy(*photon,corrPhoton) != CP::CorrectionCode::Ok)
      throw std::runtime_error("Failed to apply photon correction");
    if (_ElecIsoCorrTool->applyCorrection(*corrPhoton)  != CP::CorrectionCode::Ok)
      throw std::runtime_error("Failed to apply photon isolation correction");
    if (!corrPhoton->isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON)) {
      delete corrPhoton;
      continue;
    }
    if(!( (corrPhoton->OQ()&1073741824)!=0 ||
	  ( (corrPhoton->OQ()&134217728)!=0 &&
	    (corrPhoton->showerShapeValue(xAOD::EgammaParameters::Reta) >0.98
	     ||corrPhoton->showerShapeValue(xAOD::EgammaParameters::f1) > 0.4
	     ||(corrPhoton->OQ()&67108864)!=0)
	    ) ) ){
      //This a good photon wrt the photon cleaning
      if (_electronPhotonShowerShapeFudgeTool->applyCorrection(*corrPhoton) != CP::CorrectionCode::Ok)
	throw std::runtime_error("Failed to apply photon shower shape correction");
      int id = 0;
      if (_photonLoose->accept(corrPhoton))      id |= PhotonLoose;
      if (_photonTight->accept(corrPhoton))      id |= PhotonTight;
      if (_isoFixedCutLoose->accept(*corrPhoton)) id |= PhotonIsoFixedCutLoose;
      if (_isoFixedCutTight->accept(*corrPhoton)) id |= PhotonIsoFixedCutTight;
      tlv.SetPtEtaPhiM(corrPhoton->pt()/1000.,corrPhoton->eta(),corrPhoton->phi(),0);
      event->addPhoton(tlv,id,0,idx++);
      delete corrPhoton;
    }
  }

  idx=0;
  const xAOD::JetContainer* jets = 0;
  if ( !xaodEvent->retrieve( jets, "AntiKt4EMTopoJets").isSuccess() ) {
    throw std::runtime_error("Could not retrieve AntiKt4 jets");
  }
  for ( const auto& jet : *jets) {
    xAOD::Jet * calibJet = 0;
    _JetCalibrationTool->calibratedCopy(*jet,calibJet);

    int id=0;

    _JetVertexTaggerTool->updateJvt(*calibJet);
    if (_JetJvtEfficiencyTool->passesJvtCut(*calibJet)) id |= JVT50Jet;
    if (calibJet->pt()<20000 || (!(id&JVT50Jet))|| !_JetCleaningLooseTool->accept(*calibJet))  id |= LooseBadJet;
    if (calibJet->pt()<20000 || (!(id&JVT50Jet))|| !_JetCleaningTightTool->accept(*calibJet))  id |= TightBadJet;
    int nTracks = jet->getAttribute< std::vector<int> >(xAOD::JetAttribute::NumTrkPt500)[vtxIdx];
    if (nTracks<3) id |= LessThan3Tracks;
    if (_Btagging85Tool->accept(*calibJet)) id |= (BTag85MV2c10|BTag85MV2c20|BTag80MV2c20); //note faking old b-tagging points
    if (_Btagging77Tool->accept(*calibJet)) id |= (BTag77MV2c10|BTag77MV2c20);
    if (_Btagging70Tool->accept(*calibJet)) id |= (BTag70MV2c10|BTag70MV2c20);
    if (_Btagging60Tool->accept(*calibJet)) id |= (BTag60MV2c10);
    int flavor=0;
    if (jet->isAvailable<int>("HadronConeExclusionID"))
      flavor=jet->auxdata<int>("HadronConeExclusionID");
    else if (jet->isAvailable<int>("ConeTruthLabelID"))
      flavor=jet->auxdata<int>("ConeTruthLabelID");
    else if (jet->isAvailable<int>("PartonTruthLabelID"))
      flavor=abs(jet->auxdata<int>("PartonTruthLabelID"));
    else if (jet->isAvailable<int>("GhostBHadronsFinalCount")) {
      if (jet->auxdata<int>("GhostBHadronsFinalCount")) {
	flavor=5;
      } else if (jet->auxdata<int>("GhostCHadronsFinalCount")) {
	flavor=4;
      } else flavor=1;
    }
    if (flavor==4)       id |= TrueCJet;
    else if (flavor==5)  id |= TrueBJet;
    else if (flavor==15) id |= TrueTau;
    else                 id |= TrueLightJet; //FIXME: could be a pile-up jet as well

    tlv.SetPtEtaPhiM(calibJet->pt()/1000.,calibJet->eta(),calibJet->phi(),calibJet->m()/1000.);
    event->addJet(tlv,id,idx);
    delete calibJet;
  }

  idx=0;
  const xAOD::JetContainer* fatjets = 0;
  std::string fatjetName="AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets";

  if ( xaodEvent->contains<xAOD::JetContainer>(fatjetName)) {
    if ( !xaodEvent->retrieve( fatjets, fatjetName).isSuccess() ) {
      throw std::runtime_error("Could not retrieve fat jets with key "+fatjetName);
    }
    for ( const auto& jet : *fatjets) {
      xAOD::Jet * calibJet = 0;
      _fatJetCalibrationTool->calibratedCopy(*jet,calibJet);

      int id=0;
      tlv.SetPtEtaPhiM(calibJet->pt()/1000.,calibJet->eta(),calibJet->phi(),calibJet->m()/1000.);
      event->addFatJet(tlv,id,idx);
      delete calibJet;
    }
  }


  //Generator Filter HT (e.g. for ttbar/singleTop samples)
  float gen_ht=0.;
  if ( acc_filtHT.isAvailable(*(eventInfo)) ){
    gen_ht = eventInfo->auxdata<float>("GenFiltHT");
  }
  else{
    WARN_ONCE("Warning : No GenFiltHT decoration available. Setting HT to 0 for now...");
  }
  event->setGenHT( gen_ht/1000. );


  //Generator Filter MET (e.g. for ttbar/singleTop samples)
  float gen_met=0.;
  if ( acc_filtMET.isAvailable(*(eventInfo)) )
    gen_met = eventInfo->auxdata<float>("GenFiltMET");
  else
    WARN_ONCE("Warning : No GenFiltMETT decoration available. Setting MET to 0 for now...");
  //TODO : implement GenMET building from TruthParticle container as in xAODTruthReader //M.T.
  //...

  event->setGenMET( gen_met/1000. );


  //Get LHE3 weights
  const xAOD::TruthEventContainer* truthEvtCont;
  if( !xaodEvent->retrieve( truthEvtCont, "TruthEvents").isSuccess() )
    throw std::runtime_error("Could not retrieve truth event container with key TruthEvents");
  const xAOD::TruthEvent *truthevent = (*truthEvtCont)[0];
  const std::vector<float> weights  = truthevent->weights();
  xAOD::TruthEvent::PdfInfo pdfInfo = truthevent->pdfInfo();
  if (pdfInfo.valid()) {
    event->setPDFInfo(pdfInfo.pdgId1,pdfInfo.x1,pdfInfo.xf1,
		      pdfInfo.pdgId2,pdfInfo.x2,pdfInfo.xf2,
		      pdfInfo.Q);
  }
  event->setMCWeights(weights);

  _analysisRunner->processEvent(event,eventNumber);

  delete event;
  return true;
}

void xAODRecoReader::processFilesInternal(const std::vector<std::string>& inputFileNames, unsigned int nevents) {
  xAOD::TStore transientStorage;
  transientStorage.setActive();
  TFile *inFile=0;
  unsigned int procEvents = 0;
  for(const auto& inName : inputFileNames) {
    delete inFile;
    std::cout<<"Now reading: "<<inName<<std::endl;
    inFile=TFile::Open(inName.c_str());
    if ( ! _event->readFrom(inFile).isSuccess() ) {
      throw std::runtime_error("Could not connect TEvent to file !");
    }
    Long64_t numEntries=_event->getEntries();
    for(Long64_t index = 0; index<numEntries; index++) {
      ++procEvents;
      if (procEvents>nevents) break;
      Long64_t entry = _event->getEntry(index);
      if (entry<0) break;
      if (index%10000==0) std::cout<<"at: "<<index<<"/"<<numEntries<<std::endl;
      processEvent(_event,&transientStorage);
      transientStorage.clear();
    }
  }
}

xAODRecoReader::~xAODRecoReader() {
  return;
  delete _susytools;
  delete _event;
}
