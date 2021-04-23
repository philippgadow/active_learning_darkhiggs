#include "SimpleAnalysisFramework/TruthEvent.h"
#include "TruthSmear.h"
#include "SimpleAnalysisFramework/AnalysisClass.h"
#include <RootCore/Packages.h>

DefineSmearer(TruthSmear, "smear,s", "Smearing with UPF. Has comma-separated list smearing options (use help to see full list of options)");

#include <UpgradePerformanceFunctions/UpgradePerformanceFunctions.h>
#include "CalibrationDataInterface/CalibrationDataContainer.h"

using namespace Upgrade;

#define CHECK( EXP )                                 \
  do {                                               \
    const StatusCode sc__ = EXP;                     \
    if (!sc__.isSuccess() ) {                        \
      throw std::runtime_error("Failed in call: "#EXP);   \
    }                                                \
  } while(0)

void TruthSmear::init(std::vector<std::string>& options ) {
  smearElectrons=smearMuons=smearTaus=smearPhotons=smearJets=smearMET=true;
  addPileupJets=useHGTD0=useHGTD1=useMuonHighEta=optPhotons=false;
  
  std::string mu="None";
  std::string slayout="None";
  int seed=12345;
  UpgradePerformanceFunctions::UpgradeLayout layout = UpgradePerformanceFunctions::UpgradeLayout::run2;
  for(const auto& option : options){
    if (option.find("help")==0) {
      std::cout<<"Options for smearing:"<<std::endl;
      std::cout<<" layout=<value> - set layout (*required*, use \"run2\" or \"run4\")"<<std::endl;
      std::cout<<" noElectrons    - do not smear electrons"<<std::endl;
      std::cout<<" noMuons        - do not smear muons"<<std::endl;
      std::cout<<" noTaus         - do not smear taus"<<std::endl;
      std::cout<<" noPhotons      - do not smear photons"<<std::endl;
      std::cout<<" noJets         - do not smear jets"<<std::endl;
      std::cout<<" noMET          - do not smear MET"<<std::endl;
      std::cout<<" addPileupJets  - add pileup jets (add for HL-LHC conditions)"<<std::endl;
      std::cout<<" useHGTD0       - use HGTD v0 configuration"<<std::endl;
      std::cout<<" useHGTD1       - use HGTD v1 configuration"<<std::endl;
      std::cout<<" useMuonHighEta - Use muons out to eta of 4"<<std::endl;
      std::cout<<" optPhotons     - assume optimistic photons resolution"<<std::endl;
      //      std::cout<<" mu=<value>     - choice pile-up level (required)"<<std::endl;
      //std::cout<<"                  Only mu=200 is allowed for now"<<std::endl;
      std::cout<<" seed=<value>   - set seed value (default: "<<seed<<std::endl;
    }
    if (option.find("noElectrons")==0) smearElectrons=false;
    if (option.find("noMuons")==0)     smearMuons=false;
    if (option.find("noTaus")==0)      smearTaus=false;
    if (option.find("noPhotons")==0)   smearPhotons=false;
    if (option.find("noJets")==0)      smearJets=false;
    if (option.find("noMET")==0)       smearMET=false;
    if (option.find("addPileupJets")==0) addPileupJets=true;
    if (option.find("useHGTD0")==0)    useHGTD0=true;
    if (option.find("useHGTD1")==0)    useHGTD1=true;
    if (option.find("useMuonHighEta")==0) useMuonHighEta=true;
    if (option.find("optPhotons")==0)  optPhotons=true;
    /* FIXME: for now mu hardcoded as smearing functions only fully support mu=200
    if (option.find("mu=")==0) {
      mu=option.substr(3);
    }
    */
    if (option.find("seed=")==0) {
      seed=stoi(option.substr(5));
    }
    if (option.find("layout=")==0) {
      slayout=option.substr(7); 
      if(slayout == "run2"){
	layout = UpgradePerformanceFunctions::UpgradeLayout::run2;
	mu = "200"; //dummy value for now
      }
      else if(slayout == "run4"){
	layout = UpgradePerformanceFunctions::UpgradeLayout::Step1p6;
	mu = "200";
      }      
      else {
	throw std::runtime_error("Only \"run4\" and \"run2\" layouts supported");
      }
    }
  }

  if (slayout == "None") 
    	throw std::runtime_error("specify \"run4\" or \"run2\" layout to enable smearing");

  std::cout<<"Smearing with mu="<<mu<<" and seed="<<seed<<std::endl;

  if (mu!="200" && slayout!="run2") throw std::runtime_error("Unsupported pile-up level. Only mu=200 currently supported");

  std::string METhistfile = "UpgradePerformanceFunctions/CalibArea-00-01/sumetPU_mu200_ttbar_gold.root";
  if (slayout=="run2") 
    METhistfile = "UpgradePerformanceFunctions/CalibArea-00-01/met_resol_tst_run2.root";

  m_upgrade = new UpgradePerformanceFunctions("LooseSmear");
  m_upgradeMuonTight = new UpgradePerformanceFunctions("TightSmear");
  m_upgradeMuonHighPt = new UpgradePerformanceFunctions("HighSmear");
  CHECK(m_upgrade->setProperty("Layout",layout));
  CHECK(m_upgradeMuonTight->setProperty("Layout",layout));
  CHECK(m_upgradeMuonHighPt->setProperty("Layout",layout));
  CHECK(m_upgrade->setProperty("AvgMu",mu));
  CHECK(m_upgradeMuonTight->setProperty("AvgMu",mu));
  CHECK(m_upgradeMuonHighPt->setProperty("AvgMu",mu));
  CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::looseElectron));
  CHECK(m_upgrade->setProperty("ElectronRandomSeed",seed));
  CHECK(m_upgrade->setProperty("MuonWorkingPoint",UpgradePerformanceFunctions::looseMuon));
  CHECK(m_upgradeMuonTight->setProperty("MuonWorkingPoint",UpgradePerformanceFunctions::tightMuon));
  CHECK(m_upgradeMuonHighPt->setProperty("MuonWorkingPoint",UpgradePerformanceFunctions::highPtMuon));
  CHECK(m_upgrade->setProperty("UseMuonHighEta",useMuonHighEta));
  CHECK(m_upgradeMuonTight->setProperty("UseMuonHighEta",useMuonHighEta));
  CHECK(m_upgradeMuonHighPt->setProperty("UseMuonHighEta",useMuonHighEta));
  CHECK(m_upgrade->setProperty("PhotonWorkingPoint",UpgradePerformanceFunctions::tightPhoton));
  CHECK(m_upgrade->setProperty("PhotonNoiseScaling",1.0)); 
  if (optPhotons)
    CHECK(m_upgrade->setProperty("PhotonNoiseScaling",0.375)); 
  CHECK(m_upgrade->setProperty("TauRandomSeed",seed));
  CHECK(m_upgrade->setProperty("JetRandomSeed",seed));
  CHECK(m_upgrade->setProperty("METRandomSeed",seed));
  CHECK(m_upgrade->setProperty("METFile",METhistfile));
  CHECK(m_upgrade->loadMETHistograms());
  CHECK(m_upgrade->setProperty("PileupRandomSeed",seed));
  CHECK(m_upgrade->setProperty("UseTrackConfirmation",true));
  CHECK(m_upgrade->setProperty("PileupJetThresholdMeV",30000.));
  CHECK(m_upgrade->setProperty("PileupEfficiency",UpgradePerformanceFunctions::PU));
  CHECK(m_upgrade->setProperty("JVT_PU_Efficiency",0.02));
  //CHECK(m_upgrade->setProperty("PileupTemplatesPath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/UpgradePerformanceFunctions/"));
  //  CHECK(m_upgrade->setProperty("PhotonFakeFile","UpgradePerformanceFunctions/CalibArea-00-01/PhotonFakes.root"));
  CHECK(m_upgrade->initPhotonFakeHistograms());
  if(slayout=="run2") {
    CHECK(m_upgrade->setProperty("FlavourTaggingCalibrationFile","xAODBTaggingEfficiency/13TeV/2016-20_7-13TeV-MC15-CDI-2017-01-31_v1.root"));
    // do more complete flavour efficiency
    TFile *fh;
    std::string calibFile = PathResolverFindCalibFile("xAODBTaggingEfficiency/13TeV-Online/2020-21-13TeV-ConditionalOnlineMV2GivenOfflineMV2c10-MC16e.root");
    fh=TFile::Open(calibFile.c_str());
    setTaggerHist(fh,BTag85DL1,85,"DL1");
    setTaggerHist(fh,BTag77DL1,77,"DL1");
    setTaggerHist(fh,BTag70DL1,70,"DL1");
    setTaggerHist(fh,BTag60DL1,60,"DL1");
    setTaggerHist(fh,BTag85DL1r,85,"DL1r");
    setTaggerHist(fh,BTag77DL1r,77,"DL1r");
    setTaggerHist(fh,BTag70DL1r,70,"DL1r");
    setTaggerHist(fh,BTag60DL1r,60,"DL1r");
    setTaggerHist(fh,BTag85MV2c10,85,"MV2c10");
    setTaggerHist(fh,BTag77MV2c10,77,"MV2c10");
    setTaggerHist(fh,BTag70MV2c10,70,"MV2c10");
    setTaggerHist(fh,BTag60MV2c10,60,"MV2c10");
    setTaggerHist(fh,BTag85MV2c10|BTag60OnlineMV2,60,"ConditionalOnlineMV2GivenOfflineMV2c1085");
    setTaggerHist(fh,BTag77MV2c10|BTag60OnlineMV2,60,"ConditionalOnlineMV2GivenOfflineMV2c1077");
    setTaggerHist(fh,BTag70MV2c10|BTag60OnlineMV2,60,"ConditionalOnlineMV2GivenOfflineMV2c1070");
    setTaggerHist(fh,BTag60MV2c10|BTag60OnlineMV2,60,"ConditionalOnlineMV2GivenOfflineMV2c1060");
    setTaggerHist(fh,BTag85MV2c10|BTag40OnlineMV2,40,"ConditionalOnlineMV2GivenOfflineMV2c1085");
    setTaggerHist(fh,BTag77MV2c10|BTag40OnlineMV2,40,"ConditionalOnlineMV2GivenOfflineMV2c1077");
    setTaggerHist(fh,BTag70MV2c10|BTag40OnlineMV2,40,"ConditionalOnlineMV2GivenOfflineMV2c1070");
    setTaggerHist(fh,BTag60MV2c10|BTag40OnlineMV2,40,"ConditionalOnlineMV2GivenOfflineMV2c1060");

    setTaggerHist(fh,BTag85DL1,85,"DL1", true);
    setTaggerHist(fh,BTag77DL1,77,"DL1", true);
    setTaggerHist(fh,BTag70DL1,70,"DL1", true);
    setTaggerHist(fh,BTag60DL1,60,"DL1", true);
    setTaggerHist(fh,BTag85DL1r,85,"DL1r", true);
    setTaggerHist(fh,BTag77DL1r,77,"DL1r", true);
    setTaggerHist(fh,BTag70DL1r,70,"DL1r", true);
    setTaggerHist(fh,BTag60DL1r,60,"DL1r", true);
  }
  else
    CHECK(m_upgrade->setProperty("FlavourTaggingCalibrationFile","UpgradePerformanceFunctions/CalibArea-00-01/flavor_tags_v2.1.root"));
  CHECK(m_upgrade->setProperty("UseHGTD0",useHGTD0));
  CHECK(m_upgrade->setProperty("UseHGTD1",useHGTD1));
  CHECK( m_upgrade->initialize() );
  CHECK( m_upgradeMuonTight->initialize() );
  CHECK( m_upgradeMuonHighPt->initialize() );
  m_random.SetSeed(seed);

}

void TruthSmear::setTaggerHist(TFile *ff, int tagger, int operating_point, std::string dirName, bool trackJet) {
  // get MC efficiency from CDI files!
  static std::vector<std::string> flavours = {"B", "C", "Light", "T"};
  if (trackJet) {
    for(int flavour = 0; flavour<4; flavour++) {
      TString hfName = Form("/%s/AntiKtVR30Rmax4Rmin02TrackJets/FixedCutBEff_%d/%s/", dirName.c_str(), operating_point, flavours[flavour].c_str());
      auto dir= ff->GetDirectory(hfName,false,"cd");
      TH2D* fEff= 0;
      if (dir) {
        dir->cd();
        Analysis::CalibrationDataHistogramContainer* cHCont = (Analysis::CalibrationDataHistogramContainer*) ff->Get(hfName + "default_Eff");
        if (cHCont) {
    fEff = (TH2D*)cHCont->GetValue("result");
        } 
      }
      m_fEffs_trackjets[tagger][flavour]=fEff;
    }
  } else {
    for(int flavour = 0; flavour<4; flavour++) {
      TString hfName = Form("/%s/AntiKt4EMTopoJets_BTagging201810/FixedCutBEff_%d/%s/", dirName.c_str(), operating_point, flavours[flavour].c_str());
      auto dir= ff->GetDirectory(hfName,false,"cd");
      TH2D* fEff= 0;
      if (dir) {
        dir->cd();
        Analysis::CalibrationDataHistogramContainer* cHCont = (Analysis::CalibrationDataHistogramContainer*) ff->Get(hfName + "default_Eff");
        if (cHCont) {
    fEff = (TH2D*)cHCont->GetValue("result");
        } 
      }
      m_fEffs[tagger][flavour]=fEff;
    }
  }
}

float TruthSmear::getFlavourTagEfficiency(double ptGeV, float eta, int flavour, int tagger, bool trackJet) {
  // Modified from UpgradePerformanceFunctions 
  //
  // the function returns b/c/l-tagging efficincies obtained using ttbar samples

  // flavour is an index:
  // 0 (to get the b-tag efficiency)
  // 1 (to get the c-tag efficiency)
  // 2 (to get the mistag rate)
  // 3 (to get the PU-tag efficiency)

  // pT range protection
  //  if (ptGeV<20) ptGeV = 20;
  if (!trackJet && ptGeV < 20) return 0.;
 
  TH2D* fEffs;
  if (trackJet)
    fEffs = m_fEffs_trackjets[tagger][flavour];
  else
    fEffs = m_fEffs[tagger][flavour];

  if (fEffs == 0) return 0;
  // CDI files have switched from eta-pt binning to pt-eta binning.
  // There are more bins in pt than in eta, so check which kind of file we have based on that.
  bool invert_axes = (fEffs->GetNbinsY() > fEffs->GetNbinsX()) ? true : false;
  
  //eta&pt bins
  auto eta_bin = invert_axes ? fEffs->GetXaxis()->FindBin(eta) : fEffs->GetYaxis()->FindBin(eta);
  auto pt_bin  = invert_axes ? fEffs->GetYaxis()->FindBin(ptGeV) : fEffs->GetXaxis()->FindBin(ptGeV);
  if (eta_bin == 0) eta_bin = 1;
  if (pt_bin == 0)  pt_bin = 1;
  if (invert_axes) {
    if (eta_bin > fEffs->GetNbinsX()) eta_bin = fEffs->GetNbinsX();
    if (pt_bin  > fEffs->GetNbinsY()) pt_bin  = fEffs->GetNbinsY();
  }
  else {
    if (eta_bin > fEffs->GetNbinsY()) eta_bin = fEffs->GetNbinsY();
    if (pt_bin  > fEffs->GetNbinsX()) pt_bin  = fEffs->GetNbinsX();
  }
  
  return invert_axes ? fEffs->GetBinContent(eta_bin, pt_bin) : fEffs->GetBinContent(pt_bin, eta_bin);
}


TruthEvent *TruthSmear::smearEvent(AnalysisEvent *
event
) {
  auto electrons  = event->getElectrons(1.,4.2); // Filter on pT, eta and "ID"
  auto muons      = event->getMuons(1.,4.2);
  auto taus       = event->getTaus(20.,4.0,TauOneProng)+event->getTaus(20.,4.0,TauThreeProng);
  auto photons    = event->getPhotons(5.,4.2);
  auto jets       = event->getJets(10.,5.2);
  auto fatjets    = event->getFatJets(10.,5.2);
  auto trackjets    = event->getTrackJets(10.,2.5);
  auto met        = event->getMET();
  auto sumet      = event->getSumET();
  
  double met_x = met.Px();
  double met_y = met.Py();

  TruthEvent* smeared = new TruthEvent(sumet,met_x,met_y);
  smeared->setChannelInfo(event->getMCNumber(),event->getSUSYChannel());
  smeared->setGenMET(event->getGenMET());
  smeared->setGenHT(event->getGenHT());
  smeared->setMCWeights(event->getMCWeights());
  smeared->setTruth(event);
  smeared->setSmearing(true);
  smeared->setPDFInfo(event->getPDF_id1(),event->getPDF_x1(),event->getPDF_pdf1(),
		      event->getPDF_id2(),event->getPDF_x2(),event->getPDF_pdf2(),
		      event->getPDF_scale());
  smeared->addHSTruthList(event->getHSTruth(0,999.,0));

  const auto layout = m_upgrade->getLayout();
  
  for(const auto& electron : electrons) {
    if (smearElectrons){
      CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::looseElectron));
      CHECK(m_upgrade->setProperty("Layout",UpgradePerformanceFunctions::UpgradeLayout::Step1p6));
      float looseEff = m_upgrade->getElectronEfficiency(electron.Pt()*1000., electron.Eta());
      CHECK(m_upgrade->setProperty("Layout",layout));
      float looseFlip = m_upgrade->getElectronChargeFlipProb(electron.Pt()*1000., std::min(fabs(electron.Eta()),2.4699));
      CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::mediumElectron));
      CHECK(m_upgrade->setProperty("Layout",UpgradePerformanceFunctions::UpgradeLayout::Step1p6));
      float mediumEff = m_upgrade->getElectronEfficiency(electron.Pt()*1000., electron.Eta());
      CHECK(m_upgrade->setProperty("Layout",layout));
      float mediumFlip = m_upgrade->getElectronChargeFlipProb(electron.Pt()*1000., std::min(fabs(electron.Eta()),2.4699));
      CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::tightElectron));
      CHECK(m_upgrade->setProperty("Layout",UpgradePerformanceFunctions::UpgradeLayout::Step1p6));
      float tightEff = m_upgrade->getElectronEfficiency(electron.Pt()*1000., electron.Eta());
      CHECK(m_upgrade->setProperty("Layout",layout));
      float tightFlip = m_upgrade->getElectronChargeFlipProb(electron.Pt()*1000., std::min(fabs(electron.Eta()),2.4699));
      CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::looseElectron));
      float prob = m_random.Uniform(1.0);

      if (prob<looseEff) {
	float flipProb = looseFlip;
	int electronID=EVeryLooseLH | ELooseLH | ELooseBLLH | ED0Sigma5 | EZ05mm; //BP: not actually LH in upgrade
	//FIXME: for now we always mark electrons as isolated as there is no smearing for that
	electronID |= EIsoGood&(~EGood);
	if (prob<mediumEff) {
	  electronID |= EMediumLH;
	  flipProb = mediumFlip;  //not quite right as loose and medium flip probabilities become a bit too small
	}
	if (prob<tightEff) {
	  electronID |= ETightLH;
	  flipProb = tightFlip;
	}
	float electron_e = m_upgrade->getElectronSmearedEnergy(electron.E()*1000., electron.Eta())/1000.; 
	TLorentzVector eLV;
	eLV.SetPtEtaPhiM(electron.Pt()*electron_e/electron.E(),electron.Eta(),electron.Phi(),0.000510998910);
	int echarge = electron.charge();
	if (m_random.Uniform(1.0)< flipProb) echarge*=-1;
	smeared->addElectron(eLV,echarge,electronID,electron.motherID(),electron.index());
	if (m_random.Uniform(1.0)<m_upgrade->getElectronToPhotonFakeRate(electron.Pt()*1000., electron.Eta())) 
	  smeared->addPhoton(eLV,PhotonIsoGood,electron.motherID(),-2);
      }
    }
    else {
      smeared->addElectron(electron,electron.charge(),electron.id(),electron.motherID(),electron.index());
    }
  }
  for(const auto& muon : muons) {
    if (smearMuons){
      float looseEff  = m_upgrade->getMuonEfficiency(muon.Pt()*1000., muon.Eta(), muon.Phi());
      float tightEff  = m_upgradeMuonTight->getMuonEfficiency(muon.Pt()*1000., muon.Eta(), muon.Phi());
      float highPtEff = m_upgradeMuonHighPt->getMuonEfficiency(muon.Pt()*1000., muon.Eta(), muon.Phi());
      float prob = m_random.Uniform(1.0);
      if (prob<std::max(looseEff,highPtEff)) { //FIXME: presumably not 100% overlap between loose and highPt
	int muonID = MuVeryLoose | MuD0Sigma3 | MuZ05mm | MuNotCosmic | MuQoPSignificance;
	//FIXME: no smearing for isolation yet
	muonID |= MuIsoGood&(~MuGood);
	if (prob<looseEff) muonID |= MuLoose;
	if (prob<tightEff) muonID |= MuTight|MuMedium;
	if (prob<highPtEff) muonID |= MuHighPt;
	float muonUnsmearedPt = muon.Pt()*1000.;
	float qoverpt = muon.charge() / muonUnsmearedPt;
	float muonQOverPtResolution = m_upgrade->getMuonQOverPtResolution(muonUnsmearedPt, muon.Eta());
	qoverpt += m_random.Gaus(0., muonQOverPtResolution);
	float muonSmearedPt = fabs(1./qoverpt)/1000.;
	int muonCharge = 1;
	if (qoverpt<0) muonCharge = -1;
	TLorentzVector mLV;
	mLV.SetPtEtaPhiM(muonSmearedPt, muon.Eta(), muon.Phi(), 0.1056583715);
      	smeared->addMuon(mLV,muonCharge,muonID,muon.motherID(),muon.index());
      }
    }
    else {
      smeared->addMuon(muon,muon.charge(),muon.id(),muon.motherID(),muon.index());
    }
  }
  for(const auto& tau : taus) {
    if (smearTaus){
      short prong=1;
      if (tau.pass(TauThreeProng)) prong=3;
      float looseEff = m_upgrade->getTauEfficiency(tau.Pt()*1000., tau.Eta(), prong, 0);
      float mediumEff = m_upgrade->getTauEfficiency(tau.Pt()*1000., tau.Eta(), prong, 1);
      float tightEff = m_upgrade->getTauEfficiency(tau.Pt()*1000., tau.Eta(), prong, 2);
      float prob=m_random.Uniform(1.0);
      if (prob<looseEff) {
	int tauID = TauLoose | (tau.pass(TauThreeProng)?TauThreeProng:TauOneProng);
	if (prob<mediumEff) tauID |= TauMedium;
	if (prob<tightEff) tauID |= TauTight;
	float tau_E = m_upgrade->getTauSmearedEnergy(tau.E()*1000., tau.Eta(), prong)/1000.; 
	TLorentzVector tauLV;
	tauLV.SetPtEtaPhiM(tau.Pt()*tau_E/tau.E(),tau.Eta(),tau.Phi(),1.777682);
	smeared->addTau(tauLV,tau.charge(),tauID,tau.motherID(),tau.index());
      }
    }
    else {
      smeared->addTau(tau,tau.charge(),tau.id(),tau.motherID(),tau.index());
    }
  }
  for(const auto& photon : photons) {
    if (smearPhotons){
      float eff = m_upgrade->getPhotonEfficiency(photon.Pt()*1000./*, photon.Eta()*/);
      if (m_random.Uniform(1.0)<eff) {
	TLorentzVector pLV;
	pLV.SetPtEtaPhiM(photon.Pt()*1000.,photon.Eta(),photon.Phi(),0.);
	pLV = m_upgrade->getPhotonSmearedVector(&pLV);
	pLV.SetPtEtaPhiM(pLV.Pt()/1000.,pLV.Eta(),pLV.Phi(),0.);
	smeared->addPhoton(pLV,photon.id(),photon.motherID(),photon.index());
      }
    }
    else {
      smeared->addPhoton(photon,photon.id(),photon.motherID(),photon.index());
    }
  }
  for(const auto& jet : jets) {
    if (smearJets){
      float jetpt  = jet.Pt();
      float jetptMeV = jetpt*1000.;
      if (jetpt<1500) jetpt = m_upgrade->getJetSmearedEnergy(jetptMeV,jet.Eta(),true)/1000.; // FIXME: can only smear jets below 1500 GeV
      jetptMeV = jetpt*1000.;
      float jeteta = jet.Eta();
      float jetphi = jet.Phi();
      float jetE   = jet.E()*jetpt/jet.Pt();

      char jetType = 'L';
      int jetTypeIdx = 2;
      if (jet.pass(TrueBJet)) {
	jetType = 'B';
	jetTypeIdx = 0;
      }
      if (jet.pass(TrueCJet)) {
	jetType = 'C';
	jetTypeIdx = 1;
      }

      float tag=m_random.Uniform(1.0);
      float otag=m_random.Uniform(1.0);
      int jetid=GoodJet;

      if(m_upgrade->getLayout() != UpgradePerformanceFunctions::UpgradeLayout::run2){
	float tagEff70 = m_upgrade->getFlavourTagEfficiency(jetpt*1000., jeteta, jetType, "mv2c10", 70, m_upgrade->getPileupTrackConfSetting());
	float tagEff85 = m_upgrade->getFlavourTagEfficiency(jetpt*1000., jeteta, jetType, "mv2c10", 85, m_upgrade->getPileupTrackConfSetting());

	if (tag < tagEff70) jetid|=(BTag70MV2c10|BTag70MV2c20); 
	if (tag < tagEff85) jetid|=(BTag85MV2c10|BTag85MV2c20);

	if (tag < tagEff70) jetid=GoodBJet;
	
	if (jet.pass(TrueLightJet)) jetid|=TrueLightJet;
	if (jet.pass(TrueCJet))     jetid|=TrueCJet;
	if (jet.pass(TrueBJet))     jetid|=TrueBJet;
	if (jet.pass(TrueTau))      jetid|=TrueTau;

      }
      else{ //Run2 settings
	float tagEff60MV2c10 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag60MV2c10);
	float tagEff70MV2c10 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag70MV2c10);
	float tagEff77MV2c10 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag77MV2c10);
	float tagEff85MV2c10 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag85MV2c10);
 
	float tagEff60MV2c10Online60 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag60MV2c10|BTag60OnlineMV2);
	float tagEff70MV2c10Online60 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag70MV2c10|BTag60OnlineMV2);
	float tagEff77MV2c10Online60 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag77MV2c10|BTag60OnlineMV2);
	float tagEff85MV2c10Online60 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag85MV2c10|BTag60OnlineMV2);

	float tagEff60MV2c10Online40 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag60MV2c10|BTag40OnlineMV2);
	float tagEff70MV2c10Online40 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag70MV2c10|BTag40OnlineMV2);
	float tagEff77MV2c10Online40 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag77MV2c10|BTag40OnlineMV2);
	float tagEff85MV2c10Online40 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag85MV2c10|BTag40OnlineMV2);

	if (tag < tagEff60MV2c10) {
	  jetid|=BTag60MV2c10|BTag70MV2c10|BTag70MV2c20|BTag77MV2c10|BTag77MV2c20|BTag80MV2c20|BTag85MV2c10|BTag85MV2c20; 
	  if (otag < tagEff60MV2c10Online40) jetid|=BTag40OnlineMV2;
	  if (otag < tagEff60MV2c10Online60) jetid|=BTag60OnlineMV2;
	} else if (tag < tagEff70MV2c10) {
	  jetid|=BTag70MV2c10|BTag70MV2c20|BTag77MV2c10|BTag77MV2c20|BTag80MV2c20|BTag85MV2c10|BTag85MV2c20; 
	  if (otag < (tagEff70MV2c10*tagEff70MV2c10Online40-tagEff60MV2c10*tagEff60MV2c10Online40)/(tagEff70MV2c10-tagEff60MV2c10)) jetid|=BTag40OnlineMV2;
	  if (otag < (tagEff70MV2c10*tagEff70MV2c10Online60-tagEff60MV2c10*tagEff60MV2c10Online60)/(tagEff70MV2c10-tagEff60MV2c10)) jetid|=BTag60OnlineMV2;
	} else if (tag < tagEff77MV2c10) {
	  jetid|=BTag77MV2c10|BTag77MV2c20|BTag80MV2c20|BTag85MV2c10|BTag85MV2c20; 
	  if (otag < (tagEff77MV2c10*tagEff77MV2c10Online40-tagEff70MV2c10*tagEff70MV2c10Online40)/(tagEff77MV2c10-tagEff70MV2c10)) jetid|=BTag40OnlineMV2;
	  if (otag < (tagEff77MV2c10*tagEff77MV2c10Online60-tagEff70MV2c10*tagEff70MV2c10Online60)/(tagEff77MV2c10-tagEff70MV2c10)) jetid|=BTag60OnlineMV2;
	} else if (tag < tagEff85MV2c10) {
	  jetid|=BTag85MV2c10|BTag85MV2c20; 
	  if (otag < (tagEff85MV2c10*tagEff85MV2c10Online40-tagEff77MV2c10*tagEff77MV2c10Online40)/(tagEff85MV2c10-tagEff77MV2c10)) jetid|=BTag40OnlineMV2;
	  if (otag < (tagEff85MV2c10*tagEff85MV2c10Online60-tagEff77MV2c10*tagEff77MV2c10Online60)/(tagEff85MV2c10-tagEff77MV2c10)) jetid|=BTag60OnlineMV2;
	}

	float tagEff60DL1 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag60DL1);
	float tagEff70DL1 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag70DL1); 
	float tagEff77DL1 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag77DL1); 
	float tagEff85DL1 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag85DL1); 

	if (tag < tagEff60DL1) jetid|=BTag60DL1; 
	if (tag < tagEff70DL1) jetid|=BTag70DL1; 
	if (tag < tagEff77DL1) jetid|=BTag77DL1; 
	if (tag < tagEff85DL1) jetid|=BTag85DL1; 

	float tagEff60DL1r = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag60DL1r);
	float tagEff70DL1r = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag70DL1r); 
	float tagEff77DL1r = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag77DL1r); 
	float tagEff85DL1r = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag85DL1r); 

	if (tag < tagEff60DL1r) jetid|=BTag60DL1r; 
	if (tag < tagEff70DL1r) jetid|=BTag70DL1r; 
	if (tag < tagEff77DL1r) jetid|=BTag77DL1r; 
	if (tag < tagEff85DL1r) jetid|=BTag85DL1r; 

	if (jet.pass(TrueLightJet)) jetid|=TrueLightJet;
	if (jet.pass(TrueCJet))     jetid|=TrueCJet;
	if (jet.pass(TrueBJet))     jetid|=TrueBJet;
	if (jet.pass(TrueTau))      jetid|=TrueTau;

      }

      if (addPileupJets) {
	if ( (jetptMeV) < m_upgrade->getPileupJetPtThresholdMeV()) jetpt=0;
	else {
	  if (fabs(jet.Eta())<3.8) {
	    float trackEff = m_upgrade->getTrackJetConfirmEff(jetptMeV,jet.Eta(), "HS");
	    float hsProb = m_random.Uniform(1.0);
	    if (hsProb > trackEff) jetpt=0; // FIXME: should couple this to JVT flag
	  }
	}
      }
      if (jetpt) {
	TLorentzVector j;
	j.SetPtEtaPhiE(jetpt,jeteta,jetphi,jetE);
	smeared->addJet(j,jetid,jet.index());

	//introduce jets faking electrons, photons and taus unless jet is from one these
	if ((AnalysisClass::minDR(jet,electrons)>0.05) && (AnalysisClass::minDR(jet,photons)>0.05) && (AnalysisClass::minDR(jet,taus)>0.05)) {

	  // Add jets faking electrons    //MT : need to implement run2 version yet
	  if( m_upgrade->getLayout() != UpgradePerformanceFunctions::UpgradeLayout::run2){
	    CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::looseElectron));
	    float looseRate = m_upgrade->getElectronFakeRate(jet.Pt()*1000, jet.Eta());
	    CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::mediumElectron));
	    float mediumRate = m_upgrade->getElectronFakeRate(jet.Pt()*1000, jet.Eta());
	    CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::tightElectron));
	    float tightRate = m_upgrade->getElectronFakeRate(jet.Pt()*1000, jet.Eta());
	    CHECK(m_upgrade->setProperty("ElectronWorkingPoint",UpgradePerformanceFunctions::looseElectron));
	    float prob=m_random.Uniform(1.0);
	    if (prob<looseRate) {
	      int electronID=EVeryLooseLH | ELooseLH | ELooseBLLH | ED0Sigma5 | EZ05mm; //BP: not actually LH in upgrade
	      //FIXME: for now we always mark electrons as isolated as there is no smearing for that
	      electronID |= EIsoGradientLoose | EIsoBoosted | EIsoFixedCutTight | EIsoLooseTrack | EIsoLoose | EIsoGradient | EIsoFixedCutLoose | EIsoFixedCutTightTrackOnly;
	      if (prob<mediumRate) electronID |= EMediumLH;
	      if (prob<tightRate) electronID |= ETightLH;
	      float electron_e = m_upgrade->getElectronFakeRescaledEnergy(jet.E()*1000., jet.Eta())/1000.;
	      TLorentzVector eLV;
	      eLV.SetPtEtaPhiM(jet.Pt()*electron_e/jet.E(),jet.Eta(),jet.Phi(),0.000510998910);
	      int echarge = 1;
	      if (m_random.Uniform(1.0)<0.5) echarge = -1;
	      smeared->addElectron(eLV,echarge,EIsoGood,0,-1);
	    }
	    // Add jets faking tau
	    if (jetpt<20 || fabs(jet.Eta())>4.0) continue;
	    float looseEff  = m_upgrade->getTauFakeRate(jetpt*1000, jet.Eta(), 1, 0);
	    float mediumEff = m_upgrade->getTauFakeRate(jetpt*1000, jet.Eta(), 1, 1);
	    float tightEff  = m_upgrade->getTauFakeRate(jetpt*1000, jet.Eta(), 1, 2);
	    int tauID = TauLoose|TauOneProng;
	    prob = m_random.Uniform(1.0);
	    if (prob>=looseEff) { //if not one prong, it could be three prong
	      prob = prob - looseEff;
	      looseEff  = m_upgrade->getTauFakeRate(jetpt*1000, jet.Eta(), 3, 0);
	      mediumEff = m_upgrade->getTauFakeRate(jetpt*1000, jet.Eta(), 3, 1);
	      tightEff  = m_upgrade->getTauFakeRate(jetpt*1000, jet.Eta(), 3, 2);
	      tauID = TauLoose|TauThreeProng;
	    }
	    if (prob<looseEff) {
	      if (prob<mediumEff) tauID |= TauMedium;
	      if (prob<tightEff) tauID |= TauTight;
	      float tau_et = jetpt; // FIXME: no tau smearing exists yet
	      TLorentzVector tauLV;
	      tauLV.SetPtEtaPhiM(tau_et,jet.Eta(),jet.Phi(),1.777682);
	      int taucharge = 1;
	      if (m_random.Uniform(1.0)<0.5) taucharge = -1;
	      smeared->addTau(tauLV,taucharge,tauID,0,-1);
	    }
	    // Add jets faking photon
	    if (m_random.Uniform(1.0)<m_upgrade->getPhotonFakeRate(jet.Pt()*1000/*, jet.Eta()*/)) {
	      float photon_et = m_upgrade->getPhotonFakeRescaledET(jet.Pt()*1000/*., jet.Eta()*/)/1000.;
	      TLorentzVector pLV;
	      pLV.SetPtEtaPhiM(photon_et,jet.Eta(),jet.Phi(),0.0);
	      smeared->addPhoton(pLV,PhotonIsoGood,0,-1);
	    }
	  }
	}
      }

    }
    else {
      smeared->addJet(jet,jet.id(),jet.index());
    }
  }
  if (addPileupJets) {
    for (const auto& pujet : m_upgrade->getPileupJets()) {
      float trackEff = 1.0;
      if (fabs(pujet.Eta())<3.8) trackEff=m_upgrade->getTrackJetConfirmEff(pujet.Pt(), pujet.Eta(), "PU");
      float puProb = m_random.Uniform(1.0);
            
      if (puProb > trackEff) continue; // FIXME: should couple this to JVT flag
      float tagEff70 = m_upgrade->getFlavourTagEfficiency(pujet.Pt(), pujet.Eta(), 'P', "mv2c10", 70, m_upgrade->getPileupTrackConfSetting());
      float tagEff85 = m_upgrade->getFlavourTagEfficiency(pujet.Pt(), pujet.Eta(), 'P', "mv2c10", 85, m_upgrade->getPileupTrackConfSetting());
      float tag=m_random.Uniform(1.0);
      int jetid=GoodJet;
      if (tag<tagEff85) jetid|=BTag85MV2c20; //FIXME: check if this should set other working points too
      if (tag<tagEff70) jetid=GoodBJet;
      smeared->addJet(pujet.Px()/1000., pujet.Py()/1000., pujet.Pz()/1000., pujet.E()/1000., jetid, -1);
      // Add jets faking photon
      if (m_random.Uniform(1.0)<m_upgrade->getPhotonPileupFakeRate(pujet.Pt()/*, pujet.Eta()*/)) {
	float photon_et = m_upgrade->getPhotonPileupFakeRescaledET(pujet.Pt()/*., pujet.Eta()*/)/1000.;
	TLorentzVector pLV;
	pLV.SetPtEtaPhiM(photon_et,pujet.Eta(),pujet.Phi(),0.0);
	smeared->addPhoton(pLV,PhotonIsoGood,0,-3);
      }
    }
  }
  for(const auto& jet : fatjets) {
    smeared->addFatJet(jet,jet.id(),jet.index()); //FIXME: for now there is no smearing for fat jets
  }

  for(const auto& jet : trackjets) {
    float jetpt  = jet.Pt();
    float jeteta = jet.Eta();
    float jetphi = jet.Phi();
    float jetE   = jet.E();

    char jetType = 'L';
    int jetTypeIdx = 2;
    if (jet.pass(TrueBJet)) {
      jetType = 'B';
      jetTypeIdx = 0;
    }
    if (jet.pass(TrueCJet)) {
      jetType = 'C';
      jetTypeIdx = 1;
    }

    float tag=m_random.Uniform(1.0);
    int jetid=GoodJet;

    float tagEff60DL1 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag60DL1, true);
    float tagEff70DL1 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag70DL1, true); 
    float tagEff77DL1 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag77DL1, true); 
    float tagEff85DL1 = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag85DL1, true); 
    if (tag < tagEff60DL1) jetid|=BTag60DL1; 
    if (tag < tagEff70DL1) jetid|=BTag70DL1; 
    if (tag < tagEff77DL1) jetid|=BTag77DL1; 
    if (tag < tagEff85DL1) jetid|=BTag85DL1; 

    float tagEff60DL1r = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag60DL1r, true);
    float tagEff70DL1r = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag70DL1r, true); 
    float tagEff77DL1r = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag77DL1r, true); 
    float tagEff85DL1r = getFlavourTagEfficiency(jetpt, jeteta, jetTypeIdx, BTag85DL1r, true); 
    if (tag < tagEff60DL1r) jetid|=BTag60DL1r; 
    if (tag < tagEff70DL1r) jetid|=BTag70DL1r; 
    if (tag < tagEff77DL1r) jetid|=BTag77DL1r; 
    if (tag < tagEff85DL1r) jetid|=BTag85DL1r; 

    if (jet.pass(TrueLightJet)) jetid|=TrueLightJet;
    if (jet.pass(TrueCJet))     jetid|=TrueCJet;
    if (jet.pass(TrueBJet))     jetid|=TrueBJet;
    if (jet.pass(TrueTau))      jetid|=TrueTau;

    smeared->addTrackJet(jet,jetid,jet.index());
  }

  if (smearMET) {
    if(m_upgrade->getLayout() != UpgradePerformanceFunctions::UpgradeLayout::run2){ //upgrade settings

      auto smearedMET = m_upgrade->getMETSmeared(sumet*1000., met_x*1000., met_y*1000.);
      met_x = smearedMET.first/1000.;
      met_y = smearedMET.second/1000.;
  
    }
    else{ //Run2 settings
     
      //compute smeared MET now
      auto SDelectrons  = smeared->getElectrons(1.,4.2,EVeryLooseLH); // Filter on pT, eta and "ID"
      auto SDmuons      = smeared->getMuons(1.,4.2,MuLoose);
      auto SDjets       = smeared->getJets(10.,5.2,GoodJet);

      auto SMETjets       = AnalysisClass::overlapRemoval(SDjets,SDelectrons,0.2,NOT(BTag85MV2c10));
      auto SMETelectrons  = AnalysisClass::overlapRemoval(SDelectrons,SMETjets,0.4);
      SMETjets = AnalysisClass::overlapRemoval(SMETjets, SDmuons, 0.4, LessThan3Tracks); 
      auto SMETmuons = AnalysisClass::overlapRemoval(SDmuons, SMETjets, 0.4);

      //get vectorial sum of all objects and track the truth sum too
      AnalysisObject tsMET(0,0,0,0,0,0,MET,0,0);
      AnalysisObject sMET(0,0,0,0,0,0,MET,0,0);
      for(int ii=0; ii < (int)SMETjets.size(); ii++) {
	sMET += SMETjets[ii];
	tsMET += smeared->getTruthParticle(SMETjets[ii]);
      }
      for(int ii=0; ii < (int)SMETelectrons.size(); ii++) {
	sMET += SMETelectrons[ii];
	tsMET += smeared->getTruthParticle(SMETelectrons[ii]);
      }
      for(int ii=0; ii < (int)SMETmuons.size(); ii++) {
	sMET += SMETmuons[ii];
	tsMET += smeared->getTruthParticle(SMETmuons[ii]);
      }
      //smear TST 
      AnalysisObject ptHard = tsMET + met;
      TVector3 TST = m_upgrade->getTSTsmearing(ptHard.Vect()*1000.) * 0.001; //smearing back in GeV

      met_x = -sMET.Px()-TST.X() + ptHard.Px();
      met_y = -sMET.Py()-TST.Y() + ptHard.Py();

    }
    smeared->setMET(met_x, met_y);
  }
  
  return smeared;
}
