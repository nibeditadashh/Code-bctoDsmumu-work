

// -*- C++ -*-
//
// Package:    BcToDsMuMu
// Class:      BcToDsMuMu
// 
//=================================================
// Original by: Chandiprasad Kar                  |
//<chandiprasad.kar@cern.ch>                      |
//Major variables are added 05 Sept, 2020         |
//=================================================

// system include files
#include <memory>

// user include files
#include "bctodsmumu_analysis/BcToDsMuMuPAT/src/BcToDsMuMu.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <vector>
#include <utility>
#include <string>
#include <iostream>
//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//
const static unsigned kNBR_CLOSED_TRACKS = 20;
const double PI = 3.141592653589793;
//Define histogram

// struct HistArgs{
//   char name[128];
//   char title[128];
//   int n_bins;
//   double x_min;
//   double x_max;
// };

// enum HistName{
//   h_Dsmass,
//   h_Dspt,
  
//   h_bcmass,                                                                                                                           

//   kHistNameSize
// };

// HistArgs hist_args[kHistNameSize] = {
//   // name, title, n_bins, x_min, x_max                                                                                                              
//   {"h_Dsmass", "D_{s} mass; M(KK#pi) [GeV/^{2}]", 100, 0, 20},   
//   {"h_Dspt", "D_{s} p_{T}; p_{T} [GeV]", 100, 0, 40},   
//   {"h_bcmass", "B_{c} mass; M(D_{s}#mu#mu) [GeV]", 100, 0, 20},                                                                                           

// };

// //--------------------                                                                                                                                    
// // Define histograms                                                                                                                                       
// //--------------------                                                                                                                                    
// TH1F *histos[kHistNameSize];



//
// constructors and destructor
//
BcToDsMuMu::BcToDsMuMu(const edm::ParameterSet& iConfig):

  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  //  trakCollection_label(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  triggerPrescales_ (consumes<pat::PackedTriggerPrescales> (iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),  
  trigTable_( iConfig.getParameter<std::vector<std::string> >("TriggerNames")),   
  puToken_(consumes<std::vector< PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PuInfoTag"))), 

  v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondaryVerticesPtr"))),	       

  genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),
  dr0(0),dr1(0),
  dpt0(0),dpt1(0),
  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
  
  //************************ ****************************************************
  
  //*******************************************************
  nB(0), nMu(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0),
  
  B_Ds_mass(0), B_Ds_px(0), B_Ds_py(0), B_Ds_pz(0),
  B_Ds_pt1(0), B_Ds_px1(0), B_Ds_py1(0), B_Ds_pz1(0), 
  B_Ds_pt2(0), B_Ds_px2(0), B_Ds_py2(0), B_Ds_pz2(0), 
  B_Ds_pt3(0), B_Ds_px3(0), B_Ds_py3(0), B_Ds_pz3(0), 

  B_Ds_px1_track(0), B_Ds_py1_track(0), B_Ds_pz1_track(0), 
  B_Ds_px2_track(0), B_Ds_py2_track(0), B_Ds_pz2_track(0), 
  B_Ds_px3_track(0), B_Ds_py3_track(0), B_Ds_pz3_track(0), 
  
  pi1dxy(0), pi2dxy(0), pi1dz(0), pi2dz(0),
  pi1dxy_e(0), pi2dxy_e(0), pi1dz_e(0), pi2dz_e(0),
  B_Ds_charge1(0), B_Ds_charge2(0),B_Ds_charge3(0),

  B_mumu_mass(0), B_mumu_px(0), B_mumu_py(0), B_mumu_pz(0),
  B_mumu_pt1(0), B_mumu_px1(0), B_mumu_py1(0), B_mumu_pz1(0), 
  B_mumu_pt2(0), B_mumu_px2(0), B_mumu_py2(0), B_mumu_pz2(0), 
  B_mumu_charge1(0), B_mumu_charge2(0),

  B_Ds_chi2(0), B_mumu_chi2(0), B_chi2(0), B_chi2dof(0),
  B_Prob(0), B_mumu_Prob(0), B_Ds_Prob(0),
  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),
  
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0), VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  VDecayVtxXYE(0), VDecayVtxXZE(0), VDecayVtxYZE(0), 
  pVtxIPX(0),  pVtxIPY(0),  pVtxIPZ(0),
  pVtxIPXE(0),  pVtxIPYE(0),  pVtxIPZE(0),  pVtxIPCL(0),
  pVtxIPXYE(0),  pVtxIPXZE(0),  pVtxIPYZE(0),
  
  B_l3d(0),  B_l3dE(0),  B_lxy(0), B_lxyE(0),
  B_cosalpha(0),   B_cosalphaxy(0), alpha(0),  B_treco(0),   B_trecoe(0),  B_trecoxy(0), B_trecoxye(0),
  B_pvip(0), B_pviperr(0), B_pvips(0), B_pvlzip(0), B_pvlziperr(0), B_pvlzips(0),
  B_pv2ip(0), B_pv2iperr(0), B_pv2ips(0), B_pv2lzip(0), B_pv2lziperr(0), B_pv2lzips(0),

  B_l3d_pv2(0),  B_l3dE_pv2(0),
  B_iso(0), B_mum_iso(0), B_mup_iso(0), B_pi1_iso(0),B_pi2_iso(0),
  
  istruemum(0), istruemup(0), istruekp(0), istruekm(0), istruebs(0),
  bunchXingMC(0), numInteractionsMC(0), trueNumInteractionsMC(0),
  run(0), event(0),
  lumiblock(0)
  
{
   //now do what ever initialization is needed
}

BcToDsMuMu::~BcToDsMuMu()
{

}

// ------------ method called to for each event  ------------
void BcToDsMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> theV0PtrHandle;
  iEvent.getByToken(v0PtrCollection_,  theV0PtrHandle);
  
  
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  //edm::Handle<std::vector<pat::GenericParticle> > thePATTrackHandle;  
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  std::cout<<"size of the muon collection "<<thePATMuonHandle->size()<<" track size "<<theV0PtrHandle->size()<<" packed track "<<thePATTrackHandle->size()<<endl;
  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  fMagneticField = fieldHandle.product();

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genParticles_, pruned);

  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);
  // if(isMC_ ){
  //   edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
  //   iEvent.getByToken(puToken_, PupInfo);  
  //   for (std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); PVI++)
  //     {
  // 	bunchXingMC->push_back(PVI->getBunchCrossing());
  // 	numInteractionsMC->push_back(PVI->getPU_NumInteractions());
  // 	trueNumInteractionsMC->push_back(PVI->getTrueNumInteractions());
  //     }
  // }
  //std::cout << "Inside analyze " << std::endl;
  //***************************************
  // MC Gen information
  //***************************************
  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_ds_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_ds_vtx.SetXYZ(0.,0.,0.);
  gen_b_ct = -9999.;
 if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == 541)   && (dau->status() == 3) ) {
std::cout << "Bc Candidate "<< dau->pdgId() << std::endl; 
        foundit++;
        gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
        gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
        for (size_t k=0; k<dau->numberOfDaughters(); k++) {
          const reco::Candidate *gdau = dau->daughter(k);
          if (gdau->pdgId()==431 ) {
std::cout << "Ds Candidate "<< gdau->pdgId() << std::endl;
            foundit++;
            gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());

            int nm=0;
            for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
              const reco::Candidate *mm = gdau->daughter(l);
if (mm->pdgId()==321) { foundit++;
                if (mm->status()!=1) {
                  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
                    const reco::Candidate *mu = mm->daughter(m);
                    if (mu->pdgId()==321 ) {
                      nm++;
                      gen_pion1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
                      break;
                    }
                  }
                } else {
                 gen_pion1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
                nm++;
                }
              }
              if (mm->pdgId()==-321) { foundit++;
                if (mm->status()!=1) {
                  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
                    const reco::Candidate *mu = mm->daughter(m);
                    if (mu->pdgId()==-321 ) {
                      nm++;
                      gen_pion2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
                        break;
                    }
                  }
                } else {
                 gen_pion2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
                nm++;
                }
              }
 if (mm->pdgId()==211) { foundit++;
                if (mm->status()!=1) {
                  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
                    const reco::Candidate *mu = mm->daughter(m);
                    if (mu->pdgId()==211 ) {
                      nm++;
                      gen_pion3_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
                        break;
                    }
                  }
                } else {
                gen_pion3_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
                nm++;
                }
              }
      } // Ds Daughters
}// Ds
 	if (gdau->pdgId()==13 ) {
std::cout << "mu1 Candidate "<< gdau->pdgId() << std::endl;
            foundit++;
            gen_muon1_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
		}
	if (gdau->pdgId()==-13 ) {
std::cout << "mu2 Candidate "<< gdau->pdgId() << std::endl;
            foundit++;
            gen_muon2_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
		}
        }                                                                                                                                        
      }
 if (foundit>=7) break;
    }                                                                                                                                            
    if (foundit!=7) {
      gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_ks_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_pion1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_pion2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_pion3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_b_vtx.SetXYZ(0.,0.,0.);
      gen_jpsi_vtx.SetXYZ(0.,0.,0.);
      gen_b_ct = -9999.;
      std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
  }

  nB = 0; nMu = 0;
  if ( OnlyGen_ ) {
    tree_->Fill();
    return;
  }

  //****************************************************************
  //Triggers

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales>            triggerPrescales;
  edm::Handle<edm::TriggerResults> hltTriggerResults;

  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  iEvent.getByToken(triggerResults_Label, hltTriggerResults ); 
 
  std::string cctk_2 = "hltLowMassNonResonantTkAllConeTracksIter";
  std::stringstream myString;
  std::vector<float> obj_eta, obj_phi, obj_pt;
  std::vector<int> obj_charge;

  //bool foundOneTrig = false;
  const edm::TriggerNames &names = iEvent.triggerNames(*hltTriggerResults);
  for (unsigned int i = 0, n = hltTriggerResults->size(); i < n; ++i) {
    for (unsigned int it = 0; it < trigTable_.size(); it++){
      if (names.triggerName(i).find(trigTable_[it]) != std::string::npos && hltTriggerResults->accept(i))
	{ 
	  //triggernames->push_back(names.triggerName(i));
	  //triggerprescales->push_back(triggerPrescales->getPrescaleForIndex(i));
	  // cout<<"Trigger name "<<names.triggerName(i)<<" Prescale "<<triggerPrescales->getPrescaleForIndex(i)<<endl;
	  //foundOneTrig = true;
	}
    }
  }
  //if ( iEvent.isRealData() && !foundOneTrig) return;
     
  if (triggerObjects.isValid()) {
    //std::cout << "will try to match trigger object with track " << triggerObjects->size() << endl;    
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
      obj.unpackFilterLabels(iEvent,*hltTriggerResults);
      std::string cc1 = obj.collection();
      
      for (unsigned int i = 0; i < trigTable_.size(); i++)
        {
	  myString.clear(); myString.str(""); myString << trigTable_[i] << "*";
	  if ( obj.hasPathName(myString.str().c_str(), true, true) )
            {

	      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	      std::vector<std::string> pathNamesLast = obj.pathNames(true);
	            
	      //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
	            
	      if ( (cc1.find(cctk_2) != std::string::npos ) ){//|| (cc1.find(cctk_1) != std::string::npos) || ( cc1.find(cctk_2) != std::string::npos)) {
		
		obj_eta.push_back(obj.eta());
		obj_phi.push_back(obj.phi());
		obj_pt.push_back(obj.pt());
		obj_charge.push_back(obj.charge());
		
		// std::cout << myString.str().c_str() << std::endl;
		// std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
		// std::cout << "\t   Charge:   "<<obj.charge()<<std::endl;
		// std::cout << "\t   Collection: " << obj.collection() << std::endl;
		// std::cout << "\t   Type IDs:   ";
		// for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
		// std::cout << "\t   Filters:    ";
		// for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
		// std::cout << std::endl;
	      }
	    }
	}
    }
  }
    

  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  reco::Vertex bestVtxBS;


  // get primary vertex
  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_Label, recVtxs);
  edm::Handle<reco::VertexCollection> pvHandle_;
  iEvent.getByToken(primaryVertices_Label,pvHandle_);
  bestVtx = *(recVtxs->begin());
  
  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);
  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 

  nVtx = recVtxs->size();

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event(); 
  
  //std::cout<<"After Primary vertex check "<< std::endl;
  //Let's begin by looking for J/psi->mu+mu-
  unsigned int nMu_tmp = thePATMuonHandle->size();
 
  for(edm::View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {      
      for(edm::View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; // <-------------------------------------
	  
	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }
	  
	  if(iMuon1->track()->pt()<4.0) continue;
	  //if(iMuon1->track()->pt()<2.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;
	  //if(iMuon2->track()->pt()<2.0) continue;
	  
	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
	  
	  //Let's check the vertex and mass
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	  // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************

	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  //ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;

	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }

	  KinematicParticleVertexFitter fitter;   

	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }

	  if (!psiVertexFitTree->isValid()) 
	    {
	      std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }

	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();//masa del J/psi
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();//vertice del J/psi
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }

	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.01)
	    {
	      continue;
	    }
	  
	  //some loose cuts go here
	  
	  if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	  if(psi_vFit_noMC->currentState().mass()<1.0 || psi_vFit_noMC->currentState().mass()>4.8) continue;
	  //  ***************  		 
	  int ix = 0;
	  vector <int> trklist;
	  trklist.reserve(200);
	  for(ix = 0; (unsigned)ix<thePATTrackHandle->size(); ix++){
	    pat::PackedCandidate trackView((*thePATTrackHandle)[ix]);
	    if(trackView.charge()==0) continue;  
	    if(trackView.pt()<1.0) continue;  
	    if( !(trackView.hasTrackDetails())) continue;
	    if(!(trackView.bestTrack()->quality(reco::Track::highPurity))) continue; 
	    if ( IsTheSamePFtk(trackView,*iMuon1) || IsTheSamePFtk(trackView,*iMuon2) ) continue;
	    trklist.push_back(ix);
	  }
	  std::cout<<"Track size "<<trklist.size()<<std::endl;
	  if(trklist.size()<3)continue;
	  for(unsigned int ix1 = 0; ix1<trklist.size(); ix1++){	  
	    for(unsigned int ix2 = ix1+1; ix2<trklist.size(); ix2++){	      
	      for(unsigned int ix3 = 0; ix3<trklist.size(); ix3++){
	  	if((ix1 == ix2) || (ix1 == ix3) || (ix2 == ix3)) continue;
	  	pat::PackedCandidate iTrack1((*thePATTrackHandle)[trklist.at(ix1)]);
	  	pat::PackedCandidate iTrack2((*thePATTrackHandle)[trklist.at(ix2)]);
	  	pat::PackedCandidate iTrack3((*thePATTrackHandle)[trklist.at(ix3)]);
	  	if(iTrack1.charge()==iTrack2.charge()) continue;
		reco::TransientTrack pion1TT((*theB).build(iTrack1.bestTrack()));
		reco::TransientTrack pion2TT((*theB).build(iTrack2.bestTrack()));
		reco::TransientTrack pion3TT((*theB).build(iTrack3.bestTrack()));
	  // for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin();
	  //     iTrack1 != thePATTrackHandle->end(); ++iTrack1 )	{
	    
	  //   //cout<<"index of track "<<iTrack1<<endl;
	  //   if(iTrack1->charge()==0) continue;		  
	  //   if(iTrack1->pt()<0.8) continue;
	  //   if(!(iTrack1->trackHighPurity())) continue;
	    
	  //   for(View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1;
	  // 	iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) {

	      
	  //     if(iTrack1==iTrack2) continue;
	  //     if(iTrack2->charge()==0) continue;
	  //     if(iTrack2->pt()<1.0) continue;
	  //     if(!(iTrack2->trackHighPurity())) continue;
	      
	  //     if(iTrack1->charge() == iTrack2->charge()) continue;
	      
	      
	  //     //Now let's checks if our muons do not use the same tracks as we are using now
	  //     if ( IsTheSametk(*iTrack1,*iMuon1) || IsTheSametk(*iTrack1,*iMuon2) ) continue;
	  //     if ( IsTheSametk(*iTrack2,*iMuon1) || IsTheSametk(*iTrack2,*iMuon2) ) continue;
	      
	  //     for(View<pat::PackedCandidate>::const_iterator iTrack3 = thePATTrackHandle->begin();
	  // 	  iTrack3 != thePATTrackHandle->end(); ++iTrack3 ) 
	  // 	{
		  
	  // 	  if(iTrack3==iTrack2) continue;
	  // 	  if(iTrack3==iTrack1) continue;
	  // 	  if(iTrack3->charge()==0) continue;
	  // 	  if(iTrack3->pt()<1.0) continue;
	  // 	  if(!(iTrack3->trackHighPurity())) continue;
	  // 	  if ( IsTheSametk(*iTrack3,*iMuon1) || IsTheSametk(*iTrack3,*iMuon2) ) continue;		      
		  
		  
	  // 	  reco::TransientTrack pion1TT((*theB).build(iTrack1->pseudoTrack()));
	  // 	  reco::TransientTrack pion2TT((*theB).build(iTrack2->pseudoTrack()));
	  // 	  reco::TransientTrack pion3TT((*theB).build(iTrack3->pseudoTrack()));
		  

	  // //Now let's see if these two tracks make a vertex
		
		
		
		ParticleMass pion_mass = 0.13957018;
		ParticleMass kaon_mass = 0.493677;
		float kaon_sigma = 1.6e-5;
		ParticleMass Ds_mass = 1.9097614;
		float pion_sigma = pion_mass*1.e-6;		
		TLorentzVector pion14V,pion24V,pion34V; 
	       
		pion14V.SetXYZM(pion1TT.track().px(),pion1TT.track().py(),pion1TT.track().pz(),kaon_mass);
		pion24V.SetXYZM(pion2TT.track().px(),pion2TT.track().py(),pion2TT.track().pz(),kaon_mass);
		pion34V.SetXYZM(pion3TT.track().px(),pion3TT.track().py(),pion3TT.track().pz(),pion_mass);
		TLorentzVector Ds4V=pion14V+pion24V+pion34V;
		//histos[h_Dsmass]->Fill(Ds4V.M());
		//histos[h_Dspt]->Fill(Ds4V.Pt());
		if(Ds4V.M()<1.40 || Ds4V.M()>2.5) continue;
		if(Ds4V.Pt()<0)continue;
		
		//initial chi2 and ndf before kinematic fits.
		float chi = 0.;
		float ndf = 0.;
		vector<RefCountedKinematicParticle> pionParticles;
		// vector<RefCountedKinematicParticle> muonParticles;
		try {
		  pionParticles.push_back(pFactory.particle(pion1TT,kaon_mass,chi,ndf,kaon_sigma));
		  pionParticles.push_back(pFactory.particle(pion2TT,kaon_mass,chi,ndf,kaon_sigma));
		  pionParticles.push_back(pFactory.particle(pion3TT,pion_mass,chi,ndf,pion_sigma));
			
		}
		catch(...) {
		  std::cout<<" Exception caught ... continuing 3 "<<std::endl;
		  continue;
		}
		
		RefCountedKinematicTree DsVertexFitTree;
		try{
		  DsVertexFitTree = fitter.fit(pionParticles); 
		}
		catch(...) {
		  std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
		  continue;
		}
		if (!DsVertexFitTree->isValid()) 
		  {
		    //std::cout << "invalid vertex from the Ds vertex fit" << std::endl;
		    continue; 
		  }
		DsVertexFitTree->movePointerToTheTop();
		
		RefCountedKinematicParticle Ds_vFit_noMC = DsVertexFitTree->currentParticle();
		RefCountedKinematicVertex Ds_vFit_vertex_noMC = DsVertexFitTree->currentDecayVertex();
		
		if( Ds_vFit_vertex_noMC->chiSquared() < 0 )
		  { 
		    //std::cout << "negative chisq from ks fit" << endl;
		    continue;
		  }
		
		//some loose cuts go here
		//added for BCS 
		int v_BCS = -1;
	        double v_temp=999999999;
		
		if(Ds_vFit_vertex_noMC->chiSquared()>50) continue;
		if(Ds_vFit_noMC->currentState().mass()<1.90 || Ds_vFit_noMC->currentState().mass()>2.04) continue;
		double Dspt= sqrt(pow(Ds_vFit_noMC->currentState().globalMomentum().x(),2) +pow(Ds_vFit_noMC->currentState().globalMomentum().y(),2));
		//if(Dspt<0) continue;
		//added for BCS
		if(Dspt>v_temp)  { 
	        v_temp=Dspt;
		}

		double Ds_Prob_tmp  = TMath::Prob(Ds_vFit_vertex_noMC->chiSquared(),(int)Ds_vFit_vertex_noMC->degreesOfFreedom());
		//if(Ds_Prob_tmp<0.01)
		  //{
		  //  continue;
		 // }
		DsVertexFitTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle T1CandMC = DsVertexFitTree->currentParticle();
		
		DsVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle T2CandMC = DsVertexFitTree->currentParticle();
		
		DsVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle T3CandMC = DsVertexFitTree->currentParticle();
		
		DsVertexFitTree->movePointerToTheTop();
		RefCountedKinematicParticle ks0_vFit_withMC = DsVertexFitTree->currentParticle();
		
		//Now we are ready to combine!
		
		//cout<<"Dspt "<<Dspt<<endl;
		vector<RefCountedKinematicParticle> vFitMCParticles;
		vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		//vFitMCParticles.push_back(pFactory.particle(pion1TT,kaon_mass,chi,ndf,kaon_sigma));
		//vFitMCParticles.push_back(pFactory.particle(pion2TT,kaon_mass,chi,ndf,kaon_sigma));
		
		vFitMCParticles.push_back(ks0_vFit_withMC);
		
		KinematicConstrainedVertexFitter kcvFitter;
		RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles);
		if (!vertexFitTree->isValid()) {
		  //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		  continue;
		}
		
		vertexFitTree->movePointerToTheTop();		     
		
		RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
		if (!bDecayVertexMC->vertexIsValid()){
		  //std::cout << "B MC fit vertex is not valid" << endl;
		  continue;
		}
		
		if(bCandMC->currentState().mass()<5.7 || bCandMC->currentState().mass()>6.8) continue;
		
		if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 ) 
		  {
		    //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
		    continue;
		  }
		
		double B_Prob_tmp  = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		if(B_Prob_tmp<0.01)
		  {
		    continue;
		  }		     
		std::cout << "B mass "<<bCandMC->currentState().mass()<< std::endl;
		// get children from final B fit
		vertexFitTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
		vertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
		
		vertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle DsCandMC = vertexFitTree->currentParticle();
		
		KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
		KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
		KinematicParameters psiMupKP;
		KinematicParameters psiMumKP;
		
		if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
		if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
		if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
		if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 
		
		GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
				    mu1CandMC->currentState().globalMomentum().y(),
				    mu1CandMC->currentState().globalMomentum().z());
		
		GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
				    mu2CandMC->currentState().globalMomentum().y(),
				    mu2CandMC->currentState().globalMomentum().z());
		
		GlobalVector Dsp1vec(T1CandMC->currentState().globalMomentum().x(),
				     T1CandMC->currentState().globalMomentum().y(),
				     T1CandMC->currentState().globalMomentum().z());
		
		GlobalVector Dsp2vec(T2CandMC->currentState().globalMomentum().x(),
				     T2CandMC->currentState().globalMomentum().y(),
				     T2CandMC->currentState().globalMomentum().z());
		
		GlobalVector Dsp3vec(T3CandMC->currentState().globalMomentum().x(),
				     T3CandMC->currentState().globalMomentum().y(),
				     T3CandMC->currentState().globalMomentum().z());
		
		
		KinematicParameters DsPi1KP = T1CandMC->currentState().kinematicParameters();
		KinematicParameters DsPi2KP = T2CandMC->currentState().kinematicParameters();
		KinematicParameters DsPi3KP = T3CandMC->currentState().kinematicParameters();
		
		// fill candidate variables now
		
		if(nB==0){		    
		  nMu  = nMu_tmp;
		  //cout<< "*Number of Muons : " << nMu_tmp << endl;
		} // end nB==0		     
		
		B_mass->push_back(bCandMC->currentState().mass());
		B_px->push_back(bCandMC->currentState().globalMomentum().x());
		B_py->push_back(bCandMC->currentState().globalMomentum().y());
		B_pz->push_back(bCandMC->currentState().globalMomentum().z());
		
		B_Ds_mass->push_back( Ds_vFit_noMC->currentState().mass() );
		B_Ds_px->push_back( Ds_vFit_noMC->currentState().globalMomentum().x() );
		B_Ds_py->push_back( Ds_vFit_noMC->currentState().globalMomentum().y() );
		B_Ds_pz->push_back( Ds_vFit_noMC->currentState().globalMomentum().z() );
		
		B_mumu_mass->push_back( psi_vFit_noMC->currentState().mass() );
		B_mumu_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		B_mumu_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		B_mumu_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
		
		B_Ds_pt1->push_back(Dsp1vec.perp());
		B_Ds_px1->push_back(DsPi1KP.momentum().x());
		B_Ds_py1->push_back(DsPi1KP.momentum().y());
		B_Ds_pz1->push_back(DsPi1KP.momentum().z());
		B_Ds_px1_track->push_back(iTrack1.px());
		B_Ds_py1_track->push_back(iTrack1.py());
		B_Ds_pz1_track->push_back(iTrack1.pz());
		// B_Ds_px1_track->push_back(iTrack1->px());
		// B_Ds_py1_track->push_back(iTrack1->py());
		// B_Ds_pz1_track->push_back(iTrack1->pz());
		B_Ds_charge1->push_back(T1CandMC->currentState().particleCharge());
		
		B_Ds_pt2->push_back(Dsp2vec.perp());
		B_Ds_px2->push_back(DsPi2KP.momentum().x());
		B_Ds_py2->push_back(DsPi2KP.momentum().y());
		B_Ds_pz2->push_back(DsPi2KP.momentum().z());
		B_Ds_px2_track->push_back(iTrack2.px());
		B_Ds_py2_track->push_back(iTrack2.py());
		B_Ds_pz2_track->push_back(iTrack2.pz());
		// B_Ds_px2_track->push_back(iTrack2->px());
		// B_Ds_py2_track->push_back(iTrack2->py());
		// B_Ds_pz2_track->push_back(iTrack2->pz());
		B_Ds_charge2->push_back(T2CandMC->currentState().particleCharge());

		B_Ds_pt3->push_back(Dsp3vec.perp());
		B_Ds_px3->push_back(DsPi3KP.momentum().x());
		B_Ds_py3->push_back(DsPi3KP.momentum().y());
		B_Ds_pz3->push_back(DsPi3KP.momentum().z());
		B_Ds_px3_track->push_back(iTrack3.px());
		B_Ds_py3_track->push_back(iTrack3.py());
		B_Ds_pz3_track->push_back(iTrack3.pz());
		// B_Ds_px3_track->push_back(iTrack3->px());
		// B_Ds_py3_track->push_back(iTrack3->py());
		// B_Ds_pz3_track->push_back(iTrack3->pz());
		B_Ds_charge3->push_back(T3CandMC->currentState().particleCharge());
		
		B_mumu_pt1->push_back(Jp1vec.perp());
		B_mumu_px1->push_back(psiMu1KP.momentum().x());
		B_mumu_py1->push_back(psiMu1KP.momentum().y());
		B_mumu_pz1->push_back(psiMu1KP.momentum().z());
		B_mumu_charge1->push_back(mu1CandMC->currentState().particleCharge());
		
		B_mumu_pt2->push_back(Jp2vec.perp());
		B_mumu_px2->push_back(psiMu2KP.momentum().x());
		B_mumu_py2->push_back(psiMu2KP.momentum().y());
		B_mumu_pz2->push_back(psiMu2KP.momentum().z());
		B_mumu_charge2->push_back(mu2CandMC->currentState().particleCharge());
		
		B_Ds_chi2->push_back(Ds_vFit_vertex_noMC->chiSquared());
		B_mumu_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
		B_chi2->push_back(bDecayVertexMC->chiSquared());
		B_chi2dof->push_back(bDecayVertexMC->chiSquared()/bDecayVertexMC->degreesOfFreedom());
		
		//double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		//double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		//double Ds_Prob_tmp  = TMath::Prob(Ds_vFit_vertex_noMC->chiSquared(),(int)Ds_vFit_vertex_noMC->degreesOfFreedom());
		B_Prob    ->push_back(B_Prob_tmp);
		B_mumu_Prob  ->push_back(J_Prob_tmp);
		B_Ds_Prob ->push_back(Ds_Prob_tmp);
		
		      // ************
		bDecayVtxX->push_back((*bDecayVertexMC).position().x());
		bDecayVtxY->push_back((*bDecayVertexMC).position().y());
		bDecayVtxZ->push_back((*bDecayVertexMC).position().z());
		bDecayVtxXE->push_back(bDecayVertexMC->error().cxx());
		bDecayVtxYE->push_back(bDecayVertexMC->error().cyy());
		bDecayVtxZE->push_back(bDecayVertexMC->error().czz());
		bDecayVtxXYE->push_back(bDecayVertexMC->error().cyx());
		bDecayVtxXZE->push_back(bDecayVertexMC->error().czx());
		bDecayVtxYZE->push_back(bDecayVertexMC->error().czy());
		
		VDecayVtxX->push_back( Ds_vFit_vertex_noMC->position().x() );
		VDecayVtxY->push_back( Ds_vFit_vertex_noMC->position().y() );
		VDecayVtxZ->push_back( Ds_vFit_vertex_noMC->position().z() );
		VDecayVtxXE->push_back( Ds_vFit_vertex_noMC->error().cxx() );
		VDecayVtxYE->push_back( Ds_vFit_vertex_noMC->error().cyy() );
		VDecayVtxZE->push_back( Ds_vFit_vertex_noMC->error().czz() );
		VDecayVtxXYE->push_back( Ds_vFit_vertex_noMC->error().cyx() );
		VDecayVtxXZE->push_back( Ds_vFit_vertex_noMC->error().czx() );
		VDecayVtxYZE->push_back( Ds_vFit_vertex_noMC->error().czy() );
		
		// ********************* muon-trigger-machint**************** 
		
		const pat::Muon* muon1 = &(*iMuon1);
		const pat::Muon* muon2 = &(*iMuon2);
		
		int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;
		
		if (muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) tri_Dim25_tmp = 1;
		if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
		if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*")!=nullptr) tri_JpsiTkTk_tmp = 1;
		
		std::cout<<"Do the trigger matching "<<tri_JpsiTkTk_tmp<<std::endl;
		tri_Dim25->push_back( tri_Dim25_tmp );	       
		tri_JpsiTk->push_back( tri_JpsiTk_tmp );
		tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );
		// Match the track with the trigger object
		float dr0_t = 99999.;
		float dpt0_t = 99999.;
		for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		  float dp = pion1TT.track().phi() - obj_phi[ii];
		  float de = pion1TT.track().eta() - obj_eta[ii];      
		  if (dp>float(M_PI)) dp-=float(2*M_PI);  
		  float dr = std::sqrt(de*de + dp*dp);
		  
		  if (dr < dr0_t) dr0_t = dr;
		  float dpt = pion1TT.track().pt() - obj_pt[ii];
		  if (abs(dpt) < dpt0_t) dpt0_t = abs(dpt);
		}
		
		float dr1_t = 99999.;
		float dpt1_t = 99999.;
		for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		  float dp = pion2TT.track().phi() - obj_phi[ii];
		  float de = pion2TT.track().eta() - obj_eta[ii];    
		  if (dp>float(M_PI)) dp-=float(2*M_PI);  
		  float dr = std::sqrt(de*de + dp*dp);		     
		  if (dr < dr1_t) dr1_t = dr;
		  float dpt = pion2TT.track().pt() - obj_pt[ii];
		  if (abs(dpt) < dpt1_t) dpt1_t =abs(dpt);
		  
		}
		
		dr0->push_back(dr0_t);
		dr1->push_back(dr1_t);
		
		dpt0->push_back(dpt0_t);
		dpt1->push_back(dpt1_t);
		// ************
		
		mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
		mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
		mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
		mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
		mu1PF->push_back(iMuon1->isPFMuon());
		mu2PF->push_back(iMuon2->isPFMuon());
		mu1loose->push_back(muon::isLooseMuon(*iMuon1));
		mu2loose->push_back(muon::isLooseMuon(*iMuon2));
		
		mumC2->push_back( glbTrackM->normalizedChi2() );
		mumNHits->push_back( glbTrackM->numberOfValidHits() );
		mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
		mupC2->push_back( glbTrackP->normalizedChi2() );
		mupNHits->push_back( glbTrackP->numberOfValidHits() );
		mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
		mumdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
		mupdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
		mumdz->push_back(glbTrackM->dz(bestVtx.position()) );
		mupdz->push_back(glbTrackP->dz(bestVtx.position()) );
		muon_dca->push_back(dca);

		// pi1dxy->push_back(iTrack1->dxy());
		// pi2dxy->push_back(iTrack2->dxy());
		// pi1dz->push_back(iTrack1->dz());
		// pi2dz->push_back(iTrack2->dz());
		
		// pi1dxy_e->push_back(iTrack1->dxyError());
		// pi2dxy_e->push_back(iTrack2->dxyError());
		// pi1dz_e->push_back(iTrack1->dzError());
		// pi2dz_e->push_back(iTrack2->dzError());

		pi1dxy->push_back(iTrack1.dxy());
		pi2dxy->push_back(iTrack2.dxy());
		pi1dz->push_back(iTrack1.dz());
		pi2dz->push_back(iTrack2.dz());
		
		pi1dxy_e->push_back(iTrack1.dxyError());
		pi2dxy_e->push_back(iTrack2.dxyError());
		pi1dz_e->push_back(iTrack1.dzError());
		pi2dz_e->push_back(iTrack2.dzError());
			      
		nB++;	
		
		int pvIndex = -1;
		saveIP(vertexFitTree,*pvHandle_.product(),pvIndex);
		SaveIso(vertexFitTree, fMagneticField, theB,pvIndex,
			thePATTrackHandle, thePATMuonHandle,
			muon1TT, muon2TT,pion1TT, pion2TT, pion3TT);
		
		
		if(isMC_)saveTruthMatch(iEvent);
		
		
		pionParticles.clear();
		muonParticles.clear();
		vFitMCParticles.clear();
		
	      }
	    }
	    
	  }
	  trklist.clear();
	  //cout<<"number of event left after pt>1.0 "<<trklist.size()<<" track "<<thePATTrackHandle->size()<<endl;
	}
    }
  
  
   //fill the tree and clear the vectors
  if (nB > 0 ) 
    {
      std::cout << "filling tree " << nB<<std::endl;
      tree_->Fill();
    }
  // *********
  
  nB = 0; nMu = 0;

  B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();
  B_Ds_mass->clear(); B_Ds_px->clear(); B_Ds_py->clear(); B_Ds_pz->clear();
  
  B_mumu_mass->clear();  B_mumu_px->clear();  B_mumu_py->clear();  B_mumu_pz->clear();
  
  B_Ds_pt1->clear(); B_Ds_px1->clear(); B_Ds_py1->clear(); B_Ds_pz1->clear(); B_Ds_charge1->clear(); 
  B_Ds_pt2->clear(); B_Ds_px2->clear(); B_Ds_py2->clear(); B_Ds_pz2->clear(); B_Ds_charge2->clear(); 
  B_Ds_pt3->clear(); B_Ds_px3->clear(); B_Ds_py3->clear(); B_Ds_pz3->clear(); B_Ds_charge3->clear(); 
  
  B_Ds_px1_track->clear(); B_Ds_py1_track->clear(); B_Ds_pz1_track->clear(); 
  B_Ds_px2_track->clear(); B_Ds_py2_track->clear(); B_Ds_pz2_track->clear(); 
  B_Ds_px3_track->clear(); B_Ds_py3_track->clear(); B_Ds_pz3_track->clear(); 
  
  B_mumu_pt1->clear();  B_mumu_px1->clear();  B_mumu_py1->clear();  B_mumu_pz1->clear(), B_mumu_charge1->clear();
  B_mumu_pt2->clear();  B_mumu_px2->clear();  B_mumu_py2->clear();  B_mumu_pz2->clear(), B_mumu_charge2->clear();
  
  B_Ds_chi2->clear(); B_mumu_chi2->clear(); B_chi2->clear(); B_chi2dof->clear();
  B_Prob->clear(); B_mumu_Prob->clear(); B_Ds_Prob->clear();
  
  
  bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear(); 
  bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear(); 
  bDecayVtxXYE->clear(); bDecayVtxXZE->clear(); bDecayVtxYZE->clear();  
  
  VDecayVtxX->clear(); VDecayVtxY->clear(); VDecayVtxZ->clear();
  VDecayVtxXE->clear(); VDecayVtxYE->clear(); VDecayVtxZE->clear();
  VDecayVtxXYE->clear(); VDecayVtxXZE->clear(); VDecayVtxYZE->clear();
  
  // *********
  
  nVtx = 0;
  priVtxX = 0; priVtxY = 0; priVtxZ = 0; 
  priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0;
  priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;
  
  pi1dxy->clear(); pi2dxy->clear(); pi1dz->clear(); pi2dz->clear();
  pi1dxy_e->clear(); pi2dxy_e->clear(); pi1dz_e->clear(); pi2dz_e->clear();
  
  mumC2->clear();
  mumNHits->clear(); mumNPHits->clear();
  mupC2->clear();
  mupNHits->clear(); mupNPHits->clear();
  mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();
  
  tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear();
  dr0->clear(); dr1->clear(); dpt0->clear(); dpt1->clear();
  mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
  mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 
  
  // priRfVtxX->clear(); priRfVtxY->clear(); priRfVtxZ->clear(); priRfVtxXE->clear(); priRfVtxYE->clear(); 
  // priRfVtxZE->clear(); priRfVtxXYE->clear(); priRfVtxXZE->clear(); priRfVtxYZE->clear(); priRfVtxCL->clear(); 
  // priRfNTrkDif->clear(); 
  
  pVtxIPX->clear();  pVtxIPY->clear();  pVtxIPZ->clear();
  pVtxIPXE->clear();  pVtxIPYE->clear();  pVtxIPZE->clear();  pVtxIPCL->clear();
  pVtxIPXYE->clear();  pVtxIPXZE->clear();  pVtxIPYZE->clear();
  
  B_l3d->clear();  B_l3dE->clear();  B_lxy->clear(); B_lxyE->clear();
  B_cosalpha->clear();   B_cosalphaxy->clear(); alpha->clear();  B_treco->clear();   B_trecoe->clear();  B_trecoxy->clear(); B_trecoxye->clear();
  B_pvip->clear(); B_pviperr->clear(); B_pvips->clear(); B_pvlzip->clear(); B_pvlziperr->clear(); B_pvlzips->clear();
  B_pv2ip->clear(); B_pv2iperr->clear(); B_pv2ips->clear(); B_pv2lzip->clear(); B_pv2lziperr->clear(); B_pv2lzips->clear();
  
  B_l3d_pv2->clear();  B_l3dE_pv2->clear();
  B_iso->clear(); B_mum_iso->clear(); B_mup_iso->clear(); B_pi1_iso->clear();B_pi2_iso->clear();
  
  istruemum->clear(); istruemup->clear(); istruekp->clear(); istruekm->clear(); istruebs->clear();
  bunchXingMC->clear(); numInteractionsMC->clear(); trueNumInteractionsMC->clear();
  
  
  
}
bool BcToDsMuMu::IsTheSametk(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if(deltaR(tk.eta(),tk.phi(),mu.eta(),mu.phi()) < 0.0001 )return true;
  //if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool BcToDsMuMu::IsTheSamePFtk(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if(deltaR(tk.eta(),tk.phi(),mu.eta(),mu.phi()) < 0.0001 )return true;
  //if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool BcToDsMuMu::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool BcToDsMuMu::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
  if (ancestor == particle ) return true;
  for (size_t i=0; i< particle->numberOfMothers(); i++) {
    if (isAncestor(ancestor,particle->mother(i))) return true;
  }
  return false;
}

double BcToDsMuMu::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
  TVector3 pv_dv = decay_vtx - production_vtx;
  TVector3 b_p3  = b_p4.Vect();
  pv_dv.SetZ(0.);
  b_p3.SetZ(0.);
  Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
  return lxy*b_p4.M()/b_p3.Mag();
}

double
BcToDsMuMu::calEta (double Px, double Py, double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double
BcToDsMuMu::calPhi (double Px, double Py, double Pz)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

double
BcToDsMuMu::calEtaPhiDistance (double Px1, double Py1, double Pz1,
				double Px2, double Py2, double Pz2)
{
  double phi1 = calPhi (Px1,Py1,Pz1);
  double eta1 = calEta (Px1,Py1,Pz1);
  double phi2 = calPhi (Px2,Py2,Pz2);
  double eta2 = calEta (Px2,Py2,Pz2);
  return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}
void
BcToDsMuMu::saveTruthMatch(const edm::Event& iEvent){
  double deltaEtaPhi;

  for (vector<int>::size_type i = 0; i < B_mass->size(); i++) {//{{{
   
    //-----------------------
    // truth match with mu-
    //-----------------------
    double TruthMatchMuonMaxR_ = 0.004;
    double TruthMatchKaonMaxR_ = 0.3;
    deltaEtaPhi = calEtaPhiDistance(gen_muon1_p4.Px(), gen_muon1_p4.Py(), gen_muon1_p4.Pz(),
     				    B_mumu_px1->at(i), B_mumu_py1->at(i), B_mumu_pz1->at(i));
    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemum->push_back(true);
    } else {
      istruemum->push_back(false);
    }

    //-----------------------
    // truth match with mu+
    //-----------------------

    deltaEtaPhi = calEtaPhiDistance(gen_muon2_p4.Px(), gen_muon2_p4.Py(), gen_muon2_p4.Pz(),
     				    B_mumu_px2->at(i), B_mumu_py2->at(i), B_mumu_pz2->at(i));

    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemup->push_back(true);
    }
    else {
      istruemup->push_back(false);
    }

    //---------------------------------
    // truth match with pion+ track   
    //---------------------------------                          
    deltaEtaPhi = calEtaPhiDistance(gen_pion1_p4.Px(), gen_pion1_p4.Py(), gen_pion1_p4.Pz(),
                                    B_Ds_px1->at(i), B_Ds_py1->at(i), B_Ds_pz1->at(i));
    if (deltaEtaPhi < TruthMatchKaonMaxR_){
      istruekp->push_back(true);
    } else {
      istruekp->push_back(false);
    }

    //---------------------------------                                                                                                                           
    // truth match with pion- track                                                                                                                                 
    //---------------------------------                                                                                                                           
    deltaEtaPhi = calEtaPhiDistance(gen_pion2_p4.Px(), gen_pion2_p4.Py(), gen_pion2_p4.Pz(),
                                    B_Ds_px2->at(i), B_Ds_py2->at(i), B_Ds_pz2->at(i));
    if (deltaEtaPhi < TruthMatchKaonMaxR_){
      istruekm->push_back(true);
    } else {
      istruekm->push_back(false);
    }


    //---------------------------------------
    // truth match with Bs or Bs bar 
    //---------------------------------------                                                                                                
    if ( istruemum->back() && istruemup->back() && istruekm->back() && istruekp->back() ) {
      istruebs->push_back(true);
    } else {
      istruebs->push_back(false);
    }



    }//}}}

}
void BcToDsMuMu::savePUinMC(const edm::Event& iEvent){
  // #################################
  // # Save pileup information in MC #
  // #################################
  // edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
  // iEvent.getByToken(puToken_, PupInfo);
  // //edm::LumiReWeighting lumi_weights;
  // //lumi_weights      = edm::LumiReWeighting("PileupMC_2018.root", "DataPileupHistogram2018_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  // //lumi_weights      = edm::LumiReWeighting("PileupMC_2016.root", "DataPileupHistogram2016_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  // //lumi_weights      = edm::LumiReWeighting("PileupMC_2017.root", "DataPileupHistogram2017_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  // //float tnpv = -1;       // True number of primary vertices
  // //float wpu = 1;         // Pile-up re-weight factor
  
  // for (std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); PVI++)
  //   {
  //     bunchXingMC->push_back(PVI->getBunchCrossing());
  //     numInteractionsMC->push_back(PVI->getPU_NumInteractions());
  //     trueNumInteractionsMC->push_back(PVI->getTrueNumInteractions());
  //     //int bx = PVI->getBunchCrossing();
  //     //  if (bx == 0)tnpv = PVI->getTrueNumInteractions();        
  //   }

  //wpu = lumi_weights.weight(tnpv);
  //  fpuw8 = wpu ;
}


cov33_t BcToDsMuMu::GlobalError2SMatrix_33(GlobalError m_in) 
{
  cov33_t m_out;
  for (int i=0; i<3; i++) {
    for (int j=i; j<3; j++)  {
      m_out(i,j) = m_in.matrix()(i,j);
    }
  }
  return m_out;
}
  
cov99_t BcToDsMuMu::makeCovarianceMatrix(const cov33_t cov_vtx1,
			     const cov77_t cov_vtx2) 
{
  cov99_t cov;
  cov.Place_at(cov_vtx1,0,0);
  cov.Place_at(cov_vtx2.Sub<cov66_t>(0,0),3,3);
  return cov;
}
jac9_t BcToDsMuMu::makeJacobianVector3d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
 				     ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
				     const GlobalPoint &vtx2, const TVector3 &tv3momentum) 
{
  return makeJacobianVector3d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}

jac9_t BcToDsMuMu::makeJacobianVector3d(const AlgebraicVector3 &vtx1, 
				     const AlgebraicVector3 &vtx2, 
				     const AlgebraicVector3 &momentum) 
{
  jac9_t jac;
  const AlgebraicVector3 dist = vtx2 - vtx1;
  const double factor2 = 1. / ROOT::Math::Mag2(momentum);
  const double lifetime = ROOT::Math::Dot(dist, momentum) * factor2;
  jac.Place_at(-momentum*factor2,0);
  jac.Place_at( momentum*factor2,3);
  jac.Place_at( factor2*(dist-2*lifetime*momentum*factor2),6);
  return jac;
}
jac9_t BcToDsMuMu::makeJacobianVector2d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2,
				     const AlgebraicVector3 &momentum) {
  jac9_t jac;
  const double momentumMag = ROOT::Math::Mag(momentum);
  const AlgebraicVector3 dist = vtx2 - vtx1;
  const double distMag = ROOT::Math::Mag(dist);
  const double factorPositionComponent = 1./(distMag*momentumMag);
  const double factorMomentumComponent = 1./pow(momentumMag,3);
  jac(0)=-dist(0)*factorPositionComponent;
  jac(1)=-dist(1)*factorPositionComponent;
  jac(3)= dist(0)*factorPositionComponent;
  jac(4)= dist(1)*factorPositionComponent;
  jac(6)= momentum(0)*factorMomentumComponent;
  jac(7)= momentum(1)*factorMomentumComponent;
  return jac;
}

jac9_t BcToDsMuMu::makeJacobianVector2d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
			    ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			    const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
  return makeJacobianVector2d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}


void BcToDsMuMu::saveIP(const RefCountedKinematicTree& vertexFitTree,
			const reco::VertexCollection& vertices, int & pvIndex){

  vertexFitTree->movePointerToTheTop();

  RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
  RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
  
  std::cout<<" inside the save IP "<<std::endl;
  auto candTransientTrack = bCandMC->refittedTransientTrack();
  // find the first primary vertex
  const reco::Vertex* bestVertex_t(0);
  int bestVertexIndex(-1);
  double minDistance(999.);
  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);
    auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
    if (impactParameter3D.first and impactParameter3D.second.value() < minDistance){
      minDistance = impactParameter3D.second.value();
      bestVertex_t = &vertex;
      bestVertexIndex = i;
    }
  }
  pvIndex = bestVertexIndex;
  pVtxIPX->push_back( bestVertex_t->x());
  pVtxIPY->push_back( bestVertex_t->y());
  pVtxIPZ->push_back( bestVertex_t->z());
  pVtxIPXE->push_back( bestVertex_t->covariance(0, 0) );
  pVtxIPYE->push_back( bestVertex_t->covariance(1, 1) );
  pVtxIPZE->push_back( bestVertex_t->covariance(2, 2) );
  pVtxIPXYE->push_back(bestVertex_t->covariance(0, 1) );
  pVtxIPXZE->push_back(bestVertex_t->covariance(0, 2) );  
  pVtxIPYZE->push_back( bestVertex_t->covariance(1, 2) );
  pVtxIPCL->push_back(  ChiSquaredProbability((double)(bestVertex_t->chi2()),(double)(bestVertex_t->ndof())) );      

  // find second best vertex
  const reco::Vertex* bestVertex2(0);
  //int bestVertexIndex2(-1);
  double minDistance2(999.);
  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);    
    auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
    if (impactParameter3D.first and impactParameter3D.second.value() < minDistance2 and impactParameter3D.second.value() > minDistance){
      minDistance2 = impactParameter3D.second.value();
      bestVertex2 = &vertex;
      //bestVertexIndex2 = i;
    }
  }
  //std::cout<<" inside the save IP second vertex "<<std::endl;
  //  if (! bestVertex) continue;
  auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex_t);
  auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex_t);
  //pv = bestVertex;
  //pvIndex = bestVertexIndex;
  double longitudinalImpactParameter(0.0), longitudinalImpactParameterErr(0.0);
  double distaceOfClosestApproach(0.0), distaceOfClosestApproachErr(0.0) ;
  if (impactParameterZ.first) {
    longitudinalImpactParameter    = impactParameterZ.second.value();
    longitudinalImpactParameterErr = impactParameterZ.second.error();
  }
  if(impactParameter3D.first) {
    distaceOfClosestApproach       = impactParameter3D.second.value();
    distaceOfClosestApproachErr    = impactParameter3D.second.error();
  }
  //  std::cout<<" inside dca "<<distaceOfClosestApproach<<"\t "<<distaceOfClosestApproachErr<<"\t"<<distaceOfClosestApproachErr/distaceOfClosestApproach<<std::endl;
  //std::cout<<" inside long "<<longitudinalImpactParameter<<"\t"<<longitudinalImpactParameterErr<<"\t"<<longitudinalImpactParameterErr/longitudinalImpactParameter<<std::endl;

  B_pvip->push_back(distaceOfClosestApproach);
  B_pviperr->push_back(distaceOfClosestApproachErr);
  B_pvips->push_back(distaceOfClosestApproachErr/distaceOfClosestApproach);
  B_pvlzip->push_back(longitudinalImpactParameter);
  B_pvlziperr->push_back(longitudinalImpactParameterErr);
  B_pvlzips->push_back(longitudinalImpactParameterErr/longitudinalImpactParameter);
  
  std::cout<<" inside pvip "<<longitudinalImpactParameter<<"\t"<<longitudinalImpactParameterErr<<"\t"<<longitudinalImpactParameterErr/longitudinalImpactParameter<<std::endl;
  VertexDistance3D distance3D;
  auto dist = distance3D.distance(*bestVertex_t, bDecayVertexMC->vertexState() );
  double decayLength(-1.), decayLengthErr(0);
  decayLength    = dist.value();
  decayLengthErr = dist.error();
  
  VertexDistanceXY distanceXY;
  auto distXY = distanceXY.distance(*bestVertex_t, bDecayVertexMC->vertexState() );

  B_l3d ->push_back(decayLength);
  B_l3dE ->push_back(decayLengthErr);
  B_lxy ->push_back(distXY.value());
  B_lxyE ->push_back(distXY.error());
  
  if (bestVertex2){
    double longitudinalImpactParameter2(0.0), longitudinalImpactParameter2Err(0.0);
    double distaceOfClosestApproach2(0.0), distaceOfClosestApproach2Err(0.0) ;
    auto impactParameter3D2 = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex2);
    auto impactParameterZ2  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex2);
    //pv2 = bestVertex2;
    //pv2Index = bestVertexIndex2;
    if (impactParameterZ2.first) {
      longitudinalImpactParameter2    = impactParameterZ2.second.value();
      longitudinalImpactParameter2Err = impactParameterZ2.second.error();
    }
    if (impactParameter3D2.first) {
      distaceOfClosestApproach2       = impactParameter3D2.second.value();
      distaceOfClosestApproach2Err    = impactParameter3D2.second.error();
    }
    B_pv2ip->push_back(distaceOfClosestApproach2);
    B_pv2iperr->push_back(distaceOfClosestApproach2Err);
    B_pv2ips->push_back(distaceOfClosestApproach2Err/distaceOfClosestApproach2);
    B_pv2lzip->push_back(longitudinalImpactParameter2);
    B_pv2lziperr->push_back(longitudinalImpactParameter2Err);
    B_pv2lzips->push_back(longitudinalImpactParameter2Err/longitudinalImpactParameter2);
    
    // compute decay length
    VertexDistance3D distance3D;
    auto dist = distance3D.distance(*bestVertex2, bDecayVertexMC->vertexState() );
    B_l3d_pv2 ->push_back(dist.value());
    B_l3dE_pv2 ->push_back(dist.error());
  }
  TVector3 plab(bCandMC->currentState().globalMomentum().x(),
		bCandMC->currentState().globalMomentum().y(),
		bCandMC->currentState().globalMomentum().z());
  TVector3 p1(bestVertex_t->x(), bestVertex_t->y(), bestVertex_t->z());
  TVector3 p2(bDecayVertexMC->vertexState().position().x(), 
	      bDecayVertexMC->vertexState().position().y(), 
	      bDecayVertexMC->vertexState().position().z());
  TVector3 pDiff = p2-p1;
  TVector3 pDiffXY = TVector3(pDiff.X(), pDiff.Y(), 0.);
  TVector3 ptrans  = TVector3(plab.X(), plab.Y(), 0.);
  double cosAlpha(-999.),cosAlphaXY(-999.),decayTime(-999.),decayTimeError(-999.), decayTimeXY(-999.),decayTimeXYError(-999.);
  cosAlpha  = plab.Dot(pDiff) / (plab.Mag() * pDiff.Mag());
  cosAlphaXY  = ptrans.Dot(pDiffXY) / (ptrans.Mag() * pDiffXY.Mag());
  
  B_cosalpha -> push_back(cosAlpha);
  B_cosalphaxy -> push_back(cosAlphaXY);
  alpha->push_back(TMath::ACos(cosAlpha));
  // compute decayTime information

  const double massOverC =  bCandMC->currentState().mass()/TMath::Ccgs();

  // get covariance matrix for error propagation in decayTime calculation
  auto vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(bestVertex_t->error()),
					     bCandMC->currentState().kinematicParametersError().matrix());
  auto vtxDistanceJac3d = makeJacobianVector3d(bestVertex_t->position(), bDecayVertexMC->vertexState().position(), plab);
  auto vtxDistanceJac2d = makeJacobianVector2d(bestVertex_t->position(), bDecayVertexMC->vertexState().position(), plab);

  decayTime = dist.value() / plab.Mag() * cosAlpha * massOverC;
  decayTimeError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac3d)) * massOverC;

  decayTimeXY = distXY.value() / plab.Perp() * cosAlphaXY * massOverC;
  decayTimeXYError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac2d)) * massOverC;

  B_treco -> push_back(decayTime);
  B_trecoe -> push_back(decayTimeError);
  B_trecoxy -> push_back(decayTimeXY);
  B_trecoxye -> push_back(decayTimeXYError);
  
}
// Isolation
void BcToDsMuMu::SaveIso(const RefCountedKinematicTree& vertexFitTree, const MagneticField *fMagneticField, 			 
			 edm::ESHandle<TransientTrackBuilder> theB,
			 unsigned int pvIndex,
			 edm::Handle< edm::View<pat::PackedCandidate> > thePATTrackHandle,
			 //edm::Handle<std::vector<pat::GenericParticle> > thePATTrackHandle,
			 edm::Handle< edm::View<pat::Muon> > thePATMuonHandle,
			 const reco::TransientTrack muon1TT, 
			 const reco::TransientTrack muon2TT,
			 const reco::TransientTrack pion1TT,
			 const reco::TransientTrack pion2TT,
			 const reco::TransientTrack pion3TT		      
			 ){

  vertexFitTree->movePointerToTheTop(); // Bs --> Jpsi Kss
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex vertex = vertexFitTree->currentDecayVertex();
  
  
  ClosestApproachInRPhi ClosestApp;
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  
  double  bpx = b_KP->currentState().globalMomentum().x();
  double  bpy = b_KP->currentState().globalMomentum().y();
  double  bpz = b_KP->currentState().globalMomentum().z();
  double masstmp = b_KP->currentState().mass();
  TLorentzVector tmp_bmeson_lv;
  tmp_bmeson_lv.SetXYZM(bpx,bpy, bpz, masstmp);
  
  double sumppt(0),summpt(0),sumtrkmpt(0),sumtrkppt(0);
  double sumBpt(0);
  //std::cout<<" Inside the Isolation folder "<<std::endl;
  //std::vector<std::pair<int,std::pair<float,float> > > fNstTracks;
  for(edm::View<pat::PackedCandidate>::const_iterator iTrack = thePATTrackHandle->begin();                                                                 
      iTrack != thePATTrackHandle->end(); ++iTrack )   {
    
    if(iTrack->charge()==0) continue;                                                                                                           
    if(!iTrack->trackHighPurity()) continue;                                                                                                  
    if(!iTrack->hasTrackDetails()) continue;
    bool skip_this_track = false; 
    for(edm::View<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin(); iMuon != thePATMuonHandle->end(); ++iMuon)
      {
	if ( IsTheSametk(*iTrack,*iMuon) ) skip_this_track = true;
      }	
    if(skip_this_track)continue;
    if (iTrack->vertexRef().key()!=pvIndex) continue;//{ std::cout<<"passing our pvindex "<<std::endl;}// continue;
    
    if(skip_this_track)continue;    
    reco::TransientTrack TrackIsoTT((*theB).build(iTrack->pseudoTrack()));
    if( !(TrackIsoTT.isValid()) )continue;    
    if(iTrack->pt()<0.8) continue;
    bool skipping_track= false;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), muon1TT.track().eta(),  muon1TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==muon1TT.track().charge())) skipping_track = true;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), muon2TT.track().eta(),  muon2TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==muon2TT.track().charge())) skipping_track = true;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), pion1TT.track().eta(),  pion1TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==pion1TT.track().charge())) skipping_track = true;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), pion2TT.track().eta(),  pion2TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==pion2TT.track().charge())) skipping_track = true;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), pion3TT.track().eta(),  pion3TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==pion3TT.track().charge())) skipping_track = true;
    if(skipping_track) continue;
    
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), muon1TT.initialFreeState());
    if (ClosestApp.status() != false)
      {
	if ( ClosestApp.distance() < 0.1 ) {
	  double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), muon1TT.track().eta(),  muon1TT.track().phi());
	  if(deltaR_tmp<0.5) sumppt += iTrack->pt();
	}                 
      }
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), muon2TT.initialFreeState());
    if (ClosestApp.status() != false)
      {
	if ( ClosestApp.distance() < 0.1 ) {
	  double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), muon2TT.track().eta(),  muon2TT.track().phi());
	  if(deltaR_tmp<0.5) summpt += iTrack->pt();
	}
      }
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), pion1TT.initialFreeState());
    if (ClosestApp.status() != false)
      {
	if ( ClosestApp.distance() < 0.1 ) {
	  double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), pion1TT.track().eta(),  pion1TT.track().phi());
	    if(deltaR_tmp<0.5) sumtrkmpt += iTrack->pt();
	}                 
	}
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), pion2TT.initialFreeState());
      if (ClosestApp.status() != false)
	{
	  if ( ClosestApp.distance() < 0.1 ) {
	    double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), pion2TT.track().eta(),  pion2TT.track().phi());
	    if(deltaR_tmp<0.5) sumtrkppt += iTrack->pt();
	  }
	}  
      
      VertexDistance3D distance3D;      
      const GlobalPoint BVP = GlobalPoint( vertex->position() );
      if(vertex->vertexIsValid()){       
	TrajectoryStateOnSurface  tsos = extrapolator.extrapolate(TrackIsoTT.initialFreeState(), BVP);
	if (tsos.isValid()) {      			   
	  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex->vertexState());
	  double svDoca = doca.value();	    
	  //docatrks.push_back(svDoca);
	  // check to ensure the goodness of the track
	  if( svDoca<0.05){
	    double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), tmp_bmeson_lv.Eta(), tmp_bmeson_lv.Phi());			   
	    if(deltaR_tmp<0.7) sumBpt += iTrack->pt();			     
	  }
	  
	}
      }
    }

  float BIso_ =tmp_bmeson_lv.Pt()/(sumBpt +tmp_bmeson_lv.Pt());
  float mupIso_ = -99., mumIso_ = -99.0;
  if(muon1TT.track().charge() == 1){
    mupIso_ = muon1TT.track().pt()/(sumppt +muon1TT.track().pt() );
    mumIso_ = muon2TT.track().pt()/(summpt +muon2TT.track().pt() );
  }else {
    mupIso_ = muon2TT.track().pt()/(sumppt +muon2TT.track().pt() );
    mumIso_ = muon1TT.track().pt()/(summpt +muon1TT.track().pt() );
  }
  float trkpIso_ = pion1TT.track().pt()/(sumtrkppt +pion1TT.track().pt() );
  float trkmIso_ = pion2TT.track().pt()/(sumtrkmpt +pion2TT.track().pt() );
  B_iso->push_back(BIso_); 
  B_mum_iso->push_back(mumIso_); 
  B_mup_iso->push_back(mupIso_); 
  B_pi1_iso->push_back(trkpIso_);
  B_pi2_iso->push_back(trkmIso_);
  
}
// ------------ method called once each job just before starting event loop  ------------

void 
BcToDsMuMu::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  // for(int i=0; i<kHistNameSize; i++) {
  //   histos[i] = new TH1F(hist_args[i].name, hist_args[i].title,
  //                        hist_args[i].n_bins,
  //                        hist_args[i].x_min, hist_args[i].x_max);

  // }

  tree_ = fs->make<TTree>("ntuple","Bs->J/psi Ds ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_Ds_mass", &B_Ds_mass);
  tree_->Branch("B_Ds_px", &B_Ds_px);
  tree_->Branch("B_Ds_py", &B_Ds_py);
  tree_->Branch("B_Ds_pz", &B_Ds_pz);
 
  tree_->Branch("B_mumu_mass", &B_mumu_mass);
  tree_->Branch("B_mumu_px", &B_mumu_px);
  tree_->Branch("B_mumu_py", &B_mumu_py);
  tree_->Branch("B_mumu_pz", &B_mumu_pz);

  tree_->Branch("B_Ds_pt1", &B_Ds_pt1);
  tree_->Branch("B_Ds_px1", &B_Ds_px1);
  tree_->Branch("B_Ds_py1", &B_Ds_py1);
  tree_->Branch("B_Ds_pz1", &B_Ds_pz1);
  tree_->Branch("B_Ds_px1_track", &B_Ds_px1_track);
  tree_->Branch("B_Ds_py1_track", &B_Ds_py1_track);
  tree_->Branch("B_Ds_pz1_track", &B_Ds_pz1_track);
  tree_->Branch("B_Ds_charge1", &B_Ds_charge1); 
 
  tree_->Branch("B_Ds_pt2", &B_Ds_pt2);
  tree_->Branch("B_Ds_px2", &B_Ds_px2);
  tree_->Branch("B_Ds_py2", &B_Ds_py2);
  tree_->Branch("B_Ds_pz2", &B_Ds_pz2);
  tree_->Branch("B_Ds_px2_track", &B_Ds_px2_track);
  tree_->Branch("B_Ds_py2_track", &B_Ds_py2_track);
  tree_->Branch("B_Ds_pz2_track", &B_Ds_pz2_track);
  tree_->Branch("B_Ds_charge2", &B_Ds_charge2);

  tree_->Branch("B_Ds_pt3", &B_Ds_pt3);
  tree_->Branch("B_Ds_px3", &B_Ds_px3);
  tree_->Branch("B_Ds_py3", &B_Ds_py3);
  tree_->Branch("B_Ds_pz3", &B_Ds_pz3);
  tree_->Branch("B_Ds_px3_track", &B_Ds_px3_track);
  tree_->Branch("B_Ds_py3_track", &B_Ds_py3_track);
  tree_->Branch("B_Ds_pz3_track", &B_Ds_pz3_track);
  tree_->Branch("B_Ds_charge3", &B_Ds_charge3);

  tree_->Branch("B_mumu_pt1", &B_mumu_pt1);
  tree_->Branch("B_mumu_px1", &B_mumu_px1);
  tree_->Branch("B_mumu_py1", &B_mumu_py1);
  tree_->Branch("B_mumu_pz1", &B_mumu_pz1);
  tree_->Branch("B_mumu_charge1", &B_mumu_charge1);

  tree_->Branch("B_mumu_pt2", &B_mumu_pt2);
  tree_->Branch("B_mumu_px2", &B_mumu_px2);
  tree_->Branch("B_mumu_py2", &B_mumu_py2);
  tree_->Branch("B_mumu_pz2", &B_mumu_pz2);
  tree_->Branch("B_mumu_charge2", &B_mumu_charge2);

  tree_->Branch("B_chi2", &B_chi2);
  tree_->Branch("B_chi2dof", &B_chi2dof);
  tree_->Branch("B_Ds_chi2", &B_Ds_chi2);
  tree_->Branch("B_mumu_chi2", &B_mumu_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_Ds_Prob", &B_Ds_Prob);
  tree_->Branch("B_mumu_Prob",  &B_mumu_Prob);
       
  // *************************

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_->Branch("bDecayVtxX",&bDecayVtxX);
  tree_->Branch("bDecayVtxY",&bDecayVtxY);
  tree_->Branch("bDecayVtxZ",&bDecayVtxZ);
  tree_->Branch("bDecayVtxXE",&bDecayVtxXE);
  tree_->Branch("bDecayVtxYE",&bDecayVtxYE);
  tree_->Branch("bDecayVtxZE",&bDecayVtxZE);
  tree_->Branch("bDecayVtxXYE",&bDecayVtxXYE);
  tree_->Branch("bDecayVtxXZE",&bDecayVtxXZE);
  tree_->Branch("bDecayVtxYZE",&bDecayVtxYZE);

  tree_->Branch("VDecayVtxX",&VDecayVtxX);
  tree_->Branch("VDecayVtxY",&VDecayVtxY);
  tree_->Branch("VDecayVtxZ",&VDecayVtxZ);
  tree_->Branch("VDecayVtxXE",&VDecayVtxXE);
  tree_->Branch("VDecayVtxYE",&VDecayVtxYE);
  tree_->Branch("VDecayVtxZE",&VDecayVtxZE);
  tree_->Branch("VDecayVtxXYE",&VDecayVtxXYE);
  tree_->Branch("VDecayVtxXZE",&VDecayVtxXZE);
  tree_->Branch("VDecayVtxYZE",&VDecayVtxYZE);

  tree_->Branch("pVtxIPX",     &pVtxIPX);
  tree_->Branch("pVtxIPY",     &pVtxIPY);
  tree_->Branch("pVtxIPZ",     &pVtxIPZ);
  tree_->Branch("pVtxIPXE",     &pVtxIPXE);
  tree_->Branch("pVtxIPYE",     &pVtxIPYE);
  tree_->Branch("pVtxIPZE",     &pVtxIPZE);
  tree_->Branch("pVtxIPXYE",     &pVtxIPXYE);
  tree_->Branch("pVtxIPXZE",     &pVtxIPXZE);
  tree_->Branch("pVtxIPYZE",     &pVtxIPYZE);
  tree_->Branch("pVtxIPCL",     &pVtxIPCL);

  // tree_->Branch("priRfVtxX",&priRfVtxX);
  // tree_->Branch("priRfVtxY",&priRfVtxY);
  // tree_->Branch("priRfVtxZ",&priRfVtxZ);
  // tree_->Branch("priRfVtxXE",&priRfVtxXE);
  // tree_->Branch("priRfVtxYE",&priRfVtxYE);
  // tree_->Branch("priRfVtxZE",&priRfVtxZE);
  // tree_->Branch("priRfVtxXYE",&priRfVtxXYE);
  // tree_->Branch("priRfVtxXZE",&priRfVtxXZE);
  // tree_->Branch("priRfVtxYZE",&priRfVtxYZE);
  // tree_->Branch("priRfVtxCL",&priRfVtxCL);
  // tree_->Branch("priRfNTrkDif",&priRfNTrkDif);

  tree_->Branch("pi1dxy",&pi1dxy);
  tree_->Branch("pi2dxy",&pi2dxy);
  tree_->Branch("pi1dz",&pi1dz);
  tree_->Branch("pi2dz",&pi2dz);

  tree_->Branch("pi1dxy_e",&pi1dxy_e);
  tree_->Branch("pi2dxy_e",&pi2dxy_e);
  tree_->Branch("pi1dz_e",&pi1dz_e);
  tree_->Branch("pi2dz_e",&pi2dz_e);

  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
  tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk); 
  tree_->Branch("dr0",&dr0);
  tree_->Branch("dr1",&dr1);
  tree_->Branch("dpt0",&dpt0);
  tree_->Branch("dpt1",&dpt1);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);



  tree_->Branch("B_pvip",&B_pvip);
  tree_->Branch("B_pviperr",&B_pviperr);
  tree_->Branch("B_pvips",&B_pvips);
  tree_->Branch("B_pvlzip",&B_pvlzip);
  tree_->Branch("B_pvlziperr",&B_pvlziperr);
  tree_->Branch("B_pvlzips",&B_pvlzips);
  tree_->Branch("B_pv2ip",&B_pv2ip);
  tree_->Branch("B_pv2iperr",&B_pv2iperr);
  tree_->Branch("B_pv2ips",&B_pv2ips);
  tree_->Branch("B_pv2lzip",&B_pv2lzip);
  tree_->Branch("B_pv2lziperr",&B_pv2lziperr);
  tree_->Branch("B_pv2lzips",&B_pv2lzips);
  tree_->Branch("B_l3d_pv2",&B_l3d_pv2);
  tree_->Branch("B_l3dE_pv2",&B_l3dE_pv2);
  tree_->Branch("B_l3d",&B_l3d);
  tree_->Branch("B_l3dE",&B_l3dE);
  tree_->Branch("B_lxy", &B_lxy);
  tree_->Branch("B_lxyE",&B_lxyE);
  tree_->Branch("B_cosalpha",&B_cosalpha);  
  tree_->Branch("B_cosalphaxy",&B_cosalphaxy);
  tree_->Branch("alpha",&alpha);
  tree_->Branch("B_treco",&B_treco);
  tree_->Branch("B_trecoe",&B_trecoe);
  tree_->Branch("B_trecoxy",&B_trecoxy);
  tree_->Branch("B_trecoxye",&B_trecoxye);
  tree_->Branch("B_iso",&B_iso);
  tree_->Branch("B_mum_iso",&B_mum_iso);
  tree_->Branch("B_mup_iso",&B_mup_iso);
  tree_->Branch("B_pi1_iso",&B_pi1_iso);
  tree_->Branch("B_pi2_iso",&B_pi2_iso);
  
  
  //gen information
  if (isMC_) {
    tree_->Branch("gen_b_p4",     "TLorentzVector",  &gen_b_p4);
    tree_->Branch("gen_ks_p4",   "TLorentzVector",  &gen_ks_p4);
    tree_->Branch("gen_pion1_p4",  "TLorentzVector",  &gen_pion1_p4);
    tree_->Branch("gen_pion2_p4",  "TLorentzVector",  &gen_pion2_p4);
    tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
    tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
    tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
    tree_->Branch("gen_b_vtx",    "TVector3",        &gen_b_vtx);
    tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
    tree_->Branch("gen_b_ct",     &gen_b_ct,        "gen_b_ct/F");
    }
  tree_->Branch("istruemum",  &istruemum );
  tree_->Branch("istruemup",  &istruemup );
  tree_->Branch("istruekp",   &istruekp  );
  tree_->Branch("istruekm",   &istruekm  );
  tree_->Branch("istruebs",   &istruebs  );
  
  tree_->Branch("bunchXingMC",&bunchXingMC);
  tree_->Branch("numInteractionsMC",&numInteractionsMC);
  tree_->Branch("trueNumInteractionsMC",&trueNumInteractionsMC);
    //}


}
// ------------ method called once each job just after ending the event loop  ------------
void BcToDsMuMu::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
  // for(int i = 0; i < kHistNameSize; i++) {
  //   histos[i]->Write();
  //   histos[i]->Delete();
  // }
}

//define this as a plug-in
DEFINE_FWK_MODULE(BcToDsMuMu);

