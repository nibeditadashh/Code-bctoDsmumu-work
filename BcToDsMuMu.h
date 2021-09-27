#ifndef _BcToDsMuMu_h
#define _BcToDsMuMu_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "DataFormats/Math/interface/Vector3D.h"

//
// struct decleration
//
typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > cov33_t;
typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepSym<double,6> > cov66_t;
typedef ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > cov77_t;
typedef ROOT::Math::SMatrix<double,9,9,ROOT::Math::MatRepSym<double,9> > cov99_t;
typedef ROOT::Math::SVector<double,9> jac9_t;

//class definition
class BcToDsMuMu : public edm::EDAnalyzer {
public:
  explicit BcToDsMuMu(const edm::ParameterSet&);
  ~BcToDsMuMu();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSametk(const pat::GenericParticle& tk, const pat::Muon& mu);
  bool IsTheSamePFtk(const pat::PackedCandidate& tk, const pat::Muon& mu);
  bool IsTheSame(const reco::Track& tk, const pat::Muon& mu);
  bool   isAncestor(const reco::Candidate*, const reco::Candidate*);
  double GetLifetime(TLorentzVector, TVector3, TVector3);


  cov99_t makeCovarianceMatrix(cov33_t cov_vtx1, cov77_t cov_vtx2);
  cov33_t GlobalError2SMatrix_33(GlobalError);

  jac9_t makeJacobianVector3d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2,
			      const AlgebraicVector3 &momentum);
  //jac9_t makeJacobianVector3d(const GlobalPoint &vtx1, const GlobalPoint &vtx2, const TVector3 &tv3momentum);
  jac9_t makeJacobianVector3d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			      const GlobalPoint &vtx2, const TVector3 &tv3momentum);
  jac9_t makeJacobianVector2d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2,
			      const AlgebraicVector3 &momentum);
  //jac9_t makeJacobianVector2d(const GlobalPoint &vtx1, const GlobalPoint &vtx2, const TVector3 &tv3momentum);
  jac9_t makeJacobianVector2d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			      const GlobalPoint &vtx2, const TVector3 &tv3momentum);
  void saveIP(const RefCountedKinematicTree& vertexFitTree,
	      const reco::VertexCollection& vertices, int & pvIndex);
  void SaveIso(const RefCountedKinematicTree& vertexFitTree, const MagneticField *fMagneticField,
	       edm::ESHandle<TransientTrackBuilder> theB,
	       unsigned int  pvIndex,
	       //edm::Handle<std::vector<pat::GenericParticle> > thePATTrackHandle,
	       edm::Handle< edm::View<pat::PackedCandidate> > thePATTrackHandle,
	       edm::Handle< edm::View<pat::Muon> > thePATMuonHandle,
	       const reco::TransientTrack muon1TT,
	       const reco::TransientTrack muon2TT,
	       const reco::TransientTrack pion1TT,
	       const reco::TransientTrack pion2TT,
	       const reco::TransientTrack pion3TT);
  double calEta(double, double, double);
  double calPhi(double, double, double);
  double calEtaPhiDistance (double, double, double, double, double, double);
  void saveTruthMatch(const edm::Event& iEvent);
  void savePUinMC(const edm::Event& iEvent);

  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  //  edm::EDGetTokenT<std::vector<pat::GenericParticle>> trakCollection_label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  std::vector<std::string> trigTable_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> v0PtrCollection_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;

  const MagneticField   *fMagneticField;
  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;
  bool doMC_;
  vector<string> TriggerNames_;

  TTree*   tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;
  
  std::vector<float>       *mumC2;
  std::vector<int>         *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<float>       *muon_dca;

  std::vector<int>         *tri_Dim25, *tri_JpsiTk, *tri_JpsiTkTk;
  std::vector<float>       *dr0,*dr1, *dpt0, *dpt1;
  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  

  int                      muAcc, muTrig, weight;
 
  // vertice primario CON mayor Pt
  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;
 
  // ********************************** ************************************************************************

  unsigned int             nB;  
  unsigned int             nMu;
  
  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz;
 
  std::vector<float>       *B_Ds_mass, *B_Ds_px, *B_Ds_py, *B_Ds_pz;
  std::vector<float>       *B_Ds_pt1, *B_Ds_px1, *B_Ds_py1, *B_Ds_pz1;
  std::vector<float>       *B_Ds_pt2, *B_Ds_px2, *B_Ds_py2, *B_Ds_pz2;
  std::vector<float>       *B_Ds_pt3, *B_Ds_px3, *B_Ds_py3, *B_Ds_pz3;
  
  std::vector<float>       *B_Ds_px1_track, *B_Ds_py1_track, *B_Ds_pz1_track;
  std::vector<float>       *B_Ds_px2_track, *B_Ds_py2_track, *B_Ds_pz2_track;
  std::vector<float>       *B_Ds_px3_track, *B_Ds_py3_track, *B_Ds_pz3_track;

  std::vector<float>       *pi1dxy, *pi2dxy, *pi1dz, *pi2dz;
  std::vector<float>       *pi1dxy_e, *pi2dxy_e, *pi1dz_e, *pi2dz_e;
  std::vector<int>         *B_Ds_charge1, *B_Ds_charge2, *B_Ds_charge3;
  
  std::vector<float>       *B_mumu_mass, *B_mumu_px, *B_mumu_py, *B_mumu_pz;
  std::vector<float>       *B_mumu_pt1, *B_mumu_px1, *B_mumu_py1, *B_mumu_pz1;
  std::vector<float>       *B_mumu_pt2, *B_mumu_px2, *B_mumu_py2, *B_mumu_pz2;
  std::vector<int>         *B_mumu_charge1, *B_mumu_charge2;
  
  std::vector<float>       *B_Ds_chi2, *B_mumu_chi2, *B_chi2, *B_chi2dof;
  std::vector<float>       *B_Prob, *B_mumu_Prob, *B_Ds_Prob;


  std::vector<float>       *bDecayVtxX, *bDecayVtxY, *bDecayVtxZ;
  std::vector<double>      *bDecayVtxXE, *bDecayVtxYE, *bDecayVtxZE;
  std::vector<double>      *bDecayVtxXYE, *bDecayVtxXZE, *bDecayVtxYZE;

  std::vector<float>       *VDecayVtxX, *VDecayVtxY, *VDecayVtxZ;
  std::vector<float>       *VDecayVtxXE, *VDecayVtxYE, *VDecayVtxZE;
  std::vector<float>       *VDecayVtxXYE, *VDecayVtxXZE, *VDecayVtxYZE;

  // *************************************
  std::vector<float>       *pVtxIPX,  *pVtxIPY, *pVtxIPZ, *pVtxIPXE, *pVtxIPYE, *pVtxIPZE, *pVtxIPCL;
  std::vector<float>       *pVtxIPXYE,  *pVtxIPXZE, *pVtxIPYZE;

  // refitting the primary without the tracks in the B reco candidate
   std::vector<float>       *priRfVtxX, *priRfVtxY, *priRfVtxZ, *priRfVtxXE, *priRfVtxYE, *priRfVtxZE, *priRfVtxCL; 
   std::vector<float>       *priRfVtxXYE, *priRfVtxXZE, *priRfVtxYZE; 
   std::vector<int>         *priRfNTrkDif; 
  
  std::vector<float>       *B_l3d,  *B_l3dE,  *B_lxy, *B_lxyE ;
  std::vector<float>       *B_cosalpha ,  *B_cosalphaxy, *alpha, *B_treco ,  *B_trecoe,  *B_trecoxy, *B_trecoxye;
  std::vector<float>       *B_pvip , *B_pviperr, *B_pvips, *B_pvlzip, *B_pvlziperr,*B_pvlzips;
  std::vector<float>       *B_pv2ip , *B_pv2iperr, *B_pv2ips, *B_pv2lzip, *B_pv2lziperr,*B_pv2lzips;
  std::vector<float>       *B_l3d_pv2,  *B_l3dE_pv2;
  std::vector<float>       *B_iso, *B_mum_iso, *B_mup_iso, *B_pi1_iso, *B_pi2_iso;

  std::vector<bool>        *istruemum, *istruemup, *istruekp, *istruekm, *istruekmm,  *istruebs;
  std::vector<float>       *bunchXingMC, *numInteractionsMC, *trueNumInteractionsMC;  
  int  run, event;
  int  lumiblock;

  TLorentzVector gen_b_p4,gen_ks_p4,gen_pion1_p4,gen_pion2_p4,gen_pion3_p4,gen_jpsi_p4,gen_muon1_p4,gen_muon2_p4,gen_ds_p4;
  TVector3       gen_b_vtx,gen_jpsi_vtx,gen_ds_vtx;
  float          gen_b_ct;


};


#endif
