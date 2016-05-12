#include <map>

#include "TH3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

struct SparseHistSlice {
  std::vector< std::pair<Int_t, Float_t> > Bins;
  Float_t Max;
  Int_t NOrigBinsX;
  SparseHistSlice() : Bins(), Max(0), NOrigBinsX(0) {};
  Int_t GetNbinsX(){ return Bins.size(); }
  Float_t GetBinContent_Index(size_t i){ return Bins[i].second; }
  Float_t GetBinNumber_Index(size_t i){ return Bins[i].first; }
  Float_t GetBinContent_BinNum(Int_t i);
  Float_t GetMaximum() { return Max; }
  Int_t BuildFromTH2D(TH2D const *Hist, Int_t YBinNum) ;
  TH1D * RecreateFullSlice();
};

struct Smearception {
  TRandom3 *rnd;
  std::map<Int_t, TH2D *> SmearingMatrices;
  std::map<Int_t, TH3D *> TrueKineMaps;
  std::map<Int_t, TH3D *> ReconKineMaps;
  std::map<TH2D *, std::map<Int_t, SparseHistSlice > > ReconSliceCache;

  void AddSmearablePDG(Int_t PDG, TH2D *SmearMat, TH3D *TrueKineMap,
                       TH3D *ReconKineMap);

  Int_t GetGlobalBinNumber(TH3D *KinematicMappingMatrix,
                           Double_t (&FourMom)[4]);

  Int_t GetGlobalBinNumber(TH3D *KinematicMappingMatrix, TLorentzVector(&p4));

  SparseHistSlice & GetProj(TH2D *SmearingMatrices, Int_t TrueGlobalBinNumber);

  Int_t GetSmearedReconBinNumber(TH2D *SmearingMatrices,
                                 Int_t TrueGlobalBinNumber);

  void ThrowReconSmearedKinematics(TH3D *KinematicMappingMatrix,
                                   Int_t ReconGlobalBinNumber, Int_t StdHepPDG,
                                   Double_t (&SmearedFourMom)[4]);

  void ThrowReconSmearedKinematics(TH3D *KinematicMappingMatrix,
                                   Int_t ReconGlobalBinNumber, Int_t StdHepPDG,
                                   TLorentzVector(&SmearedFourMom));

  bool SmearTrueParticle(Int_t StdHepPDG, Double_t (&StdHepP4)[4],
                         Double_t (&SmearedFourMom)[4]);

  bool SmearTrueParticle(Int_t StdHepPDG, TLorentzVector(&p4),
                         TLorentzVector(&Smearedp4));

  Smearception(Int_t RndSeed);

  ~Smearception() { delete rnd; }
};
