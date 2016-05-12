#include <algorithm>
#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TParticle.h>
#include <TVector3.h>

#include "Smearception.hxx"

#define SMEARCEPTION_DEBUG

std::ostream &operator<<(std::ostream &os, TLorentzVector const &tl) {
  return os << "[ |q3|: " << tl.Vect().Mag() << ", theta: " << tl.Vect().Theta()
            << ", phi: " << tl.Vect().Phi() << ", Mass: " << tl.M() << " ]";
}

Int_t SparseHistSlice::BuildFromTH2D(TH2D const *Hist, Int_t YBinNum) {
  for (Int_t i = 0; i < (Hist->GetXaxis()->GetNbins() + 2); ++i) {
    Int_t MatrixBin = Hist->GetBin(i, YBinNum);
    Float_t cont = Hist->GetBinContent(MatrixBin);
    if (cont <= 0) {
      continue;
    }
    Bins.push_back(std::make_pair(i, cont));
    Max = std::max(Max, Bins.back().second);
  }

  NOrigBinsX = Hist->GetXaxis()->GetNbins();
  return GetNbinsX();
}

Float_t SparseHistSlice::GetBinContent_BinNum(Int_t binnum) {
  for (size_t i = 0; i < Bins.size(); ++i) {
    if (Bins[i].first == binnum) {
      return Bins[i].second;
    }
  }
  return 0;
}

TH1D * SparseHistSlice::RecreateFullSlice(){
  TH1D * rec = new TH1D("recSlice","",(NOrigBinsX-2),0,(NOrigBinsX-2));
  for (size_t i = 0; i < NOrigBinsX; ++i){
    rec->SetBinContent(i,GetBinContent_BinNum(i));
  }
  return rec;
}

//################################################################
//#               Getting the Global Bin Number                  #
//################################################################

Int_t Smearception::GetGlobalBinNumber(TH3D *KinematicMappingMatrix,
                                       Double_t (&FourMom)[4]) {
  TLorentzVector p4(FourMom[0], FourMom[1], FourMom[2], FourMom[3]);
  return GetGlobalBinNumber(KinematicMappingMatrix, p4);
}

Int_t Smearception::GetGlobalBinNumber(TH3D *KinematicMappingMatrix,
                                       TLorentzVector(&p4)) {
  KinematicMappingMatrix->FindBin(p4.Vect().Mag(), p4.Theta(), p4.Phi());

  Int_t gbin =
      KinematicMappingMatrix->FindBin(p4.Vect().Mag(), p4.Theta(), p4.Phi());

#ifdef SMEARCEPTION_DEBUG
  std::cout << "[INFO]: True particle kinematics: " << p4
            << " , global bin: " << gbin << std::endl;
#endif
  return gbin;
}

//#################################################################
//#      Get Projection Slice and Apply Rejection Sampling        #
//#################################################################

SparseHistSlice &Smearception::GetProj(TH2D *SmearMatrix,
                                       Int_t TrueGlobalBinNumber) {
  if (ReconSliceCache[SmearMatrix].count(TrueGlobalBinNumber)) {
    return ReconSliceCache[SmearMatrix][TrueGlobalBinNumber];
  }

  ReconSliceCache[SmearMatrix][TrueGlobalBinNumber].BuildFromTH2D(
      SmearMatrix, TrueGlobalBinNumber);

#ifdef SMEARCEPTION_DEBUG
  Int_t smearMatPdg = 0;
  for (std::map<Int_t, TH2D *>::iterator sm_it = SmearingMatrices.begin();
       sm_it != SmearingMatrices.end(); ++sm_it) {
    if (sm_it->second == SmearMatrix) {
      smearMatPdg = sm_it->first;
      break;
    }
  }
  if (smearMatPdg) {
    std::cout << "[INFO]: Calculated new slice for PDG " << smearMatPdg << "("
              << SmearMatrix << ")"
              << " -- Bin:" << TrueGlobalBinNumber << ", contains: "
              << ReconSliceCache[SmearMatrix][TrueGlobalBinNumber].GetNbinsX()
              << " entries." << std::endl;
  }
#endif
  return ReconSliceCache[SmearMatrix][TrueGlobalBinNumber];
}

Int_t Smearception::GetSmearedReconBinNumber(TH2D *SmearMatrix,
                                             Int_t TrueGlobalBinNumber) {
  SparseHistSlice &myproj = GetProj(SmearMatrix, TrueGlobalBinNumber);

  if (myproj.GetNbinsX() < 1) {
    return -1;
  }

  // Accept-Reject Method
  bool Accepted = false;
  Int_t xbin, BinContentX, yrndm, ymax, xrndm;

  while (!Accepted) {
    xbin = myproj.GetNbinsX();
    xrndm = rnd->Integer(xbin);
    ymax = myproj.GetMaximum();
    yrndm = rnd->Integer(ymax);

    //+1 because GetBinContent 0 is the underflow bin
    BinContentX = myproj.GetBinContent_Index(xrndm);
    Accepted = (BinContentX > yrndm);
  }
  return myproj.GetBinNumber_Index(xrndm);
}

//################################################################
//#                 Get Fully Smeared Kinematics                 #
//################################################################

void Smearception::ThrowReconSmearedKinematics(TH3D *KinematicMappingMatrix,
                                               Int_t ReconGlobalBinNumber,
                                               Int_t StdHepPDG,
                                               Double_t (&SmearedFourMom)[4]) {
  TLorentzVector p4(SmearedFourMom[0], SmearedFourMom[1], SmearedFourMom[2],
                    SmearedFourMom[3]);
  ThrowReconSmearedKinematics(KinematicMappingMatrix, ReconGlobalBinNumber,
                              StdHepPDG, p4);
  SmearedFourMom[0] = p4.X();
  SmearedFourMom[1] = p4.Y();
  SmearedFourMom[2] = p4.Z();
  SmearedFourMom[3] = p4.E();
}

void Smearception::ThrowReconSmearedKinematics(
    TH3D *KinematicMappingMatrix, Int_t ReconGlobalBinNumber, Int_t StdHepPDG,
    TLorentzVector(&SmearedFourMom)) {
  Int_t ReconMuonBin, ReconTBin, ReconPhiBin;

  KinematicMappingMatrix->GetBinXYZ(ReconGlobalBinNumber, ReconMuonBin,
                                    ReconTBin, ReconPhiBin);

  Double_t LowEdgeBinX, LowEdgeBinY, LowEdgeBinZ;
  Double_t UpEdgeBinX, UpEdgeBinY, UpEdgeBinZ;
  Double_t ReconMuonMomX, ReconThetaY, ReconPhiZ;

  LowEdgeBinX = KinematicMappingMatrix->GetXaxis()->GetBinLowEdge(ReconMuonBin);
  UpEdgeBinX = KinematicMappingMatrix->GetXaxis()->GetBinUpEdge(ReconMuonBin);

  LowEdgeBinY = KinematicMappingMatrix->GetYaxis()->GetBinLowEdge(ReconTBin);
  UpEdgeBinY = KinematicMappingMatrix->GetYaxis()->GetBinUpEdge(ReconTBin);

  LowEdgeBinZ = KinematicMappingMatrix->GetZaxis()->GetBinLowEdge(ReconPhiBin);
  UpEdgeBinZ = KinematicMappingMatrix->GetZaxis()->GetBinUpEdge(ReconPhiBin);

  // Smeared kinematics
  // assume 3mom mag
  ReconMuonMomX = LowEdgeBinX + rnd->Uniform(UpEdgeBinX - LowEdgeBinX);
  // 3mom theta
  ReconThetaY = LowEdgeBinY + rnd->Uniform(UpEdgeBinY - LowEdgeBinY);
  // 3mom phi
  ReconPhiZ = LowEdgeBinZ + rnd->Uniform(UpEdgeBinZ - LowEdgeBinZ);

  TVector3 p3;
  p3.SetMagThetaPhi(ReconMuonMomX, ReconThetaY, ReconPhiZ);

  TParticle particle;
  particle.SetPdgCode(StdHepPDG);

  SmearedFourMom.SetXYZM(p3.X(), p3.Y(), p3.Z(), particle.GetMass() * 1000.0);

}

void Smearception::AddSmearablePDG(Int_t PDG, TH2D *SmearMat, TH3D *TrueKineMap,
                                   TH3D *ReconKineMap) {
  SmearingMatrices[PDG] = (TH2D *)SmearMat->Clone();
  TrueKineMaps[PDG] = (TH3D *)TrueKineMap->Clone();
  ReconKineMaps[PDG] = (TH3D *)ReconKineMap->Clone();

  std::cout << "[INFO]: Added smearing matrix definitions for "
               "PDG code: "
            << PDG << std::endl;
}

bool Smearception::SmearTrueParticle(Int_t StdHepPDG, Double_t (&StdHepP4)[4],
                                     Double_t (&SmearedFourMom)[4]) {
  TLorentzVector p4(StdHepP4[0], StdHepP4[1], StdHepP4[2], StdHepP4[3]);
  TLorentzVector Sp4(0, 0, 0, 0);
  bool smtpr = SmearTrueParticle(StdHepPDG, p4, Sp4);
  SmearedFourMom[0] = Sp4.X();  // momentum X
  SmearedFourMom[1] = Sp4.Y();  // momentum Y
  SmearedFourMom[2] = Sp4.Z();  // momentum Z
  SmearedFourMom[3] = Sp4.E();  // Energy
  return smtpr;
}

// A copy function with TLorentzVector instead of Doubles
bool Smearception::SmearTrueParticle(Int_t StdHepPDG, TLorentzVector(&p4),
                                     TLorentzVector(&Smearedp4)) {
  if (!SmearingMatrices.count(StdHepPDG)) {
    return false;
  }

  Int_t globbinnum = GetGlobalBinNumber(TrueKineMaps[StdHepPDG], p4);
  Int_t reconglobbinnum = -1;

  Int_t MomBin, dumdum, dumdumdumdum;
  //So that we don't return an overflow bin.
  do {
    reconglobbinnum =
        GetSmearedReconBinNumber(SmearingMatrices[StdHepPDG], globbinnum);

    ReconKineMaps[StdHepPDG]->GetBinXYZ(reconglobbinnum, MomBin, dumdum,
                                        dumdumdumdum);
  } while ((ReconKineMaps[StdHepPDG]->GetXaxis()->GetNbins() + 1) == MomBin);

  if (reconglobbinnum == -1) {
    return false;
  }

  if (!reconglobbinnum) {
    return false;
  }
  ThrowReconSmearedKinematics(ReconKineMaps[StdHepPDG], reconglobbinnum,
                              StdHepPDG, Smearedp4);
  return true;
}

Smearception::Smearception(Int_t RndSeed) { rnd = new TRandom3(RndSeed); }
