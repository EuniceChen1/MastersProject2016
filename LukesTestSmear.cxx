#include <algorithm>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

// nuwro reweight
#include "NuwroReWeightSimple.h"

// nuwro vali
#include "SimpleAnalysisFormat.h"

#include "Smearception.hxx"

bool IsSignal(PODSimpleAnalysisFormat const& ev) {
  // Define Signal cut here
  return true;
}

int main(int argc, char const* argv[]) {
  if (argc != 3) {
    std::cout << "[ERROR]: Expected 2 CLI args. <Input SAF data> <Muon "
                 "Smearing Matrix File>"
              << std::endl;
    return 1;
  }

  std::vector<PODSimpleAnalysisFormat> SAFData;
  if (!SRW::LoadSignalEventsIntoSAFVector(argv[1], SAFData, &IsSignal)) {
    return 1;
  }

  Smearception smearer(0);

  TFile* InputFile = TFile::Open(argv[2]);
  TH2D* SmearMat = dynamic_cast<TH2D*>(InputFile->Get("smear3D"));
  TH3D* TrueKineMap = dynamic_cast<TH3D*>(InputFile->Get("pthetaphiTH3_true"));
  TH3D* ReconKineMap = dynamic_cast<TH3D*>(InputFile->Get("pthetaphiTH3_rec"));

  smearer.AddSmearablePDG(13, SmearMat, TrueKineMap, ReconKineMap);

  // Create output File and output Tree
  TFile* outputFile = TFile::Open("SmearedOutputTree.root", "RECREATE");
  TTree* outputTree = new TTree("SmearedData", "");

  TLorentzVector* Smeared = 0;
  Bool_t DidSmear = false;

  outputTree->Branch("Smeared_4Mom_GeV", &Smeared);
  outputTree->Branch("DidSmear", &DidSmear, "DidSmear/O");

  TH1::SetDefaultSumw2(true);

  for (size_t i = 0; i < SAFData.size(); ++i) {

    TLorentzVector ScaledUp(SAFData[i].HMFSLepton_4Mom);
    ScaledUp *= 1000.0;
    Smeared->SetXYZT(0, 0, 0, 0);
    DidSmear = smearer.SmearTrueParticle(SAFData[i].HMFSLepton_PDG, ScaledUp,
                                         (*Smeared));

    if (DidSmear) {
      (*Smeared) *= (1.0 / 1000.0);
    } else {
      (*Smeared) = SAFData[i].HMFSLepton_4Mom;
    }
    outputTree->Fill();
  }

  outputFile->Close();
}
