#include <algorithm>
#include <iostream>
#include <vector>

// POSIX
#include <sys/time.h>

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

double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main(int argc, char const* argv[]) {
  if (argc != 2) {
    std::cout << "[ERROR]: Expected 2 CLI args. <Muon Smearing Matrix File>"
              << std::endl;
    return 1;
  }

  Smearception smearer(0);

  TFile* InputFile = TFile::Open(argv[1]);
  TH2D* SmearMat = dynamic_cast<TH2D*>(InputFile->Get("smear3D")->Clone());
  TH3D* TrueKineMap =
      dynamic_cast<TH3D*>(InputFile->Get("pthetaphiTH3_true")->Clone());
  TH3D* ReconKineMap =
      dynamic_cast<TH3D*>(InputFile->Get("pthetaphiTH3_rec")->Clone());
  SmearMat->SetDirectory(0);
  TrueKineMap->SetDirectory(0);
  ReconKineMap->SetDirectory(0);
  InputFile->Close();
  smearer.AddSmearablePDG(13, SmearMat, TrueKineMap, ReconKineMap);

  Int_t GlobalBinNum = 0;
  for (Int_t i = 1; i < (SmearMat->GetYaxis()->GetNbins() + 1); ++i) {
    if (smearer.GetProj(SmearMat, i).GetNbinsX() > 20) {
      GlobalBinNum = i;
      break;
    }
  }

  std::cout << "[INFO]: Chose Global bin number: " << GlobalBinNum << " to test."
            << std::endl;

  // Create output File and output Tree
  TFile* outputFile = TFile::Open("TestRejSampling.root", "RECREATE");

  TH1::SetDefaultSumw2(true);

  TH1D* Proj = smearer.GetProj(SmearMat, GlobalBinNum).RecreateFullSlice();
  Proj->SetNameTitle("Projected", "");
  TH1D* Thrown = (TH1D*)Proj->Clone();

  Proj->SetDirectory(outputFile);
  Thrown->SetDirectory(outputFile);

  for (size_t i = 0; i < (Proj->GetNbinsX() + 2); ++i) {
    Thrown->SetBinContent(i, 0);
    Thrown->SetBinError(i, 0);
  }
  Thrown->Sumw2(false);
  Thrown->Sumw2(true);
  Thrown->SetNameTitle("Thrown", "");

  size_t NThrows = 1000000;
  std::cout << "[INFO]: Testing with " << NThrows << " throws..." << std::endl;
  double start = get_wall_time();

  for (size_t i = 0; i < NThrows; ++i) {
    Thrown->Fill(smearer.GetSmearedReconBinNumber(SmearMat, GlobalBinNum) - 1);
    std::cout << "\r[THROWN]: " << int((i + 1) * 100 / NThrows) << "%"
              << std::flush;
  }
  std::cout << std::endl;
  std::cout << NThrows << " throws of the same distribution took "
            << get_wall_time() - start << "s." << std::endl;

  Proj->SetLineColor(kRed);
  Proj->SetLineWidth(2);
  Proj->SetLineStyle(1);
  Thrown->SetLineColor(kBlack);
  Thrown->SetLineWidth(2);
  Thrown->SetLineStyle(2);

  Proj->Scale(1.0 / Proj->Integral());
  Thrown->Scale(1.0 / Thrown->Integral());

  outputFile->Write();

  start = get_wall_time();
  size_t nthr = 0;
  while (nthr < NThrows) {
    Int_t i = GlobalBinNum;
    for (; i < (SmearMat->GetYaxis()->GetNbins() + 1); ++i) {
      if (smearer.GetProj(SmearMat, i).GetNbinsX() > 20) {
        GlobalBinNum = i;
        break;
      }
    }
    if (i == (SmearMat->GetYaxis()->GetNbins() + 1)) {
      GlobalBinNum = 1;
      continue;
    }
    smearer.GetSmearedReconBinNumber(SmearMat, GlobalBinNum);
    nthr++;
  }
  std::cout << NThrows
            << " throws of random distributions (with more than 20 bins) took "
            << get_wall_time() - start << "s." << std::endl;
}
