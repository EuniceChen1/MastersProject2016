#include <iostream>
#include <vector>
#include <algorithm>

#include "TH1.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

//nuwro reweight
#include "NuwroReWeightSimple.h"

//nuwro vali
#include "SimpleAnalysisFormat.h"

#include "Smearception.hxx"

bool IsSignal(PODSimpleAnalysisFormat const & ev){
  //Define Signal cut here
  // return ev.NuWroCC;
   return true;
}

int main(int argc, char const * argv[]){

  std::vector<PODSimpleAnalysisFormat> SAFData;
  if(!SRW::LoadSignalEventsIntoSAFVector(argv[1],SAFData,&IsSignal)){
    return NULL;
  }
  
  Smearception smearer(0);

  TFile *InputFile = TFile::Open(argv[2]);
  TH2D* SmearMat = dynamic_cast<TH2D*>(InputFile->Get("smear3D"));
  TH3D* TrueKineMap = dynamic_cast<TH3D*>(InputFile->Get("pthetaphiTH3_true"));
  TH3D* ReconKineMap = dynamic_cast<TH3D*>(InputFile->Get("pthetaphiTH3_rec"));
  
  smearer.AddSmearablePDG(13,SmearMat,TrueKineMap,ReconKineMap);

  //Create output File and output Tree
  TFile *outputFile = TFile::Open("SmearedOutputTree.root","RECREATE");
  TTree *outputTree = new TTree("SmearedData","");

  TLorentzVector *Smeared = 0;
  Bool_t DidSmear = false;

  outputTree->Branch("Smeared_4Mom_GeV",&Smeared);
  outputTree->Branch("DidSmear",&DidSmear,"DidSmear/O");  

  TH1::SetDefaultSumw2(true);

  //Creating new histograms for each graph
  TH1D* DatapHist = 
    new TH1D("DatapHisto",";#mu Momentum (GeV);#sigma",100,0,2);
  TH1D* DatathetaHist = 
    new TH1D("DatathetaHisto",";#theta;#sigma",100,0,2);
  TH1D* DataphiHist = 
    new TH1D("DataphiHisto",";#phi;#sigma",100,0,2);  

  TH1D* TrueAcceptedHist = 
    new TH1D("TrueAcceptedHisto",";#mu Momentum (GeV);#sigma",100,0,2);

  TH1D* Efficiency = 
    new TH1D("Efficiency", ";#mu Momentum (GeV);Efficiency",100,0,2);

  TH1D* SmearpHist = 
    new TH1D("SmearedpHisto",";#mu Momentum (GeV);#sigma",100,0,2);
  TH1D* SmearthetaHist = 
    new TH1D("SmearedthetaHisto",";#theta;#sigma",100,0,2);
  TH1D* SmearphiHist = 
    new TH1D("SmearedphiHisto",";#phi;#sigma",100,0,2);

  TH1D* UnSmearpHist = 
    new TH1D("UnSmearedHisto",";#mu Momentum (GeV);#sigma",100,0,2);
  TH1D* UnSmearthetaHist = 
    new TH1D("UnSmearedthetaHisto",";#theta;#sigma",100,0,2);
  TH1D* UnSmearphiHist = 
    new TH1D("UnSmearedphiHisto",";#phi;#sigma",100,0,2);

  TH1D* pDiffHist = 
    new TH1D("pDiffHist",";#mu Momentum (GeV);#sigma",100,0,2);



  for (size_t i = 0; i < std::min(size_t(-1),SAFData.size()); ++i){

    DatapHist->Fill(SAFData[i].HMFSLepton_4Mom.Vect().Mag(), 
		   SAFData[i].EvtWght*1E38);
    DatathetaHist->Fill(SAFData[i].HMFSLepton_4Mom.Vect().Theta(),
		     SAFData[i].EvtWght*1E38);
    DataphiHist->Fill(SAFData[i].HMFSLepton_4Mom.Vect().Phi(),
		      SAFData[i].EvtWght*1E38);

    TLorentzVector ScaledUp(SAFData[i].HMFSLepton_4Mom);
    ScaledUp *= 1000.0;
    Smeared->SetXYZT(0,0,0,0);
    DidSmear = smearer.SmearTrueParticle(SAFData[i].HMFSLepton_PDG,
					 ScaledUp,
					 (*Smeared));

     if(DidSmear){
        (*Smeared) *= (1.0/1000.0);
        SmearpHist->Fill(Smeared->Vect().Mag(), 
			SAFData[i].EvtWght*1E38);
        SmearthetaHist->Fill(Smeared->Vect().Theta(),
			     SAFData[i].EvtWght*1E38);
	SmearphiHist->Fill(Smeared->Vect().Phi(),
			   SAFData[i].EvtWght*1E38);

	//True mumom that was accepted through rejection sampling
	TrueAcceptedHist->Fill(SAFData[i].HMFSLepton_4Mom.Vect().Mag(), 
			       SAFData[i].EvtWght*1E38);
	
	Efficiency->Divide(TrueAcceptedHist,DatapHist);
        
	//Obtaining the difference in histograms
        pDiffHist->Fill((Smeared->Vect().Mag()-SAFData[i].HMFSLepton_4Mom.Vect().Mag()), SAFData[i].EvtWght*1E38);


     }
     else {
      ScaledUp *= (1.0/1000.0);
      UnSmearpHist->Fill(ScaledUp.Vect().Mag(),
			SAFData[i].EvtWght*1E38);
      UnSmearthetaHist->Fill(ScaledUp.Vect().Theta(),
			     SAFData[i].EvtWght*1E38);
      UnSmearphiHist->Fill(ScaledUp.Vect().Phi(),
			   SAFData[i].EvtWght*1E38);

     }
     outputTree->Fill();
  }

  //Set Line Color and Style for plots
  SmearpHist->SetLineColor(kRed);
  SmearpHist->SetLineStyle(2);
  SmearthetaHist->SetLineColor(kRed);
  SmearthetaHist->SetLineStyle(2);
  SmearphiHist->SetLineColor(kRed);
  SmearphiHist->SetLineStyle(2);

  TrueAcceptedHist->SetLineColor(kRed);
  TrueAcceptedHist->SetLineStyle(2);

  UnSmearpHist->SetLineColor(kBlack);
  UnSmearpHist->SetLineStyle(3);
  UnSmearthetaHist->SetLineColor(kBlack);
  UnSmearthetaHist->SetLineStyle(3);
  UnSmearphiHist->SetLineColor(kBlack);
  UnSmearphiHist->SetLineStyle(3);

  //Plots and saving canvases
  TCanvas* c1 = new TCanvas("canv","");
  DatapHist->Draw("EHIST");
  UnSmearpHist->Draw("EHIST SAME");
  SmearpHist->Draw("EHIST SAME");
  c1->SaveAs("muSmearing500k.png");

  TCanvas *c2 = new TCanvas("c2","");
  DatathetaHist->Draw("EHIST");
  UnSmearthetaHist->Draw("EHIST SAME");
  SmearthetaHist->Draw("EHIST SAME");
   c2->SaveAs("thetaSmearing500k.png"); 

  TCanvas *c3 = new TCanvas("c3","");
  DataphiHist->Draw("EHIST");
  UnSmearphiHist->Draw("EHIST SAME");
  SmearphiHist->Draw("EHIST SAME");
  c3->SaveAs("phiSmearing500k.png"); 

  TCanvas *c4 = new TCanvas("c4","");
  DatapHist->Draw("EHIST");
  TrueAcceptedHist->Draw("EHIST SAME");
  c4->SaveAs("TrueAccepted.png");

  TCanvas *c5 = new TCanvas("eff","");
  Efficiency->Draw("HIST");
  c5->SaveAs("Efficiency.png");
 
  TCanvas *c6 = new TCanvas("diff","");
  pDiffHist->Draw("HIST");
  c6->SaveAs("pDiffHist.png");
 

  outputTree->Write();
  outputFile->Write();
  outputFile->Close();
}
