#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
//#include "tree2DataSet.C"

// New function: build invariant mass histogram for a given pT bin
TH1D* tree2MassHist(string inputName, string TreeName,
                    double ptLow, double ptHigh,
                    int nMassBins = 50, double massMin = 2.5, double massMax = 4.0, double ptCut = 0.0) 
{
  // Read list of ROOT files
  vector<string> InputFileNames;
  readList(inputName, InputFileNames);

  // Create TChain
  TChain *theTree = new TChain(TreeName.c_str(), "");
  for (auto &fname : InputFileNames) theTree->Add(fname.c_str());

  // Branches we need
  float fMass, fPt, fEta1, fEta2, fPt1, fPt2, fChi2MatchMCHMFT1, fChi2MatchMCHMFT2;
  int fSign;
  theTree->SetBranchAddress("fMass", &fMass);
  theTree->SetBranchAddress("fPt", &fPt);
  theTree->SetBranchAddress("fEta1", &fEta1);
  theTree->SetBranchAddress("fEta2", &fEta2);
  theTree->SetBranchAddress("fPt1", &fPt1);
  theTree->SetBranchAddress("fPt2", &fPt2);
  theTree->SetBranchAddress("fSign", &fSign);
  theTree->SetBranchAddress("fChi2MatchMCHMFT1", &fChi2MatchMCHMFT1);
  theTree->SetBranchAddress("fChi2MatchMCHMFT2", &fChi2MatchMCHMFT2);

  // Create histogram
  TString hname = Form("hMass_pt%.1f_%.1f", ptLow, ptHigh);
  TH1D *hMass = new TH1D(hname, hname, nMassBins, massMin, massMax);
  hMass->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/c^{2})");
  hMass->GetYaxis()->SetTitle("Counts");

  Long64_t nentries = theTree->GetEntries();
  std::cout << "[INFO] Filling mass histogram for "
            << ptLow << " < pT < " << ptHigh
            << " with " << nentries << " entries" << std::endl;

  // Loop over entries
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    theTree->GetEntry(jentry);

    // Apply same event selection as in tree2DataSet
    if (fEta1 < -2.5 && fEta1 > -4.0 && fEta2 < -2.5 && fEta2 > -4.0) {
      if (fPt1 > ptCut && fPt2 > ptCut) {
        if (fSign == 2 || fSign == -2) {
          if (fPt > ptLow && fPt <= ptHigh) {
          if (fIsAmbig1 == 0 && fIsAmbig2 == 0) {
            if (fMass > massMin && fMass < massMax) {
              hMass->Fill(fMass);
            }
	  }
          }
        }
      }
    }
  }

  theTree->Reset();
  delete theTree;
  return hMass;
}

