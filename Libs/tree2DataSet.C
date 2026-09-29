//#include "Common/CCDB/EventSelectionParams.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "initOniaTree.C"
#include <TH3.h>

void readList(string filename, std::vector<string> &dataMap) {
  std::ifstream inputFile(filename);
  if (!inputFile) {
    std::cerr << "Error opening file." << std::endl;
    return;
  }
  std::string line;
  while (std::getline(inputFile, line)) {
    dataMap.push_back(line);
  }
  inputFile.close();
}

RooDataSet *tree2DataSet(string inputName, string TreeName,
                         RooRealVar *mass, RooRealVar *pt, RooRealVar *y, double ptCut = 0.7, bool useLS = 0) {
  vector<string> InputFileNames;
  readList(inputName, InputFileNames);

  TChain *theTree = new TChain(TreeName.c_str(), "");
  cout << "[INFO] Creating TChain" << endl;
  getTChain(theTree, InputFileNames, TreeName); // Import files to TChain
  initOniaTree(theTree, TreeName);              // Initialize the Onia Tree
  iniBranch(theTree, TreeName);                 // Initialize the Branches

  RooArgSet *cols = new RooArgSet(*mass, *pt, *y); // PENSER A AJOUTER MASSE

  RooDataSet *data = new RooDataSet("dOS", "dOS", *cols);

  Long64_t nentries = theTree->GetEntries();

  cout << "[INFO] Starting to process " << nentries << " nentries" << endl;
  
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    theTree->GetEntry(jentry);

    mass->setVal(fMass);

    pt->setVal(fPt);

    float pz = fPt * sinh(fEta);
    float E  = sqrt(fMass*fMass + fPt*fPt + pz*pz);
    float rap = - 0.5 * log((E + pz)/(E - pz));
    y->setVal(rap);

//    if (fEta1 < -2.5 && fEta1 > -4.0 && fEta2 < -2.5 && fEta2 > -4.0) {
//        if (fPt1 > ptCut && fPt2 > ptCut) {
//          if (fIsAmbig1 == 0 && fIsAmbig2 == 0) {
            if ((fSign == 0 && !useLS) || (useLS && (fSign == 2 || fSign == -2))) {
              if (rap > 2.5 && rap < 3.6 && fChi2MatchMCHMFT1 < 40 && fChi2MatchMCHMFT2 < 40) {
                      data->add(*cols, 1.0); // Signal and background dimuons
              } // rap
            } // sign
//          } // ambig
//        } // single pT
//    } // eta
  } // for loop
  theTree->Reset();
  delete theTree;
  return data;
}



// New function: build invariant mass histogram for a given pT bin
TH1D* tree2MassHist(string inputName, string TreeName,
                    double ptLow, double ptHigh,
                    int nMassBins = 50, double massMin = 2.5, double massMax = 4.0, double ptCut = 0.0)
{
  // Read list of ROOT files
  vector<string> InputFileNames;
  readList(inputName, InputFileNames);

  TChain *theTree = new TChain(TreeName.c_str(), "");
  cout << "[INFO] Creating TChain" << endl;
  getTChain(theTree, InputFileNames, TreeName); // Import files to TChain
  initOniaTree(theTree, TreeName);              // Initialize the Onia Tree
  iniBranch(theTree, TreeName);                 // Initialize the Branches

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

