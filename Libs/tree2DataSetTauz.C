//#include "Common/CCDB/EventSelectionParams.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "initOniaTree.C"
#include <TH3.h>

void readListTauz(string filename, std::vector<string> &dataMap) {
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

RooDataSet *tree2DataSetTauz(string inputName, string TreeName,
                             RooRealVar *mass, RooRealVar *pt, RooRealVar *y,
                             RooRealVar *tauz, RooRealVar *chi2mchmft1,
                             RooRealVar *chi2mchmft2, double ptCut = 0.7,
                             bool useLS = 0, double rapMax = 3.6,
                             bool applyChi2Cut = true) {
  vector<string> InputFileNames;
  readListTauz(inputName, InputFileNames);

  TChain *theTree = new TChain(TreeName.c_str(), "");
  cout << "[INFO] Creating TChain" << endl;
  getTChain(theTree, InputFileNames, TreeName);
  initOniaTree(theTree, TreeName);
  iniBranch(theTree, TreeName);

  RooArgSet *cols =
      new RooArgSet(*mass, *pt, *y, *tauz, *chi2mchmft1, *chi2mchmft2);
  RooDataSet *data = new RooDataSet("dOS", "dOS", *cols);

  Long64_t nentries = theTree->GetEntries();
  cout << "[INFO] Starting to process " << nentries << " nentries" << endl;

  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    theTree->GetEntry(jentry);
    float ct = fTauz * 299792458.e-7 * 10.;

    mass->setVal(fMass);
    pt->setVal(fPt);
    tauz->setVal(ct);
    chi2mchmft1->setVal(fChi2MatchMCHMFT1);
    chi2mchmft2->setVal(fChi2MatchMCHMFT2);

    float pz = fPt * sinh(fEta);
    float E = sqrt(fMass * fMass + fPt * fPt + pz * pz);
    float rap = -0.5 * log((E + pz) / (E - pz));
    y->setVal(rap);

    bool signOK = (fSign == 0 && !useLS) ||
                  (useLS && (fSign == 2 || fSign == -2));
    bool rapOK = rap > 2.5 && rap < rapMax;
    bool chi2OK = !applyChi2Cut ||
                  (fChi2MatchMCHMFT1 < 40 && fChi2MatchMCHMFT2 < 40);
    if (signOK && rapOK && chi2OK) {
      data->add(*cols, 1.0);
    }
  }

  theTree->Reset();
  delete theTree;
  return data;
}
