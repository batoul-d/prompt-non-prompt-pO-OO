#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TSystem.h>
#include <TClass.h>
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TMath.h"
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>

#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooStats/SPlot.h"

#include "RooExtCBShape.h"
#include "TSystem.h"
extern TSystem *gSystem;

using namespace std;
using namespace  RooFit;

TChain* getInputTree(string inputTreeFile, string treeName);
void findAndAddTrees(TDirectory* dir, const std::string& treeName, TChain* chain, const std::string& fileName);
std::array<RooDataSet*, 4> createDataset(bool ispO = true, bool isMC=false, const char *methodName = "sPlot");
double getMax(RooHist* hist);
void fixPullStyle(RooPlot* pullFrame); //move from here maybe
bool parseFile(string FileName, vector< map<string, string> >& data);
bool parseString(string input, string delimiter, vector<double>& output);
bool readFile(string FileName, vector< vector<string> >& content, const int nCol=-1, int nRow=-1);
bool setParameters(map<string, string> row, struct KinCuts& cut, map<string, string>& parIni);
bool addParameters(string InputFile,  vector< struct KinCuts >& cutVector, vector< map<string, string> >&  parIniVector);


std::array<RooDataSet*, 4> createDataset(bool ispO, bool isMC, const char *methodName) {
  bool doSplot = false; bool doSideBands = false;
  if (strcmp(methodName, "sPlot") == 0) doSplot = true;
  if (strcmp(methodName, "sideBands") == 0) doSideBands = true;

  TChain *fChain = getInputTree(Form("inputFiles/input_%s_%s.txt", isMC ? "MC" : "data", ispO ? "pO" : "OO"), "O2rtdimuonall");

  Float_t fMass; fChain->SetBranchAddress("fMass", &fMass);
  Float_t fPt; fChain->SetBranchAddress("fPt", &fPt);
  Float_t fEta; fChain->SetBranchAddress("fEta", &fEta);
  Float_t fEta1; fChain->SetBranchAddress("fEta1", &fEta1);
  Float_t fEta2; fChain->SetBranchAddress("fEta2", &fEta2);
  Float_t fPhi; fChain->SetBranchAddress("fPhi", &fPhi);
  Int_t fSign; fChain->SetBranchAddress("fSign", &fSign);
  Float_t fChi2pca; fChain->SetBranchAddress("fChi2pca", &fChi2pca);
  Float_t fTauz; fChain->SetBranchAddress("fTauz", &fTauz);
  Float_t fTauzErr; fChain->SetBranchAddress("fTauzErr", &fTauzErr);
  Int_t fIsAmbig1; fChain->SetBranchAddress("fIsAmbig1", &fIsAmbig1);
  Int_t fIsAmbig2; fChain->SetBranchAddress("fIsAmbig2", &fIsAmbig2);
  Float_t fChi2MatchMCHMFT1; fChain->SetBranchAddress("fChi2MatchMCHMFT1", &fChi2MatchMCHMFT1);
  Float_t fChi2MatchMCHMFT2; fChain->SetBranchAddress("fChi2MatchMCHMFT2", &fChi2MatchMCHMFT2);
  // Do the bottom part only for MC
  auto fMcDecision = static_cast<uint32_t>(0);
  UShort_t fMcMask1;
  UShort_t fMcMask2;
  if (isMC) {
    fChain->SetBranchAddress("fMcDecision", &fMcDecision);
    fChain->SetBranchAddress("fMcMask1", &fMcMask1);
    fChain->SetBranchAddress("fMcMask2", &fMcMask2);
  }

  fChain->SetBranchStatus("*", 0);
  fChain->SetBranchStatus("fMass", 1);
  fChain->SetBranchStatus("fPt", 1);
  fChain->SetBranchStatus("fEta", 1);
  fChain->SetBranchStatus("fEta1", 1);
  fChain->SetBranchStatus("fEta2", 1);
  fChain->SetBranchStatus("fPhi", 1);
  fChain->SetBranchStatus("fSign", 1);
  fChain->SetBranchStatus("fChi2pca", 1);
  fChain->SetBranchStatus("fTauz", 1);
  fChain->SetBranchStatus("fTauzErr", 1);
  fChain->SetBranchStatus("fIsAmbig1", 1);
  fChain->SetBranchStatus("fIsAmbig2", 1);
  fChain->SetBranchStatus("fChi2MatchMCHMFT1", 1);
  fChain->SetBranchStatus("fChi2MatchMCHMFT2", 1);
  if (isMC) {
    fChain->SetBranchStatus("fMcDecision", 1);
    fChain->SetBranchStatus("fMcMask1", 1);
    fChain->SetBranchStatus("fMcMask2", 1);
  }

  RooRealVar* mass = new RooRealVar("mass", "Mass_{#mu^{+}#mu^{-}}", 2.4, 4.0, "GeV/c^{2}");
  RooRealVar* pt = new RooRealVar("pt", "p_{T, #mu^{+}#mu^{-}}", 0, 20, "GeV/c");
  RooRealVar* y = new RooRealVar("y", "y_{#mu^{+}#mu^{-}}", -5, -2);
  RooRealVar* sign = new RooRealVar("sign", "dimuon sign", -3, 3);
  RooRealVar* tauz = new RooRealVar("tauz", "#tau_{z,#mu^{+}#mu^{-}}", -0.07, 0.07);
  RooRealVar* chi2_1 = new RooRealVar("chi2_1", "#chi^{2}_{1, #mu^{+}#mu^{-}}", 0, 1000);
  RooRealVar* chi2_2 = new RooRealVar("chi2_2", "#chi^{2}_{2, #mu^{+}#mu^{-}}", 0, 1000);

  // Normal dataset
  RooArgSet* varSet = new RooArgSet(*mass, *pt, *y, *sign, *tauz, *chi2_1, *chi2_2);
  RooDataSet* data = new RooDataSet("data", "data for dimuon pairs", *varSet);
  RooRealVar* massPeak = new RooRealVar("massPeak", "Mass_{#mu^{+}#mu^{-}} peak", 2.0, 4.5, "GeV/c^{2}");
  RooRealVar* massSB = new RooRealVar("massSB", "Mass_{#mu^{+}#mu^{-}} sidebands", 2.0, 4.5, "GeV/c^{2}");

  // Sideband datasets contain the same variables needed for the
  // tauz/background study, but use massSB rather than mass.
  RooArgSet* sidebandVarSet = new RooArgSet(*massSB, *pt, *y, *sign, *tauz, *chi2_1, *chi2_2);
  RooDataSet* dataPeak = new RooDataSet("dataPeak", "peak region", *sidebandVarSet);
  RooDataSet* dataSidebandLow = new RooDataSet("dataSidebandLow", "low mass sideband", *sidebandVarSet);
  RooDataSet* dataSidebandHigh = new RooDataSet("dataSidebandHigh", "high mass sideband", *sidebandVarSet);

  // Investigate the matching distributions
  // Only in MC
  // TODO: initialise them only in MC
  TH1D *hMatchesTrueTrue = new TH1D("hMatchesTrueTrue", "True matches mass_{#mu^{+}#mu^{-}}", 100, 2.4, 4.0);
  TH1D *hMatchesTrueFake = new TH1D("hMatchesTrueFake", "True-fake or fake-true mass_{#mu^{+}#mu^{-}}", 100, 2.4, 4.0);
  TH1D *hMatchesFakeFake = new TH1D("hMatchesFakeFake", "Fake-fake mass_{#mu^{+}#mu^{-}}", 100, 2.4, 4.0);

  int n_entries = fChain->GetEntries();
  for (int nEv = 0; nEv < n_entries; nEv++) {
    if (nEv % 100000 == 0) cout << "processing evt " << nEv << "/" << n_entries << endl;
    fChain->GetEntry(nEv);

    float rap = TMath::Log((TMath::Sqrt(fMass * fMass + fPt * fPt * TMath::CosH(fEta) * TMath::CosH(fEta)) + fPt * TMath::SinH(fEta)) / (TMath::Sqrt(fMass * fMass + fPt * fPt)));

    if (fMass < 2.0 || fMass > 4.5) continue;
    if (fTauz < -0.07 || fTauz > 0.07) continue;
    if (fEta1 < -3.6 || fEta1 > -2.5) continue;
    if (fEta2 < -3.6 || fEta2 > -2.5) continue;
    if (fChi2MatchMCHMFT1 > 500 || fChi2MatchMCHMFT2 > 500) continue;
    if (fIsAmbig1 || fIsAmbig2) continue;

    // For MC we select only the MC truth signal
    // 4 --> Jpsi + psi2S
    // 5 --> Only Jpsi
    // We also plot the MFT-MCH matching distributions
    // McMask == 0 --> true global match
    if (isMC) {
      if (fMcMask1 == 0 && fMcMask2 == 0) hMatchesTrueTrue->Fill(fMass);
      if ((fMcMask1 == 0 && fMcMask2 != 0) || (fMcMask1 != 0 && fMcMask2 == 0)) hMatchesTrueFake->Fill(fMass);
      if (fMcMask1 != 0 && fMcMask2 != 0) hMatchesFakeFake->Fill(fMass);
      // if (fMcDecision == 0) continue;  
    }

    pt->setVal(fPt);
    y->setVal(rap);
    sign->setVal(fSign);
    tauz->setVal(fTauz);
    chi2_1->setVal(fChi2MatchMCHMFT1);
    chi2_2->setVal(fChi2MatchMCHMFT2);

    // TODO: remove the sidebands magic numbers...
    if (fMass >= 2.4 && fMass <= 4.0) { mass->setVal(fMass); data->add(*varSet); }
    if (fMass >= 2.9 && fMass <= 3.2) { massPeak->setVal(fMass); dataPeak->add(*sidebandVarSet); }
    else if (fMass >= 2.4 && fMass < 2.8) { massSB->setVal(fMass); dataSidebandLow->add(*sidebandVarSet); }
    else if (fMass > 3.4 && fMass <= 4.0) { massSB->setVal(fMass); dataSidebandHigh->add(*sidebandVarSet); }
  }

  if (doSideBands) {
    cout << endl;
    cout << "Sideband datasets created:" << endl;
    cout << "  Peak :  " << dataPeak->numEntries() << " events (2.9 < m < 3.2 GeV/c^2)" << endl;
    cout << "  Low :  " << dataSidebandLow->numEntries() << " events (2.4 < m < 2.8 GeV/c^2)" << endl;
    cout << "  High:  " << dataSidebandHigh->numEntries() << " events (3.4 < m < 4.0 GeV/c^2)" << endl;
  }
  cout << "Main dataset:" << endl;
  cout << "  Full mass region: " << data->numEntries()  << " events (2.4 < m < 4.0 GeV/c^2)" << endl;

  if (isMC) {
    TCanvas *cMatches = new TCanvas("cMatches", "MFT-MCH matches mass");
    cMatches->cd();
    gPad->SetLogy();
    hMatchesTrueTrue->GetYaxis()->SetRangeUser(1e0, 1e6);
    hMatchesTrueTrue->SetLineColor(kRed);
    hMatchesTrueFake->SetLineColor(kGreen);
    hMatchesFakeFake->SetLineColor(kBlue);
    hMatchesTrueTrue->Draw("hist");
    hMatchesTrueFake->Draw("same hist");
    hMatchesFakeFake->Draw("same hist");
    cMatches->SaveAs("output/mass_distribution_MFT_MCH_matches.pdf");
  }

  return {data, dataPeak, dataSidebandLow, dataSidebandHigh};
}

TChain* getInputTree(string inputTreeFile, string treeName){
  std::ifstream infile(inputTreeFile);
  std::string line;
  std::vector<std::string> rootFiles;
  
  while (std::getline(infile, line)) {
    if (!line.empty()) rootFiles.push_back(line);
  }
  
  if (rootFiles.empty()) {
    std::cerr << "No input ROOT files found in "<< inputTreeFile << std::endl;
  }
  
  TChain* mergedChain = new TChain(treeName.c_str());

  for (const auto& file : rootFiles) {
    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) {
      std::cerr << "Failed to open: " << file << std::endl;
      continue;
    }
    findAndAddTrees(f, treeName.c_str(), mergedChain, file);
    f->Close();
  }
  std::cout << "Total entries in merged chain: " << mergedChain->GetEntries() << std::endl;
  return mergedChain;
}

void findAndAddTrees(TDirectory* dir, const std::string& treeName, TChain* chain, const std::string& fileName) {
    TIter nextKey(dir->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)nextKey())) {
        TObject* obj = key->ReadObj();

        if (obj->InheritsFrom(TDirectory::Class())) {
            findAndAddTrees((TDirectory*)obj, treeName, chain, fileName); // Recurse into subdirectories
        }
        else if (obj->InheritsFrom(TTree::Class())) {
            TTree* tree = (TTree*)obj;
            std::string fullPath = std::string(dir->GetPath()) + "/" + tree->GetName();

            if (treeName.empty() || tree->GetName() == treeName) {
                std::string fullTreePath = fileName + "/" + fullPath.substr(fullPath.find(':') + 2);
                chain->Add(fullTreePath.c_str());
                std::cout << "Added tree: " << fullTreePath << std::endl;
            }
        }
    }
}



typedef struct StartEnd {
  int Start, End;
} StartEnd;

typedef struct MinMax {
  double Min, Max;
} MinMax;

typedef struct KinCuts {
  StartEnd cent;
  MinMax pt;
  MinMax rap;
  MinMax chi2;
} KinCuts;

bool isEqualKinCuts(struct KinCuts cutA, struct KinCuts cutB, bool ispO) {
  bool cond = true;
  if (!ispO) {
    cond = cond && (cutA.cent.Start    == cutB.cent.Start);
    cond = cond && (cutA.cent.End      == cutB.cent.End);
  }
  cond = cond && (cutA.pt.Min    == cutB.pt.Min);
  cond = cond && (cutA.pt.Max      == cutB.pt.Max);

  cond = cond && (cutA.rap.Min    == cutB.rap.Min);
  cond = cond && (cutA.rap.Max      == cutB.rap.Max);
  
  cond = cond && (cutA.chi2.Min    == cutB.chi2.Min);
  cond = cond && (cutA.chi2.Max      == cutB.chi2.Max);

  return cond;
    
}

std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        result.push_back(token);
    }
    return result;
}

std::string removeSurroundingBrackets(const std::string& str) {
    if (str.size() >= 2 && str.front() == '[' && str.back() == ']') {
        return str.substr(1, str.size() - 2);
    }
    return str;
}

bool addParameters(string InputFile,  vector< struct KinCuts >& cutVector, vector< map<string, string> >&  parIniVector)
{
  vector< map<string, string> >  data;
  if(!parseFile(InputFile, data)) { return false; }
  
    for(vector< map<string, string> >::iterator row=data.begin(); row!=data.end(); ++row) {
      struct KinCuts cut; map<string, string> parIni;
      if(!setParameters(*row, cut, parIni)) { return false; }
      cutVector.push_back(cut);  parIniVector.push_back(parIni);
    }
  return true;
};

bool setParameters(map<string, string> row, struct KinCuts& cut, map<string, string>& parIni)
{
  cut.pt.Min = 0.;
  cut.pt.Max = 150.;

  cut.rap.Min = -10.;
  cut.rap.Max = 10.;

  cut.chi2.Min = 0.;
  cut.chi2.Max = 150.;

  cut.cent.Start = 0;
  cut.cent.End = 100;
  
  // set parameters from file
  for(map<string, string>::iterator col=row.begin(); col!=row.end(); ++col) {

    string label = col->first;
    if (label=="pt") {
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'pt' has invalid value: " << col->second << endl; return false;
      }  
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'pt' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.pt.Min = v.at(0); 
      cut.pt.Max = v.at(1);
    } 
    else if (label=="y") {
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'y' has invalid value: " << col->second << endl; return false;
      }  
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'y' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.rap.Min = v.at(0); 
      cut.rap.Max = v.at(1);
    } 
    else if (label=="chi2"){
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'chi2' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'chi2' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.chi2.Min = v.at(0); 
      cut.chi2.Max = v.at(1);
    }
    else if (label=="cent"){
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'cent' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v;
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'cent' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }
      cut.cent.Start = (int) (v.at(0));
      cut.cent.End = (int) (v.at(1));
    }
    else if (label == "doIterativeFit") {
      string val = col->second;
      // convert to 0/1 (if not already)
      if (val == "true") val = "1";
      else if (val == "false") val = "0";
      parIni[col->first] = val;
    }
    else if (label.find("model") != std::string::npos || label == "fitStat") {
      if (col->second=="") {
        cout << "[ERROR] Input column "<<label<<" has empty value" << endl; return false;
      }
      parIni[col->first] = col->second;
    }
    else {
      if (col->second != "") {
	string value = col->second;
	// check that initial parameters format is correct: [ num, num, num ]
	if ((value.find("[")==std::string::npos)||(value.find("]")==std::string::npos)) {
	  // Special cases like parameter constrains could be set here but for now, let's keep it simple
	  cout << "[ERROR] Either ']' or '[' are missing in the initial parameter values for " <<label<<"!"<< endl; return false;
	} else {
	  value.erase(value.find("["), string("[").length());
	  value.erase(value.find("]"), string("]").length());
	}
	std::vector<double> v; 
	if(!parseString(value, ",", v)) { return false; }
        if (v.size()>3 || v.size()<1) {
          cout << "[ERROR] Initial parameter " << col->first << " has incorrect number of values, it has: " << v.size() << endl; return false;
        }
	// everything seems alright, then proceed to save the values
	if (v.size()==1){
	  // if only one value is given i.e. [ num ], consider it a constant value
	  parIni[col->first] = Form("[%.6f]", v.at(0));
	}
	else
	  {
	    parIni[col->first] = col->second;
	  }
      } else {
        parIni[col->first] = "";
      }
    }
  }
  return true;
};

bool parseString(string input, string delimiter, vector<double>& output)
{
  // remove spaces from input string 
  input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
  // proceed to parse input string
  char *end;
  while(input!="") {
    double d = strtod(input.c_str(), &end);
    if (end != input) {
      output.push_back(d);
    } else {
      cout << "[ERROR] The conversion from string to double failed!"; return false;
    }
    input = end; 
    if(input.find(delimiter.c_str())!= std::string::npos){ input.erase(input.find(delimiter.c_str()), delimiter.length()); }
  }
  return true;
};

bool parseFile(string FileName, vector< map<string, string> >& data) {
  vector< vector<string> > content, tmp; 
  if(!readFile(FileName, tmp, -1, 1)){ return false; }
  vector<string> header = tmp.at(0);
  if (header.size()==0) { cout << "[ERROR] The header is null!" << endl; return false; }
  if(!readFile(FileName, content, header.size())){ return false; }
  for(vector<string>::iterator rHeader=header.begin(); rHeader!=header.end(); ++rHeader) {
    if (*rHeader=="") { cout << "[ERROR] A column has no label!" << endl; return false; }
  }

  for(vector< vector<string> >::iterator row=content.begin()+1; row!=content.end(); ++row) {
    map<string, string> col;
    for (unsigned int i=0; i<header.size(); i++) {
      if (i<row->size()) {
	col[header.at(i)] = row->at(i);
      } else {
	col[header.at(i)] = "";
      }
    }
    data.push_back(col);
  }

  return true;
};

bool readFile(string FileName, vector< vector<string> >& content, const int nCol, int nRow) {
  if (nCol==0 || nRow==0) { 
    cout << "[WARNING] Ignoring content of File: " << FileName << endl; return true; 
  }
  if (nRow!=1) { cout << "[INFO] Reading file: " << FileName << endl; }
  ifstream myfile(FileName.c_str());
  if (myfile.is_open()){ 
    string line;
    while ( getline(myfile, line) ){
      if (nRow==0) break; else {nRow=nRow-1;}
      stringstream row(line);
      vector<string> cols; int i=0;
      while (true){
	string col; getline(row, col, ';');
	if ( (nCol>=0) ? (i>=nCol) : (col=="") ){ break; }
	cols.push_back(col);
	i++;
      }
      content.push_back(cols);
    }
  } else {
    cout << "[ERROR] File: " << FileName << " was not found!" << endl; return false;
  }
  return true;
};

/*
string cuts(float minPt, float maxPt, float minRap, float maxRap, float maxChi2, bool isOS) {
  return Form("pt > %f && pt < %f && y > %f && y< %f && chi2_1< %f && chi2_2< %f && %s", minPt, maxPt, minRap, maxRap, minChi2, maxChi2, isOS?"sign == 0":"sign != 0");
}

string label(float minPt, float maxPt, float minRap, float maxRap, float maxChi2, bool isOS) {
  return Form("_pt_%d_%d_y_%d_%d_chi2_%d_%s", (int) minPt, (int) maxPt, (int) (minRap*-10), (int) (maxRap*-10), (int) maxChi2, isOS?"OS":"SS");
}
*/

