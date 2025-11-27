// Author: Batoul Diab
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <TLatex.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TKey.h>

#include "SignalExtraction.C"

// TODO: move to some header
struct PtBinningScanResults {
  std::vector<int> vPtMinValues;
  std::vector<int> vPtMaxValues;
  std::vector<double> vJpsiPeaks;
  std::vector<double> vJpsiPeaksErrors;
  std::vector<double> vJpsiWidths;
  std::vector<double> vJpsiWidthsErrors;
  std::vector<double> vNJpsi;
  std::vector<double> vNJpsiErrors;
  std::vector<double> vNBkg;
  std::vector<double> vNBkgErrors;
};

void signalExtraction(PtBinningScanResults &ptBinningScanResult, bool ispO=true, bool isMC =false, const char *caseName = "nominal", bool remakeDS =false, bool fitMass=true, bool fitTauz=false);
void plotResult(const char *caseName = "nominal", string axisName = "pt", int incMinCent=0, int incMaxCent=100, float incMinPt=0., float incMaxPt=50., float incMinRap=-3.5, float incMaxRap=-2.5, float incMinChi2=0, float incMaxChi2=50);

void StudyPtBinning(bool ispO=false, bool isMC=false, const char *caseName = "nominal", bool remakeDS = false, bool fitMass1D=false, bool fitTauz1D=false, bool fit2D=true, bool plotResults = false) {
  gSystem->Load("RooExtCBShape.cxx+");
  //cout<<"plotResults = "<<plotResults<<endl;
  PtBinningScanResults ptBinningScanResult;
  if (fitMass1D)
    signalExtraction(ptBinningScanResult, ispO, isMC, caseName, remakeDS, true, false);
  if (fitTauz1D)
    signalExtraction(ptBinningScanResult, ispO, isMC, caseName, remakeDS, false, true);
  if (fit2D)
    signalExtraction(ptBinningScanResult, ispO, isMC, caseName, remakeDS, true, true);
  

  if (plotResults) {
    int minCent = 0;
    int maxCent = 100;
    float minPt = 0.;
    float maxPt = 20.;
    float minRap = -3.5;
    float maxRap = -2.6;
    float minChi2 = 0;
    float maxChi2 = 50;
    string axisName = "pt";
    plotResult(caseName, axisName.c_str(), minCent, maxCent, minPt, maxPt, minRap, maxRap, minChi2, maxChi2);
  }

  // Assuming all bins tested have the global maximum value
  // e.g. bins [0, 1], [1, 3], [3, 5] have global ptMax = 5
  // This can be used to see how many bins we take in one comparison
  // As all the bins will be plotted in one bin of the template histogram
  int ptGlobalMax = 15;
  int nPtBins = 0;

  std::vector<int> vPtMaxValues = ptBinningScanResult.vPtMaxValues;
  for (int iPtBin = 0; iPtBin < vPtMaxValues.size(); iPtBin++) {
    int ptBinMax = vPtMaxValues[iPtBin];
    if (ptBinMax == ptGlobalMax) { nPtBins++; }
  }

  // TODO: this can be done with really just one canvas and histogram
  TCanvas *cJpsiPeaks = new TCanvas("cJpsiPeaks", "cJpsiPeaks", 800, 600);
  TCanvas *cJpsiWidths = new TCanvas("cJpsiWidths", "cJpsiWidths", 800, 600);
  TCanvas *c_nJpsi = new TCanvas("c_nJpsi", "c_nJpsi", 800, 600);
  TCanvas *c_s_over_b = new TCanvas("c_s_over_b", "c_s_over_b", 800, 600);
  TH1D *hTemplateJpsiPeaks = new TH1D("hTemplateJpsiPeaks", "hTemplateJpsiPeaks", nPtBins, 0, nPtBins);
  TH1D *hTemplateJpsiWidths = new TH1D("hTemplateJpsiWidths", "hTemplateJpsiWidths", nPtBins, 0, nPtBins);
  TH1D *hTemplate_nJpsi = new TH1D("hTemplate_nJpsi", "hTemplate_nJpsi", nPtBins, 0, nPtBins);
  TH1D *hTemplate_s_over_b = new TH1D("hTemplate_s_over_b", "hTemplate_s_over_b", nPtBins, 0, nPtBins);
  TH1D *hJpsiPeaks = new TH1D("hJpsiPeaks", "hJpsiPeaks", nPtBins, 0, nPtBins);
  TH1D *hJpsiWidths = new TH1D("hJpsiWidths", "hJpsiWidths", nPtBins, 0, nPtBins);
  TH1D *h_nJpsi = new TH1D("h_nJpsi", "h_nJpsi", nPtBins, 0, nPtBins);
  TH1D *h_s_over_b = new TH1D("h_s_over_b", "h_s_over_b", nPtBins, 0, nPtBins);

  // Add matching chi2 cut to title
  if (ispO) { hTemplateJpsiPeaks->SetTitle(Form("J/psi mass for pO")); }
  else { hTemplateJpsiPeaks->SetTitle(Form("J/psi mass for OO")); }
  hTemplateJpsiPeaks->GetXaxis()->SetTitle("pT bin collection");
  hTemplateJpsiPeaks->GetYaxis()->SetRangeUser(2.5, 4.0);
  if (ispO) { hTemplateJpsiWidths->SetTitle(Form("J/psi sigma for pO")); }
  else { hTemplateJpsiWidths->SetTitle(Form("J/psi sigma for OO")); }
  hTemplateJpsiWidths->GetXaxis()->SetTitle("pT bin collection");
  hTemplateJpsiWidths->GetYaxis()->SetRangeUser(0, 0.1);
  if (ispO) { 
    hTemplate_nJpsi->SetTitle(Form("N(J/psi) for pO"));
    hTemplate_nJpsi->GetYaxis()->SetRangeUser(0, 4000); 
  }
  else { 
    hTemplate_nJpsi->SetTitle(Form("N(J/psi) for OO"));
    hTemplate_nJpsi->GetYaxis()->SetRangeUser(0, 20000); 
  }
  hTemplate_nJpsi->GetXaxis()->SetTitle("pT bin collection");
  if (ispO) { hTemplate_s_over_b->SetTitle(Form("signal/background for pO")); }
  else { hTemplate_s_over_b->SetTitle(Form("signal/background for OO")); }
  hTemplate_s_over_b->GetXaxis()->SetTitle("pT bin collection");
  hTemplate_s_over_b->GetYaxis()->SetRangeUser(0, 1.5);

  cJpsiPeaks->cd(); hTemplateJpsiPeaks->Draw("HIST");
  cJpsiWidths->cd(); hTemplateJpsiWidths->Draw("HIST");
  c_nJpsi->cd(); hTemplate_nJpsi->Draw("HIST");
  c_s_over_b->cd(); hTemplate_s_over_b->Draw("HIST");

  std::vector<int> vLineColours = {1, 2, 6, 7, 8, 9};
  std::vector<std::string> vBinLabels;
  std::vector<TH1D*> vJpsiPeaksClones;
  std::vector<TH1D*> vJpsiWidthsClones;
  std::vector<TH1D*> v_nJpsiClones;
  std::vector<TH1D*> v_s_over_bClones;

  int iPtBin = 0;
  int maxColourIndex = 0;

  for (int iTemplateBin = 0; iTemplateBin < nPtBins; iTemplateBin++) {

    vBinLabels.clear();
    int colourIndex = 0; // just for consistent drawing
    
    // Keep going until the global max pT is found
    while (1) {
      int ptBinMin = ptBinningScanResult.vPtMinValues[iPtBin];
      int ptBinMax = ptBinningScanResult.vPtMaxValues[iPtBin];
      std::string binLabel = Form("[%i, %i]", ptBinMin, ptBinMax);
      vBinLabels.push_back(binLabel);
      double jpsiPeak = ptBinningScanResult.vJpsiPeaks[iPtBin];
      double jpsiPeakError = ptBinningScanResult.vJpsiPeaksErrors[iPtBin];
      double jpsiWidth = ptBinningScanResult.vJpsiWidths[iPtBin];
      double jpsiWidthError = ptBinningScanResult.vJpsiWidthsErrors[iPtBin];
      double nJpsi = ptBinningScanResult.vNJpsi[iPtBin];
      double nJpsiError = ptBinningScanResult.vNJpsiErrors[iPtBin];
      double nBkg = ptBinningScanResult.vNBkg[iPtBin];
      double nBkgError = ptBinningScanResult.vNBkgErrors[iPtBin];
      double s_over_b = nJpsi / nBkg;
      double s_over_bError = s_over_b * std::sqrt((nJpsiError/nJpsi)*(nJpsiError/nJpsi) + (nBkgError/nBkg)*(nBkgError/nBkg));
      // clone the histograms for every ptBin considered (only necessary once)
      TH1D *hJpsiPeaksClone = (TH1D*)hJpsiPeaks->Clone(Form("hJpsiPeaks_%i_%i", iTemplateBin, iPtBin));
      TH1D *hJpsiWidthsClone = (TH1D*)hJpsiWidths->Clone(Form("hJpsiWidths_%i_%i", iTemplateBin, iPtBin));
      TH1D *h_nJpsiClone = (TH1D*)h_nJpsi->Clone(Form("h_nJpsi_%i_%i", iTemplateBin, iPtBin));
      TH1D *h_s_over_bClone = (TH1D*)h_s_over_b->Clone(Form("h_s_over_b_%i_%i", iTemplateBin, iPtBin));
      hJpsiPeaksClone->SetLineColor(vLineColours[colourIndex]);
      hJpsiWidthsClone->SetLineColor(vLineColours[colourIndex]);
      h_nJpsiClone->SetLineColor(vLineColours[colourIndex]);
      h_s_over_bClone->SetLineColor(vLineColours[colourIndex]);
      std::cout << "setting Jpsi peak = " << jpsiPeak << std::endl;
      std::cout << "setting Jpsi width = " << jpsiWidth << std::endl;
      std::cout << "setting nJpsi = " << nJpsi << std::endl;
      std::cout << "setting s/b = " << s_over_b << std::endl;
      hJpsiPeaksClone->SetBinContent(iTemplateBin + 1, jpsiPeak);
      hJpsiPeaksClone->SetBinError(iTemplateBin + 1, jpsiPeakError);
      hJpsiWidthsClone->SetBinContent(iTemplateBin + 1, jpsiWidth);
      hJpsiWidthsClone->SetBinError(iTemplateBin + 1, jpsiWidthError);
      h_nJpsiClone->SetBinContent(iTemplateBin + 1, nJpsi);
      h_nJpsiClone->SetBinError(iTemplateBin + 1, nJpsiError);
      h_s_over_bClone->SetBinContent(iTemplateBin + 1, s_over_b);
      h_s_over_bClone->SetBinError(iTemplateBin + 1, s_over_bError);
      vJpsiPeaksClones.push_back(hJpsiPeaksClone);
      vJpsiWidthsClones.push_back(hJpsiWidthsClone);
      v_nJpsiClones.push_back(h_nJpsiClone);
      v_s_over_bClones.push_back(h_s_over_bClone);
      if (ptBinMax == ptGlobalMax) {
        std::string mergedBinLabel = "{";
        for (size_t i = 0; i < vBinLabels.size(); i++) {
          if (i > 0) mergedBinLabel += ", ";
          mergedBinLabel += vBinLabels[i];
        }
        mergedBinLabel += "}";
        hTemplateJpsiPeaks->GetXaxis()->SetBinLabel(iTemplateBin + 1, Form("%s", mergedBinLabel.c_str()));
        hTemplateJpsiWidths->GetXaxis()->SetBinLabel(iTemplateBin + 1, Form("%s", mergedBinLabel.c_str()));
        hTemplate_nJpsi->GetXaxis()->SetBinLabel(iTemplateBin + 1, Form("%s", mergedBinLabel.c_str()));
        hTemplate_s_over_b->GetXaxis()->SetBinLabel(iTemplateBin + 1, Form("%s", mergedBinLabel.c_str()));
        std::cout << "--- break ---" << std::endl;
        iPtBin++;
        if (colourIndex + 1 > maxColourIndex) { maxColourIndex = colourIndex + 1; }
        break; 
      }
      iPtBin++;
      colourIndex++;
    }
    std::cout << ">> go to the next template bin" << std::endl;
  }

  TLegend *leg = new TLegend(0.70, 0.70, 0.90, 0.90);
  for (int i = 0; (i < vLineColours.size() && i < maxColourIndex); i++) {
      std::string entryLabel = Form("bin %d", i+1);
      // create a dummy object just to show the colour
      TH1D *hDummy = new TH1D(Form("hDummy_%d", i), "", 1, 0, 1);
      hDummy->SetLineColor(vLineColours[i]);
      hDummy->SetLineWidth(3);
      leg->AddEntry(hDummy, entryLabel.c_str(), "l");
  }

  cJpsiPeaks->cd(); for (auto h : vJpsiPeaksClones) h->Draw("APE SAME"); leg->Draw();
  cJpsiWidths->cd(); for (auto h : vJpsiWidthsClones) h->Draw("APE SAME"); leg->Draw();
  c_nJpsi->cd(); for (auto h : v_nJpsiClones) h->Draw("APE SAME"); leg->Draw();
  c_s_over_b->cd(); for (auto h : v_s_over_bClones) h->Draw("APE SAME"); leg->Draw();

  // TODO: add pT binning here too in name
  std::string outputChi2MatchingName;
  if (ispO) { outputChi2MatchingName = "output/output_pTBinningScan_pO.pdf"; }
  else { outputChi2MatchingName = "output/output_pTBinningScan_OO.pdf"; }
  TCanvas *cdummy = new TCanvas("cdummy", "cdummy", 800, 600);
  cdummy->SaveAs(Form("%s(", outputChi2MatchingName.c_str()));
  cJpsiPeaks->SaveAs(outputChi2MatchingName.c_str());
  cJpsiWidths->SaveAs(outputChi2MatchingName.c_str());
  c_nJpsi->SaveAs(outputChi2MatchingName.c_str());
  c_s_over_b->SaveAs(outputChi2MatchingName.c_str());
  cdummy->SaveAs(Form("%s)", outputChi2MatchingName.c_str()));
}

void signalExtraction(PtBinningScanResults &ptBinningScanResult, bool ispO, bool isMC, const char *caseName, bool remakeDS, bool fitMass, bool fitTauz) {
  
  //do the fits or at least some of them
  vector< struct KinCuts >       cutVector;
  vector< map<string, string> >  parIniVector;
  vector< map<string, double> >  allResults;
  
  if (!addParameters(Form("inputFiles/initialPars_%s_%s_%s.txt", fitMass?"mass":"tauz", ispO?"pO":"OO", caseName), cutVector, parIniVector)) { return; } //if 1D fit on mass or tauz it reads the corresponding file, if it's 2D it reads the mass but the other variables get added from the default values and then they are fixed
  
  string outputName = Form("output/output_fit%s%s_%s%s_%s.root", fitMass?"Mass":"", fitTauz?"Tauz":"", ispO?"pO":"OO", isMC?"_MC":"", caseName);
  if (gSystem->AccessPathName("output")) gSystem->mkdir("output", true);
  std::ifstream fileCheck(outputName.c_str());
  bool fileExists = fileCheck.good();
  fileCheck.close();
  TFile* fSave;
  
  if (fileExists) 
    fSave = TFile::Open(outputName.c_str(),"UPDATE");
  else
    fSave = TFile::Open(outputName.c_str(),"RECREATE");
  if (!fSave || fSave->IsZombie()) {
    cout<<"[ERROR] Problem with the output file "<<outputName<<endl;
    return;
  }
  
  string dsFileName = Form("output/dataset_%s%s.root", ispO?"pO":"OO", isMC?"_MC":"");
  fileExists = std::filesystem::exists(dsFileName);
  bool fileIsGood = false;
  
  TFile *dsInputFile = NULL;
  RooDataSet* data = NULL;
  
  if (fileExists) {
    cout <<"the DS file exists"<<endl;
    dsInputFile = TFile::Open(dsFileName.c_str(), "READ");
    if (dsInputFile && !dsInputFile->IsZombie() && dsInputFile->IsOpen()) {
      fileIsGood = true;
      cout <<"the DS file is good"<<endl;
      data = (RooDataSet*) dsInputFile->Get("data");
    }
  }
  
  if (remakeDS || !fileExists || !fileIsGood) {
    gSystem->mkdir("output");
    dsInputFile = TFile::Open(dsFileName.c_str(), "RECREATE");
    data = createDataset(ispO, isMC);
    dsInputFile->cd();
    cout <<"saving the dataset"<<endl;
    data->Write("data");
  }

  for (uint j = 0; j < cutVector.size(); j++) {
    if (parIniVector[j]["fitStat"] != "todo") {
      cout<<"[INFO] This fit does not have the 'todo' status, so it will be skipped!"<<endl;
      continue;
    }
    
    std::string rangeLabel = Form("pt_%d_%d_rap_%d_%d_cent_%d_%d_chi2_%d_%d",
				  (int) cutVector[j].pt.Min,
				  (int) cutVector[j].pt.Max,
				  (int) cutVector[j].rap.Min,
				  (int) cutVector[j].rap.Max,			    
				  (int) cutVector[j].cent.Start,
				  (int) cutVector[j].cent.End,
				  (int) cutVector[j].chi2.Min,
				  (int) cutVector[j].chi2.Max);


    string dsCuts = Form("pt > %f && pt < %f && y > %f && y < %f && chi2_1 > %f && chi2_1 < %f && chi2_2 > %f && chi2_2 < %f && sign == 0",
			 cutVector[j].pt.Min,
			 cutVector[j].pt.Max,
			 -1*cutVector[j].rap.Max,
			 -1*cutVector[j].rap.Min,
			 cutVector[j].chi2.Min,
			 cutVector[j].chi2.Max,
			 cutVector[j].chi2.Min,
			 cutVector[j].chi2.Max);
    cout<<"cutting on "<<dsCuts<<endl;
    
    RooWorkspace* ws = new RooWorkspace(Form("ws_fit%s%s_%s", fitMass?"Mass":"", fitTauz?"Tauz":"", rangeLabel.c_str()));
    map<string, double> resultsFit;

    cout<<"[INFO] number of entries in the dataset before reduction = "<<data->numEntries()<<endl;
    RooDataSet* cutDataset = (RooDataSet*) data->reduce(Form("%s", dsCuts.c_str()));
    cout<<"[INFO] number of entries in the dataset after reduction = "<<cutDataset->numEntries()<<endl;
    ws->import(*cutDataset);

    if (fitTauz && !fitMass && !isMC) { //import the sPlot datasets
      string sPlotFileName = Form("output/output_fitMass_%s_%s.root", ispO?"pO":"OO", caseName);
      TFile* sPlotFile = TFile::Open(sPlotFileName.c_str());
      RooDataSet* sPlotDs = (RooDataSet*) sPlotFile->Get(Form("sPlotDS_%s", rangeLabel.c_str()));
      sPlotDs = (RooDataSet*) sPlotDs->reduce(Form("%s", dsCuts.c_str()));
      RooDataSet* sPlotDsS = (RooDataSet*) sPlotDs->reduce("fJpsi_mass_sw > 0 && fJpsi_mass_sw < 50");
      RooDataSet* sPlotDsB = (RooDataSet*) sPlotDs->reduce("fBkg_mass_sw > 0 && fBkg_mass_sw < 50");
      
      cout<<"[INFO] found and reduced the sPlot"<<endl;
      
      const RooArgSet* varSet = sPlotDs->get();
      RooDataSet* sPlotDsSig = new RooDataSet("sPlotDsSig", "Signal-weighted dataset", sPlotDsS, *varSet, 0, "fJpsi_mass_sw");
      RooDataSet* sPlotDsBkg = new RooDataSet("sPlotDsBkg", "Background-weighted dataset", sPlotDsB, *varSet, 0, "fBkg_mass_sw");
      ws->import(*sPlotDsSig);
      ws->import(*sPlotDsBkg);
    }
    
    resultsFit.clear();
    resultsFit = SignalExtraction1Fit(parIniVector[j], ws, caseName, rangeLabel, ispO, isMC, fitMass, fitTauz, cutVector[j]);
    resultsFit["centMin"] = cutVector[j].cent.Start;
    resultsFit["centMax"] = cutVector[j].cent.End;
    resultsFit["ptMin"] = cutVector[j].pt.Min;
    resultsFit["ptMax"] = cutVector[j].pt.Max;
    resultsFit["rapMin"] = cutVector[j].rap.Min;
    resultsFit["rapMax"] = cutVector[j].rap.Max;
    resultsFit["chi2Min"] = cutVector[j].chi2.Min;
    resultsFit["chi2Max"] = cutVector[j].chi2.Max;

    // Study certain parameters as a function of matching chi2
    ptBinningScanResult.vPtMinValues.push_back(resultsFit["ptMin"]);
    ptBinningScanResult.vPtMaxValues.push_back(resultsFit["ptMax"]);
    ptBinningScanResult.vJpsiPeaks.push_back(resultsFit["mean_mass"]);
    ptBinningScanResult.vJpsiPeaksErrors.push_back(resultsFit["mean_mass_err"]);
    ptBinningScanResult.vJpsiWidths.push_back(resultsFit["sigma_mass"]);
    ptBinningScanResult.vJpsiWidthsErrors.push_back(resultsFit["sigma_mass_err"]);
    ptBinningScanResult.vNJpsi.push_back(resultsFit["fJpsi_mass"]);
    ptBinningScanResult.vNJpsiErrors.push_back(resultsFit["fJpsi_mass_err"]);
    ptBinningScanResult.vNBkg.push_back(resultsFit["fBkg_mass"]);
    ptBinningScanResult.vNBkgErrors.push_back(resultsFit["fBkg_mass_err"]);
    
    allResults.push_back(resultsFit);

    //transform the results into trees
    fSave->cd();
    TTree* resTree = new TTree(Form("tree_%s_new", rangeLabel.c_str()), "Tree of results");
    
    std::map<std::string, double*> branches;
    
    // Create branches for each map entry
    for (const auto& entry : resultsFit) {
        branches[entry.first] = new double;
        resTree->Branch(entry.first.c_str(), branches[entry.first]);
    }

    // Fill the resTree with data from the map
    for (const auto& entry : resultsFit) {
        *(branches[entry.first]) = entry.second;
	cout <<"filling tree with "<<entry.first<<" = "<<entry.second<<endl;
    }
    resTree->Fill(); //just fill one entry in the tree (make a different tree for each fit to allow updating)
    
    //check if the ws exists, if yes delete it and save it again
    TObject* obj = fSave->Get(Form("tree_%s", rangeLabel.c_str()));
    if (obj) {
      fSave->Delete(Form("tree_%s;*",rangeLabel.c_str()));
    }
    if (fitMass && !fitTauz && fSave->Get(Form("sPlotDS_%s", rangeLabel.c_str()))) {
      fSave->Delete(Form("sPlotDS_%s;*",rangeLabel.c_str()));
    }
    
    fSave->cd();
    if (fitMass && !fitTauz) ws->data("data_sPlot")->Write(Form("sPlotDS_%s", rangeLabel.c_str()));
    resTree->Write(Form("tree_%s", rangeLabel.c_str()));//, TObject::kOverwrite);
  }
  fSave->Close();
}


void plotResult(const char *caseName, string axisName, int incMinCent, int incMaxCent, float incMinPt, float incMaxPt, float incMinRap, float incMaxRap, float incMinChi2, float incMaxChi2){
  cout <<"[INFO] The plotting function is yet to be done"<<endl;
  gStyle->SetOptStat(0);
  
  map<string, double> resultsFit;
  double *binEdges = new double[100];
  int nbins;
  double *resVal_fb = new double[100];
  double *resVal_pr = new double[100];
  double *resVal_npr = new double[100];
  double *resErr_fb = new double[100];
  double *resErr_pr = new double[100];
  double *resErr_npr = new double[100];
  
  string fileName = Form("output/output_%s.root",caseName);
  TFile* fHist = TFile::Open(fileName.c_str(),"READ");
  if (!fHist || fHist->IsZombie()) {
    cout<<"[ERROR] Problem with the result file that contains all the histograms"<<fileName<<endl;
    return;
  }
    
  cout<<"[INFO] Opening the result input file"<<endl;
  TIter next(fHist->GetListOfKeys());
  TKey *key; int iBin = 0; double lastEdge = 0.;
  while ((key = (TKey*)next())) {
    // Get the class name of the object
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl) continue;
    if (cl->InheritsFrom("TTree")) {
      TTree *resTree = (TTree*) key->ReadObj();//fHist->Get(key->GetName());
      
      cout<<"[INFO] Reading the input tree"<<resTree->GetName()<<endl;

      double centMin = 0.; resTree->SetBranchAddress("centMin", &centMin);
      double centMax = 0.; resTree->SetBranchAddress("centMax", &centMax);
      double ptMin = 0.; resTree->SetBranchAddress("ptMin", &ptMin);
      double ptMax = 0.; resTree->SetBranchAddress("ptMax", &ptMax);
      double rapMin = 0.; resTree->SetBranchAddress("rapMin", &rapMin);
      double rapMax = 0.; resTree->SetBranchAddress("rapMax", &rapMax);
      double chi2Min = 0.; resTree->SetBranchAddress("chi2Min", &chi2Min);
      double chi2Max = 0.; resTree->SetBranchAddress("chi2Max", &chi2Max);
      
      double fb_jpsi = 0.; resTree->SetBranchAddress("b_jpsi_tauzMass", &fb_jpsi);
      double fb_jpsi_err = 0.; resTree->SetBranchAddress("b_jpsi_tauzMass_err", &fb_jpsi_err);
      double N_jpsi = 0.; resTree->SetBranchAddress("fJpsi_tauzMass_jpsi", &N_jpsi);
      double N_jpsi_err = 0.; resTree->SetBranchAddress("fJpsi_tauzMass_err", &N_jpsi_err);
      /*
      double fb_psi2s = 0.; resTree->SetBranchAddress("b_psi2s_tauzMass", &fb_psi2s);
      double fb_psi2s_err = 0.; resTree->SetBranchAddress("b_psi2s_tauzMass_err", &fb_psi2s_err);
      double N_psi2s = 0.; resTree->SetBranchAddress("fPsi2s_tauzMass", &N_Psi2s);
      double N_psi2s_err = 0.; resTree->SetBranchAddress("fPsi2s_tauzMass_err", &N_Psi2s_err);
      */
      
      resTree->GetEntry(0);
      
      if (axisName.find("centrality")!=std::string::npos) { //if plotting as function of centrality, the centrality of the result needs to be in the range but the pt need to be the exact edges
	
        if (centMin < incMinCent || centMax > incMaxCent)
          continue;
        if (fabs(ptMin - incMinPt)>0.00001 || fabs(ptMax - incMaxPt)>0.00001)
          continue;
        if (fabs(rapMin - incMinRap)>0.00001 || fabs(rapMax - incMaxRap)>0.00001)
          continue;
	if (fabs(chi2Min - incMinChi2)>0.00001 || fabs(chi2Max - incMaxChi2)>0.00001)
	  continue;
      
	cout<<"[INFO] This tree passed the bin selection, getting the results for pt ["<<ptMin<<"-"<<ptMax<<"], ptAssoc ["<<rapMin<<"-"<<rapMax<<"], cent ["<<centMin<<"-"<<centMax<<"]"<<endl;
	binEdges[iBin] = centMin;
	lastEdge =centMax;
      }//end of if centrality
      else if (axisName.find("pt")!=std::string::npos) {
        if (ptMin < incMinPt || ptMax > incMaxPt)
          continue;
        if (fabs(centMin - incMinCent)>0.00001 || fabs(centMax - incMaxCent)>0.00001)
          continue;
        if (fabs(rapMin - incMinRap)>0.00001 || fabs(rapMax - incMaxRap)>0.00001)
          continue;
	if (fabs(chi2Min - incMinChi2)>0.00001 || fabs(chi2Max - incMaxChi2)>0.00001)
	  continue;
	binEdges[iBin] = ptMin;
	lastEdge = ptMax;
      }//end of if pt
      else if (axisName.find("rap")!=std::string::npos) {
	if (rapMin < incMinRap || rapMax > incMaxRap)
          continue;
        if (fabs(centMin - incMinCent)>0.00001 || fabs(centMax - incMaxCent)>0.00001)
          continue;
        if (fabs(ptMin - incMinPt)>0.00001 || fabs(ptMax - incMaxPt)>0.00001)
          continue;
	if (fabs(chi2Min - incMinChi2)>0.00001 || fabs(chi2Max - incMaxChi2)>0.00001)
	  continue;
	binEdges[iBin] = rapMin;
	lastEdge = rapMax;
      }//end of rap
      else if (axisName.find("chi2")!=std::string::npos) {
	if (chi2Min < incMinChi2 || chi2Max > incMaxChi2)
          continue;
        if (fabs(centMin - incMinCent)>0.00001 || fabs(centMax - incMaxCent)>0.00001)
          continue;
        if (fabs(rapMin - incMinRap)>0.00001 || fabs(rapMax - incMaxRap)>0.00001)
          continue;
	if (fabs(ptMin - incMinPt)>0.00001 || fabs(ptMax - incMaxPt)>0.00001)
	  continue;
	binEdges[iBin] = chi2Min;
	lastEdge = chi2Max;
      }//end of chi2
	
      resVal_pr[iBin] = N_jpsi*(1-fb_jpsi);
      resErr_pr[iBin] = resVal_pr[iBin]*sqrt(pow((fb_jpsi_err/fb_jpsi),2)+pow((N_jpsi_err/N_jpsi),2)); //correlation needs to be taken into account
      resVal_npr[iBin] = N_jpsi*fb_jpsi;
      resErr_npr[iBin] = resVal_npr[iBin]*sqrt(pow((fb_jpsi_err/fb_jpsi),2)+pow((N_jpsi_err/N_jpsi),2)); //correlation needs to be taken into account
      resVal_fb[iBin] = fb_jpsi;
      resErr_fb[iBin] = fb_jpsi_err;

      iBin++;
    }//end of if TTree
  }//end of while loop on file components

  //cout<<"iBin"
  binEdges[iBin] = lastEdge;
  
  TCanvas *c1 = new TCanvas("c1", "", 800, 800);
  c1->cd();
    
  nbins = iBin;
  TH1D* results_fb = new TH1D(Form("results_fb_%s", axisName.c_str()), "", nbins, binEdges);

  results_fb->SetLineColor(kBlue+2);
  results_fb->SetLineWidth(2);
  results_fb->SetMarkerStyle(kFullCircle);
  results_fb->SetMarkerColor(kBlue+1);
  
  results_fb->GetYaxis()->SetLabelSize(0.04);
  results_fb->GetYaxis()->SetTitleSize(0.04);
  results_fb->GetYaxis()->SetTitleOffset(1.2);
  results_fb->GetYaxis()->SetTitleFont(42);
  results_fb->GetYaxis()->CenterTitle(kTRUE);
  
  results_fb->GetXaxis()->SetLabelSize(0.04);
  results_fb->GetXaxis()->CenterTitle(kTRUE);
  results_fb->GetXaxis()->SetTitleSize(0.04);
  results_fb->GetXaxis()->SetTitleFont(42);
  results_fb->GetXaxis()->SetTitleOffset(1.);
  if (axisName.find("centrality")!=std::string::npos) {
    results_fb->GetXaxis()->SetTitle("Centrality"); 
  }
  else if (axisName.find("pt")!=std::string::npos) {
    results_fb->GetXaxis()->SetTitle("p_{T}"); 
  }
  else if (axisName.find("rap")!=std::string::npos) {
    results_fb->GetXaxis()->SetTitle("y"); 
  }
  else if (axisName.find("chi2")!=std::string::npos) {
    results_fb->GetXaxis()->SetTitle("#chi^{2}_{matching}"); 
  }

  TH1D* results_pr = (TH1D*) results_fb->Clone(Form("results_pr_%s", axisName.c_str()));
  TH1D* results_npr = (TH1D*) results_fb->Clone(Form("results_npr_%s", axisName.c_str()));

  results_fb->GetYaxis()->SetTitle("b fraction");
  results_fb->GetYaxis()->SetRangeUser(0,*std::max_element(resVal_fb, resVal_fb+nbins)*1.2);
  results_pr->GetYaxis()->SetTitle("N_{prompt J/#psi}");
  results_pr->GetYaxis()->SetRangeUser(0,*std::max_element(resVal_pr, resVal_pr+nbins)*1.2);
  results_npr->GetYaxis()->SetTitle("N_{non-prompt J/#psi}");
  results_npr->GetYaxis()->SetRangeUser(0,*std::max_element(resVal_npr, resVal_npr+nbins)*1.2);
  //results_fb->Draw();
  
  for (int j=0; j<nbins; j++) {
    int jBin = results_fb->FindFixBin((binEdges[j]+binEdges[j+1])/2);
    results_fb->SetBinContent(jBin, resVal_fb[j]);
    results_fb->SetBinError(jBin, resErr_fb[j]);
    results_pr->SetBinContent(jBin, resVal_pr[j]);
    results_pr->SetBinError(jBin, resErr_pr[j]);
    results_npr->SetBinContent(jBin, resVal_npr[j]);
    results_npr->SetBinError(jBin, resErr_npr[j]);
  }

  string rangeName = Form("_pt_%d_%d_rap_%d_%d_cent_%d_%d_chi2_%d_%d",
			  (int) incMinPt,
			  (int) incMaxPt,
			  (int) incMinRap,
			  (int) incMaxRap,
			  (int) incMinCent,
			  (int) incMaxCent,
			  (int) incMinChi2,
			  (int) incMaxChi2);
  
  results_fb->Draw("E1");
  c1->SaveAs(Form("output/results_FB_%s_%s.pdf",caseName,rangeName.c_str()));
  results_pr->Draw("E1");
  c1->SaveAs(Form("output/results_PR_%s_%s.pdf",caseName,rangeName.c_str()));
  results_npr->Draw("E1");
  c1->SaveAs(Form("output/results_NPR_%s_%s.pdf",caseName,rangeName.c_str()));
  fHist->Close();
}
