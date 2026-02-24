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
#include <TLatex.h>

#include "SignalExtraction.C"

void signalExtraction(bool ispO=true, bool isMC =false, const char *caseName = "nominal", bool remakeDS =false, bool fitMass=true, bool fitTauz=false);
void plotResult(bool ispO=true, const char *caseName = "nominal", string axisName = "pt", int incMinCent=0, int incMaxCent=100, float incMinPt=0., float incMaxPt=50., float incMinRap=-3.5, float incMaxRap=-2.5, float incMinChi2=0, float incMaxChi2=50, bool diffChi2=false);

void InputToResults(bool ispO=true, bool isMC=false, const char *caseName = "nominal", bool remakeDS = false, bool fitMass1D=true, bool fitTauz1D=false, bool fit2D=false, bool plotResults = false) {
  gSystem->Load("RooExtCBShape.cxx+");
  //cout<<"plotResults = "<<plotResults<<endl;
  if (fitMass1D) {
    signalExtraction(ispO, isMC, caseName, remakeDS, true, false);
    remakeDS=false;
  }
  if (fitTauz1D) {
    signalExtraction(ispO, isMC, caseName, remakeDS, false, true);
    remakeDS=false;
  }
  if (fit2D)
    signalExtraction(ispO, isMC, caseName, remakeDS, true, true);
  

  if (plotResults) {
    int minCent = 0;
    int maxCent = 100;
    float minPt = 0.;
    float maxPt = 20.;
    float minRap = 2.5;//-3.5;
    float maxRap = 3.6;//-2.6;
    float minChi2 = 0;
    float maxChi2 = 50;
    string axisName = "pt";
    bool diffChi2 = false;
    if (strstr(caseName, "diffChi2") != NULL) diffChi2 = true; 
    plotResult(ispO, caseName, axisName.c_str(), minCent, maxCent, minPt, maxPt, minRap, maxRap, minChi2, maxChi2, diffChi2);
  }
  
}

void signalExtraction(bool ispO, bool isMC, const char *caseName, bool remakeDS, bool fitMass, bool fitTauz) {
  
  //do the fits or at least some of them
  vector< struct KinCuts >       cutVector;
  vector< map<string, string> >  parIniVector;
  vector< map<string, double> >  allResults;
  
  if (!addParameters(Form("inputFiles/initialPars_%s_%s_%s.txt", fitMass?(fitTauz?"tauz":"mass"):"tauz", ispO?"pO":"OO", caseName), cutVector, parIniVector)) { return; } //if 1D fit on mass or tauz it reads the corresponding file, if it's 2D it reads the mass but the other variables get added from the default values and then they are fixed
  
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
				  (int) (10*cutVector[j].pt.Min),
				  (int) (10*cutVector[j].pt.Max),
				  (int) (10*cutVector[j].rap.Min),
				  (int) (10*cutVector[j].rap.Max),			    
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


void plotResult(bool ispO, const char *caseName, string axisName, int incMinCent, int incMaxCent, float incMinPt, float incMaxPt, float incMinRap, float incMaxRap, float incMinChi2, float incMaxChi2, bool diffChi2){
  //cout <<"[INFO] The plotting function is yet to be done"<<endl;
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
  
  //string fileName = Form("output/output_%s.root",caseName);
  string fileName = Form("output/output_fitMassTauz_%s_%s.root", ispO?"pO":"OO", caseName);
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
      
      cout<<"[INFO] Reading the input tree "<<resTree->GetName()<<endl;

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
      double N_jpsi = 0.; resTree->SetBranchAddress("fJpsi_tauzMass", &N_jpsi);
      double N_jpsi_err = 0.; resTree->SetBranchAddress("fJpsi_tauzMass_err", &N_jpsi_err);
      /*
      double fb_psi2s = 0.; resTree->SetBranchAddress("b_psi2s_tauzMass", &fb_psi2s);
      double fb_psi2s_err = 0.; resTree->SetBranchAddress("b_psi2s_tauzMass_err", &fb_psi2s_err);
      double N_psi2s = 0.; resTree->SetBranchAddress("fPsi2s_tauzMass", &N_Psi2s);
      double N_psi2s_err = 0.; resTree->SetBranchAddress("fPsi2s_tauzMass_err", &N_Psi2s_err);
      */
      
      resTree->GetEntry(0);

      cout<<"[INFO] This tree has the following bin selection, getting the results for pt ["<<ptMin<<"-"<<ptMax<<"], rap ["<<rapMin<<"-"<<rapMax<<"], cent ["<<centMin<<"-"<<centMax<<"], chi2 ["<<chi2Min<<"-"<<chi2Max<<"]"<<endl;
      if (axisName.find("centrality")!=std::string::npos) { //if plotting as function of centrality, the centrality of the result needs to be in the range but the pt need to be the exact edges
	
        if (centMin < incMinCent || centMax > incMaxCent)
          continue;
	if (fabs(centMin - incMinCent)<0.00001 && fabs(centMax - incMaxCent)<0.00001)
	  continue;
        if (fabs(ptMin - incMinPt)>0.00001 || fabs(ptMax - incMaxPt)>0.00001)
          continue;
        if (fabs(rapMin - incMinRap)>0.00001 || fabs(rapMax - incMaxRap)>0.00001)
          continue;
	if (!diffChi2 && (fabs(chi2Min - incMinChi2)>0.00001 || fabs(chi2Max - incMaxChi2)>0.00001))
	  continue;
      
	//cout<<"[INFO] This tree passed the bin selection, getting the results for pt ["<<ptMin<<"-"<<ptMax<<"], rap ["<<rapMin<<"-"<<rapMax<<"], cent ["<<centMin<<"-"<<centMax<<"]"<<endl;
	binEdges[iBin] = centMin;
	lastEdge = centMax;
      }//end of if centrality
      else if (axisName.find("pt")!=std::string::npos) {
        if (ptMin < incMinPt || ptMax > incMaxPt)
          continue;
	if (fabs(ptMin - incMinPt) < 0.00001 && fabs(ptMax - incMaxPt) <0.00001)
	  continue;
        if (fabs(centMin - incMinCent)>0.00001 || fabs(centMax - incMaxCent)>0.00001)
          continue;
        if (fabs(rapMin - incMinRap)>0.00001 || fabs(rapMax - incMaxRap)>0.00001)
          continue;
	if (!diffChi2 && (fabs(chi2Min - incMinChi2)>0.00001 || fabs(chi2Max - incMaxChi2)>0.00001))
	  continue;

	//cout<<"[INFO] This tree passed the bin selection, getting the results for pt ["<<ptMin<<"-"<<ptMax<<"], rap ["<<rapMin<<"-"<<rapMax<<"], cent ["<<centMin<<"-"<<centMax<<"]"<<endl;
		
	binEdges[iBin] = ptMin;
	lastEdge = ptMax;
      }//end of if pt
      else if (axisName.find("rap")!=std::string::npos) {
	if (rapMin < incMinRap || rapMax > incMaxRap)
          continue;
	if (fabs(rapMin - incMinRap)<0.00001 && fabs(rapMax - incMaxRap)<0.00001)
	  continue;
        if (fabs(centMin - incMinCent)>0.00001 || fabs(centMax - incMaxCent)>0.00001)
          continue;
        if (fabs(ptMin - incMinPt)>0.00001 || fabs(ptMax - incMaxPt)>0.00001)
          continue;
	if (!diffChi2 && (fabs(chi2Min - incMinChi2)>0.00001 || fabs(chi2Max - incMaxChi2)>0.00001))
	  continue;

	//cout<<"[INFO] This tree passed the bin selection, getting the results for pt ["<<ptMin<<"-"<<ptMax<<"], rap ["<<rapMin<<"-"<<rapMax<<"], cent ["<<centMin<<"-"<<centMax<<"]"<<endl;
	
	binEdges[iBin] = rapMin;
	lastEdge = rapMax;
      }//end of rap
      else if (axisName.find("chi2")!=std::string::npos) {
	if (chi2Min < incMinChi2 || chi2Max > incMaxChi2)
          continue;
	if (fabs(chi2Min - incMinChi2)<0.00001 && fabs(chi2Max - incMaxChi2)<0.00001)
	  continue;
        if (fabs(centMin - incMinCent)>0.00001 || fabs(centMax - incMaxCent)>0.00001)
          continue;
        if (fabs(rapMin - incMinRap)>0.00001 || fabs(rapMax - incMaxRap)>0.00001)
          continue;
	if (fabs(ptMin - incMinPt)>0.00001 || fabs(ptMax - incMaxPt)>0.00001)
	  continue;
	//cout<<"[INFO] This tree passed the bin selection, getting the results for pt ["<<ptMin<<"-"<<ptMax<<"], rap ["<<rapMin<<"-"<<rapMax<<"], cent ["<<centMin<<"-"<<centMax<<"]"<<endl;
	
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
      
  nbins = iBin;
  TH1D* results_fb = new TH1D(Form("results_fb_%s", axisName.c_str()), "", nbins, binEdges);

  results_fb->SetLineColor(kBlue+2);
  results_fb->SetLineWidth(2);
  results_fb->SetMarkerStyle(kFullCircle);
  results_fb->SetMarkerColor(kBlue+1);
  
  results_fb->GetYaxis()->SetLabelSize(0.03);
  results_fb->GetYaxis()->SetTitleSize(0.04);
  results_fb->GetYaxis()->SetTitleOffset(1.3);
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
  results_fb->GetYaxis()->SetRangeUser(0,0.5);//*std::max_element(resVal_fb, resVal_fb+nbins)*1.2);
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
			  (int) (10*incMinPt),
			  (int) (10*incMaxPt),
			  (int) (10*incMinRap),
			  (int) (10*incMaxRap),
			  (int) incMinCent,
			  (int) incMaxCent,
			  (int) incMinChi2,
			  (int) incMaxChi2);


  TCanvas *cfb = new TCanvas("cfb", "", 800, 800);
  cfb->cd();
  results_fb->Draw("E1");
  
  float xText=0.5; float yText=0.5;
  
  TLatex* textAlice_fb = AliceText(ispO);
  TLatex* textCut_fb = cutTextResult(ispO, axisName,xText, yText, incMinCent, incMaxCent, incMinPt, incMaxPt, incMaxRap, incMinRap, incMinChi2, incMaxChi2, diffChi2);
  textAlice_fb->Draw();
  textCut_fb->Draw();
  cfb->SaveAs(Form("output/results_%s_FB_vs%s_%s%s.pdf", ispO?"pO":"OO", axisName.c_str(), caseName, rangeName.c_str()));

  
  TCanvas *cpr = new TCanvas("cpr", "", 800, 800);
  cpr->cd();
  results_pr->Draw("E1");
  TLatex* textAlice_pr = AliceText(ispO);
  TLatex* textCut_pr = cutTextResult(ispO, axisName,xText, yText, incMinCent, incMaxCent, incMinPt, incMaxPt, incMaxRap, incMinRap, incMinChi2, incMaxChi2, diffChi2);
  textAlice_pr->Draw();
  textCut_pr->Draw("same");
  cpr->SaveAs(Form("output/results_%s_PR_vs%s_%s%s.pdf", ispO?"pO":"OO", axisName.c_str(), caseName, rangeName.c_str()));

  TCanvas *cnpr = new TCanvas("cnpr", "", 800, 800);
  cnpr->cd();
  results_npr->Draw("E1");
  TLatex* textAlice_npr = AliceText(ispO);
  TLatex* textCut_npr = cutTextResult(ispO, axisName,xText, yText, incMinCent, incMaxCent, incMinPt, incMaxPt, incMaxRap, incMinRap, incMinChi2, incMaxChi2, diffChi2);
  textAlice_npr->Draw("same");
  textCut_npr->Draw("same");
  cnpr->SaveAs(Form("output/results_%s_NPR_vs%s_%s%s.pdf", ispO?"pO":"OO", axisName.c_str(), caseName, rangeName.c_str()));
  
  fHist->Close();
}
