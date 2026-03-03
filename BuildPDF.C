// Author: Batoul Diab
#include <vector>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "TF1.h"

#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooPoisson.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooBifurGauss.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPolynomial.h"
#include "RooFormulaVar.h"
#include "RooAddModel.h"

#include "inputUtils.h"

//using namespace RooFit;
//using namespace std;

void buildPDF_mass(RooWorkspace* ws, map<string, string> parIni, bool isMC);
void buildPDF_tauzRes(RooWorkspace* ws, map<string, string> parIni);
void buildPDF_tauzSig(RooWorkspace* ws, map<string, string> parIni);
void buildPDF_tauzBkg(RooWorkspace* ws, map<string, string> parIni);
void buildPDF_tauzTrue(RooWorkspace* ws, map<string, string> parIni);
void buildPDF_2D(RooWorkspace* ws, bool isMC);
void setDefaultParameters(map<string, string>& parIni, double nEntriesDS);

void BuildPDF(RooWorkspace* ws, map<string, string>& parIni, bool isMC, bool fitMass, bool fitTauz) {
  setDefaultParameters(parIni, ws->data("data")->numEntries());

  cout<<"[INFO] number of entries in the dataset = "<<ws->data("data")->numEntries()<<endl;
  
  for(auto it = parIni.cbegin(); it != parIni.cend(); ++it) {
    if(it->first.find("fit")!=std::string::npos || it->first.find("model")!=std::string::npos) continue;
    
    RooRealVar* var (0x0);
    // read the initial parameter string and get the initial values
    
    char delimiter = ',';
    std::string processedInput = removeSurroundingBrackets(it->second);
    std::vector<std::string> substrings = splitString(processedInput, delimiter);
    if (substrings.size() == 3) {
      cout <<"[INFO] Let's treat variable "<<it->first<<" with values "<<substrings[0]<<", "<<substrings[1]<<", "<<substrings[2]<<endl;
      var = new RooRealVar(it->first.c_str(),it->first.c_str(), std::stof(substrings[0]), std::stof(substrings[1]), std::stof(substrings[2]));
    }
    else if (substrings.size() == 1)
      var = new RooRealVar(it->first.c_str(),it->first.c_str(), std::stof(substrings[0]));
    else { cout <<"[ERROR] This initial parameter does not have the correct number of values."<<endl; continue;}

    ws->import(*var);
  }
  

  if (fitMass) buildPDF_mass(ws, parIni, isMC);
  if (fitTauz) {
    buildPDF_tauzRes(ws, parIni);
    buildPDF_tauzSig(ws, parIni);
    if (!isMC) buildPDF_tauzBkg(ws, parIni);
    //else buildPDF_tauzTrue(ws, parIni);
  }
  if (fitMass && fitTauz) {
    buildPDF_2D(ws, isMC);
    
  }
}

void buildPDF_mass(RooWorkspace* ws, map<string, string> parIni, bool isMC){
    ws->factory("expr::mean_psi2s('mean_mass + 0.589188', mean_mass)"); // 0.589188 is the m_psi2 - m_jpsi according to PDG
   
  if (parIni["modelBkg_mass"]=="Uniform") {
    ws->factory("Uniform::bkgUni_mass(mass)");
    ws->factory("RooExtendPdf::bkgPDF_mass(bkgUni_mass, fBkg_mass)");
  }
  else if (parIni["modelBkg_mass"]=="ChevPol1") {
    ws->factory("RooChebychev::bkgChev_mass(mass, {c0_mass, c1_mass})");
    ws->factory("RooExtendPdf::bkgPDF_mass(bkgChev_mass, fBkg_mass)");
  }
  else if (parIni["modelBkg_mass"]=="ChevPol2") {
    ws->factory("RooChebychev::bkgChev_mass(mass, {c0_mass, c1_mass, c2_mass})");
    ws->factory("RooExtendPdf::bkgPDF_mass(bkgChev_mass, fBkg_mass)");
  }
  else if (parIni["modelBkg_mass"]=="ChevPol3") {
    ws->factory("RooChebychev::bkgChev_mass(mass, {c0_mass, c1_mass, c2_mass, c3_mass})");
    ws->factory("RooExtendPdf::bkgPDF_mass(bkgChev_mass, fBkg_mass)");
  }
  else if (parIni["modelBkg_mass"]=="ChevPol4") {
    ws->factory("RooChebychev::bkgChev_mass(mass, {c0_mass, c1_mass, c2_mass, c3_mass, c4_mass})");
    ws->factory("RooExtendPdf::bkgPDF_mass(bkgChev_mass, fBkg_mass)");
  }
  else if (parIni["modelBkg_mass"]=="ChevVarPol3") { // 'var' means we can use it for iterative fitting (parameters constant)
    ws->Print();
    ws->factory("RooChebychev::bkgChev_mass(mass, {c0_mass, c1_mass, c2_mass, c3_mass})");
    ws->factory("RooExtendPdf::bkgPDF_mass(bkgChev_mass, fBkg_mass)");
    // Set coefficients constant, except first (to be released order-by-order later depending on the fit chi2)
    for (int i = 1; i <= 3; ++i) {
      RooRealVar* coeff = ws->var(Form("c%d_mass", i));
      if (coeff) coeff->setConstant(kTRUE);
    }
    ws->Print();
  }
  else if (parIni["modelBkg_mass"]=="ChevVarPol6") { // 'var' means we can use it for iterative fitting (parameters constant)
    ws->Print();
    ws->factory("RooChebychev::bkgChev_mass(mass, {c0_mass, c1_mass, c2_mass, c3_mass, c4_mass, c5_mass, c6_mass})");
    ws->factory("RooExtendPdf::bkgPDF_mass(bkgChev_mass, fBkg_mass)");
    // Set coefficients constant, except first (to be released order-by-order later depending on the fit chi2)
    for (int i = 1; i <= 6; ++i) {
      RooRealVar* coeff = ws->var(Form("c%d_mass", i));
      if (coeff) coeff->setConstant(kTRUE);
    }
    ws->Print();
  }

  if (parIni["modelSig_mass"]=="Gauss") {
    ws->factory("RooGaussian::jpsiGaus_mass(mass, mean_mass, sigma_mass)");
    ws->factory("RooExtendPdf::jpsiPDF_mass(jpsiGaus_mass, fJpsi_mass)");
    ws->factory("RooGaussian::psi2sGaus_mass(mass, mean_psi2s, sigma_mass)");
    ws->factory("RooExtendPdf::psi2sPDF_mass(psi2sGaus_mass, fPsi2s_mass)");
  }
  else if (parIni["modelSig_mass"]=="extCB") {
    ws->factory("RooExtCBShape::jpsiExtCB_mass(mass, mean_mass, sigma_mass, alpha0_mass, n0_mass, alpha1_mass, n1_mass)");
    ws->factory("RooExtendPdf::jpsiPDF_mass(jpsiExtCB_mass, fJpsi_mass)");
    ws->factory("RooExtCBShape::psi2sExtCB_mass(mass, mean_psi2s, sigma_mass, alpha0_mass, n0_mass, alpha1_mass, n1_mass)");
    ws->factory("RooExtendPdf::psi2sPDF_mass(psi2sExtCB_mass, fPsi2s_mass)");
  }
  
  // Combine the components into a composite model
  RooAddPdf* totPDF_mass = new RooAddPdf ("totPDF_mass", "model for mass fit", RooArgList(*ws->pdf("jpsiPDF_mass"), *ws->pdf("psi2sPDF_mass"), *ws->pdf("bkgPDF_mass")));
  if (isMC) totPDF_mass = new RooAddPdf ("totPDF_mass", "model for mass fit", RooArgList(*ws->pdf("jpsiPDF_mass"), *ws->pdf("psi2sPDF_mass")));
  ws->import(*totPDF_mass);

}

void buildPDF_tauzRes(RooWorkspace* ws, map<string, string> parIni){
  if (parIni["model_tauzRes"]=="extCB") {
    ws->factory("RooExtCBShape::tauzResPDF(tauz, mean_tauzRes, sigma_tauzRes, alpha_tauzRes, n_tauzRes, alpha_tauzRes, n_tauzRes)");
  }
  else if (parIni["model_tauzRes"]=="CB") {
    ws->factory("RooCBShape::tauzResPDF(tauz, mean_tauzRes, sigma_tauzRes, alpha_tauzRes, n_tauzRes)");
  }
  else if (parIni["model_tauzRes"]=="Gauss3") {
    ws->factory("Gaussian::gauss0_tauzRes(tauz, mean_tauzRes, sigma0_tauzRes)");
    ws->factory("Gaussian::gauss1_tauzRes(tauz, mean_tauzRes, sigma1_tauzRes)");
    ws->factory("Gaussian::gauss2_tauzRes(tauz, mean_tauzRes, sigma2_tauzRes)");
    ws->factory("expr::fGaus1m01_tauzRes('1 - (fGaus0_tauzRes + fGaus1_tauzRes)', {fGaus0_tauzRes, fGaus1_tauzRes})");
    ws->factory("SUM::tauzResPDF(fGaus0_tauzRes*gauss0_tauzRes, fGaus1_tauzRes*gauss1_tauzRes, fGaus1m01_tauzRes*gauss2_tauzRes)");
    
    ws->factory("RooGaussModel::gauss0(tauz, mean_tauzRes, sigma0_tauzRes)");
    ws->factory("RooGaussModel::gauss1(tauz, mean_tauzRes, sigma1_tauzRes)");
    ws->factory("RooGaussModel::gauss2(tauz, mean_tauzRes, sigma2_tauzRes)");

    ws->factory("RooAddModel::tauzResModel({gauss0, gauss1, gauss2},{fGaus0_tauzRes, fGaus1_tauzRes})");
  }

    else if (parIni["model_tauzRes"]=="Gauss4") {
    ws->factory("Gaussian::gauss0_tauzRes(tauz, mean_tauzRes, sigma0_tauzRes)");
    ws->factory("Gaussian::gauss1_tauzRes(tauz, mean_tauzRes, sigma1_tauzRes)");
    ws->factory("Gaussian::gauss2_tauzRes(tauz, mean_tauzRes, sigma2_tauzRes)");
    ws->factory("Gaussian::gauss3_tauzRes(tauz, mean_tauzRes, sigma3_tauzRes)");
    ws->factory("expr::fGaus1m012_tauzRes('1 - (fGaus0_tauzRes + fGaus1_tauzRes + fGaus2_tauzRes)', {fGaus0_tauzRes, fGaus1_tauzRes, fGaus2_tauzRes})");
    ws->factory("SUM::tauzResPDF(fGaus0_tauzRes*gauss0_tauzRes, fGaus1_tauzRes*gauss1_tauzRes, fGaus2_tauzRes*gauss2_tauzRes, fGaus1m012_tauzRes*gauss3_tauzRes)");

    
    ws->factory("RooGaussModel::gauss0(tauz, mean_tauzRes, sigma0_tauzRes)");
    ws->factory("RooGaussModel::gauss1(tauz, mean_tauzRes, sigma1_tauzRes)");
    ws->factory("RooGaussModel::gauss2(tauz, mean_tauzRes, sigma2_tauzRes)");
    ws->factory("RooGaussModel::gauss3(tauz, mean_tauzRes, sigma3_tauzRes)");

    ws->factory("RooAddModel::tauzResModel({gauss0, gauss1, gauss2, gauss3},{fGaus0_tauzRes, fGaus1_tauzRes, fGaus2_tauzRes})");
  }
}
void buildPDF_tauzBkg(RooWorkspace* ws, map<string, string> parIni){
  if (parIni["modelNpr_tauzBkg"]=="TripleDecay"){
    ws->factory("Decay::ssdNpr_tauzBkg(tauz, lambdaDssNpr_tauzBkg, tauzResModel, RooDecay::SingleSided)");
    ws->factory("Decay::dfNpr_tauzBkg(tauz, lambdaDfNpr_tauzBkg, tauzResModel, RooDecay::Flipped)");
    ws->factory("Decay::dsdNpr_tauzBkg(tauz, lambdaDdsNpr_tauzBkg, tauzResModel, RooDecay::DoubleSided)");
    ws->factory("SUM::tauzNprBkgPDF1(fDfssNpr_tauzBkg*ssdNpr_tauzBkg, dfNpr_tauzBkg)");
    ws->factory("SUM::tauzNprBkgPDF(fDNpr_tauzBkg*tauzNprBkgPDF1, dsdNpr_tauzBkg)");
  }
  ws->factory("SUM::tauzPrBkgPDF(tauzResPDF)");
  //RooAbsPdf* tauzResPDF = ws->pdf("tauzResPDF");
  //ws->import(*tauzResPDF, Rename("tauzPrBkgPDF"));
  ws->factory("SUM::tauzBkgPDF(fb_tauzBkg*tauzNprBkgPDF, tauzPrBkgPDF)");
  
}
void buildPDF_tauzSig(RooWorkspace* ws, map<string, string> parIni){
  if (parIni["modelNpr_tauzSig"]=="SingleSidedDecay"){
    ws->factory("Decay::tauzNprSigPDF(tauz, lambdaDssNpr_tauzSig, tauzResModel, RooDecay::SingleSided)");
  }
  ws->factory("SUM::tauzPrSigPDF(tauzResPDF)");
  ws->factory("SUM::tauzSigPDF(fb_tauzSig*tauzNprSigPDF, tauzPrSigPDF)");
}

void buildPDF_2D(RooWorkspace* ws, bool isMC){  
  ws->factory("PROD::tauzMassPrJpsiPDF(tauzPrSigPDF, jpsiPDF_mass)");
  ws->factory("PROD::tauzMassNprJpsiPDF(tauzNprSigPDF, jpsiPDF_mass)");
  ws->factory("SUM::tauzMassTotJpsiPDF(b_jpsi_tauzMass*tauzMassNprJpsiPDF, tauzMassPrJpsiPDF)");

  ws->factory("PROD::tauzMassPrPsi2sPDF(tauzPrSigPDF, psi2sPDF_mass)");
  ws->factory("PROD::tauzMassNprPsi2sPDF(tauzNprSigPDF, psi2sPDF_mass)");
  ws->factory("SUM::tauzMassTotPsi2sPDF(b_psi2s_tauzMass*tauzMassNprPsi2sPDF, tauzMassPrPsi2sPDF)");

  if (!isMC) {
    ws->factory("PROD::tauzMassPrBkgPDF(tauzPrBkgPDF, bkgPDF_mass)");
    ws->factory("PROD::tauzMassNprBkgPDF(tauzNprBkgPDF, bkgPDF_mass)");
    ws->factory("SUM::tauzMassTotBkgPDF(b_bkg_tauzMass*tauzMassNprBkgPDF, tauzMassPrBkgPDF)");
    
    ws->factory("SUM::totPDF_2D(fJpsi_tauzMass*tauzMassTotJpsiPDF, fPsi2s_tauzMass*tauzMassTotPsi2sPDF, fBkg_tauzMass*tauzMassTotBkgPDF)");
  }
  else {
    ws->factory("SUM::totPDF_2D(fJpsi_tauzMass*tauzMassTotJpsiPDF, fPsi2s_tauzMass*tauzMassTotPsi2sPDF");
  }
}

void setDefaultParameters(map<string, string>& parIni, double nEntriesDS){
  std::map<std::string, std::string> modelMap;
  modelMap["modelBkg_mass"] = "ChevPol2";
  modelMap["modelSig_mass"] = "extCB";
  modelMap["model_tauzRes"] = "Gauss3";
  modelMap["modelNpr_tauzBkg"] = "TripleDecay";
  modelMap["modelNpr_tauzSig"] = "SingleSidedDecay";
  
  std::map<std::string, std::vector<double>> varMap;
  varMap["b_jpsi_tauzMass"] = {0.2, 0.01, 0.5};
  varMap["b_psi2s_tauzMass"] = {0.2, 0, 1};
  varMap["b_bkg_tauzMass"] = {0.2, 0, 1};
  varMap["fJpsi_tauzMass"] = {0.01*nEntriesDS, 0, 2.0*nEntriesDS};
  varMap["fPsi2s_tauzMass"] = {0.01*nEntriesDS, 0, 2.0*nEntriesDS};
  varMap["fBkg_tauzMass"] = {0.1*nEntriesDS, 0, 2.0*nEntriesDS};
  
  varMap["fJpsi_mass"] = {0.1*nEntriesDS, 0, 2.0*nEntriesDS};
  varMap["fPsi2s_mass"] = {0.001*nEntriesDS, 0, 0.2*nEntriesDS};
  varMap["fBkg_mass"] = {0.1*nEntriesDS, 0, 2.0*nEntriesDS};
  
  varMap["mean_mass"] = {3.096,2.9,3.3};
  varMap["sigma_mass"] = {0.02,0,0.1};
  varMap["alpha0_mass"] = {1.,0,5};
  varMap["n0_mass"] = {3.,0,5};
  varMap["alpha1_mass"] = {1.,0,5};
  varMap["n1_mass"] = {3.,0,5};
  
  varMap["c0_mass"] = {0,-2,2};
  varMap["c1_mass"] = {0,-2,2};
  varMap["c2_mass"] = {0,-2,2};
  varMap["c3_mass"] = {0,-2,2};
  varMap["c4_mass"] = {0,-2,2};
  
  varMap["mean_tauzRes"] = {0,-0.1,0.1};
  varMap["sigma_tauzRes"] = {0.00045,0.0002,0.0010};
  varMap["alpha_tauzRes"] = {1,0.,3};
  varMap["n_tauzRes"] = {1.5,0.,10};
  varMap["sigma0_tauzRes"] = {0.00045,0.0002,0.0010};
  varMap["sigma1_tauzRes"] = {0.00045,0.0002,0.0010};
  varMap["sigma2_tauzRes"] = {0.00045,0.0002,0.0010};
  varMap["sigma3_tauzRes"] = {0.00045,0.0002,0.0030};
  varMap["fGaus0_tauzRes"] = {0.3, 0, 1};
  varMap["fGaus1_tauzRes"] = {0.06, 0, 1};
  varMap["fGaus2_tauzRes"] = {0.05, 0, 1};
  varMap["fGaus3_tauzRes"] = {0.05, 0, 1};
    
  varMap["fb_tauzSig"] = {0.2, 0.01, 0.5};
  varMap["lambdaDssNpr_tauzSig"] = {1., 0.01, 2.0};
  
  varMap["fDfssNpr_tauzBkg"] = {0.8, 0.2, 0.9};
  varMap["fDNpr_tauzBkg"] = {0.9, 0.2, 0.95};
  varMap["fb_tauzBkg"] = {0.8, 0., 1.};
  varMap["lambdaDssNpr_tauzBkg"] = {0.055, 0.0001, 0.02};
  varMap["lambdaDfNpr_tauzBkg"] = {0.030, 0.00001, 0.02};
  varMap["lambdaDdsNpr_tauzBkg"] = {0.045, 0.0001, 0.02};
  
  for (auto it = varMap.begin(); it != varMap.end(); ++it) {
    cout<<"[INFO] adding the default parameters for "<<it->first<<endl;
    if (parIni.count(it->first)==0 || parIni[it->first]=="") {
      cout<<"[INFO] this parameter does not exist in the input files so adding it now"<<endl;
      parIni[it->first] = Form("[%f,%f,%f]", it->second[0], it->second[1], it->second[2]);
      cout<<"[Info] added "<<it->first<<" = "<< parIni[it->first]<<endl;
    }
  }

  for (auto it = modelMap.begin(); it != modelMap.end(); ++it) {
    if (parIni.count(it->first)==0 || parIni[it->first]=="") {
      cout<<"[INFO] this fit model does not exist in the input files so adding it now"<<endl;
      parIni[it->first] = it->second;
    }
  }
}

void fixParPDF(RooWorkspace* ws, RooFitResult* fitResult, map<string, string> &parIni, bool ispO, string rangeLabel, const char *caseName, bool fromMassPDF, bool fromTauzResPDF, bool fromTauzBkgPDF) {
  std::vector<std::string> fixedPars;
  if (fromMassPDF) {
    if (parIni["modelSig_mass"]=="Gauss") {
      fixedPars.push_back("mean_mass");
      fixedPars.push_back("sigma_mass");
    }
    else if (parIni["modelSig_mass"]=="extCB") {
      fixedPars.push_back("mean_mass");
      fixedPars.push_back("sigma_mass");
      fixedPars.push_back("alpha0_mass");
      fixedPars.push_back("n0_mass");
      fixedPars.push_back("alpha1_mass");
      fixedPars.push_back("n1_mass");
    }
    if (parIni["modelBkg_mass"]=="ChevPol1") {
      fixedPars.push_back("c0_mass");
      fixedPars.push_back("c1_mass");
    }
    else if (parIni["modelBkg_mass"]=="ChevPol2") {
      fixedPars.push_back("c0_mass");
      fixedPars.push_back("c1_mass");
      fixedPars.push_back("c2_mass");
    }
    else if (parIni["modelBkg_mass"]=="ChevPol3") {
      fixedPars.push_back("c0_mass");
      fixedPars.push_back("c1_mass");
      fixedPars.push_back("c2_mass");
      fixedPars.push_back("c3_mass");
    }
    else if (parIni["modelBkg_mass"]=="ChevPol4") {
      fixedPars.push_back("c0_mass");
      fixedPars.push_back("c1_mass");
      fixedPars.push_back("c2_mass");
      fixedPars.push_back("c3_mass");
      fixedPars.push_back("c4_mass");
    }
  }
  else if (fromTauzResPDF) {
    if (parIni["model_tauzRes"]=="extCB" || parIni["model_tauzRes"]=="CB") {
      //fixedPars.push_back("mean_tauzRes");
      fixedPars.push_back("sigma_tauzRes");
      fixedPars.push_back("alpha_tauzRes");
      fixedPars.push_back("n_tauzRes");
    }
    else if (parIni["model_tauzRes"]=="Gauss3" || parIni["model_tauzRes"]=="Gauss4") {
      //RooRealVar* sigma0_free_tauzRes = (RooRealVar*) ws->var("sigma0_tauzRes");
      //ws->import(*sigma0_free_tauzRes, Rename("sigma0_free_tauzRes"));
      
      fixedPars.push_back("sigma0_tauzRes");
      fixedPars.push_back("sigma1_tauzRes");
      fixedPars.push_back("sigma2_tauzRes");
      fixedPars.push_back("fGaus0_tauzRes");
      fixedPars.push_back("fGaus1_tauzRes");
      //fixedPars.push_back("fGaus2_tauzRes");
      if (parIni["model_tauzRes"]=="Gauss4") {
	fixedPars.push_back("sigma3_tauzRes");
	fixedPars.push_back("fGaus2_tauzRes");
      }
      for (const auto& par : fixedPars) {
	RooRealVar* freePar = (RooRealVar*) ws->var(par);
	if (!freePar) {
	  std::cerr << "Warning: parameter " << par << " not found in workspace!\n";
	  continue;
	}
	RooRealVar clonePar(*freePar, Form("free_%s", par.c_str()));
	ws->import(clonePar);

	string parFree = "free_" + par;
	if (parIni.count(parFree)==0 || parIni[parFree]=="") {
	  parIni[parFree] = Form("[%f,%f,%f]", freePar->getVal(), freePar->getVal()-freePar->getError(), freePar->getVal()+freePar->getError());
	}
	
      }
    }
  }
  else if (fromTauzBkgPDF) {
    if (parIni["modelNpr_tauzBkg"]=="TripleDecay") {
      fixedPars.push_back("fDfssNpr_tauzBkg");
      fixedPars.push_back("fDNpr_tauzBkg");
      fixedPars.push_back("lambdaDssNpr_tauzBkg");
      fixedPars.push_back("lambdaDdsNpr_tauzBkg");
      fixedPars.push_back("lambdaDfNpr_tauzBkg");
    }
  }

  
  if (!fitResult) {
    cout<<"[INFO] fixing parameters from previous fits"<<endl;
    string resFileName = Form("output/output_fitMass_%s_%s.root", ispO?"pO":"OO", /*isMC?"_MC":"",*/ caseName);
    TFile* resFile =  TFile::Open(resFileName.c_str(),"READ");
    TTree* resTree = (TTree*) resFile->Get(Form("tree_%s", rangeLabel.c_str()));
    std::map<std::string, double> fixValues;
    for (const auto& par : fixedPars) {
      fixValues[par] = 0.0;
      resTree->SetBranchAddress(par.c_str(), &fixValues[par]);
    }
    resTree->GetEntry(0);
    for (const auto& [name, val] : fixValues) {
      ws->var(name.c_str())->setVal(val);
      ws->var(name.c_str())->setConstant(kTRUE);
    }
  }
  else if (fitResult) {
    fitResult->Print();
    for (const auto& par : fixedPars) {
      cout<<"[INFO] let's fix "<<par<<endl;
      RooRealVar* var = (RooRealVar*) fitResult->floatParsFinal().find(par.c_str());
      ws->var(par.c_str())->setVal(var->getVal());
      ws->var(par.c_str())->setConstant(kTRUE);
    }
  }
}
