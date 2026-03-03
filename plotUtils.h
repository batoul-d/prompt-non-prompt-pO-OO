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

string varFancyLabel(string varLabel="") {
  std::map<std::string, std::string> fancyLabels;

  fancyLabels["b_jpsi_tauzMass"] = "f_{b}^{J/#psi}";
  fancyLabels["b_psi2s_tauzMass"] = "f_{b}^{#psi(2S)}";
  fancyLabels["b_bkg_tauzMass"] = "f_{b}^{Bkg}";
  fancyLabels["fJpsi_tauzMass"] = "N^{J/#psi}";
  fancyLabels["fPsi2s_tauzMass"] = "N^{#psi(2S)}";
  fancyLabels["fBkg_tauzMass"] = "N^{Bkg}";
  
  fancyLabels["fJpsi_mass"] = "N^{J/#psi}";
  fancyLabels["fPsi2s_mass"] = "N^{#psi(2S)}";
  fancyLabels["fBkg_mass"] = "N^{Bkg}";
  
  fancyLabels["mean_mass"] = "m^{J/#psi}";
  fancyLabels["mean_psi2s"] = "m^{#psi(2S)}";
  fancyLabels["sigma_mass"] = "#sigma^{J/#psi}";
  fancyLabels["alpha0_mass"] = "#alpha_{0}^{J/#psi}";
  fancyLabels["n0_mass"] = "n_{0}^{J/#psi}";
  fancyLabels["alpha1_mass"] = "#alpha_{1}^{J/#psi}";
  fancyLabels["n1_mass"] = "n_{1}^{J/#psi}";
  fancyLabels["c0_mass"] = "c^{Bkg}_{0}";
  fancyLabels["c1_mass"] = "c^{Bkg}_{1}";
  fancyLabels["c2_mass"] = "c^{Bkg}_{2}";
  fancyLabels["c3_mass"] = "c^{Bkg}_{3}";
  fancyLabels["c4_mass"] = "c^{Bkg}_{4}";
  
  fancyLabels["mean_tauzRes"] = "#tau_{z0}";
  fancyLabels["sigma_tauzRes"] = "#sigma^{#tau_{z}}";
  fancyLabels["alpha_tauzRes"] = "#alpha^{#tau_{z}}";
  fancyLabels["n_tauzRes"] = "n^{#tau_{z}}";
  fancyLabels["sigma0_tauzRes"] = "#sigma_{0}^{#tau_{z}}";
  fancyLabels["sigma1_tauzRes"] = "#sigma_{1}^{#tau_{z}}";
  fancyLabels["sigma2_tauzRes"] = "#sigma_{2}^{#tau_{z}}";
  fancyLabels["sigma3_tauzRes"] = "#sigma_{3}^{#tau_{z}}";
  fancyLabels["free_sigma0_tauzRes"] = "#sigma_{0}^{#tau_{z}}";
  fancyLabels["free_sigma1_tauzRes"] = "#sigma_{1}^{#tau_{z}}";
  fancyLabels["free_sigma2_tauzRes"] = "#sigma_{2}^{#tau_{z}}";
  fancyLabels["free_sigma3_tauzRes"] = "#sigma_{3}^{#tau_{z}}";
  fancyLabels["fGaus0_tauzRes"] = "f^{gaus}_{0}";
  fancyLabels["fGaus1_tauzRes"] = "f^{gaus}_{1}";
  fancyLabels["fGaus2_tauzRes"] = "f^{gaus}_{2}";
  fancyLabels["fGaus3_tauzRes"] = "f^{gaus}_{3}";
  fancyLabels["free_fGaus0_tauzRes"] = "f^{gaus}_{0}";
  fancyLabels["free_fGaus1_tauzRes"] = "f^{gaus}_{1}";
  fancyLabels["free_fGaus2_tauzRes"] = "f^{gaus}_{2}";
  fancyLabels["free_fGaus3_tauzRes"] = "f^{gaus}_{3}";
  fancyLabels["fb_tauzSig"] = "f_{b}^{Sig}";
  fancyLabels["fDfssNpr_tauzBkg"] = "f_{dfss}";
  fancyLabels["fDNpr_tauzBkg"] = "f_{d}";
  fancyLabels["fb_tauzBkg"] = "f_{b}^{Bkg}";
  fancyLabels["lambdaDssNpr_tauzSig"] = "#lambda^{SSD}";
  fancyLabels["lambdaDssNpr_tauzBkg"] = "#lambda^{SSD}";
  fancyLabels["lambdaDfNpr_tauzBkg"] = "#lambda^{FD}";
  fancyLabels["lambdaDdsNpr_tauzBkg"] = "#lambda^{DSD}";
  
  if (fancyLabels.find(varLabel) != fancyLabels.end())
    return fancyLabels[varLabel]; 
  else return varLabel;
}

double getMax(RooHist* hist){
  if (hist) {
    int n = hist->GetN();
    double* x = hist->GetX();
    double* y = hist->GetY();
    
    double ymax = y[0];
    double xmax = x[0];
    for (int i = 1; i < n; ++i) {
      if (y[i] > ymax) {
	ymax = y[i];
	xmax = x[i];
      }
    }
    
    return xmax;
  }
  else return -999;
}

void fixFrameStyle(RooPlot* distFrame, bool setLog){
  distFrame->SetTitle(" ");
  distFrame->GetXaxis()->SetTitleSize(0.035);
  distFrame->GetYaxis()->SetTitleSize(0.035);
  distFrame->GetXaxis()->SetLabelSize(0.025);
  distFrame->GetYaxis()->SetLabelSize(0.025);
  distFrame->GetYaxis()->SetNdivisions(505);
  distFrame->GetYaxis()->SetTitleOffset(1.2);
  distFrame->GetXaxis()->SetTitleOffset(0.9);
  if (setLog)
    distFrame->GetYaxis()->SetRangeUser(0.1, 10*distFrame->GetMaximum());
  else 
    distFrame->GetYaxis()->SetRangeUser(0, 1.8*distFrame->GetMaximum());
}
void fixPullStyle(RooPlot* pullFrame){
  pullFrame->SetTitle(" ");
  pullFrame->GetYaxis()->SetTitle("(Data-Fit)/#sigma");
  pullFrame->SetMinimum(-5);
  pullFrame->SetMaximum(5);
  pullFrame->GetXaxis()->SetTitleSize(0.12);
  pullFrame->GetYaxis()->SetTitleSize(0.1);
  pullFrame->GetXaxis()->SetLabelSize(0.10);
  pullFrame->GetYaxis()->SetLabelSize(0.08);
  pullFrame->GetYaxis()->SetNdivisions(505);
  pullFrame->GetYaxis()->SetTitleOffset(0.4);
  pullFrame->GetXaxis()->SetTitleOffset(1.);
  
  /*
   // make the histograms look better
  padPull->cd();
  pullFrame->SetTitle("");
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetRangeUser(-5,5);
  pullFrame->GetYaxis()->SetLabelSize(0.07);
  pullFrame->GetYaxis()->SetTitleSize(0.1);
  pullFrame->GetYaxis()->SetTitleOffset(0.2);
  pullFrame->GetYaxis()->SetTitleFont(42);
  pullFrame->GetYaxis()->CenterTitle(kTRUE);
      
  pullFrame->GetXaxis()->SetLabelSize(0.1);
  pullFrame->GetXaxis()->CenterTitle(kTRUE);

  pullFrame->GetXaxis()->SetTitleSize(0.15);
  pullFrame->GetXaxis()->SetTitleFont(42);
  pullFrame->GetXaxis()->SetTitleOffset(0.7);
 */
}

TLatex *AliceText (bool ispO) {
    TLatex *textAlice = new TLatex();
    textAlice->SetNDC();
    textAlice->SetTextAlign(12);
    textAlice->SetTextFont(43);
    textAlice->SetTextSize(25); // Size in pixel height
    textAlice->DrawLatex(0.55, 0.92, Form("ALICE %s @#sqrt{s} = %s TeV", ispO?"pO":"OO", ispO?"9.62":"5.36"));
    return textAlice;
}

TLatex *varLatex (RooWorkspace *ws, map<string, string> parIni, double chi2ndf, bool fitMass, bool fitTauz, bool fitTauzRes, bool fitTauzBkg, float xText, float yText) {
  TLatex *textVar = new TLatex();
  textVar->SetNDC();
  textVar->SetTextAlign(12);
  textVar->SetTextFont(43);
  textVar->SetTextSize(17); // Size in pixel height
  textVar->DrawLatex(xText, yText, Form("#chi^{2}/ndf = %.2f", chi2ndf));
  yText = yText-0.04;
  
  for (auto it = parIni.cbegin(); it != parIni.cend(); ++it) {
    cout<<"[INFO] checking variable "<<it->first.c_str()<<" = "<<endl;
    if(it->first.find("fit")!=std::string::npos || it->first.find("model")!=std::string::npos) continue;
    if (ws->var(it->first.c_str())->isConstant()) continue;
    if (fitMass && fitTauz) {
      if (!(it->first.find("_tauzMass")!=std::string::npos || it->first.find("_tauzSig")!=std::string::npos)) //include the lambda for the tauz nonprompt
	continue;
      
    }
    else if (fitMass && !fitTauz) {
      if (!(it->first.find("_mass")!=std::string::npos)) continue;
    }
    else if (!fitMass && fitTauz) {
      if (fitTauzRes) {
	if (!(it->first.find("_tauzRes")!=std::string::npos)) continue; 
      }
      else if(fitTauzBkg) {
	if (!(it->first.find("_tauzBkg")!=std::string::npos)) continue; 
      }
      
      else if (!(it->first.find("_tauzSig")!=std::string::npos)) continue; 
    }
    if (it->first.find("si2s")!=std::string::npos) continue;
    
    if (ws->var(it->first.c_str())->getError()==0) continue;
    cout<<"[INFO] writing variable "<<it->first.c_str()<<" = ";
    cout<<ws->var(it->first.c_str())->getValV();
    cout<<" err = "<<ws->var(it->first.c_str())->getError()<<endl;
    
    textVar->DrawLatex(xText, yText, Form("%s = %g #pm  %g", varFancyLabel(it->first.c_str()).c_str(), ws->var(it->first.c_str())->getValV(), ws->var(it->first.c_str())->getError()));
    yText = yText-0.04;
  }
  return textVar;
}

TLatex *cutLatex (bool ispO, struct KinCuts cutVector, float xText, float yText) {
  TLatex *textCut = new TLatex();
  textCut->SetNDC();
  textCut->SetTextAlign(12);
  textCut->SetTextFont(43);
  textCut->SetTextSize(20); // Size in pixel height
  textCut->DrawLatex(xText, yText, Form("%g < p_{T} < %g GeV", cutVector.pt.Min,cutVector.pt.Max)); yText = yText-0.04;
  textCut->DrawLatex(xText, yText, Form("%g < y < %g", -1*cutVector.rap.Max,-1*cutVector.rap.Min)); yText = yText-0.04;
  textCut->DrawLatex(xText, yText, Form("%g < #chi^{2}_{matching} < %g", cutVector.chi2.Min,cutVector.chi2.Max)); yText = yText-0.04;
  if (!ispO)
    textCut->DrawLatex(xText, yText, Form("%d < cent < %d %%", cutVector.cent.Start ,cutVector.cent.End)); //yText = yText-0.04;
  
    return textCut;
}

TLatex* cutTextResult(bool ispO, string axisName, float xText, float yText, float incMinCent, float incMaxCent, float incMinPt, float incMaxPt, float incMaxRap, float incMinRap, float incMinChi2, float incMaxChi2, bool diffChi2) {
    TLatex* textCut = new TLatex();
  textCut->SetNDC();
  textCut->SetTextAlign(12);
  textCut->SetTextFont(43);
  textCut->SetTextSize(25); // Size in pixel height

  if (!(axisName.find("centrality")!=std::string::npos) && !ispO) {
    textCut->DrawLatex(xText, yText, Form("%d < cent < %d %%", (int) incMinCent, (int) incMaxCent));
    yText = yText-0.04;
  }
  if (!(axisName.find("pt")!=std::string::npos)){
    if (incMinPt==0)
      textCut->DrawLatex(xText, yText, Form("p_{T} < %g GeV", incMaxPt));
    else
      textCut->DrawLatex(xText, yText, Form("%g < p_{T} < %g GeV", incMinPt, incMaxPt));
    yText = yText-0.04;
  }
  if (!(axisName.find("rap")!=std::string::npos)) {
    textCut->DrawLatex(xText, yText, Form("%g < y < %g", -1*incMaxRap, -1*incMinRap));
    yText = yText-0.04;
  }
  if (!(axisName.find("chi2")!=std::string::npos) && !diffChi2) {
    textCut->DrawLatex(xText, yText, Form("%g < #chi^{2}_{matching} < %g", incMinChi2, incMaxChi2));
    yText = yText-0.04;
  }
  return textCut;
}
TLegend* makePlotLegend(RooPlot* frame, map<string, vector<string>> legendEntries, double xmin, double ymin, double xmax, double ymax) {
  TLegend *legend_lines = new TLegend(xmin, ymin, xmax, ymax); // Define the legend position and size
  legend_lines->SetBorderSize(0);
  legend_lines->SetFillStyle(0);
  legend_lines->SetTextAlign(12);
  legend_lines->SetTextFont(43);
  legend_lines->SetTextSize(17);
  
  for (const auto& [key, vec] : legendEntries) {
    legend_lines->AddEntry(frame->RooPlot::findObject(key.c_str()), vec[0].c_str(), vec[1].c_str());
    
  }
  return legend_lines;
}
