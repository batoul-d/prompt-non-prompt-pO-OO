// Author: Batoul Diab
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <TStyle.h>
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TFitResult.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <fstream>
#include <string>
#include <sstream>

#include "BuildPDF.C"
#include "plotUtils.h"


//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
map<string, double> SignalExtraction1Fit(map<string, string>& parIni, RooWorkspace* ws, const char *caseName, string rangeLabel, bool ispO, bool isMC, bool fitMass, bool fitTauz, struct KinCuts cutVector) {
  gStyle->SetOptStat(0);
  // make output directory
  string outDirName = Form("output/outputFits_%s%s_%s", ispO?"pO":"OO", isMC?"_MC":"", caseName);
  gSystem->mkdir(outDirName.c_str(), kTRUE);
  
  // bulding the model from input
  BuildPDF(ws, parIni, isMC, fitMass, fitTauz);
  
  if (fitMass && !fitTauz) {
    if (parIni["doIterativeFit"] == "1") {
      std::cout << "[INFO] Applying iterative fitting on background function in 1D mass fits" << std::endl;
      RooFitResult* fitResult_mass = ws->pdf("totPDF_mass")->fitTo(*ws->data("data"), Extended(kTRUE), SumW2Error(true), RooFit::Save());
      RooDataSet* sPlotDs = (RooDataSet*) ws->data("data")->Clone("data_sPlot");
      RooStats::SPlot* sData = new RooStats::SPlot("sData", "sPlot", *sPlotDs, ws->pdf("totPDF_mass"), RooArgList(*ws->var("fJpsi_mass"), *ws->var("fBkg_mass")));
      ws->import(*sPlotDs);
    }
    else {
      RooFitResult* fitResult_mass = ws->pdf("totPDF_mass")->fitTo(*ws->data("data"), Extended(kTRUE), SumW2Error(true), RooFit::Save());
      RooDataSet* sPlotDs = (RooDataSet*) ws->data("data")->Clone("data_sPlot");
      RooStats::SPlot* sData = new RooStats::SPlot("sData", "sPlot", *sPlotDs, ws->pdf("totPDF_mass"), RooArgList(*ws->var("fJpsi_mass"), *ws->var("fBkg_mass")));
      ws->import(*sPlotDs);
    }
  }
  else if (!fitMass && fitTauz) {
    RooPlot* tauzResFrame = ws->var("tauz")->frame();
    ws->data("sPlotDsSig")->plotOn(tauzResFrame, DataError(RooAbsData::SumW2));
    RooHist* hist = (RooHist*) tauzResFrame->getObject(0);  // dataset is usually the first
    double xMaxRes = getMax(hist); //to get the mean of the resolution function
    ws->var("mean_tauzRes")->setVal(xMaxRes);
    ws->var("mean_tauzRes")->setConstant(kTRUE);
    ws->var("tauz")->setRange("neg", ws->var("tauz")->getMin(), xMaxRes);
    RooDataSet* sPlotDsNeg = (RooDataSet*) ws->data("sPlotDsSig")->reduce(Form("tauz < %f", xMaxRes));
    ws->import(*sPlotDsNeg, RooFit::Rename("sPlotDsSigNeg"));
    
    RooFitResult* fitResult_tauzRes = ws->pdf("tauzResPDF")->fitTo(*ws->data("sPlotDsSigNeg"), Extended(kFALSE), SumW2Error(true), RooFit::Save(), Range("neg"));
    cout<<"[INFO] done with tauz resolution and now let's fix the parpameters to fit the bkg"<<endl;
    fixParPDF(ws, fitResult_tauzRes, parIni, ispO, rangeLabel, caseName, 0, 1, 0);
    //cout<<"[INFO] done with fixing the params"<<endl;
   RooFitResult* fitResult_tauzBkg = ws->pdf("tauzBkgPDF")->fitTo(*ws->data("sPlotDsBkg"), Extended(kFALSE), SumW2Error(true), RooFit::Save());
   //cout<<"[INFO] done with tauz bkg fit"<<endl;  
   //RooFitResult* fitResult_tauz = ws->pdf("tauzSigPDF")->fitTo(*ws->data("sPlotDsSig"), Extended(kTRUE), SumW2Error(true), RooFit::Save());
   //cout<<"[INFO] done with the tauz Sig fit"<<endl;
  }
  else if (fitMass && fitTauz) {
    fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, 1, 1, 1);
    RooFitResult* fitResult_tauzMass = ws->pdf("totPDF_2D")->fitTo(*ws->data("data"), Extended(kTRUE), SumW2Error(true), RooFit::Save());
  }

  TCanvas* can = new TCanvas("can","",800,800);
  can->cd();
  can->Divide(1, 2);
  
  TPad *padDist = (TPad*)can->cd(1);//new TPad("padDist","",0,.23,1,1);
  padDist->SetPad(0.0, 0.24, 1.0, 1.0);
  padDist->SetBottomMargin(0.02);

  //can->cd();
  TPad *padPull = (TPad*)can->cd(2);//new TPad("padPull","",0,0,1,.228);
  padPull->SetPad(0.0, 0.0, 1.0, 0.25);
  padPull->SetTopMargin(0.016);
  padPull->SetBottomMargin(0.3);

  map<string, vector<string>> legendEntries;
  
  RooPlot* massFrame = ws->var("mass")->frame();
  RooPlot* tauzFrame = ws->var("tauz")->frame();
  RooPlot* tauzResFrame = ws->var("tauz")->frame();
  RooPlot* tauzBkgFrame = ws->var("tauz")->frame();

  int nPar;
  double chi2ndf = -999;
  
  if (fitMass && !fitTauz) {
    padDist->cd();
    TCanvas* canFitChi2 = new TCanvas("canFitChi2", "Fit chi2 evolution", 800, 600);
    if (parIni["doIterativeFit"] == "1") {
      // plot the fit chi2 as a function of the iteration steps
      std::vector<double> vChi2ndfResults;

        // prepare iterative fitting
        Double_t chi2ndf = -1.;
        nPar = ws->pdf("totPDF_mass")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
        int nCoeffs = 0;
        for (int i = 0; ; ++i) {
            RooRealVar* coeff = ws->var(Form("c%d_mass", i));
            if (!coeff) break;
            nCoeffs++;
        }

        // start itterative fitting
        // define the max fitting chi2 you accept
        double chi2ndfMax = 2.0;
        for (Int_t i = 0; i < nCoeffs; i++) {
          if (chi2ndf < 0. || chi2ndf > chi2ndfMax || std::isnan(chi2ndf)) {
            RooRealVar* coeff = ws->var(Form("c%d_mass", i));
            if (coeff) coeff->setConstant(kFALSE);
            ws->pdf("totPDF_mass")->fitTo(*ws->data("data"),Extended(kTRUE),SumW2Error(true),RooFit::Save());

            // avoid drawing on top of each other for every iteration
            RooPlot* tmpFrame = ws->var("mass")->frame();
            ws->data("data")->plotOn(tmpFrame, Name("data"));
            ws->pdf("totPDF_mass")->plotOn(tmpFrame, Name("background_mass"), Components(RooArgSet(*ws->pdf("bkgPDF_mass"))),DrawOption("F"), FillColor(kGray), LineColor(kGray));
            ws->pdf("totPDF_mass")->plotOn(tmpFrame, Name("signalPsi2s_mass"), Components(RooArgSet(*ws->pdf("psi2sPDF_mass"))),DrawOption("L"), LineColor(kGreen+4));
            ws->pdf("totPDF_mass")->plotOn(tmpFrame, Name("signalJpsi_mass"), Components(RooArgSet(*ws->pdf("jpsiPDF_mass"))),DrawOption("L"), LineColor(kGreen+2));
            ws->pdf("totPDF_mass")->plotOn(tmpFrame, Name("total_mass"), LineColor(kRed));

            nPar = ws->pdf("totPDF_mass")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
            chi2ndf = tmpFrame->chiSquare(nPar);
            std::cout << "chi2ndf = " << chi2ndf << std::endl;
            vChi2ndfResults.push_back(chi2ndf);
            // delete tmpFrame; ?
          }
          else break;
          if (i == nCoeffs - 1 && (chi2ndf < 0. || chi2ndf > chi2ndfMax)) {
            std::cout << "WARNING. Fitting failed. No convergence" << std::endl;
          }
        }

        // Visualise the fit chi2 evolution for the iterative fitting
        int n = vChi2ndfResults.size();
        std::vector<double> xVals(n);
        for (int i = 0; i < n; ++i) { xVals[i] = i; }
        TGraph* gr = new TGraph(n, xVals.data(), vChi2ndfResults.data());
        gr->SetTitle("Fit chi2 vs step index;step index (i);fit chi2");
        gr->SetMarkerStyle(20);
        canFitChi2->cd();
        gr->Draw("APL");

      padDist->cd();
      massFrame = ws->var("mass")->frame();
      ws->data("data")->plotOn(massFrame, Name("data")); legendEntries["data"] = {"data","P"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("background_mass"), Components(RooArgSet(*ws->pdf("bkgPDF_mass"))),DrawOption("F"), FillColor(kGray), LineColor(kGray)); legendEntries["background_mass"] = {"Background", "F"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalPsi2s_mass"), Components(RooArgSet(*ws->pdf("psi2sPDF_mass"))),DrawOption("L"), LineColor(kGreen+4)); //legendEntries["signalPsi2s_mass"] = {"#psi(2S) signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalJpsi_mass"), Components(RooArgSet(*ws->pdf("jpsiPDF_mass"))),DrawOption("L"), LineColor(kGreen+2)); legendEntries["signalJpsi_mass"] = {"J/#psi signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("total_mass"), LineColor(kRed)); legendEntries["total_mass"] = {"total fit","L"};
    }
    else {
      ws->data("data")->plotOn(massFrame, Name("data")); legendEntries["data"] = {"data","P"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("background_mass"), Components(RooArgSet(*ws->pdf("bkgPDF_mass"))),DrawOption("F"), FillColor(kGray), LineColor(kGray)); legendEntries["background_mass"] = {"Background", "F"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalPsi2s_mass"), Components(RooArgSet(*ws->pdf("psi2sPDF_mass"))),DrawOption("L"), LineColor(kGreen+4)); //legendEntries["signalPsi2s_mass"] = {"#psi(2S) signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalJpsi_mass"), Components(RooArgSet(*ws->pdf("jpsiPDF_mass"))),DrawOption("L"), LineColor(kGreen+2)); legendEntries["signalJpsi_mass"] = {"J/#psi signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("total_mass"), LineColor(kRed)); legendEntries["total_mass"] = {"total fit","L"};
    }

    RooHist *hpull = massFrame->pullHist();
    RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution"));
    pullFrame->addPlotable(hpull,"P");

    nPar = ws->pdf("totPDF_mass")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = massFrame->chiSquare(nPar);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;
    
    ws->data("data")->plotOn(massFrame);
    padDist->cd();
    fixFrameStyle(massFrame, false);
    massFrame->Draw();
    TLatex* textVar = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 0, 0.57, 0.8); textVar->Draw("same");
    TLatex* textCut = cutLatex(ispO, cutVector, 0.25,0.6); textCut->Draw("same");
    TLegend* leg = makePlotLegend(massFrame, legendEntries, 0.25, 0.7, 0.5, 0.85); leg->Draw("same");
    TLatex *textAlice = AliceText(ispO); textAlice->Draw("same");
    padPull->cd();
    fixPullStyle(pullFrame);
    pullFrame->Draw();
    TLine* linePull = new TLine(ws->var("mass")->getMin(), 0, ws->var("mass")->getMax(),0);
    linePull->SetLineColor(kRed); linePull->SetLineStyle(2); linePull->Draw("same");
    //std::string pdfPath = "";
    can->SaveAs(Form("%s/massFit1D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
    if (parIni["doIterativeFit"] == "1") {
      canFitChi2->SaveAs(Form("%s/massFit1D_%s_fitChi2Trend.pdf", outDirName.c_str(), rangeLabel.c_str()));
    }
  }
  else if (fitTauz && !fitMass) {
    //tauz resolution first
    padDist->cd();
    ws->Print(); 
    ws->data("sPlotDsSigNeg")->plotOn(tauzResFrame, Name("sPlotDsSigNeg"), DataError(RooAbsData::SumW2)); legendEntries["sPlotDsSigNeg"] = {"sPlot signal-like data","P"};
    int nGauss = 0;
    //int gausColor[] = {kBlue, kGreen+2, kPink+4, kPink+4};
    if (parIni["model_tauzRes"]=="Gauss3") nGauss = 3;
    else if (parIni["model_tauzRes"]=="Gauss4") nGauss = 4;
    for (int i=0; i<nGauss; i++) {
      ws->pdf("tauzResPDF")->plotOn(tauzResFrame, Components(RooArgSet(*ws->pdf(Form("gauss%d_tauzRes", i)))), Name(Form("gauss%d_tauzRes", i)), LineColor(kBlue+i), LineStyle(2+i)); legendEntries[Form("gauss%d_tauzRes", i)] = {Form("gaussian %d", i),"L"};
    } 
    
    ws->pdf("tauzResPDF")->plotOn(tauzResFrame, Name("tauzResPDF"), LineColor(kRed)); legendEntries["tauzResPDF"] = {"total fit","L"};
    
    RooHist* hpull = tauzResFrame->pullHist();
    RooPlot* pullFrame = ws->var("tauz")->frame(Title("Pull Distribution"));
    pullFrame->addPlotable(hpull,"P");

    nPar = ws->pdf("tauzResPDF")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = tauzResFrame->chiSquare(nPar);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;
    
    ws->data("sPlotDsSigNeg")->plotOn(tauzResFrame, DataError(RooAbsData::SumW2));
    //padDist->cd();
    fixFrameStyle(tauzResFrame, true);
    tauzResFrame->Draw();
    TLatex* textVarRes = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 1, 0, 0.57, 0.85); textVarRes->Draw("same");
    TLatex* textCutRes = cutLatex(ispO, cutVector, 0.2, 0.8); textCutRes->Draw("same");
    TLegend* legRes = makePlotLegend(tauzResFrame, legendEntries, 0.15, 0.5, 0.3, 0.65); legRes->Draw("same");

    TLatex *textAlice = AliceText(ispO); //textAlice->Draw("same");
    
    padPull->cd();
    fixPullStyle(pullFrame);
    pullFrame->Draw();
    TLine* linePull = new TLine(ws->var("tauz")->getMin(), 0, ws->var("tauz")->getMax(),0);
    linePull->SetLineColor(kRed); linePull->SetLineStyle(2); linePull->Draw("same");
    //std::string pdfPath = "";
    padDist->SetLogy();
    can->SaveAs(Form("%s/tauzResFit1D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));

    //////// then tauz Bkg
    padDist->cd();
    legendEntries.clear();
    ws->data("sPlotDsBkg")->plotOn(tauzBkgFrame, Name("sPlotDsBkg"), DataError(RooAbsData::SumW2)); legendEntries["sPlotDsBkg"] = {"sPlot background-like data","P"};
    ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzNprBkgPDF"), Components(RooArgSet(*ws->pdf("tauzNprBkgPDF"))),DrawOption("L"), LineColor(kBlue)); legendEntries["tauzNprBkgPDF"] = {"non-prompt-like", "L"};
    ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzPrBkgPDF"), Components(RooArgSet(*ws->pdf("tauzPrBkgPDF"))),DrawOption("L"), LineColor(kGreen+2)); legendEntries["tauzPrBkgPDF"] = {"prompt-like", "L"};
    ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzBkgPDF"), LineColor(kRed)); legendEntries["tauzBkgPDF"] = {"total fit","L"};
    
    RooHist* hpullBkg = tauzBkgFrame->pullHist();
    RooPlot* pullBkgFrame = ws->var("tauz")->frame(Title("Pull Distribution"));
    pullBkgFrame->addPlotable(hpullBkg,"P");

    nPar = ws->pdf("tauzBkgPDF")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = tauzBkgFrame->chiSquare(nPar);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;
    
    ws->data("sPlotDsBkg")->plotOn(tauzBkgFrame, DataError(RooAbsData::SumW2));
    //padDist->cd();
    fixFrameStyle(tauzBkgFrame, true);
    tauzBkgFrame->Draw();
    TLatex* textVarBkg = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 1, 0.57, 0.85); textVarBkg->Draw("same");
    TLatex* textCutBkg = cutLatex(ispO, cutVector, 0.2,0.8); textCutBkg->Draw("same");
    TLegend* legBkg = makePlotLegend(tauzBkgFrame, legendEntries, 0.15, 0.5, 0.3, 0.65); legBkg->Draw("same");
    TLatex *textAliceBkg = AliceText(ispO); //textAlice->Draw("same");
    
    
    padPull->cd();
    fixPullStyle(pullBkgFrame);
    pullBkgFrame->Draw();
    TLine* linePullBkg = new TLine(ws->var("tauz")->getMin(), 0, ws->var("tauz")->getMax(),0);
    linePullBkg->SetLineColor(kRed); linePullBkg->SetLineStyle(2); linePullBkg->Draw("same");
    //std::string pdfPath = "";
    //can->SetLogy();
    can->SaveAs(Form("%s/tauzBkgFit1D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
    
  }
  else if (fitTauz && fitMass) {
    padDist->cd();
    ws->data("data")->plotOn(massFrame, Name("data"), DataError(RooAbsData::SumW2)); legendEntries["data"] = {"data","P"};
    ws->pdf("totPDF_2D")->plotOn(massFrame, Name("background"), Components(RooArgSet(*ws->pdf("tauzMassTotBkgPDF"))),DrawOption("F"), FillColor(kGray), LineColor(kGray)); legendEntries["background"] = {"Background", "F"};
    ws->pdf("totPDF_2D")->plotOn(massFrame, Name("tauzMassPrJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassPrJpsiPDF"))),DrawOption("L"), LineColor(kGreen+2)); legendEntries["tauzMassPrJpsiPDF"] = {"prompt J/#psi signal","L"};
    ws->pdf("totPDF_2D")->plotOn(massFrame, Name("tauzMassNprJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassNprJpsiPDF"))),DrawOption("L"), LineColor(kOrange+2)); legendEntries["tauzMassNprJpsiPDF"] = {"non-prompt J/#psi signal","L"};
    //ws->pdf("totPDF_2D")->plotOn(massFrame, Name("tauzMassTotJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassTotJpsiPDF"))),DrawOption("L"), LineColor(kBlack)); legendEntries["tauzMassTotJpsiPDF"] = {"J/#psi signal","L"};
    //ws->pdf("totPDF_mass")->plotOn(massFrame, Name("tauzMassTotPsi2sPDF"), Components(RooArgSet(*ws->pdf("tauzMassTotPsi2sPDF"))),DrawOption("L"), LineColor(kBlue+4)); legendEntries["tauzMassTotPsi2sPDF"] = {"#psi(2S) signal","L"};
    ws->pdf("totPDF_2D")->plotOn(massFrame, Name("totPDF_2D"), LineColor(kRed)); legendEntries["totPDF_2D"] = {"total fit","L"};
    
    RooHist* hpull = massFrame->pullHist();
    RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution"));
    pullFrame->addPlotable(hpull,"P");

    nPar = ws->pdf("totPDF_2D")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = massFrame->chiSquare(nPar);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;
    
    ws->data("data")->plotOn(massFrame, DataError(RooAbsData::SumW2));
    //padDist->cd();
    fixFrameStyle(massFrame, false);
    massFrame->Draw();
    TLatex* textVar = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 0, 0.57, 0.8); textVar->Draw("same");
    TLatex* textCut = cutLatex(ispO, cutVector, 0.15,0.6); textCut->Draw("same");
    TLegend* leg = makePlotLegend(massFrame, legendEntries, 0.15, 0.7, 0.5, 0.85); leg->Draw("same");

    TLatex *textAlice = AliceText(ispO); textAlice->Draw("same");
    padPull->cd();
    fixPullStyle(pullFrame);
    pullFrame->Draw();
    TLine* linePull = new TLine(ws->var("mass")->getMin(), 0, ws->var("mass")->getMax(),0);
    linePull->SetLineColor(kRed); linePull->SetLineStyle(2); linePull->Draw("same");
    //padDist->SetLogy();
    //std::string pdfPath = "";
    can->SaveAs(Form("%s/massFitProj2D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));

    ///////tauz frame
    padDist->cd();
    legendEntries.clear();
    ws->data("data")->plotOn(tauzFrame, Name("data"), DataError(RooAbsData::SumW2)); legendEntries["data"] = {"data","P"};
    ws->pdf("totPDF_2D")->plotOn(tauzFrame, Name("background"), Components(RooArgSet(*ws->pdf("tauzMassTotBkgPDF"))),DrawOption("F"), FillColor(kGray), LineColor(kGray)); legendEntries["background"] = {"Background", "F"};
    ws->pdf("totPDF_2D")->plotOn(tauzFrame, Name("tauzMassPrJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassPrJpsiPDF"))),DrawOption("L"), LineColor(kGreen+2)); legendEntries["tauzMassPrJpsiPDF"] = {"prompt J/#psi signal","L"};
    ws->pdf("totPDF_2D")->plotOn(tauzFrame, Name("tauzMassNprJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassNprJpsiPDF"))),DrawOption("L"), LineColor(kOrange+2)); legendEntries["tauzMassNprJpsiPDF"] = {"non-prompt J/#psi signal","L"};
    //ws->pdf("totPDF_2D")->plotOn(tauzFrame, Name("tauzMassTotJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassTotJpsiPDF"))),DrawOption("L"), LineColor(kBlack)); legendEntries["tauzMassTotJpsiPDF"] = {"J/#psi signal","L"};
    //ws->pdf("totPDF_mass")->plotOn(tauzFrame, Name("tauzMassTotPsi2sPDF"), Components(RooArgSet(*ws->pdf("tauzMassTotPsi2sPDF"))),DrawOption("L"), LineColor(kBlue+4)); legendEntries["tauzMassTotPsi2sPDF"] = {"#psi(2S) signal","L"};
    ws->pdf("totPDF_2D")->plotOn(tauzFrame, Name("totPDF_2D"), LineColor(kRed)); legendEntries["totPDF_2D"] = {"total fit","L"};
    
    RooHist* hpullTauz = tauzFrame->pullHist();
    RooPlot* pullTauzFrame = ws->var("tauz")->frame(Title("Pull Distribution"));
    pullTauzFrame->addPlotable(hpullTauz,"P");

    nPar = ws->pdf("totPDF_2D")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = massFrame->chiSquare(nPar);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;
    
    ws->data("data")->plotOn(tauzFrame, DataError(RooAbsData::SumW2));
    //padDist->cd();
    fixFrameStyle(tauzFrame, true);
    tauzFrame->Draw();
    TLatex* textVarTauz = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 0, 0.57, 0.8); textVarTauz->Draw("same");
    TLatex* textCutTauz = cutLatex(ispO, cutVector, 0.15,0.6); textCutTauz->Draw("same");
    TLegend* legTauz = makePlotLegend(tauzFrame, legendEntries, 0.15, 0.7, 0.5, 0.85); legTauz->Draw("same");
    TLatex *textAliceTauz = AliceText(ispO); textAliceTauz->Draw("same");
    
    padPull->cd();
    fixPullStyle(pullTauzFrame);
    pullTauzFrame->Draw();
    TLine* linePullTauz = new TLine(ws->var("tauz")->getMin(), 0, ws->var("tauz")->getMax(),0);
    linePullTauz->SetLineColor(kRed); linePullTauz->SetLineStyle(2); linePullTauz->Draw("same");
    //std::string pdfPath = "";
    padDist->SetLogy();
    can->SaveAs(Form("%s/tauzFitProj2D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
  }
  
  map<string, double> results;
  for (auto it = parIni.cbegin(); it != parIni.cend(); ++it) {
    if(it->first.find("fit")!=std::string::npos || it->first.find("model")!=std::string::npos) continue;
    results[it->first.c_str()] = ws->var(it->first.c_str())->getValV();
    results[Form("%s_err",it->first.c_str())] = ws->var(it->first.c_str())->getError();
  }
      
  results["chi2ndf"] = chi2ndf;
  //results["ndf"] = nPar;//fitResult->ndf();
      
  return results;
      
} //end of SignalExtraction Function

//___________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
