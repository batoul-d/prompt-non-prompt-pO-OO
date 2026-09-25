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
map<string, double> SignalExtraction1Fit(map<string, string>& parIni, RooWorkspace* ws, const char *caseName, string rangeLabel, bool ispO, bool isMC, bool fitMass, bool fitTauz, bool fitNpTauz, bool fitTauzBkg, struct KinCuts cutVector) {
  gStyle->SetOptStat(0);
  // make output directory
  string outDirName = Form("output/outputFits_%s%s_%s", ispO?"pO":"OO", isMC?"_MC":"", caseName);
  gSystem->mkdir(outDirName.c_str(), kTRUE);
  
  // bulding the model from input
  BuildPDF(ws, parIni, isMC, fitMass, fitTauz);
  cout << "PDFs correctly built" << endl;

  // Perform the tauz background fit?
  // NOTE: this is also set to 'false' and 'true' as part of the 2-step resolution fit. Might have to be changed below too.
  if (!fitTauzBkg) cout << "[INFO] NOTE: RUNNING WITHOUT TAUZ BACKGROUND FIT -------------------------------------------" << endl;
  
  // Number of bins to be drawn (does not affect fitting)
  int nBins = 200;

  // Set range for plotting and fitting the mass
  double massMin = 2.4;
  double massMax = 4.0;

  // Set range for plotting the tau_z
  double tauzMin = -0.02;
  double tauzMax = 0.02;

  int nParTauzRes = -1;

  // Set range for fitting the tau_z
  // TODO: Remove this again!!!
  // Double_t tauzFittingMin = -0.05;
  // Double_t tauzFittingMax = 0.05;

  // TODO: remove this
  // Define the named RooFit range
  // ws->var("tauz")->setRange("tauzFit", tauzFittingMin, tauzFittingMax);

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
      if (!isMC) { RooDataSet* sPlotDs = (RooDataSet*) ws->data("data")->Clone("data_sPlot");
        RooStats::SPlot* sData = new RooStats::SPlot("sData", "sPlot", *sPlotDs, ws->pdf("totPDF_mass"), RooArgList(*ws->var("fJpsi_mass"), *ws->var("fBkg_mass")));
        ws->import(*sPlotDs);
      }
    }
  }
  else if (!fitMass && fitTauz && !isMC) {
    RooPlot* tauzResFrame = ws->var("tauz")->frame(Range(tauzMin, 0), Bins(nBins));
    ws->data("sPlotDsSig")->plotOn(tauzResFrame, DataError(RooAbsData::SumW2));
    RooHist* hist = (RooHist*) tauzResFrame->getObject(0);
    ws->var("xMaxRes")->setVal(getMax(hist));
    double xMaxRes = ws->var("xMaxRes")->getVal(); // to get the mean of the resolution function
    // double xMaxRes = 0.;
    cout << " [INFO] xMaxRes = " << xMaxRes << " ns" << endl;
    ws->var("mean_tauzRes")->setVal(xMaxRes);
    double bias = 0.3;
    ws->var("mean_tauzRes")->setMin(ws->var("mean_tauzRes")->getVal() - bias * std::abs(ws->var("mean_tauzRes")->getVal()));
    ws->var("mean_tauzRes")->setMax(ws->var("mean_tauzRes")->getVal() + bias * std::abs(ws->var("mean_tauzRes")->getVal()));
    // ws->var("mean_tauzRes")->setConstant(kTRUE);

    // Step 1: resolution fit on negative tauz range
    if (!fitNpTauz) {
      ws->var("tauz")->setRange("neg", ws->var("tauz")->getMin(), xMaxRes);
      RooDataSet* sPlotDsNeg = (RooDataSet*) ws->data("sPlotDsSig")->reduce(Form("tauz < %f", xMaxRes));
      ws->import(*sPlotDsNeg, RooFit::Rename("sPlotDsSigNeg"));
      RooFitResult* fitResult_tauzRes = ws->pdf("tauzResPDF")->fitTo(*ws->data("sPlotDsSigNeg"), Extended(kFALSE), SumW2Error(true), RooFit::Save(), Range("neg"));
      cout << "[INFO] I found this value for the ws->var(mean_tauzRes)->getVal()" << endl;
      nParTauzRes = ws->pdf("tauzResPDF")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
      fixParPDF(ws, fitResult_tauzRes, parIni, ispO, rangeLabel, caseName, 0, 1, 1, 0, 0);
      // Save params, as is done for the 2D fit with function
      cout << "[INFO] done with fixing the params for tauz resolution 1st step" << endl;
    }
    // Step 2: include signal fit of non-prompt decay (initialised from MC)
    else { 
      // Take parameters from 1st step
      fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, 0, 1, 1, 0, 0);
      RooFitResult* fitResult_tauz = ws->pdf("tauzSigPDF")->fitTo(*ws->data("sPlotDsSig"), Extended(kFALSE), SumW2Error(true), RooFit::Save(), Range(tauzMin, tauzMax));
      cout << "[INFO] done with the tauz Sig fit" << endl;
      // Before fixing the parameters, we should get the number of free parameters used in the fit for the chi2 later
      nParTauzRes = ws->pdf("tauzSigPDF")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
      // TODO: fix parameters and bias the lambda (with special flag)
      fixParPDF(ws, fitResult_tauz, parIni, ispO, rangeLabel, caseName, 0, 1, 0, 0, 1);
      cout << "[INFO] done with fixing the params for tauz resolution 2nd step" << endl;
    }
    
    // Do the background fit using the resolution
    if (fitTauzBkg) { 
      cout << "[INFO] done with tauz resolution and now let's fix the parpameters to fit the bkg" << endl;
      RooFitResult* fitResult_tauzBkg = ws->pdf("tauzBkgPDF")->fitTo(*ws->data("sPlotDsBkg"), Extended(kFALSE), SumW2Error(true), RooFit::Save()); }
      cout << "[INFO] done with tauz bkg fit" << endl;
  }
  else if (!fitMass && fitTauz && isMC) {
    // Fit the non-prompt signal: exponential decay convoluted with the detector resolution (MC)
    RooFitResult* fitResult_tauzNpr = ws->pdf("tauzNprSigPDF")->fitTo(*ws->data("data"), Extended(kFALSE), SumW2Error(true), RooFit::Save());
    cout << "[INFO] Done with MC non-prompt tauz fit" << endl;
    // fitResult_tauzNpr->Print();
  }
  else if (fitMass && fitTauz && !isMC) {

    // Fix mass-shape parameters from the 1D mass fit
    fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, true, false, false, false, false);
    // Fix tauz-resolution parameters from the 1D tauz fit
    fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, false, true, false, false, false);
    // Fix tauz-background parameters from the 1D tauz fit
    if (!isMC && fitTauzBkg) { fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, false, false, false, true, false); }

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
  
  RooPlot* massFrame = ws->var("mass")->frame(Range(massMin, massMax), Bins(nBins));
  RooPlot* tauzFrame = ws->var("tauz")->frame(Range(tauzMin, tauzMax), Bins(nBins));
  RooPlot* tauzResFrame = ws->var("tauz")->frame(Range(tauzMin, tauzMax), Bins(nBins));
  RooPlot* tauzBkgFrame = ws->var("tauz")->frame(Range(tauzMin, tauzMax), Bins(nBins));

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
      massFrame = ws->var("mass")->frame(Range(massMin, massMax), Bins(nBins));
      ws->data("data")->plotOn(massFrame, Name("data")); legendEntries["data"] = {"data","P"};
      if (!isMC) { ws->pdf("totPDF_mass")->plotOn(massFrame, Name("background_mass"), Components(RooArgSet(*ws->pdf("bkgPDF_mass"))),DrawOption("F"), FillColor(kGray), LineColor(kGray)); legendEntries["background_mass"] = {"Background", "F"}; }
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalPsi2s_mass"), Components(RooArgSet(*ws->pdf("psi2sPDF_mass"))),DrawOption("L"), LineColor(kGreen+4)); //legendEntries["signalPsi2s_mass"] = {"#psi(2S) signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalJpsi_mass"), Components(RooArgSet(*ws->pdf("jpsiPDF_mass"))),DrawOption("L"), LineColor(kGreen+2)); legendEntries["signalJpsi_mass"] = {"J/#psi signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("total_mass"), LineColor(kRed)); legendEntries["total_mass"] = {"total fit","L"};
    }
    else {
      ws->data("data")->plotOn(massFrame, Name("data")); legendEntries["data"] = {"data","P"};
      if (!isMC) { ws->pdf("totPDF_mass")->plotOn(massFrame, Name("background_mass"), Components(RooArgSet(*ws->pdf("bkgPDF_mass"))),DrawOption("F"), FillColor(kGray), LineColor(kGray)); legendEntries["background_mass"] = {"Background", "F"}; }
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalPsi2s_mass"), Components(RooArgSet(*ws->pdf("psi2sPDF_mass"))),DrawOption("L"), LineColor(kGreen+4)); //legendEntries["signalPsi2s_mass"] = {"#psi(2S) signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalJpsi_mass"), Components(RooArgSet(*ws->pdf("jpsiPDF_mass"))),DrawOption("L"), LineColor(kGreen+2)); legendEntries["signalJpsi_mass"] = {"J/#psi signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("total_mass"), LineColor(kRed)); legendEntries["total_mass"] = {"total fit","L"};
    }

    RooHist *hpull = massFrame->pullHist();
    RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution"), Range(massMin, massMax));
    pullFrame->addPlotable(hpull,"P");

    nPar = ws->pdf("totPDF_mass")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = massFrame->chiSquare(nPar);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;
    
    ws->data("data")->plotOn(massFrame);
    padDist->cd();
    fixFrameStyle(massFrame, false);
    massFrame->Draw();
    TLatex* textVar = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 0, 0.57, 0.8);
    TLatex* textCut = cutLatex(ispO, cutVector, 0.25,0.6);
    TLegend* leg = makePlotLegend(massFrame, legendEntries, 0.25, 0.7, 0.5, 0.85); leg->Draw("same");
    TLatex *textAlice = AliceText(ispO);
    padPull->cd();
    fixPullStyle(pullFrame);
    pullFrame->Draw();
    TLine* linePull = new TLine(massMin, 0, massMax, 0);
    linePull->SetLineColor(kRed); linePull->SetLineStyle(2); linePull->Draw("same");
    //std::string pdfPath = "";
    can->SaveAs(Form("%s/massFit1D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
    if (parIni["doIterativeFit"] == "1") {
      canFitChi2->SaveAs(Form("%s/massFit1D_%s_fitChi2Trend.pdf", outDirName.c_str(), rangeLabel.c_str()));
    }
  }
  else if (fitTauz && !fitMass && !isMC) {
    // Draw tauz resolution first
    padDist->cd();
    ws->Print(); 
    // change this line to display only the negative data points (or all data points)
    ws->data("sPlotDsSig")->plotOn(tauzResFrame, Name("sPlotDsSig"), DataError(RooAbsData::SumW2)); legendEntries["sPlotDsSigNeg"] = {"sPlot signal-like data","P"};
    int nGauss = 0;
    if (parIni["model_tauzRes"]=="Gauss2") nGauss = 2;
    else if (parIni["model_tauzRes"]=="Gauss3" || parIni["model_tauzRes"]=="Gauss3VarMeans") nGauss = 3;
    else if (parIni["model_tauzRes"]=="Gauss4") nGauss = 4;
    if (!fitNpTauz) {
      for (int i=0; i<nGauss; i++) {
        ws->pdf("tauzResPDF")->plotOn(tauzResFrame, Components(RooArgSet(*ws->pdf(Form("gauss%d_tauzRes", i)))), Name(Form("gauss%d_tauzRes", i)), LineColor(kBlue+i), LineStyle(2+i)); legendEntries[Form("gauss%d_tauzRes", i)] = {Form("gaussian %d", i),"L"};
      } 
      ws->pdf("tauzResPDF")->plotOn(tauzResFrame, Name("tauzResPDF"), LineColor(kRed)); legendEntries["tauzResPDF"] = {"total fit","L"};
    }
    // Draw non-prompt decay x resolution fit and components
    else {
      // ws->pdf("tauzResPDF")->plotOn(tauzResFrame, Name("tauzResPDF"), LineColor(kBlue), LineStyle(kDashed)); legendEntries["tauzResPDF"] = {"resolution fit","L"};
      ws->pdf("tauzSigPDF")->plotOn(tauzResFrame, Components("tauzNprSigPDF"), Name("tauzNprSigPDF"), LineColor(kOrange+2), DrawOption("L")); legendEntries["tauzNprSigPDF"] = {"Non-prompt", "L"};
      ws->pdf("tauzSigPDF")->plotOn(tauzResFrame, Components("tauzPrSigPDF"), Name("tauzPrSigPDF"), LineColor(kGreen+2), DrawOption("L")); legendEntries["tauzPrSigPDF"] = {"Prompt", "L"};
      for (int i=0; i<nGauss; i++) {
        ws->pdf("tauzSigPDF")->plotOn(tauzResFrame, Components(RooArgSet(*ws->pdf(Form("gauss%d_tauzRes", i)))), Name(Form("gauss%d_tauzRes", i)), LineColor(kBlue+i), LineStyle(2+i)); legendEntries[Form("gauss%d_tauzRes", i)] = {Form("exp x Gaussian %d", i),"L"};
      }
      ws->pdf("tauzSigPDF")->plotOn(tauzResFrame, Name("tauzSigPDF"), LineColor(kRed), DrawOption("L")); legendEntries["tauzSigPDF"] = {"Pr + Np signal fit", "L"};
    }

    RooHist* hpull = tauzResFrame->pullHist();
    RooPlot* pullFrame = ws->var("tauz")->frame(Title("Pull Distribution"), Range(tauzMin, tauzMax));
    pullFrame->addPlotable(hpull,"P");

    // Calculate chi2 with number of parameters from before fixing them
    /*
    RooPlot* tauzResChi2Frame = ws->var("tauz")->frame(Range(tauzMin, 0), Bins(nBins/2));
    ws->data("sPlotDsSigNeg")->plotOn(tauzResChi2Frame, DataError(RooAbsData::SumW2));
    ws->pdf("tauzResPDF")->plotOn(tauzResChi2Frame);
    chi2ndf = tauzResChi2Frame->chiSquare(nParTauzRes);
    std::cout << "tauzRes chi2/ndf old = " << tauzResChi2Frame->chiSquare(nParTauzRes/6) << std::endl;
    std::cout << "[INFO] tauzRes floating parameters = " << nParTauzRes << std::endl;
    std::cout << "[INFO] tauzRes chi2/ndf = " << chi2ndf << std::endl;
    */

    // previous calculation (doesn't take into account the range?)
    // nPar = ws->pdf("tauzResPDF")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = tauzResFrame->chiSquare(nParTauzRes);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;

    // ws->data("sPlotDsSigNeg")->plotOn(tauzResFrame, DataError(RooAbsData::SumW2));
    //padDist->cd();
    fixFrameStyle(tauzResFrame, true);
    tauzResFrame->Draw();
    TLatex* textVarRes = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 1, 0, 0.57, 0.85);
    TLatex* textCutRes = cutLatex(ispO, cutVector, 0.2, 0.8);
    TLegend* legRes = makePlotLegend(tauzResFrame, legendEntries, 0.15, 0.5, 0.3, 0.65); legRes->Draw("same");
    TLatex *textAlice = AliceText(ispO);
    
    padPull->cd();
    fixPullStyle(pullFrame);
    pullFrame->Draw();
    // TLine* linePull = new TLine(ws->var("tauz")->getMin(), 0, ws->var("tauz")->getMax(),0);
    TLine* linePull = new TLine(tauzMin, 0, tauzMax, 0);
    linePull->SetLineColor(kRed); linePull->SetLineStyle(2); linePull->Draw("same");
    padDist->SetLogy();
    if (fitNpTauz) cout << "HIJ IS AAN" << endl;
    can->SaveAs(Form("%s/tauzResFit1D%s_%s.pdf", outDirName.c_str(), fitNpTauz?"_npTauz":"", rangeLabel.c_str()));
    // ROOT file containing the editable objects
    /*
    TFile *fOut = new TFile(Form("%s/tauzResFit1D%s_%s.root", outDirName.c_str(), fitNpTauz?"_npTauz":"", rangeLabel.c_str()), "RECREATE");
    can->Write("can");
    tauzResFrame->Write("tauFrame");
    pullFrame->Write("pullFrame");
    ws->Write("ws");
    fOut->Close();
    */

    //////// then tauz Bkg
    if (fitTauzBkg) {
      padDist->cd();
      legendEntries.clear();
      ws->data("sPlotDsBkg")->plotOn(tauzBkgFrame, Name("sPlotDsBkg"), DataError(RooAbsData::SumW2)); legendEntries["sPlotDsBkg"] = {"sPlot background-like data","P"};
      ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzNprBkgPDF"), Components(RooArgSet(*ws->pdf("tauzNprBkgPDF"))), DrawOption("L"), LineColor(kBlue)); legendEntries["tauzNprBkgPDF"] = {"non-prompt-like", "L"};
      ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzPrBkgPDF"), Components(RooArgSet(*ws->pdf("tauzPrBkgPDF"))), DrawOption("L"), LineColor(kGreen+2)); legendEntries["tauzPrBkgPDF"] = {"prompt-like", "L"};
      ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzBkgPDF"), LineColor(kRed)); legendEntries["tauzBkgPDF"] = {"total fit","L"};
      
      RooHist* hpullBkg = tauzBkgFrame->pullHist();
      RooPlot* pullBkgFrame = ws->var("tauz")->frame(Title("Pull Distribution"), Range(tauzMin, tauzMax));
      pullBkgFrame->addPlotable(hpullBkg,"P");

      nPar = ws->pdf("tauzBkgPDF")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
      chi2ndf = tauzBkgFrame->chiSquare(nPar);
      std::cout << "fit chi2 = " << chi2ndf << std::endl;
      
      ws->data("sPlotDsBkg")->plotOn(tauzBkgFrame, DataError(RooAbsData::SumW2));
      //padDist->cd();
      fixFrameStyle(tauzBkgFrame, true);
      tauzBkgFrame->Draw();
      TLatex* textVarBkg = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 1, 0.57, 0.85); // textVarBkg->Draw("same");
      TLatex* textCutBkg = cutLatex(ispO, cutVector, 0.2,0.8); // textCutBkg->Draw("same");
      TLegend* legBkg = makePlotLegend(tauzBkgFrame, legendEntries, 0.15, 0.5, 0.3, 0.65); legBkg->Draw("same");
      TLatex *textAliceBkg = AliceText(ispO); // textAlice->Draw("same");
      
      
      padPull->cd();
      fixPullStyle(pullBkgFrame);
      pullBkgFrame->Draw();
      // TLine* linePullBkg = new TLine(ws->var("tauz")->getMin(), 0, ws->var("tauz")->getMax(),0);
      TLine* linePullBkg = new TLine(tauzMin, 0, tauzMax, 0);
      linePullBkg->SetLineColor(kRed); linePullBkg->SetLineStyle(2); linePullBkg->Draw("same");
      //std::string pdfPath = "";
      //can->SetLogy();
      can->SaveAs(Form("%s/tauzBkgFit1D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
    }
  }
  else if (fitTauz && !fitMass && isMC) {
    padDist->cd();
    ws->data("data")->plotOn(tauzFrame, Name("data"), DataError(RooAbsData::SumW2));
    legendEntries["data"] = {"MC non-prompt signal", "P"};
    ws->pdf("tauzNprSigPDF")->plotOn(tauzFrame, Name("tauzNprSigPDF"), LineColor(kRed));
    legendEntries["tauzNprSigPDF"] = {"non-prompt signal fit", "L"};

    // Plot the resolution-convolved signal
    RooHist* hpull = tauzFrame->pullHist();
    RooPlot* pullFrame = ws->var("tauz")->frame(Title("Pull Distribution"), Range(tauzMin, tauzMax));
    pullFrame->addPlotable(hpull, "P");

    nPar = ws->pdf("tauzNprSigPDF")->getParameters(*ws->data("data"))->selectByAttrib("Constant", kFALSE)->getSize();
    chi2ndf = tauzFrame->chiSquare(nPar);
    cout << "[INFO] MC non-prompt tauz fit chi2 = " << chi2ndf << endl;

    ws->data("data")->plotOn(tauzFrame, DataError(RooAbsData::SumW2));
    fixFrameStyle(tauzFrame, true); tauzFrame->Draw();
    TLatex* textVar = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 0, 0.57, 0.8); // textVar->Draw("same");

    TLegend* leg = makePlotLegend(tauzFrame, legendEntries, 0.15, 0.7, 0.3, 0.85); leg->Draw("same");
    padPull->cd();
    fixPullStyle(pullFrame);
    pullFrame->Draw();
    TLine* linePull = new TLine(tauzMin, 0, tauzMax, 0);
    linePull->SetLineColor(kRed);
    linePull->SetLineStyle(2);
    linePull->Draw("same");
    padDist->SetLogy();

    can->SaveAs(Form("%s/tauzNprMCFit1D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
  }
  else if (fitTauz && fitMass && !isMC) {
    padDist->cd();
    ws->data("data")->plotOn(massFrame, Name("data"), DataError(RooAbsData::SumW2)); legendEntries["data"] = {"data","P"};
    ws->pdf("totPDF_2D")->plotOn(massFrame, Name("background"), Components(RooArgSet(*ws->pdf("tauzMassTotBkgPDF"))),DrawOption("F"), FillColor(kGray), LineColor(kGray)); legendEntries["background"] = {"Background", "F"};
    ws->pdf("totPDF_2D")->plotOn(massFrame, Name("tauzMassPrJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassPrJpsiPDF"))),DrawOption("L"), LineColor(kGreen+2)); legendEntries["tauzMassPrJpsiPDF"] = {"prompt J/#psi signal","L"};
    ws->pdf("totPDF_2D")->plotOn(massFrame, Name("tauzMassNprJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassNprJpsiPDF"))),DrawOption("L"), LineColor(kOrange+2)); legendEntries["tauzMassNprJpsiPDF"] = {"non-prompt J/#psi signal","L"};
    //ws->pdf("totPDF_2D")->plotOn(massFrame, Name("tauzMassTotJpsiPDF"), Components(RooArgSet(*ws->pdf("tauzMassTotJpsiPDF"))),DrawOption("L"), LineColor(kBlack)); legendEntries["tauzMassTotJpsiPDF"] = {"J/#psi signal","L"};
    //ws->pdf("totPDF_mass")->plotOn(massFrame, Name("tauzMassTotPsi2sPDF"), Components(RooArgSet(*ws->pdf("tauzMassTotPsi2sPDF"))),DrawOption("L"), LineColor(kBlue+4)); legendEntries["tauzMassTotPsi2sPDF"] = {"#psi(2S) signal","L"};
    ws->pdf("totPDF_2D")->plotOn(massFrame, Name("totPDF_2D"), LineColor(kRed)); legendEntries["totPDF_2D"] = {"total fit","L"};
    
    RooHist* hpull = massFrame->pullHist();
    RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution"), Range(massMin, massMax));
    pullFrame->addPlotable(hpull,"P");

    nPar = ws->pdf("totPDF_2D")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = massFrame->chiSquare(nPar);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;
    
    ws->data("data")->plotOn(massFrame, DataError(RooAbsData::SumW2));
    //padDist->cd();
    fixFrameStyle(massFrame, false);
    massFrame->Draw();
    TLatex* textVar = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 0, 0.57, 0.8); // textVar->Draw("same");
    TLatex* textCut = cutLatex(ispO, cutVector, 0.15,0.6); // textCut->Draw("same");
    TLegend* leg = makePlotLegend(massFrame, legendEntries, 0.15, 0.7, 0.5, 0.85); leg->Draw("same");
    TLatex *textAlice = AliceText(ispO); // textAlice->Draw("same");
    padPull->cd();
    fixPullStyle(pullFrame);
    pullFrame->Draw();
    TLine* linePull = new TLine(massMin, 0, massMax, 0);
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
    RooPlot* pullTauzFrame = ws->var("tauz")->frame(Title("Pull Distribution"), Range(tauzMin, tauzMax));
    pullTauzFrame->addPlotable(hpullTauz,"P");

    nPar = ws->pdf("totPDF_2D")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
    chi2ndf = tauzFrame->chiSquare(nPar);
    std::cout << "fit chi2 = " << chi2ndf << std::endl;
    
    ws->data("data")->plotOn(tauzFrame, DataError(RooAbsData::SumW2));
    //padDist->cd();
    fixFrameStyle(tauzFrame, true);
    tauzFrame->Draw();
    TLatex* textVarTauz = varLatex(ws, parIni, chi2ndf, fitMass, fitTauz, 0, 0, 0.57, 0.8); // textVarTauz->Draw("same");
    TLatex* textCutTauz = cutLatex(ispO, cutVector, 0.15,0.6); // textCutTauz->Draw("same");
    TLegend* legTauz = makePlotLegend(tauzFrame, legendEntries, 0.15, 0.7, 0.5, 0.85); legTauz->Draw("same");
    TLatex *textAliceTauz = AliceText(ispO); // textAliceTauz->Draw("same");
    
    padPull->cd();
    fixPullStyle(pullTauzFrame);
    pullTauzFrame->Draw();
    // TLine* linePullTauz = new TLine(ws->var("tauz")->getMin(), 0, ws->var("tauz")->getMax(),0);
    TLine* linePullTauz = new TLine(tauzMin, 0, tauzMax, 0);
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
