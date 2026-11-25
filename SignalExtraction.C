// Author: Batoul Diab, Paul Veen
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
map<string, double> SignalExtraction1Fit(map<string, string>& parIni, RooWorkspace* ws, const char *caseName, const char *methodName, string rangeLabel, bool ispO, bool isMC, bool fitMass, bool fitTauz, bool fitNpTauz, bool fitTauzBkg, struct KinCuts cutVector, int idCall) {
  gStyle->SetOptStat(0);
  // make output directory
  string outDirName = Form("output/outputFits_%s%s_%s", ispO?"pO":"OO", isMC?"_MC":"", caseName);
  gSystem->mkdir(outDirName.c_str(), kTRUE);
  
  bool doSplot = false; bool doSideBands = false; bool doSplot_binned = false;
  if (strcmp(methodName, "sPlot") == 0) doSplot = true;
  if (strcmp(methodName, "sideBands") == 0) doSideBands = true;
  if (strcmp(methodName, "sPlot_binned") == 0) doSplot_binned = true;

  const char *sigDataSet; const char *bkgDataSet;
  if (doSplot || doSplot_binned) {
    if (!isMC) { 
      sigDataSet = "sPlotDsSig";
      bkgDataSet = "sPlotDsBkg";
    }
    else { sigDataSet = "data"; }
  }
  if (doSideBands) {
    sigDataSet = "tauzSignalData";
    bkgDataSet = "tauzSBData";
  }

  // Number of bins to be drawn (does not affect fitting)
  // For now: hardcoded bins in the BINNED fit case
  int nBinsMass = 80; // used for mass only
  int nBins = 100; // used for tauz only

  // Set range for plotting and fitting the mass
  // Make sure to also modify this accordingly in the inputUtils.h !!
  double massMin = 2.4;
  double massMax = 4.0;

  // Set range for plotting the tau_z
  // pO zoom: 0.005
  // OO zoom: 0.002
  // nominal fitting range: 0.02 (both cases)
  double tauzMin = -0.02;
  double tauzMax = 0.02;
  double tauzMinFitting = -0.02;
  double tauzMaxFitting = 0.02;

  int nParTauzRes = -1;
  RooDataHist* binnedData = nullptr;

  // bulding the model from input
  // Two options: either reduce the data here or make several datasets already in inputUtils
  // Make it spit out an array of data sets and take the 0th one for doSplot, define that as the previous data
  // For sidebands also take sideBand1 and sideBand2
  BuildPDF(ws, parIni, isMC, fitMass, fitTauz);
  cout << "PDFs correctly built" << endl;

  // Perform the tauz background fit?
  // NOTE: this is also set to 'false' and 'true' as part of the 2-step resolution fit. Might have to be changed below too
  if (!fitTauzBkg) cout << "[INFO] NOTE: RUNNING WITHOUT TAUZ BACKGROUND FIT -------------------------------------------" << endl;

  if (fitMass && !fitTauz) {
    if (parIni["doIterativeFit"] == "1") {
      std::cout << "[INFO] Applying iterative fitting on background function in 1D mass fits" << std::endl;
      RooFitResult* fitResult_mass = ws->pdf("totPDF_mass")->fitTo(*ws->data("data"), Extended(kTRUE), SumW2Error(true), RooFit::Save());
      RooDataSet* sPlotDs = (RooDataSet*) ws->data("data")->Clone("data_sPlot");
      RooStats::SPlot* sData = new RooStats::SPlot("sData", "sPlot", *sPlotDs, ws->pdf("totPDF_mass"), RooArgList(*ws->var("fJpsi_mass"), *ws->var("fBkg_mass")));
      ws->import(*sPlotDs);
    }
    else {
      std::cout << "before mass fit sigma = " << ws->var("sigma_mass")->getVal() << " +/- " << ws->var("sigma_mass")->getError() << std::endl;
      RooFitResult* fitResult_mass = ws->pdf("totPDF_mass")->fitTo(*ws->data("data"), Extended(true), Range(massMin, massMax), SumW2Error(true), RooFit::Save());
      std::cout << "after mass fit sigma = " << ws->var("sigma_mass")->getVal() << " +/- " << ws->var("sigma_mass")->getError() << std::endl;
      if (!isMC) {
        if (doSplot || doSplot_binned) {
          RooDataSet* sPlotDs = (RooDataSet*) ws->data("data")->Clone("data_sPlot");
          RooStats::SPlot* sData = new RooStats::SPlot("sData", "sPlot", *sPlotDs, ws->pdf("totPDF_mass"), RooArgList(*ws->var("fJpsi_mass"), *ws->var("fBkg_mass")));
          ws->import(*sPlotDs);

          // Investigate the sWeights
          // TODO: add boolean to do this or not
          RooAbsData* dataSplotMass = ws->data("data_sPlot");

          if (!dataSplotMass) {
            cout << "[ERROR] data_sPlot not found in workspace!" << endl;
          }

          dataSplotMass->get(0)->Print("v");

          cout << "[INFO] Mass sPlot dataset: " << dataSplotMass->GetName() << endl;
          cout << "[INFO] Number of entries: " << dataSplotMass->numEntries() << endl;

          // ************************************************************************************
          // Mass sWeight QA

          if (!dataSplotMass) {
              cout << "[ERROR] Mass sPlot dataset not found!" << endl;
          }

          const char* weightName = "fJpsi_mass_sw";

          double minWeightMass = 1e10;
          double maxWeightMass = -1e10;
          double sumWeightMass = 0.;
          double sumWeight2Mass = 0.;
          int nNegativeMass = 0;

          for (int i = 0; i < dataSplotMass->numEntries(); ++i) {
              const RooArgSet* row = dataSplotMass->get(i);
              double w = row->getRealValue(weightName);

              minWeightMass = std::min(minWeightMass, w);
              maxWeightMass = std::max(maxWeightMass, w);
              sumWeightMass += w;
              sumWeight2Mass += w * w;

              if (w < 0.) nNegativeMass++;
          }

          cout << "\n[INFO] Mass sWeight diagnostics:" << endl;
          cout << "[INFO] Minimum sWeight = " << minWeightMass << endl;
          cout << "[INFO] Maximum sWeight = " << maxWeightMass << endl;
          cout << "[INFO] Sum of sWeights = " << sumWeightMass << endl;
          cout << "[INFO] Sum of squared weights = " << sumWeight2Mass << endl;
          cout << "[INFO] Negative weights = " << nNegativeMass << endl;
          cout << "[INFO] Effective entries = " << sumWeightMass * sumWeightMass / sumWeight2Mass << endl;
          cout << "[INFO] Fitted signal yield = " << ws->var("fJpsi_mass")->getVal() << endl;


          // -----------------------------------------------------------------------------
          // 1. Individual sWeight distribution

          TH1D* hWeightsMass = new TH1D(
              Form("hWeightsMass_%s", rangeLabel.c_str()),
              "Signal sWeight distribution; Signal sWeight; Candidates",
              100, minWeightMass, maxWeightMass
          );

          for (int i = 0; i < dataSplotMass->numEntries(); ++i) {
              const RooArgSet* row = dataSplotMass->get(i);
              double w = row->getRealValue(weightName);
              hWeightsMass->Fill(w);
          }

          TCanvas* cWeightsMass = new TCanvas("cWeightsMass", "Mass sWeight distribution", 900, 700);
          hWeightsMass->SetLineColor(kBlue + 1);
          hWeightsMass->Draw("HIST");

          TLine* lineZeroMass = new TLine(0., 0., 0., hWeightsMass->GetMaximum());
          lineZeroMass->SetLineColor(kRed);
          lineZeroMass->SetLineStyle(2);
          lineZeroMass->Draw("SAME");

          cWeightsMass->SaveAs(Form("%s/cWeights_Mass_Distribution_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));


          // -----------------------------------------------------------------------------
          // 2. Signal sWeight versus invariant mass

          TH2D* hWeightVsMass = new TH2D(
              Form("hWeightVsMass_%s", rangeLabel.c_str()),
              "Signal sWeight versus invariant mass; m_{#mu#mu} (GeV/c^{2}); Signal sWeight",
              100, massMin, massMax,
              100, minWeightMass, maxWeightMass
          );

          for (int i = 0; i < dataSplotMass->numEntries(); ++i) {
              const RooArgSet* row = dataSplotMass->get(i);
              double mass = row->getRealValue("mass");
              double w = row->getRealValue(weightName);

              hWeightVsMass->Fill(mass, w);
          }

          TCanvas* cWeightVsMass = new TCanvas("cWeightVsMass", "sWeight versus mass", 900, 700);
          hWeightVsMass->Draw("COLZ");
          cWeightVsMass->SaveAs(Form("%s/cWeights_Mass_VsMass_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));


          // -----------------------------------------------------------------------------
          // 3. Weighted versus unweighted invariant-mass distributions

          TH1D* hMassUnweighted = new TH1D(
              Form("hMassUnweighted_%s", rangeLabel.c_str()),
              ";m_{#mu#mu} (GeV/c^{2}); Candidates",
              nBins, massMin, massMax
          );

          TH1D* hMassWeighted = new TH1D(
              Form("hMassWeighted_%s", rangeLabel.c_str()),
              ";m_{#mu#mu} (GeV/c^{2}); Weighted candidates",
              nBins, massMin, massMax
          );

          hMassWeighted->Sumw2();

          for (int i = 0; i < dataSplotMass->numEntries(); ++i) {
              const RooArgSet* row = dataSplotMass->get(i);
              double mass = row->getRealValue("mass");
              double w = row->getRealValue(weightName);

              hMassUnweighted->Fill(mass);
              hMassWeighted->Fill(mass, w);
          }

          TCanvas* cMassComparison = new TCanvas("cMassComparison", "Weighted versus unweighted mass", 900, 700);

          hMassUnweighted->SetLineColor(kBlue);
          hMassWeighted->SetLineColor(kRed);

          hMassWeighted->Draw("E");
          hMassUnweighted->Draw("HIST SAME");

          TLegend* legMass = new TLegend(0.60, 0.75, 0.88, 0.88);
          legMass->AddEntry(hMassUnweighted, "Unweighted", "l");
          legMass->AddEntry(hMassWeighted, "Signal sWeighted", "lep");
          legMass->Draw();

          cMassComparison->SaveAs(Form("%s/cWeights_Mass_WeightedVsUnweighted_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));

        }
        // ****************************************************************************
        else if (doSideBands) {
          cout << "[INFO] Constructing sideband-subtracted tauz distribution" << endl;

          RooDataSet* dataPeak = dynamic_cast<RooDataSet*>(ws->data("dataPeak"));
          RooDataSet* dataSBLow = dynamic_cast<RooDataSet*>(ws->data("dataSBLow"));
          RooDataSet* dataSBHigh = dynamic_cast<RooDataSet*>(ws->data("dataSBHigh"));

          double NJpsi = ws->var("fJpsi_mass")->getVal();
          double Nbkg  = ws->var("fBkg_mass")->getVal();
          cout << "[INFO] NJpsi = " << NJpsi << endl;
          cout << "[INFO] Nbkg  = " << Nbkg  << endl;

          RooAbsPdf* bkgPDF = ws->pdf("bkgPDF_mass");
          RooRealVar* mass = ws->var("mass");
          if (!bkgPDF || !mass) { cout << "[ERROR] Could not find mass or bkgPDF_mass!" << endl; }

          RooArgSet massSet(*mass);
          mass->setRange("peak", 2.9, 3.2);
          mass->setRange("lowSB", 2.4, 2.8);
          mass->setRange("highSB", 3.4, 4.0);

          // Integral of the background PDF in the peak region
          double bkgPeakIntegral = bkgPDF->createIntegral(*mass, RooFit::NormSet(*mass), RooFit::Range("peak"))->getVal();
          // Integral in lower sideband
          double bkgLowIntegral = bkgPDF->createIntegral(*mass, RooFit::NormSet(*mass), RooFit::Range("lowSB"))->getVal();
          // Integral in upper sideband
          double bkgHighIntegral = bkgPDF->createIntegral(*mass, RooFit::NormSet(*mass), RooFit::Range("highSB"))->getVal();
          double bkgSBIntegral = bkgLowIntegral + bkgHighIntegral;
          double alpha = bkgPeakIntegral / bkgSBIntegral;
          cout << "[INFO] Background integral in peak = " << bkgPeakIntegral << endl;
          cout << "[INFO] Background integral in sidebands = " << bkgSBIntegral << endl;
          cout << "[INFO] Sideband scale alpha = " << alpha << endl;

          // Now create the tauz histograms from this info
          // Tauz fit will be done later in the tauz fitting section...

          const int nTauzBins = nBins;
          const double tauzMinSB = tauzMin;
          const double tauzMaxSB = tauzMax;
          TH1D* hTauzPeak = new TH1D("hTauzPeak", ";#tau_{z};Entries", nTauzBins, tauzMinSB, tauzMaxSB);
          TH1D* hTauzSB = new TH1D("hTauzSB", ";#tau_{z};Entries", nTauzBins, tauzMinSB, tauzMaxSB);

          hTauzPeak->Sumw2();
          hTauzSB->Sumw2();

          for (int i = 0; i < dataPeak->numEntries(); ++i) {
            const RooArgSet* row = dataPeak->get(i);
            double tauz = row->getRealValue("tauz");
            hTauzPeak->Fill(tauz);
          }
          for (int i = 0; i < dataSBLow->numEntries(); ++i) {
            const RooArgSet* row = dataSBLow->get(i);
            double tauz = row->getRealValue("tauz");
            hTauzSB->Fill(tauz);
          }
          for (int i = 0; i < dataSBHigh->numEntries(); ++i) {
            const RooArgSet* row = dataSBHigh->get(i);
            double tauz = row->getRealValue("tauz");
            hTauzSB->Fill(tauz);
          }

          // Scale and subtract background under peak
          TH1D* hTauzBkg = (TH1D*)hTauzSB->Clone("hTauzBkg");
          hTauzBkg->Scale(alpha);
          TH1D* hTauzSignal = (TH1D*)hTauzPeak->Clone("hTauzSignal");
          hTauzSignal->Add(hTauzBkg, -1.0);

          cout << "[INFO] Peak tauz entries = " << hTauzPeak->Integral() << endl;
          cout << "[INFO] Scaled background = " << hTauzBkg->Integral() << endl;
          cout << "[INFO] Unscaled background tauz entries = " << hTauzSB->Integral() << endl;
          cout << "[INFO] Signal - background in peak = " << hTauzSignal->Integral() << endl;
          cout << "[INFO] Mass-fit Jpsi yield = " << NJpsi << endl;

          ws->import(*hTauzPeak);
          ws->import(*hTauzBkg);
          ws->import(*hTauzSB);
          ws->import(*hTauzSignal);
          cout << "[INFO] Sideband-subtracted tauz distribution created" << endl;

          TCanvas *canTauz = new TCanvas("canTauz", "comparison tauz distributions");
          canTauz->cd();
          // gPad->SetLogy();
          hTauzPeak->SetLineColor(kRed);
          // hTauzPeak->GetYaxis()->SetRangeUser(1e-1, 3e2);
          hTauzBkg->SetLineColor(kBlue);
          hTauzSignal->SetLineColor(kGreen+4);
          hTauzPeak->Draw("hist");
          hTauzBkg->Draw("hist same");
          hTauzSignal->Draw("hist same");
          TLegend *legTauz = new TLegend(0.65, 0.7, 0.88, 0.88);
          legTauz->AddEntry(hTauzPeak, "Peak", "l");
          legTauz->AddEntry(hTauzBkg, "Scaled background", "l");
          legTauz->AddEntry(hTauzSignal, "Signal = Peak - scaled background", "l");
          legTauz->Draw("same");
          canTauz->SaveAs(Form("%s/comparisonTauz_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
        }
      }
    }
  }
  else if (!fitMass && fitTauz && !isMC) {
    // Step 1: resolution fit on negative tauz range
    if (!fitNpTauz) {
      RooPlot* tauzResFrame = ws->var("tauz")->frame(Range(tauzMinFitting, 0), Bins(nBins));
      ws->data(sigDataSet)->plotOn(tauzResFrame, DataError(RooAbsData::SumW2));

      // ************************************************************************************
      // We want to investigate the splot weights
      // TODO: add a bool to do this or not
      RooAbsData* dataSplot = ws->data(sigDataSet);

      double minWeight = 1e10;
      double maxWeight = -1e10;

      for (int i = 0; i < dataSplot->numEntries(); ++i) {
        dataSplot->get(i);
        double w = dataSplot->weight();
        minWeight = std::min(minWeight, w);
        maxWeight = std::max(maxWeight, w);
      }

      cout << "[INFO] Minimum sWeight = " << minWeight << endl;
      cout << "[INFO] Maximum sWeight = " << maxWeight << endl;

      TH2D* hWeightVsTauz = new TH2D(
          "hWeightVsTauz",
          "sWeight versus #tau_{z}; #tau_{z} (ns); sWeight",
          100, tauzMinFitting, 0.02,
          100, minWeight, maxWeight
      );

      for (int i = 0; i < dataSplot->numEntries(); ++i) {
          const RooArgSet* row = dataSplot->get(i);
          double tauz = row->getRealValue("tauz");
          double w = dataSplot->weight();

          hWeightVsTauz->Fill(tauz, w);
      }

      TCanvas* cWeightVsTauz = new TCanvas("cWeight_VsTauz", "sWeight versus tauz", 900, 700);
      hWeightVsTauz->Draw("COLZ");
      cWeightVsTauz->SaveAs(Form("%s/cWeights_Tauz_VsTauz_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));

      TH1D* hTauzUnweighted = new TH1D("hTauzUnweighted", ";#tau_{z} (ns); Candidates", nBins, tauzMinFitting, tauzMaxFitting);
      TH1D* hTauzWeighted = new TH1D("hTauzWeighted", ";#tau_{z} (ns); Weighted candidates", nBins, tauzMinFitting, tauzMaxFitting);

      hTauzWeighted->Sumw2();

      for (int i = 0; i < dataSplot->numEntries(); ++i) {
          const RooArgSet* row = dataSplot->get(i);
          double tauz = row->getRealValue("tauz");
          double w = dataSplot->weight();

          hTauzUnweighted->Fill(tauz);
          hTauzWeighted->Fill(tauz, w);
      }

      TCanvas* cTauzComparison = new TCanvas("cTauzComparison", "Weighted versus unweighted tauz", 900, 700);

      hTauzUnweighted->SetLineColor(kBlue);
      hTauzWeighted->SetLineColor(kRed);

      hTauzUnweighted->Draw("HIST");
      hTauzWeighted->Draw("E SAME");

      TLegend* leg = new TLegend(0.60, 0.75, 0.88, 0.88);
      leg->AddEntry(hTauzUnweighted, "Unweighted", "l");
      leg->AddEntry(hTauzWeighted, "sWeighted", "lep");
      leg->Draw();

      cTauzComparison->cd();
      gPad->SetLogy();
      cTauzComparison->SaveAs(Form("%s/cWeights_Tauz_WeightedVsUnweighted_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));

      // ************************************************************************************

      // We need a very fine binning to estimate the xMaxRes
      RooPlot* tauzResFrameXmaxRes = ws->var("tauz")->frame(Range(-0.000005, 0), Bins(nBins));
      ws->data(sigDataSet)->plotOn(tauzResFrameXmaxRes, DataError(RooAbsData::SumW2));
      RooHist* hist = (RooHist*) tauzResFrameXmaxRes->getObject(0);
      ws->var("xMaxRes")->setVal(getMax(hist));
      double xMaxRes = ws->var("xMaxRes")->getVal(); // to get the mean of the resolution function
      // xMaxRes = 0.;
      cout << "[INFO] xMaxRes = " << xMaxRes << " ns" << endl;
      ws->var("mean_tauzRes")->setVal(xMaxRes);
      ws->var("mean_tauzRes")->setConstant(kTRUE);
      // double bias = 0;
      // ws->var("mean_tauzRes")->setMin(ws->var("mean_tauzRes")->getVal() - bias * std::abs(ws->var("mean_tauzRes")->getVal()));
      // ws->var("mean_tauzRes")->setMax(ws->var("mean_tauzRes")->getVal() + bias * std::abs(ws->var("mean_tauzRes")->getVal()));
      // std::cout << "[INFO] mean_tauzRes_min = " << ws->var("mean_tauzRes")->getMin() << std::endl;
      // std::cout << "[INFO] mean_tauzRes_max = " << ws->var("mean_tauzRes")->getMax() << std::endl;

      ws->var("tauz")->setRange("neg", tauzMinFitting, xMaxRes);
      RooDataSet* sPlotDsNeg = (RooDataSet*) ws->data(sigDataSet)->reduce(Form("tauz < %f", xMaxRes));
      ws->import(*sPlotDsNeg, RooFit::Rename("sPlotDsSigNeg"));
      RooFitResult* fitResult_tauzRes = ws->pdf("tauzResPDF")->fitTo(*ws->data("sPlotDsSigNeg"), Extended(kFALSE), SumW2Error(true), RooFit::Save(), Range("neg"));
      nParTauzRes = ws->pdf("tauzResPDF")->getParameters(*ws->data(sigDataSet))->selectByAttrib("Constant",kFALSE)->getSize();
      fixParPDF(ws, fitResult_tauzRes, parIni, ispO, rangeLabel, caseName, 0, 1, 1, 0, 0);
      // Save params, as is done for the 2D fit with function
      cout << "[INFO] done with fixing the params for tauz resolution 1st step" << endl;
    }
    
    // Step 2: include signal fit of non-prompt decay (initialised from MC)
    else {
      // Take parameters from 1st step
      fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, 0, 1, 1, 0, 0);
      RooAbsData* dataToFit = ws->data(sigDataSet);

      if (doSplot_binned) {
        RooRealVar* tauz = ws->var("tauz");
        if (!tauz || !dataToFit) { cout << "[ERROR] Missing tauz variable or signal dataset!" << endl; }
        cout << "[INFO] Performing BINNED non-prompt tauz fit" << endl;
        // cout << "[INFO] Inheriting tauz binning: " << tauz->getBinning().numBins() << " bins" << endl;
        cout << "[INFO] Using " << nBins << " tauz bins" << endl;
        tauz->setBins(nBins);
        binnedData = new RooDataHist("sPlotDsSig_binned", "Binned non-prompt signal dataset", RooArgList(*tauz), *dataToFit);
        dataToFit = binnedData;
      }
      else if (doSplot) { cout << "[INFO] Performing unbinned non-prompt tauz fit" << endl; }

      // Debug
      cout << "\n[DEBUG] Dataset passed to fitTo():"
      << "\n  Name:       " << dataToFit->GetName()
      << "\n  Class:      " << dataToFit->ClassName()
      << "\n  Entries:    " << dataToFit->numEntries()
      << "\n  Sum weights:" << dataToFit->sumEntries()
      << "\n  Weighted:   " << dataToFit->isWeighted()
      << endl;

      // Perform fit
      // Prevent Roofit from taking giant initial step in minimisation procedure
      ws->var("fb_tauzSig")->setError(0.01);
      // TODO: ADD BACK SumW2Error(kTRUE) !!! or AsymptoticError(kTRUE)
      // TODO: REMOVE BINNING
      cout << "Default tauz range: "
      << ws->var("tauz")->getMin() << " to "
      << ws->var("tauz")->getMax() << endl;
      ws->var("tauz")->setRange(tauzMinFitting, tauzMaxFitting);
      // RooFitResult* fitResult_tauz = ws->pdf("tauzSigPDF")->fitTo(*dataToFit, Extended(kFALSE), AsymptoticError(kTRUE), RooFit::Save(), Range(tauzMinFitting, tauzMaxFitting));
      RooFitResult* fitResult_tauz = ws->pdf("tauzSigPDF")->fitTo(*dataToFit, Extended(kFALSE), AsymptoticError(kTRUE), RooFit::Save());

      // Create the NLL to see what is going on with the minimisation procedure
      // RooAbsReal* nll = ws->pdf("tauzSigPDF")->createNLL(*dataToFit);
      // RooAbsReal* nll = ws->pdf("tauzSigPDF")->createNLL(*dataToFit, Extended(kFALSE), Range(tauzMinFitting, tauzMaxFitting));
      RooAbsReal* nll = ws->pdf("tauzSigPDF")->createNLL(*dataToFit, Extended(kFALSE));
      // RooAbsReal* nll = ws->pdf("tauzSigPDF")->createNLL(*dataToFit);
      // Plot the NLL as a function of fb_tauzSig
      RooPlot* nllFrame = ws->var("fb_tauzSig")->frame(RooFit::Range(0.0, 1.0));
      nll->plotOn(nllFrame, RooFit::LineColor(kBlue));
      TCanvas *cNNL = new TCanvas("cNNL", "NNL scan fb", 900, 700);
      nllFrame->Draw();
      cNNL->SaveAs(Form("%s/c0NLL_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));

      cout << "[INFO] done with the tauz Sig fit" << endl;

      // Before fixing the parameters, get the number of free parameters
      nParTauzRes = ws->pdf("tauzSigPDF")->getParameters(*dataToFit)->selectByAttrib("Constant", kFALSE)->getSize();
      fixParPDF(ws, fitResult_tauz, parIni, ispO, rangeLabel, caseName, 0, 1, 0, 0, 1);
      cout << "[INFO] done with fixing the params for tauz resolution 2nd step" << endl;

      cout << "\n[INFO] Non-prompt tauz fit diagnostics:" << "\n  status   = " << fitResult_tauz->status() << "\n  covQual  = " << fitResult_tauz->covQual() << "\n  EDM      = " << fitResult_tauz->edm() << "\n  minNLL   = " << fitResult_tauz->minNll() << endl;
      fitResult_tauz->Print("v");
      // delete binnedData;
    }

    // Do the background fit using the resolution
    if (fitTauzBkg) {
      /*
      TH1* hTauzBkg = (TH1*)ws->data(bkgDataSet)->createHistogram("h1", *ws->var("tauz"));
      // Below is wrong, because histogram is not filled with ->Fill()
      std::cout << "number of entries = " << hTauzBkg->GetEntries() << std::endl;
      TH1* hNoNegativeEntriesTauzBkg = (TH1*)hTauzBkg->Clone("hNoNegativeEntriesTauzBkg");
      for (int i_bin = 1; i_bin < hTauzBkg->GetNbinsX(); i_bin++) {
        double binCenter = hTauzBkg->GetBinCenter(i_bin);
        double content   = hTauzBkg->GetBinContent(i_bin);
        double error   = hTauzBkg->GetBinError(i_bin);
        std::cout << "Bin " << i_bin  << " | tauz val = " << binCenter  << " | content = " << content << std::endl;
        // Remove negative entries
        if (content < 0) {
          std::cout << "NEGATIVE entry = " << content << std::endl;
          hNoNegativeEntriesTauzBkg->SetBinContent(i_bin, 0);
          hNoNegativeEntriesTauzBkg->SetBinError(i_bin, 0);
        } else {
          hNoNegativeEntriesTauzBkg->SetBinContent(i_bin, content);
          hNoNegativeEntriesTauzBkg->SetBinError(i_bin, error);
        }
      }
      for (int i_bin = 1; i_bin < hNoNegativeEntriesTauzBkg->GetNbinsX(); i_bin++) {
        double binCenter = hNoNegativeEntriesTauzBkg->GetBinCenter(i_bin);
        double content = hNoNegativeEntriesTauzBkg->GetBinContent(i_bin);
        double error = hNoNegativeEntriesTauzBkg->GetBinError(i_bin);
        std::cout << "Bin after removing negative/0 " << i_bin  << " | tauz val = " << binCenter  << " | content = " << content << std::endl;
      }
      // transform again into roofit dataset
      // perform fit again
      // RooDataHist cleanDataHist("cleanDataHist", "Binned dataset without negatieve bins", *ws->var("tauz"), hNoNegativeEntriesTauzBkg);
      */
      // RooFitResult* fitResult_tauzBkg = ws->pdf("tauzBkgPDF")->fitTo(cleanDataHist, Extended(false), SumW2Error(true), RooFit::Save());
      RooFitResult* fitResult_tauzBkg = ws->pdf("tauzBkgPDF")->fitTo(*ws->data(bkgDataSet), Extended(false), SumW2Error(true), RooFit::Save());
      cout << "[INFO] done with tauz bkg fit" << endl;
      fixParPDF(ws, fitResult_tauzBkg, parIni, ispO, rangeLabel, caseName, 0, 0, 0, 1, 0);
      cout << "[INFO] done with fixing the params for tauz background fit" << endl;
    }
  }
  else if (!fitMass && fitTauz && isMC) {
    // Option 1. Fit the non-prompt signal: exponential decay convoluted with the detector resolution (MC)
    // RooFitResult* fitResult_tauzNpr = ws->pdf("tauzNprSigPDF")->fitTo(*ws->data("data"), Extended(false), SumW2Error(true), RooFit::Save());
    // Option 2. Fit the MC in the same way as is done for data (Gaussians for resolution, no NP component ek proteron)
    RooPlot* tauzResFrame = ws->var("tauz")->frame(Range(tauzMinFitting, 0), Bins(nBins));
    ws->data(sigDataSet)->plotOn(tauzResFrame, DataError(RooAbsData::SumW2));

    // We need a very fine binning to estimate the xMaxRes
    RooPlot* tauzResFrameXmaxRes = ws->var("tauz")->frame(Range(-0.000005, 0), Bins(nBins));
    ws->data(sigDataSet)->plotOn(tauzResFrameXmaxRes, DataError(RooAbsData::SumW2));
    RooHist* hist = (RooHist*) tauzResFrameXmaxRes->getObject(0);
    ws->var("xMaxRes")->setVal(getMax(hist));
    double xMaxRes = ws->var("xMaxRes")->getVal(); // to get the mean of the resolution function
    // xMaxRes = 0.;
    cout << "[INFO] xMaxRes = " << xMaxRes << " ns" << endl;
    ws->var("mean_tauzRes")->setVal(xMaxRes);
    // double bias = 0.2;
    // ws->var("mean_tauzRes")->setMin(ws->var("mean_tauzRes")->getVal() - bias * std::abs(ws->var("mean_tauzRes")->getVal()));
    // ws->var("mean_tauzRes")->setMax(ws->var("mean_tauzRes")->getVal() + bias * std::abs(ws->var("mean_tauzRes")->getVal()));
    // std::cout << "[INFO] mean_tauzRes_min = " << ws->var("mean_tauzRes")->getMin() << std::endl;
    // std::cout << "[INFO] mean_tauzRes_max = " << ws->var("mean_tauzRes")->getMax() << std::endl;
    // ws->var("mean_tauzRes")->setConstant(kTRUE);

    ws->var("tauz")->setRange("neg", tauzMinFitting, xMaxRes);
    RooDataSet* sPlotDsNeg = (RooDataSet*) ws->data(sigDataSet)->reduce(Form("tauz < %f", xMaxRes));
    ws->import(*sPlotDsNeg, RooFit::Rename("sPlotDsSigNeg"));
    RooFitResult* fitResult_tauzRes = ws->pdf("tauzResPDF")->fitTo(*ws->data("sPlotDsSigNeg"), Extended(kFALSE), SumW2Error(true), RooFit::Save(), Range("neg"));
    nParTauzRes = ws->pdf("tauzResPDF")->getParameters(*ws->data(sigDataSet))->selectByAttrib("Constant",kFALSE)->getSize();
    fixParPDF(ws, fitResult_tauzRes, parIni, ispO, rangeLabel, caseName, 0, 1, 1, 0, 0);
    // Save params, as is done for the 2D fit with function
    cout << "[INFO] Done with MC non-prompt tauz fit" << endl;
    fitResult_tauzRes->Print();
  }
  else if (fitMass && fitTauz && !isMC) {
    // Fix mass-shape parameters from the 1D mass fit
    fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, true, false, false, false, false);
    // Fix tauz-resolution parameters from the 1D tauz fit
    // Also get the lambda parameter
    // TODO: actually don't do this when you don't do the 2nd step resolution fit (check if true?)
    fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, false, true, false, false, true);
    // Fix tauz-background parameters from the 1D tauz fit
    /*
    if (!isMC) { 
      fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, false, false, false, true, false); 
    }
    */
   if (!isMC) {
    /*
    std::cout << "\n[INFO] Parameters BEFORE fixParPDF:" << std::endl;
    for (const auto& par : parIni) {
      if (auto* var = ws->var(par.first.c_str())) {
        std::cout << "  " << par.first << " = " << var->getVal() << " ± " << var->getError() << ", constant = " << (var->isConstant() ? "true" : "false") << std::endl; 
      }
    }
    */
    fixParPDF(ws, NULL, parIni, ispO, rangeLabel, caseName, false, false, false, true, false);
    /*
    std::cout << "\n[INFO] Parameters AFTER fixParPDF:" << std::endl;
    for (const auto& par : parIni) {
      if (auto* var = ws->var(par.first.c_str())) { 
        std::cout << "  " << par.first << " = " << var->getVal() << " ± " << var->getError() << ", constant = " << (var->isConstant() ? "true" : "false") << std::endl;
      }
    */
    }
    RooFitResult* fitResult_tauzMass = ws->pdf("totPDF_2D")->fitTo(*ws->data("data"), Extended(kTRUE), SumW2Error(true), RooFit::Save());
    std::cout << "\n[INFO] Parameters AFTER fitting 2D data:" << std::endl;
    for (const auto& par : parIni) {
      if (auto* var = ws->var(par.first.c_str())) { 
        std::cout << "  " << par.first << " = " << var->getVal() << " ± " << var->getError() << ", constant = " << (var->isConstant() ? "true" : "false") << std::endl;
      }
    }
    fitResult_tauzMass->Print("v");
    std::cout << "status = " << fitResult_tauzMass->status() << std::endl;
    std::cout << "covQual = " << fitResult_tauzMass->covQual() << std::endl;
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
  
  RooPlot* massFrame = ws->var("mass")->frame(Range(massMin, massMax), Bins(nBinsMass));
  RooPlot* tauzFrame = ws->var("tauz")->frame(Range(tauzMin, tauzMax), Bins(nBins));
  RooPlot* tauzResFrame = ws->var("tauz")->frame(Range(tauzMin, tauzMax), Bins(nBins));
  RooPlot* tauzBkgFrame = ws->var("tauz")->frame(Range(tauzMin, tauzMax), Bins(nBins));

  int nPar;
  double chi2ndf = -999;
  double chi2ndfTauzNp = -999;
  double chi2ndfTauzBkg = -999;
  
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
            ws->pdf("totPDF_mass")->fitTo(*ws->data("data"), Extended(kTRUE), SumW2Error(true), RooFit::Save());

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
      massFrame = ws->var("mass")->frame(Range(massMin, massMax), Bins(nBinsMass));
      ws->data("data")->plotOn(massFrame, Name("data")); legendEntries["data"] = {"data","P"};
      if (!isMC) { ws->pdf("totPDF_mass")->plotOn(massFrame, Name("background_mass"), Components(RooArgSet(*ws->pdf("bkgPDF_mass"))), DrawOption("F"), FillColor(kGray), LineColor(kGray), VLines()); legendEntries["background_mass"] = {"Background", "F"}; }
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalPsi2s_mass"), Components(RooArgSet(*ws->pdf("psi2sPDF_mass"))), DrawOption("L"), LineColor(kGreen+4)); //legendEntries["signalPsi2s_mass"] = {"#psi(2S) signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalJpsi_mass"), Components(RooArgSet(*ws->pdf("jpsiPDF_mass"))), DrawOption("L"), LineColor(kGreen+2)); legendEntries["signalJpsi_mass"] = {"J/#psi signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("total_mass"), LineColor(kRed)); legendEntries["total_mass"] = {"total fit","L"};
    }
    else {
      ws->data("data")->plotOn(massFrame, Name("data")); legendEntries["data"] = {"data","P"};
      if (!isMC) { ws->pdf("totPDF_mass")->plotOn(massFrame, Name("background_mass"), Components(RooArgSet(*ws->pdf("bkgPDF_mass"))), DrawOption("F"), FillColor(kGray), LineColor(kGray), VLines()); legendEntries["background_mass"] = {"Background", "F"}; }
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalPsi2s_mass"), Components(RooArgSet(*ws->pdf("psi2sPDF_mass"))), DrawOption("L"), LineColor(kBlue+4)); legendEntries["signalPsi2s_mass"] = {"#psi(2S) signal","L"};
      ws->pdf("totPDF_mass")->plotOn(massFrame, Name("signalJpsi_mass"), Components(RooArgSet(*ws->pdf("jpsiPDF_mass"))), DrawOption("L"), LineColor(kGreen+2)); legendEntries["signalJpsi_mass"] = {"J/#psi signal","L"};
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
    // std::string pdfPath = "";
    can->SaveAs(Form("%s/massFit1D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
    // can->SaveAs(Form("%s/massFit1D_%s.root", outDirName.c_str(), rangeLabel.c_str()));
    if (parIni["doIterativeFit"] == "1") {
      canFitChi2->SaveAs(Form("%s/massFit1D_%s_fitChi2Trend.pdf", outDirName.c_str(), rangeLabel.c_str()));
    }
  }
  else if (fitTauz && !fitMass && !isMC) {
    // Draw tauz resolution first
    padDist->cd();
    ws->Print(); 
    // change this line to display only the negative data points (or all data points)
    // ws->data(sigDataSet)->plotOn(tauzResFrame, Name(sigDataSet), DataError(RooAbsData::SumW2)); legendEntries["sPlotDsSigNeg"] = {"sPlot signal-like data","P"};
    RooAbsData* tauzPlotData = ws->data(sigDataSet);

    if (fitNpTauz && doSplot_binned) {
      if (!binnedData) { cout << "[ERROR] Binned dataset is missing for the tauz plot!" << endl; }
      tauzPlotData = binnedData;
      cout << "[INFO] Plotting binned non-prompt tauz dataset: " << tauzPlotData->GetName() << endl; }
    else if (fitNpTauz && doSplot) { cout << "[INFO] Plotting unbinned non-prompt tauz dataset: " << tauzPlotData->GetName() << endl; }
    tauzPlotData->plotOn(tauzResFrame, Name("tauzData"), DataError(RooAbsData::SumW2)); legendEntries["tauzData"] = {"sPlot signal-like data", "P"};

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
    can->SaveAs(Form("%s/tauzResFit1D%s_%s_%i.pdf", outDirName.c_str(), fitNpTauz?"_npTauz":"", rangeLabel.c_str(), idCall));
    // can->SaveAs(Form("%s/tauzResFit1D%s_%s_%i.png", outDirName.c_str(), fitNpTauz?"_npTauz":"", rangeLabel.c_str(), idCall));
    // ROOT file containing the editable objects
    TFile *fOut = new TFile(Form("%s/tauzResFit1D%s_%s.root", outDirName.c_str(), fitNpTauz?"_npTauz":"", rangeLabel.c_str()), "RECREATE");
    can->Write("can");
    tauzResFrame->Write("tauFrame");
    pullFrame->Write("pullFrame");
    ws->Write("ws");
    fOut->Close();

    //////// then tauz Bkg
    if (fitTauzBkg) {
      padDist->cd();
      legendEntries.clear();
      ws->data(bkgDataSet)->plotOn(tauzBkgFrame, Name(bkgDataSet), DataError(RooAbsData::SumW2)); legendEntries[bkgDataSet] = {"sPlot background-like data","P"};
      ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzNprBkgPDF"), Components(RooArgSet(*ws->pdf("tauzNprBkgPDF"))), DrawOption("L"), LineColor(kBlue)); legendEntries["tauzNprBkgPDF"] = {"non-prompt-like", "L"};
      ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzPrBkgPDF"), Components(RooArgSet(*ws->pdf("tauzPrBkgPDF"))), DrawOption("L"), LineColor(kGreen+2)); legendEntries["tauzPrBkgPDF"] = {"prompt-like", "L"};
      ws->pdf("tauzBkgPDF")->plotOn(tauzBkgFrame, Name("tauzBkgPDF"), LineColor(kRed)); legendEntries["tauzBkgPDF"] = {"total fit","L"};
      
      RooHist* hpullBkg = tauzBkgFrame->pullHist();
      RooPlot* pullBkgFrame = ws->var("tauz")->frame(Title("Pull Distribution"), Range(tauzMin, tauzMax));
      pullBkgFrame->addPlotable(hpullBkg,"P");

      nPar = ws->pdf("tauzBkgPDF")->getParameters(*ws->data("data"))->selectByAttrib("Constant",kFALSE)->getSize();
      chi2ndfTauzBkg = tauzBkgFrame->chiSquare(nPar);
      std::cout << "fit chi2 = " << chi2ndfTauzBkg << std::endl;
      
      ws->data(bkgDataSet)->plotOn(tauzBkgFrame, DataError(RooAbsData::SumW2));
      // padDist->cd();
      fixFrameStyle(tauzBkgFrame, true);
      tauzBkgFrame->Draw();
      TLatex* textVarBkg = varLatex(ws, parIni, chi2ndfTauzBkg, fitMass, fitTauz, 0, 1, 0.57, 0.85); // textVarBkg->Draw("same");
      TLatex* textCutBkg = cutLatex(ispO, cutVector, 0.2,0.8); // textCutBkg->Draw("same");
      TLegend* legBkg = makePlotLegend(tauzBkgFrame, legendEntries, 0.15, 0.5, 0.3, 0.65); legBkg->Draw("same");
      TLatex *textAliceBkg = AliceText(ispO); // textAlice->Draw("same");
      
      padPull->cd();
      fixPullStyle(pullBkgFrame);
      pullBkgFrame->Draw();
      // TLine* linePullBkg = new TLine(ws->var("tauz")->getMin(), 0, ws->var("tauz")->getMax(),0);
      TLine* linePullBkg = new TLine(tauzMin, 0, tauzMax, 0);
      linePullBkg->SetLineColor(kRed); linePullBkg->SetLineStyle(2); linePullBkg->Draw("same");
      // std::string pdfPath = "";
      // can->SetLogy();
      can->SaveAs(Form("%s/tauzBkgFit1D_%s.pdf", outDirName.c_str(), rangeLabel.c_str()));
    }
  }
  else if (fitTauz && !fitMass && isMC) {
    // Option 1.
    /*
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
    */

    // Option 2.
    // Draw tauz resolution first
    padDist->cd();
    ws->Print(); 
    // change this line to display only the negative data points (or all data points)
    ws->data(sigDataSet)->plotOn(tauzResFrame, Name(sigDataSet), DataError(RooAbsData::SumW2)); legendEntries["sPlotDsSigNeg"] = {"sPlot signal-like data","P"};
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
    // TODO: this is used for the bkg fit chi2. Change every chi2ndf to a unique one and put it in the results map in the bottom of the code
    chi2ndfTauzBkg = tauzResFrame->chiSquare(nParTauzRes);
    std::cout << "fit chi2 = " << chi2ndfTauzBkg << std::endl;

    // ws->data("sPlotDsSigNeg")->plotOn(tauzResFrame, DataError(RooAbsData::SumW2));
    //padDist->cd();
    fixFrameStyle(tauzResFrame, true);
    tauzResFrame->Draw();
    TLatex* textVarRes = varLatex(ws, parIni, chi2ndfTauzBkg, fitMass, fitTauz, 1, 0, 0.57, 0.85);
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
    can->SaveAs(Form("%s/tauzResMCFit1D%s_%s.pdf", outDirName.c_str(), fitNpTauz?"_npTauz":"", rangeLabel.c_str()));
    // ROOT file containing the editable objects
    /*
    TFile *fOut = new TFile(Form("%s/tauzResFit1D%s_%s.root", outDirName.c_str(), fitNpTauz?"_npTauz":"", rangeLabel.c_str()), "RECREATE");
    can->Write("can");
    tauzResFrame->Write("tauFrame");
    pullFrame->Write("pullFrame");
    ws->Write("ws");
    fOut->Close();
    */
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
