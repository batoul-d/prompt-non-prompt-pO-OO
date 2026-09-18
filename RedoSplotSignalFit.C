// Sandbox macro to redo the signal ctau fit from SPlot output.
//
// The input ROOT file must contain a RooDataSet named "data" with the
// variable "ctau" and sWeights attached (produced by SUBARooFitter::PlotSPlot).

#include <RooAddModel.h>
#include <RooAddPdf.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooDecay.h>
#include <RooFFTConvPdf.h>
#include <RooFitResult.h>
#include <RooGaussModel.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooHist.h>
#include <RooHistPdf.h>
#include <RooKeysPdf.h>
#include <RooNumIntConfig.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooTruthModel.h>
#include <RooWorkspace.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

void RedoSplotSignalFit(const char* inputFile, const char* configFile,
                        const char* inputFolder = "plots",
                        const char* dataLabel = "#tau_{z} : pT integrated") {
  gStyle->SetOptStat(0);

  // ========================================================================
  // 1. Load the SPlot signal dataset
  // ========================================================================
  TFile* f = TFile::Open(inputFile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << inputFile << std::endl;
    return;
  }

  RooDataSet* data = dynamic_cast<RooDataSet*>(f->Get("data"));
  if (!data) {
    std::cerr << "ERROR: No RooDataSet named 'data' in file" << std::endl;
    return;
  }

  std::cout << "Loaded dataset with " << data->numEntries() << " entries, "
            << "sum of weights = " << data->sumEntries() << std::endl;

  // ========================================================================
  // 1b. Load JSON config
  // ========================================================================
  std::ifstream cfgStream(configFile);
  if (!cfgStream.is_open()) {
    std::cerr << "ERROR: Cannot open config file: " << configFile << std::endl;
    return;
  }
  nlohmann::json cfg;
  try {
    cfg = nlohmann::json::parse(cfgStream);
  } catch (const nlohmann::json::exception& e) {
    std::cerr << "ERROR: Failed to parse config JSON: " << e.what()
              << std::endl;
    return;
  }

  // --- Extract top-level settings (needed before any model construction) ---
  double ctauMin = cfg["fit_range"]["ctauMin"].get<double>();
  double ctauMax = cfg["fit_range"]["ctauMax"].get<double>();
  std::string mainModeStr = cfg["mode"]["mainMode"].get<std::string>();
  std::string subModeStr = cfg["mode"]["subMode"].get<std::string>();
  const char* mainMode = mainModeStr.c_str();
  const char* subMode = subModeStr.c_str();
  bool useBinnedData = cfg["mode"]["useBinnedData"].get<bool>();
  bool includeG4 = cfg["final_model"]["includeG4"].get<bool>();
  bool useCb2OwnMean = cfg["reso_G1_CB2"]["cb2_use_own_mean"].get<bool>();
  bool useCbRatio = cfg["reso_G1_CB2"]["use_cb_ratio"].get<bool>();
  double cbRalpha = cfg["reso_G1_CB2"]["cb_Ralpha"].get<double>();
  double cbRn = cfg["reso_G1_CB2"]["cb_Rn"].get<double>();

  std::cout << "Config: " << configFile << std::endl;
  std::cout << "  pT bin    : " << cfg.value("pT_bin", "?") << std::endl;
  std::cout << "  mainMode  : " << mainMode << std::endl;
  std::cout << "  subMode   : " << subMode << std::endl;
  std::cout << "  includeG4 : " << std::boolalpha << includeG4 << std::endl;
  std::cout << "  ctau range: [" << ctauMin << ", " << ctauMax << "]"
            << std::endl;

  // ========================================================================
  // 2. Define the ctau observable
  // ========================================================================
  // Detect the variable name from the dataset (could be "ctau", "ctau_diff",
  // "ctau_res1", etc. from QualityControl.C) and rename it to "ctau" so the
  // rest of this macro works unchanged.
  const RooArgSet* dsVars = data->get();
  if (!dsVars || dsVars->empty()) {
    std::cerr << "ERROR: Dataset has no variables" << std::endl;
    return;
  }
  std::string origVarName = dsVars->first()->GetName();
  std::cout << "Dataset variable name: " << origVarName << std::endl;

  // If the variable is not already called "ctau", rename it
  if (origVarName != "ctau") {
    std::cout << "Renaming variable '" << origVarName << "' -> 'ctau'"
              << std::endl;
    RooRealVar* origVar =
        dynamic_cast<RooRealVar*>(dsVars->find(origVarName.c_str()));
    if (origVar) origVar->SetName("ctau");
  }

  RooRealVar ctau("ctau", "c#tau (mm)", 0.5 * (ctauMin + ctauMax), ctauMin,
                  ctauMax);

  // Reduce data to fit range
  std::string sel = Form("ctau >= %.6f && ctau <= %.6f", ctauMin, ctauMax);
  RooDataSet* dataReduced =
      dynamic_cast<RooDataSet*>(data->reduce(sel.c_str()));
  std::cout << "After range cut: " << dataReduced->numEntries() << " entries"
            << std::endl;

  // ========================================================================
  // 3. Build Resolution Model: 4-Gaussian with shared mean
  // ========================================================================
  // Toggle which resolution config to use by commenting/uncommenting blocks.
  // The fractions use recursive interpretation:
  //   w1 = f1, w2 = (1-f1)*f2, w3 = (1-f1)*(1-f2)*f3, w4 = remainder

  // --- Shared mean (like exp_conv_gauss_CB2.C) ---
  // Use ONE shared mean for both Gaussian and CB (symmetric range!)
  RooRealVar reso_mean("reso_mean", "mean", 0.001, -0.05, 0.05);
  // reso_mean.setConstant(true);  // FIX initially (like exp_conv_gauss_CB2.C)

  // Keep these for multi-Gaussian models but we won't use them for CB2+Gauss
  RooRealVar reso_mean1("reso_mean_1", "mean1", 0.0007217, -0.05, 0.05);
  // RooRealVar reso_mean1("reso_mean_shared1", "shared mean1", 0.0, -0.5, 0.5);
  // reso_mean1.setConstant(true);  // FIX initially (like exp_conv_gauss_CB2.C)
  RooRealVar reso_mean2("reso_mean_2", "mean2", 0.0, -0.3, 0.3);
  RooRealVar cb2_mean("cb2_mean", "cb2_mean", 0.0, -0.15, 0.05);
  // RooRealVar cb2_mean("cb2_mean", "cb2_mean", 0.0, -0.1, 0.0);
  RooRealVar reso_mean3("reso_mean_3", "mean3", 0.0, -2.0, 0.0);
  RooRealVar reso_mean5("reso_mean_5", "mean5", 0.0, -2.0, 0.005);
  RooRealVar reso_mean4("reso_mean_4", "mean4", -1.2337, -3.0, 0.2);
  // reso_mean1.setConstant(true);
  // cb2_mean.setConstant(true);
  // reso_mean2.setConstant(true);
  // reso_mean3.setConstant(true);
  // reso_mean4.setConstant(true);
  // --- Gaussian sigmas (values from exp_conv_gauss_CB2.C) ---
  RooRealVar sigma1("reso_gauss1_sigma", "#sigma_{G1}", 0.04, 0.005, 0.05);
  // RooRealVar sigma1("reso_gauss1_sigma", "#sigma_{G1}", 0.05881, 0.04, 0.07);
  RooRealVar sigma2("reso_gauss2_sigma", "#sigma_{CB}", 0.10, 0.06, 0.18);
  // RooRealVar sigma2("reso_gauss2_sigma", "#sigma_{G2}", 0.15, 0.1, 0.2);
  RooRealVar cb2_sigma("cb2_sigma", "#sigma_{CB}", 0.06, 0.05, 0.15);
  // RooRealVar cb2_sigma("cb2_sigma", "#sigma_{CB}", 0.1068, 0.07, 0.15);
  RooRealVar sigma3("reso_gauss3_sigma", "#sigma_{G3}", 0.5, 0.1, 1.5);
  RooRealVar sigma4("reso_gauss4_sigma", "#sigma_{G4}", 6.29, 4.5, 8.0);
  RooRealVar sigma5("reso_gauss5_sigma", "#sigma_{G5}", 6.0, 1.2, 7.0);
  // RooRealVar sigma4("reso_gauss4_sigma", "#sigma_{G4}", 6.0, 4.0, 6.7);
  // FIX resolution sigmas initially (like exp_conv_gauss_CB2.C)
  // sigma1.setConstant(true);
  // cb2_sigma.setConstant(true);
  // sigma2.setConstant(true);
  // sigma3.setConstant(true);
  // sigma4.setConstant(true);

  // --- Fractions (values from exp_conv_gauss_CB2.C) ---
  RooRealVar frac1("reso_fraction1", "f_{1}", 0.6472, 0.0, 1.0);
  // RooRealVar frac1("reso_fraction1", "f_{G}", 0.53, 0.3, 0.8);
  RooRealVar frac2("reso_fraction2", "f_{2}", 0.9552, 0.0, 1.0);
  RooRealVar frac3("reso_fraction3", "f_{3}", 0.661007, 0.0, 1.0);
  RooRealVar frac4("reso_fraction4", "f_{4}", 0.661007, 0.0, 1.0);
  // RooRealVar frac2("reso_fraction2", "f_{2}", 0.4056, 0.1, 0.8);
  // RooRealVar frac3("reso_fraction3", "f_{3}", 0.761007, 0.3, 0.85);
  // FIX fraction initially
  // frac1.setConstant(true);
  // frac2.setConstant(true);
  // frac3.setConstant(true);

  // --- Fraction h for R+T model ---
  // h = (N_np + N_p) / (N_np + N_p + N_t)
  // Represents fraction of R-like behavior vs T (tails)
  RooRealVar frac_h("frac_h", "h (R vs T fraction)", 0.85, 0.65, 1.0);
  // fraction test : a1 * Prompt + a2 * NPrompt + (1-a1-a2) * Tails (fB is
  // calculated after)
  RooRealVar a("a", "a", 0.8156, 0.0, 1.0);
  // RooRealVar a1("a1", "a_{1}", 0.1056, 0.02, 0.17);
  // RooRealVar a2("a2", "a_{2}", 0.761007, 0.3, 0.85);
  RooRealVar a1("a1", "a_{1}", 0.1856, 0.0, 1.0);
  RooRealVar a2("a2", "a_{2}", 0.831, 0.0, 1.0);
  //   frac_h.setConstant(true);  // uncomment to fix
  // -----------------------------------------------------------------------
  // 3a. Standalone prompt PDF (RooGaussian + RooAddPdf)
  //     Used as the prompt component in the signal model.
  //     Always use RooGaussian here (not RooGaussModel) to avoid
  //     normalization issues when used standalone.
  // -----------------------------------------------------------------------
  //   RooGaussian g1("reso_gauss1", "G1", ctau, reso_mean, sigma1);
  //   RooGaussian g2("reso_gauss2", "G2", ctau, reso_mean, sigma2);
  //   RooGaussian g3("reso_gauss3", "G3", ctau, reso_mean, sigma3);
  //   RooGaussian g4("reso_gauss4", "G4", ctau, reso_mean, sigma4);

  RooGaussian g1("reso_gauss1", "G1", ctau, reso_mean1, sigma1);
  RooGaussian g2("reso_gauss2", "G2", ctau, reso_mean1, sigma2);
  // RooGaussian g2("reso_gauss2", "G2", ctau, reso_mean2, sigma2);
  RooGaussian g3("reso_gauss3", "G3", ctau, reso_mean3, sigma3);
  RooGaussian g4("reso_gauss4", "G4", ctau, reso_mean4, sigma4);
  RooGaussian g5("reso_gauss5", "G5", ctau, reso_mean5, sigma5);

  //   RooGaussian g1("reso_gauss1", "G1", ctau, reso_mean, sigma1);
  //   RooGaussian g2("reso_gauss2", "G2", ctau, reso_mean, sigma2);
  //   RooGaussian g3("reso_gauss3", "G3", ctau, reso_mean, sigma3);
  // RooGaussian g3("reso_gauss3", "G3", ctau, reso_mean4, sigma3);
  //  RooGaussian g4("reso_gauss4", "G4", ctau, reso_mean4, sigma4);

  // // pt0-1
  // double Ralpha = 1.0 / 1.053;
  // double Rn = 1.0 / 0.791;

  // // pt1-2
  // double Ralpha = 1.0 / 1.031;
  // double Rn = 1.0 / 0.833;

  // // pt2-3
  // double Ralpha = 1.0 / 1.015;
  // double Rn = 1.0 / 0.878;

  // // pt3-4
  // double Ralpha = 1.0 / 0.992;
  // double Rn = 1.0 / 0.959;

  // // pt4-5
  // double Ralpha = 1.0 / 1.001;
  // double Rn = 1.0 / 0.941;

  // // pt5-6
  // double Ralpha = 1.0 / 0.994;
  // double Rn = 1.0 / 0.976;

  // // pt6-8
  // double Ralpha = 1.0 / 0.995;
  // double Rn = 1.0 / 0.980;

  // // pt8-10
  // double Ralpha = 1.0 / 0.970;
  // double Rn = 1.0 / 1.023;

  // // pt10-15
  // double Ralpha = 1.0 / 0.963;
  // double Rn = 1.0 / 1.008;

  // // pt15-30
  // double Ralpha = 1.0 / 1.065;
  // double Rn = 1.0 / 0.884;

  RooRealVar cbAlpha1("cb_alpha1", "#alpha_{1}", 1.617, 1.3, 2.4);
  RooRealVar cbN1("cb_n1", "n_{1}", 2.464, 1.1, 4.0);
  RooRealVar cbN2("cb_n2", "n_{2}", 2.812, 1.1, 4.5);
  RooRealVar cbAlpha2("cb_alpha2", "#alpha_{2}", 1.976, 1.3, 2.6);

  // cbAlpha1.setConstant(true);
  // cbN1.setConstant(true);
  // cbAlpha2.setConstant(true);
  // cbN2.setConstant(true);

  // When use_cb_ratio=true: cbAlpha2 and cbN2 are linked to cbAlpha1/cbN1
  // via ratios from JSON (RooFormulaVar). Otherwise they float independently.
  std::unique_ptr<RooFormulaVar> cbAlpha2Fml, cbN2Fml;
  if (useCbRatio) {
    cbAlpha2Fml = std::make_unique<RooFormulaVar>("cb_alpha2", "#alpha_{2}",
                                                  Form("%f*@0", cbRalpha),
                                                  RooArgList(cbAlpha1));
    cbN2Fml = std::make_unique<RooFormulaVar>(
        "cb_n2", "n_{2}", Form("%f*@0", cbRn), RooArgList(cbN1));
  }
  RooAbsReal& cbAlpha2Ref =
      useCbRatio ? static_cast<RooAbsReal&>(*cbAlpha2Fml) : cbAlpha2;
  RooAbsReal& cbN2Ref = useCbRatio ? static_cast<RooAbsReal&>(*cbN2Fml) : cbN2;

  // Select CB2 mean: shared with G1 (reso_mean1) or its own (cb2_mean)
  RooRealVar& cb2MeanRef = useCb2OwnMean ? cb2_mean : reso_mean1;

  RooCrystalBall cb2("cb2", "symmetric CB2", ctau,
                     cb2MeanRef,  //
                     //  sigma1,      // CB width
                     //  cb2_mean,   //
                     cb2_sigma,    // CB width
                                   //  reso_mean3,  //
                                   //  sigma3,      // CB width
                                   //  cbAlpha1,   // alphaL (left tail)
                                   //  cbN1,       // nL (left tail)
                                   //  cbAlpha1,   // alphaR (righst tail)
                                   //  cbN1);      // nR (right tail)
                     cbAlpha1,     // alphaL (left tail)
                     cbN1,         // nL (left tail)
                                   //  cbAlpha1,   // alphaR (right tail)
                                   //  cbN1);      // nR (right tail)
                     cbAlpha2Ref,  // alphaR (right tail)
                     cbN2Ref);     // nR (right tail)

  // -----------------------------------------------------------------------
  // Custom exponential PDF models (shifted, one-sided, double-sided)
  // -----------------------------------------------------------------------
  // Shared rate parameter (lambda = 1/tau, in mm^-1) and per-model shift means
  RooRealVar lambda_exp1("lambda_exp1", "#lambda_{exp}", 1.0, 0.2, 20.0);
  RooRealVar lambda_exp2("lambda_exp2", "#lambda_{exp}", 1.0, 0.2, 20.0);
  RooRealVar lambda_exp3("lambda_exp3", "#lambda_{exp}", 0.4194, 0.31, 0.75);
  RooRealVar mean_exp1("mean_exp1", "#mu_{exp}^{R}", 0.0, 0.0, 1.0);
  RooRealVar mean_exp2("mean_exp2", "#mu_{exp}^{L}", 0.0, -1.0, 0.0);
  RooRealVar mean_exp3("mean_exp3", "#mu_{exp}^{DS}", -0.589, -1.5, -0.0);
  // mean_exp1.setConstant(true);
  // mean_exp2.setConstant(true);
  // mean_exp3.setConstant(true);
  // lambda_exp1.setConstant(true);
  // lambda_exp2.setConstant(true);
  // lambda_exp3.setConstant(true);
  // 1. Right-side exponential: exp(-(ctau - mean1) * lambda)  for ctau > mean1
  //    Indicator (@0 > @1) ensures it vanishes on the left side.
  RooGenericPdf expRight(
      "expRight", "Right-side exp: exp(-(ctau-#mu_{1})#lambda)",
      "(@0 > @1) * exp(-(@0-@1)*@2)", RooArgList(ctau, mean_exp1, lambda_exp1));

  // 2. Left-side exponential: exp((ctau - mean2) * lambda)  for ctau < mean2
  //    Indicator (@0 < @1) ensures it vanishes on the right side.
  RooGenericPdf expLeft("expLeft", "Left-side exp: exp((ctau-#mu_{2})#lambda)",
                        "(@0 < @1) * exp((@0-@1)*@2)",
                        RooArgList(ctau, mean_exp2, lambda_exp2));

  // 3. Double-sided exponential: exp(-|ctau - mean3| * lambda)
  //    Symmetric around mean3, continuous at ctau = mean3.
  RooGenericPdf expDoubleSided(
      "expDoubleSided", "Double-sided exp: exp(-|ctau-#mu_{3}|#lambda)",
      "exp(-abs(@0-@2)*@1)", RooArgList(ctau, lambda_exp3, mean_exp3));

  // ------------------------------------------------s-----------------------
  // Same shapes via RooDecay + RooTruthModel (analytically exact, no smearing)
  // RooTruthModel is a delta-function resolution -> pure exponential.
  // RooDecay uses tau = 1/lambda (lifetime in mm). No shift parameter.
  // -----------------------------------------------------------------------
  RooRealVar tau_exp1("tau_exp1", "#tau_{exp} (mm)", 1.5198, 0.05, 3.0);
  RooRealVar tau_exp2("tau_exp2", "#tau_{exp} (mm)", 1.6393, 0.05, 3.0);
  RooRealVar tau_exp3("tau_exp3", "#tau_{exp} (mm)", 2.1, 0.05, 3.0);
  // tau_exp1.setConstant(true);
  // tau_exp2.setConstant(true);
  // tau_exp3.setConstant(true);

  RooTruthModel truthModelExp("truthModelExp", "Delta function (no smearing)",
                              ctau);

  // 1. Right-side: theta(ctau) * exp(-ctau / tau)
  RooDecay decayRight("decayRight", "Right-side exp (SingleSided)", ctau,
                      tau_exp1, truthModelExp, RooDecay::SingleSided);

  // 2. Left-side: theta(-ctau) * exp(+ctau / tau)
  RooDecay decayLeft("decayLeft", "Left-side exp (Flipped)", ctau, tau_exp2,
                     truthModelExp, RooDecay::Flipped);

  // 3. Double-sided: exp(-|ctau| / tau)
  RooDecay decayDoubleSided("decayDoubleSided",
                            "Double-sided exp (DoubleSided)", ctau, tau_exp3,
                            truthModelExp, RooDecay::DoubleSided);

  // RooAddPdf resoModel("resoModel", "4-Gauss resolution (prompt)",
  //                     RooArgList(g1, g2, g3, g4),
  //                     RooArgList(frac1, frac2, frac3),
  //                     true);  // recursive fractions
  // RooAddPdf resoModel("resoModel", "3-Gauss resolution (prompt)",
  //                     RooArgList(g1, g2, g3), RooArgList(frac1, frac2),
  //                     true);  // recursive fractions

  // RooAddPdf resoModel("resoModel", "2-Gauss + CB2 resolution (prompt)",
  //                     RooArgList(g1, g2, cb2), RooArgList(frac1, frac2),
  //                     true);  // recursive fractions
  RooAddPdf resoModel("resoModel", "1-Gauss + CB2", RooArgList(g1, cb2),
                      RooArgList(frac1));  // recursive fractions
  // When using a single PDF, alias it directly (no RooAddPdf wrapper needed)
  // RooAbsPdf& resoModel = cb2;
  // RooAddPdf resoModel("resoModel", "2-Gauss resolution (prompt)",
  //                     RooArgList(g1, g2),
  //                     RooArgList(frac1));  // recursive fractions
  // RooAddPdf resoModel("resoModel", "1-Gauss + 2exp resolution",
  //                     RooArgList(g1, decayLeft, decayRight),
  //                     RooArgList(frac1, frac2),
  //                     true);  // recursive fractions
  // RooAddPdf resoModel("resoModel", "cb2 + 2exp resolution",
  //                     RooArgList(cb2, expLeft, expRight),
  //                     RooArgList(frac1, frac2),
  //                     true);  // recursive fractions
  // RooAddPdf resoModel("resoModel", "cb2 + exp double sided",
  //                     RooArgList(cb2, expDoubleSided), RooArgList(frac1));
  // RooAddPdf resoModel("resoModel", "G1 + CB2 + Eds + G4",
  //                     RooArgList(g1, cb2, expDoubleSided, g4),
  //                     RooArgList(frac1, frac2, frac3), true);
  // RooAddPdf resoModel("resoModel", "G1 + CB2 + G4", RooArgList(g1, cb2, g4),
  //                     RooArgList(frac1, frac2), true);
  // RooAddPdf resoModel("resoModel", "gauss + cb2 + exp double sided",
  //                     RooArgList(g1, cb2, expDoubleSided),
  //                     RooArgList(frac1, frac2), true);
  // RooAddPdf resoModel("resoModel", "gauss + cb2", RooArgList(g1, cb2),
  //                     RooArgList(frac1));
  // RooAddPdf resoModel("resoModel", "2gauss + exp double sided",
  //                     RooArgList(g1, g2, expDoubleSided),
  //                     RooArgList(frac1, frac2), true);
  // RooAddPdf resoModel("resoModel", "3gauss + exp double sided",
  //                     RooArgList(g1, g2, g3, expDoubleSided),
  //                     RooArgList(frac1, frac2, frac3), true);
  // RooAddPdf resoModel("resoModel", "3gauss", RooArgList(g1, g2, g3),
  //                     RooArgList(frac1, frac2), true);
  // RooAddPdf resoModel("resoModel", "4gauss", RooArgList(g1, g2, g3, g4),
  //                     RooArgList(frac1, frac2, frac3), true);
  // RooAddPdf resoModel("resoModel", "5gauss", RooArgList(g1, g2, g3, g4, g5),
  //                     RooArgList(frac1, frac2, frac3, frac4), true);
  // RooAddPdf resoModel("resoModel", "cb2 + exp double sided",
  //                     RooArgList(cb2, decayDoubleSided), RooArgList(frac1));
  // RooAddPdf resoModel("resoModel", "cb2 + 2exp resolution",
  //                     RooArgList(cb2, decayLeft, decayRight),
  //                     RooArgList(frac1, frac2),
  //                     true);  // recursive fractions
  // RooAddPdf resoModel("resoModel", "1-Gauss + 2exp resolution",
  //                     RooArgList(g1, cb2, decayLeft, decayRight),
  //                     RooArgList(frac1, frac2, frac3),
  //                     true);  // recursive fractions
  // RooAbsPdf& resoModel = g4;
  RooAbsPdf& tails = g4;
  // RooAddPdf tails("tails", "3&4-Gauss", RooArgList(g3, g4),
  //                 RooArgList(frac3));  // recursive fractions

  // -----------------------------------------------------------------------
  // Alternative tails: Double-flipped exponential (symmetric or asymmetric)
  // Uncomment one of the options below to use instead of Gaussian tails
  // -----------------------------------------------------------------------
  // --- Option 1: Symmetric double-sided exponential ---
  // exp(-|ctau|/tau_tail)
  //   RooRealVar tau_tail("tau_tail", "#tau_{tail}", 2.5, 0.6, 5.5);
  //   RooGenericPdf tailsExp("tails", "Double-sided exponential tails",
  //                          "exp(-abs(@0-@1)/@2)",
  //                          RooArgList(ctau, reso_mean4, tau_tail));

  //   RooAddPdf tails("tails", "3&4-Gauss", RooArgList(g3, tailsExp),
  //                   RooArgList(frac3));  // recursive fractions

  // --- Option 2: Asymmetric double-flipped exponential ---
  // exp(-ctau/tau_pos) for ctau >= 0, exp(ctau/tau_neg) for ctau < 0
  // RooRealVar tau_tail_pos("tau_tail_pos", "#tau_{tail}^{+}", 0.5, 0.1, 2.0);
  // RooRealVar tau_tail_neg("tau_tail_neg", "#tau_{tail}^{-}", 0.3, 0.05, 1.5);
  // RooGenericPdf tailsExp("tailsExp", "Asymmetric double-flipped exponential",
  //                        "(@0>=0)*exp(-@0/@1) + (@0<0)*exp(@0/@2)",
  //                        RooArgList(ctau, tau_tail_pos, tau_tail_neg));

  // --- Option 3: Asymmetric with resolution convolution (RooDecay) ---
  // Uses RooDecay::DoubleSided for analytical convolution with resolution
  //   RooRealVar tau_tail_decay("tau_tail_decay", "#tau_{tail}", 2.5,
  //   0.5, 4.0); RooGaussModel* gm_tail =
  //       new RooGaussModel("gm_tail", "GM tail", ctau, reso_mean4, sigma3);
  //   RooDecay* tails_decay =
  //       new RooDecay("tails_decay", "Double-sided decay tail", ctau,
  //                    tau_tail_decay, *gm_tail, RooDecay::DoubleSided);

  // --- Option 4: Pure double-sided exponential (RooDecay + RooTruthModel) ---
  // Uses RooTruthModel (delta function) - no resolution smearing, pure
  // exp(-|t|/tau) RooRealVar tau_tail_pure("tau_tail_pure", "#tau_{tail}", 0.5,
  // 0.1, 2.0); RooTruthModel truthModel("truthModel", "Delta function", ctau);
  // RooDecay tailsDecayPure("tailsDecayPure", "Pure double-sided exponential",
  //                         ctau, tau_tail_pure, truthModel,
  //                         RooDecay::DoubleSided);
  // -----------------------------------------------------------------------

  //   RooAddPdf resoModel(
  //       "resoModel", "4-Gauss resolution (prompt)", RooArgList(g1, g2, g3,
  //       g4), RooArgList(frac1, frac2, frac3));  // non recursive fractions
  //   RooAddPdf resoModel("resoModel", "4-Gauss resolution (prompt)",
  //                       RooArgList(g1, g2, g3),
  //                       RooArgList(frac1, frac2));  // non recursive

  // ========================================================================
  // 4. Prompt and Non-prompt parameters (like exp_conv_gauss_CB2.C)
  // ========================================================================
  // f_B = non-prompt fraction (f_sig in exp_conv_gauss_CB2.C is similar
  // concept)
  RooRealVar f_B("f_B", "f_{B}", 0.09, 0.0, 1.0);  // wider range
  // f_B.setConstant(true);
  RooRealVar tauB("tauB", "#tau_{B}", 0.475, 0.35, 0.6);
  // tauB.setConstant(true);
  // RooGenericPdf expoNP("expo_NP", "NP exponential", "(@0>0.0)*exp(-@0/@1)",
  //                      RooArgList(ctau, tauB));
  // Use RooDecay::SingleSided with RooTruthModel for a proper single-sided
  // exponential decay: theta(t)*exp(-t/tauB). No decomposition of the
  // resolution model is needed when using RooFFTConvPdf (FFT handles the
  // full RooAddPdf at once).
  RooDecay expoNP("expo_NP", "NP exponential (SingleSided)", ctau, tauB,
                  truthModelExp, RooDecay::SingleSided);

  double totalEvents = dataReduced->sumEntries();
  RooRealVar Nprompt("Nprompt", "Prompt yield", 0.8 * totalEvents, 0.0,
                     totalEvents);
  RooRealVar Nnonprompt("Nnonprompt", "Non-Prompt yield", 0.1 * totalEvents,
                        0.0, totalEvents);
  RooRealVar Ntails("Ntails", "tails yield", 0.1 * totalEvents, 0.0,
                    totalEvents);
  // ========================================================================
  // 5. Choose convolution method
  //    Comment/uncomment ONE of the two blocks below.
  // ========================================================================

  // --- Main mode (mainMode parameter) ------------------------------------
  // "HybridConv"    → fittedPdf = sngModel (Hybrid Decay+FFT convolution)
  // "ResoConvWEds"  → fittedPdf = resoConv: G1 + CB2 + Eds
  // "ResoConvWOEds" → fittedPdf = resoConv: G1 + CB2  (no Eds)
  // "RooHistPdf"    → non-parametric histogram PDF from sWeighted data
  //                   subMode = "order0" | "order1" | "order2" (interpolation)
  // "RooKeysPdf"    → non-parametric KDE PDF from sWeighted data
  //                   subMode = "NoMirror" | "MirrorLeft" | "MirrorRight" |
  //                             "MirrorBoth" | "MirrorAsymLeft" |
  //                             "MirrorAsymRight" | "MirrorAsymBoth"
  bool useHybridConv = (mainModeStr == "HybridConv");
  bool resoWithEds = (mainModeStr == "ResoConvWEds");
  bool useHistPdf = (mainModeStr == "RooHistPdf");
  bool useKeysPdf = (mainModeStr == "RooKeysPdf");
  // "SingleGaussian" → fit a single RooGaussian (G1: reso_mean1 + sigma_G1)
  // mean and sigma_G1 free (ranges from JSON), normalization not extended
  bool useSingleGauss = (mainModeStr == "SingleGaussian");

  // --- Sub-mode (subMode parameter) --------------------------------------
  // "hybridWOEds"    → G1 + CB2 + G4 (no exponential double-sided)
  // "hybridWEds"     → G1 + CB2 + Eds + G4 (full hybrid with Eds)
  // "hybridAsymDecay"→ asymmetric decay convolution (R ⊗ AsymExp)
  bool hybridNoEds = (std::string(subMode) == "hybridWOEds");
  bool hybridAsymDecay = (std::string(subMode) == "hybridAsymDecay");
  bool hybridAsymExpLeft = false;  // set from tauLambda.active in JSON

  // Hybrid Decay+FFT components
  RooAddPdf* hybridReso = nullptr;       // Mixed resolution (prompt)
  RooGaussModel* hgm1 = nullptr;         // G1 GaussModel for RooDecay
  RooGaussModel* hgm2 = nullptr;         // G2 GaussModel for RooDecay
  RooGaussModel* hgm4 = nullptr;         // G4 GaussModel for RooDecay
  RooDecay* hdec1 = nullptr;             // Ess ⊗ G1 (analytical)
  RooDecay* hdec2 = nullptr;             // Ess ⊗ G2 (analytical)
  RooDecay* hdec4 = nullptr;             // Ess ⊗ G4 (analytical)
  RooGenericPdf* hexpoNP_cb = nullptr;   // Ess copy for CB2 FFT
  RooGenericPdf* hexpoNP_eds = nullptr;  // Ess copy for Eds FFT
  RooCrystalBall* hcb2_conv = nullptr;   // CB2 copy for FFT kernel
  RooGenericPdf* heds_conv = nullptr;    // Eds copy for FFT kernel
  RooFFTConvPdf* hfft_cb2 = nullptr;     // Ess ⊗ CB2 (FFT)
  RooFFTConvPdf* hfft_g1 = nullptr;      // Ess ⊗ G1 (FFT)
  RooFFTConvPdf* hfft_g4 = nullptr;      // Ess ⊗ G4 (FFT)
  RooFFTConvPdf* hfft_eds = nullptr;     // Ess ⊗ Eds (FFT)
  RooAddPdf* hNPconv = nullptr;          // Combined NP (hybrid)

  // --- Asymmetric decay (hybridAsymDecay) pointers ---
  // Physics kernel: f_pr·δ + (1-f_pr)·f_R·ExpRight + (1-f_pr)·(1-f_R)·ExpLeft
  // R = G1 + CB2 + Eds [+ G4]  (Eds IS inside resolution)
  // Model: f_pr·R + (1-f_pr)·f_R·(R⊗ExpRight) + (1-f_pr)·(1-f_R)·(R⊗ExpLeft)
  RooRealVar* tauD_var = nullptr;       // τ_D: right-side decay lifetime
  RooRealVar* tauLambda_var = nullptr;  // τ_λ: left-side decay lifetime
  RooRealVar* f_pr_var = nullptr;       // f_pr: prompt fraction (recursive)
  RooRealVar* f_R_var =
      nullptr;  // f_R: NP-right fraction (recursive within non-prompt)
  RooDecay* asymDecR_G1 = nullptr;         // G1 ⊗ ExpRight (analytical)
  RooDecay* asymDecL_G1 = nullptr;         // G1 ⊗ ExpLeft  (analytical)
  RooDecay* asymDecR_G4 = nullptr;         // G4 ⊗ ExpRight (analytical)
  RooDecay* asymDecL_G4 = nullptr;         // G4 ⊗ ExpLeft  (analytical)
  RooGenericPdf* asymExpR_cb = nullptr;    // ExpRight kernel for CB2 FFT
  RooGenericPdf* asymExpL_cb = nullptr;    // ExpLeft  kernel for CB2 FFT
  RooGenericPdf* asymExpR_eds = nullptr;   // ExpRight kernel for Eds FFT
  RooGenericPdf* asymExpL_eds = nullptr;   // ExpLeft  kernel for Eds FFT
  RooGenericPdf* asymEds_conv = nullptr;   // Eds PDF copy for FFT
  RooFFTConvPdf* asymFFT_cb2_R = nullptr;  // CB2 ⊗ ExpRight (FFT)
  RooFFTConvPdf* asymFFT_cb2_L = nullptr;  // CB2 ⊗ ExpLeft  (FFT)
  RooFFTConvPdf* asymFFT_eds_R = nullptr;  // Eds ⊗ ExpRight (FFT)
  RooFFTConvPdf* asymFFT_eds_L = nullptr;  // Eds ⊗ ExpLeft  (FFT)
  RooAddPdf* asymConvRight = nullptr;      // R ⊗ ExpRight (combined)
  RooAddPdf* asymConvLeft = nullptr;       // R ⊗ ExpLeft  (combined)
  RooCrystalBall* asymCB2_conv = nullptr;  // CB2 copy for asym FFT

  RooAddPdf* sngModel = nullptr;
  RooAbsPdf* resoConv = nullptr;  // fitted PDF for "ResoConv" mainMode

  // Non-parametric PDF pointers (RooHistPdf / RooKeysPdf modes)
  RooDataHist* histPdfData = nullptr;
  RooHistPdf* histPdf = nullptr;
  RooKeysPdf* keysPdf = nullptr;
  std::string keysMirrorStr = "NoMirror";  // resolved mirror name (for legend)
  double keysRhoFinal = 1.0;               // resolved rho (for legend)

  // ========================================================================
  // Build resoConv: the resolution model used directly as fittedPdf when
  // mainMode is "ResoConvWEds", "ResoConvWOEds", or "SingleGaussian".
  //   "SingleGaussian" → pure G1 (RooGaussian, reso_mean1 + sigma_G1 free)
  //   "ResoConvWEds"   → G1 + CB2 + Eds
  //   "ResoConvWOEds"  → G1 + CB2  (no Eds)
  // ========================================================================
  if (useSingleGauss) {
    // Assign g1 directly — no allocation, no ownership transfer
    resoConv = &g1;
    std::cout << "\n===== SINGLE GAUSSIAN MODE =====" << std::endl;
    std::cout << "  Fitting: G1 = RooGaussian(reso_mean1, sigma_G1)"
              << std::endl;
  } else if (!useHybridConv && !useHistPdf && !useKeysPdf) {
    if (resoWithEds) {
      resoConv = new RooAddPdf("resoConv", "G1 + CB2 + Eds (ResoConv)",
                               RooArgList(g1, cb2, expDoubleSided),
                               RooArgList(frac1, frac2), true);
    } else {
      resoConv = new RooAddPdf("resoConv", "G1 + CB2 (ResoConv)",
                               RooArgList(g1, cb2), RooArgList(frac1));
    }
  }

  // ========================================================================
  // RooHistPdf mode: non-parametric histogram PDF from sWeighted data
  // ========================================================================
  // sWeights can be negative (SPlot artefact). Strategy:
  //   1. Fill a TH1D from dataReduced using event weights via FillHistogram.
  //   2. Clamp negative bins to 0 (they are statistical noise, not physics).
  //   3. Wrap in RooDataHist → RooHistPdf.
  // Auto-binning via Scott's rule on the effective sample size (ESS):
  //   ESS = (ΣW)² / ΣW²   (accounts for weight variance from SPlot)
  //   nBins = max(10, round(ESS^(1/3) * (ctauMax-ctauMin) / (3.5*σ_data)))
  // subMode controls the interpolation order passed to RooHistPdf:
  //   "order0" → 0 (nearest-neighbour step function)
  //   "order1" → 1 (linear, default, smooth)
  //   "order2" → 2 (quadratic, smoothest)
  // ========================================================================
  if (useHistPdf) {
    std::cout << "\n===== RooHistPdf (non-parametric histogram) ====="
              << std::endl;

    // --- Interpolation order from subMode ---
    int interpOrder = 1;
    if (subModeStr == "order0")
      interpOrder = 0;
    else if (subModeStr == "order2")
      interpOrder = 2;
    std::cout << "  Interpolation order: " << interpOrder
              << " (subMode=" << subModeStr << ")" << std::endl;

    // --- Compute ESS = (ΣW)² / ΣW² ---
    double sumW = 0.0, sumW2 = 0.0;
    for (int iev = 0; iev < dataReduced->numEntries(); ++iev) {
      dataReduced->get(iev);
      double w = dataReduced->weight();
      sumW += w;
      sumW2 += w * w;
    }
    double ess =
        (sumW2 > 0.0) ? (sumW * sumW) / sumW2 : dataReduced->numEntries();
    std::cout << "  sumW=" << sumW << "  sumW2=" << sumW2 << "  ESS=" << ess
              << std::endl;

    // --- Determine nBins: JSON "nBins" overrides Scott's rule ---
    double range = ctauMax - ctauMin;
    int nBinsHist = 0;
    bool customBins = cfg.contains("nBins") && !cfg["nBins"].is_null();
    if (customBins) {
      nBinsHist = cfg["nBins"].get<int>();
      std::cout << "  nBinsHist=" << nBinsHist << " (from JSON nBins field)"
                << std::endl;
    } else {
      // Scott's rule: h = 3.5*σ*ESS^(-1/3), nBins = range/h
      double mean_w = 0.0;
      for (int iev = 0; iev < dataReduced->numEntries(); ++iev) {
        const RooArgSet* row = dataReduced->get(iev);
        double val = static_cast<RooRealVar*>(row->find("ctau"))->getVal();
        mean_w += dataReduced->weight() * val;
      }
      mean_w /= sumW;
      double var_w = 0.0;
      for (int iev = 0; iev < dataReduced->numEntries(); ++iev) {
        const RooArgSet* row = dataReduced->get(iev);
        double val = static_cast<RooRealVar*>(row->find("ctau"))->getVal();
        double d = val - mean_w;
        var_w += dataReduced->weight() * d * d;
      }
      var_w /= sumW;
      double sigma_w = std::sqrt(std::max(var_w, 1e-6));
      double h_scott = 3.5 * sigma_w * std::pow(ess, -1.0 / 3.0);
      nBinsHist = std::max(10, static_cast<int>(std::round(range / h_scott)));
      // Clamp to ESS/bin in [10, 30] — robust against asymmetric distributions
      int nBinsMin = std::max(10, static_cast<int>(ess / 30.0));
      int nBinsMax = std::max(10, static_cast<int>(ess / 10.0));
      nBinsHist = std::clamp(nBinsHist, nBinsMin, nBinsMax);
      std::cout << "  sigma_w=" << sigma_w << "  h_Scott=" << h_scott
                << "  nBinsHist=" << nBinsHist
                << "  ESS/bin=" << ess / nBinsHist << std::endl;
    }

    // --- Fill TH1D with sWeights, then clamp negatives ---
    TH1D* hCtau =
        new TH1D("hCtau_hist", "sWeighted ctau", nBinsHist, ctauMin, ctauMax);
    hCtau->SetDirectory(nullptr);  // prevent ROOT global directory registration
    hCtau->Sumw2();
    for (int iev = 0; iev < dataReduced->numEntries(); ++iev) {
      const RooArgSet* row = dataReduced->get(iev);
      double val = static_cast<RooRealVar*>(row->find("ctau"))->getVal();
      hCtau->Fill(val, dataReduced->weight());
    }
    int nNeg = 0;
    for (int ib = 1; ib <= hCtau->GetNbinsX(); ++ib) {
      if (hCtau->GetBinContent(ib) < 0.0) {
        hCtau->SetBinContent(ib, 0.0);
        hCtau->SetBinError(ib, 0.0);
        ++nNeg;
      }
    }
    if (nNeg > 0)
      std::cout << "  WARNING: clamped " << nNeg
                << " negative bins to 0 (sWeight artefact)" << std::endl;

    histPdfData = new RooDataHist("histPdfData", "sWeighted ctau histogram",
                                  RooArgList(ctau), hCtau);
    delete hCtau;  // RooDataHist copies the bin contents; safe to delete now
    histPdf = new RooHistPdf("histPdf", "RooHistPdf (non-param)",
                             RooArgSet(ctau), *histPdfData, interpOrder);
    std::cout << "  RooHistPdf built: " << nBinsHist << " bins, "
              << "interpOrder=" << interpOrder << std::endl;
  }

  // ========================================================================
  // RooKeysPdf mode: non-parametric KDE from sWeighted data
  // ========================================================================
  // RooKeysPdf uses the Cranmer adaptive KDE, fully weight-aware — it handles
  // negative sWeights correctly by including them in the kernel sum.
  //
  // subMode encodes both the boundary mirror option and an optional rho
  // override in a single compact string:
  //
  //   "<MirrorOption>"        → use keys_rho from JSON (e.g. "MirrorBoth")
  //   "<MirrorOption>,<rho>"  → override rho inline  (e.g. "MirrorBoth,2")
  //
  // Mirror options (boundary treatment):
  //   "NoMirror"        → RooKeysPdf::NoMirror       (default, no correction)
  //   "MirrorLeft"      → RooKeysPdf::MirrorLeft      (reflect left edge)
  //   "MirrorRight"     → RooKeysPdf::MirrorRight     (reflect right edge)
  //   "MirrorBoth"      → RooKeysPdf::MirrorBoth      (reflect both —
  //   recommended) "MirrorAsymLeft"  → RooKeysPdf::MirrorAsymLeft  (asymmetric
  //   left) "MirrorAsymRight" → RooKeysPdf::MirrorAsymRight (asymmetric right)
  //   "MirrorAsymBoth"  → RooKeysPdf::MirrorAsymBoth  (asymmetric both)
  //
  // rho is the bandwidth scale factor (the only smoothness knob in RooKeysPdf):
  //   rho = 1.0  → standard KEYS bandwidth (default)
  //   rho < 1.0  → narrower kernel, more peaked, risk of overfitting
  //   rho > 1.0  → wider kernel, smoother, risk of over-smoothing
  //   (what ROOT examples show as integer "2" is rho=2.0 as a double literal)
  //
  // Resolution order of rho:
  //   1. Comma suffix in subMode  (e.g. "MirrorBoth,2"  → rho=2.0)
  //   2. keys_rho field in JSON   (e.g. "keys_rho": 1.5 → rho=1.5)
  //   3. Hard-coded default 1.0
  // ========================================================================
  if (useKeysPdf) {
    std::cout << "\n===== RooKeysPdf (non-parametric KDE) =====" << std::endl;

    // --- Parse subMode: split on first comma into mirror name + optional rho
    // ---
    std::string mirrorStr = subModeStr;
    double keysRhoInline = -1.0;  // sentinel: not set inline
    {
      auto commaPos = subModeStr.find(',');
      if (commaPos != std::string::npos) {
        mirrorStr = subModeStr.substr(0, commaPos);
        std::string rhoStr = subModeStr.substr(commaPos + 1);
        try {
          keysRhoInline = std::stod(rhoStr);
        } catch (...) {
          std::cerr << "WARNING: Cannot parse rho from subMode suffix '"
                    << rhoStr << "', falling back to keys_rho field."
                    << std::endl;
        }
      }
    }

    // --- Resolve rho: inline > JSON field > default 1.0 ---
    double keysRho =
        (keysRhoInline > 0.0) ? keysRhoInline : cfg.value("keys_rho", 1.0);

    // --- Map mirror string to RooKeysPdf::Mirror enum ---
    RooKeysPdf::Mirror mirrorOpt = RooKeysPdf::NoMirror;
    if (mirrorStr == "MirrorLeft")
      mirrorOpt = RooKeysPdf::MirrorLeft;
    else if (mirrorStr == "MirrorRight")
      mirrorOpt = RooKeysPdf::MirrorRight;
    else if (mirrorStr == "MirrorBoth")
      mirrorOpt = RooKeysPdf::MirrorBoth;
    else if (mirrorStr == "MirrorAsymLeft")
      mirrorOpt = RooKeysPdf::MirrorAsymLeft;
    else if (mirrorStr == "MirrorAsymRight")
      mirrorOpt = RooKeysPdf::MirrorAsymRight;
    else if (mirrorStr == "MirrorAsymBoth")
      mirrorOpt = RooKeysPdf::MirrorAsymBoth;
    else if (mirrorStr != "NoMirror")
      std::cerr << "WARNING: Unknown mirror option '" << mirrorStr
                << "', defaulting to NoMirror." << std::endl;

    keysMirrorStr = mirrorStr;
    keysRhoFinal = keysRho;

    std::cout << "  Mirror option: " << mirrorStr << "  rho=" << keysRho
              << (keysRhoInline > 0.0 ? " (from subMode inline)"
                                      : " (from keys_rho field)")
              << std::endl;

    // RooKeysPdf normalisation = Σ wᵢ · kernel_integral_i.  Negative sWeights
    // (SPlot artefact) can drive this sum to zero or negative, producing NaN
    // everywhere.  Clamp negative weights to 0 before building the KDE,
    // exactly as RooHistPdf does with its negative-bin clamping.
    RooDataSet* dataForKeys =
        new RooDataSet("dataForKeys", "sWeighted ctau (w>=0)", RooArgSet(ctau),
                       RooFit::WeightVar());
    int nClamped = 0;
    for (int iev = 0; iev < dataReduced->numEntries(); ++iev) {
      const RooArgSet* row = dataReduced->get(iev);
      double w = dataReduced->weight();
      double val = static_cast<RooRealVar*>(row->find("ctau"))->getVal();
      if (w <= 0.0) {
        ++nClamped;
        continue;
      }  // skip negative-weight events
      ctau.setVal(val);
      dataForKeys->add(RooArgSet(ctau), w);
    }
    if (nClamped > 0)
      std::cout << "  WARNING: skipped " << nClamped
                << " events with w<=0 for RooKeysPdf (sWeight artefact)"
                << std::endl;
    std::cout << "  Building RooKeysPdf from " << dataForKeys->numEntries()
              << " entries (sum of weights = " << dataForKeys->sumEntries()
              << ")..." << std::flush;
    auto t0 = std::chrono::steady_clock::now();

    RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG);
    keysPdf = new RooKeysPdf("keysPdf", "RooKeysPdf (non-param)", ctau,
                             *dataForKeys, mirrorOpt, keysRho);
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    delete dataForKeys;

    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << " done (" << std::fixed << std::setprecision(1) << elapsed
              << " s)" << std::endl;
  }

  // R+T model components
  RooAbsPdf* promptRT = nullptr;     // P = h*R + (1-h)*T
  RooAbsPdf* nonpromptRT = nullptr;  // NP = h*conv(Exp,R) + (1-h)*T

  if (useHybridConv) {
    // =================================================================
    // Option E: Hybrid Decay+FFT convolution
    // =================================================================
    // Exploits linearity of convolution:
    //   M(x) = [f_B * Ess(x) + (1-f_B) * delta(x)] ⊗ R(x)
    //        = f_B * (Ess ⊗ R) + (1-f_B) * R
    //
    // Expanding Ess ⊗ R by distributing across resolution components:
    //   Ess ⊗ R = c1*(Ess⊗G1) + c2*(Ess⊗CB2) + c3*(Ess⊗Eds) + c4*(Ess⊗G2)
    //
    // Route each convolution to the optimal method:
    //   Ess ⊗ G1  → RooDecay (analytical via complex error functions)
    //   Ess ⊗ G2  → RooDecay (analytical via complex error functions)
    //   Ess ⊗ CB2 → RooFFTConvPdf (no closed-form solution)
    //   Ess ⊗ Eds → RooFFTConvPdf (no closed-form solution)
    // =================================================================
    std::cout << "\n===== HYBRID DECAY+FFT CONVOLUTION =====" << std::endl;

    int fftBins = 50000;
    double fftBuffer = 0.65;
    ctau.setBins(fftBins, "cache");

    // --- Analytical convolutions: Ess ⊗ Gaussian ---
    // RooDecay REQUIRES RooGaussModel (inherits from RooResolutionModel),
    // NOT RooGaussian. Parameters are the SAME RooRealVar pointers as
    // the standalone g1/g4, ensuring simultaneous fitting.
    hgm1 = new RooGaussModel("hgm1", "GM1 (hybrid)", ctau, reso_mean1, sigma1);
    hgm2 = new RooGaussModel("hgm4", "GM4 (hybrid)", ctau, reso_mean2, sigma2);
    hgm4 = new RooGaussModel("hgm4", "GM4 (hybrid)", ctau, reso_mean4, sigma4);

    hdec1 = new RooDecay("hdec1", "Ess #otimes G1 (analytical)", ctau, tauB,
                         *hgm1, RooDecay::SingleSided);
    hdec2 = new RooDecay("hdec2", "Ess #otimes G2 (analytical)", ctau, tauB,
                         *hgm2, RooDecay::SingleSided);
    hdec4 = new RooDecay("hdec4", "Ess #otimes G4 (analytical)", ctau, tauB,
                         *hgm4, RooDecay::SingleSided);

    // --- Numerical convolutions: Ess ⊗ non-Gaussian ---
    // Separate copies of the physics PDF (Ess) for each FFT to avoid
    // normalization cross-talk between convolution products.
    // Using RooGenericPdf with Heaviside step: theta(ctau)*exp(-ctau/tauB)
    hexpoNP_cb =
        new RooGenericPdf("hexpoNP_cb", "NP exp (hybrid, for CB2)",
                          "(@0>=0.0)*exp(-@0/@1)", RooArgList(ctau, tauB));

    // Separate resolution component copies for FFT kernels (avoid
    // normalization cross-talk with standalone hybridReso).
    // CB2: same parameters as cb2, different PDF object.
    // hcb2_conv =
    //     new RooCrystalBall("hcb2_conv", "CB2 conv (hybrid)", ctau, cb2_mean,
    //                        cb2_sigma, cbAlpha1, cbN1, cbAlpha2, cbN2);
    hcb2_conv =
        new RooCrystalBall("hcb2_conv", "CB2 conv (hybrid)", ctau, cb2MeanRef,
                           cb2_sigma, cbAlpha1, cbN1, cbAlpha2Ref, cbN2Ref);
    // hcb2_conv =
    //     new RooCrystalBall("hcb2_conv", "CB2 conv (hybrid)", ctau, cb2_mean,
    //                        cb2_sigma, cbAlpha1, cbN1, cbAlpha1, cbN1);
    // hcb2_conv =
    //     new RooCrystalBall("hcb2_conv", "CB2 conv (hybrid)", ctau,
    //     reso_mean1,
    //                        cb2_sigma, cbAlpha1, cbN1, cbAlpha1, cbN1);
    // hcb2_conv =
    //     new RooCrystalBall("hcb2_conv", "CB2 conv (hybrid)", ctau,
    //     reso_mean1,
    //                        cb2_sigma, cbAlpha1, cbN1, cbAlpha2, cbN2);

    hfft_cb2 = new RooFFTConvPdf("hfft_cb2", "Ess #otimes CB2 (FFT)", ctau,
                                 *hexpoNP_cb, *hcb2_conv);
    // hfft_cb2 = new RooFFTConvPdf("hfft_cb2", "Ess #otimes CB2 (FFT)", ctau,
    //                              *hexpoNP_cb, cb2);
    hfft_g1 = new RooFFTConvPdf("hfft_g1", "Ess #otimes G1 (FFT)", ctau,
                                *hexpoNP_cb, g1);
    hfft_g4 = new RooFFTConvPdf("hfft_g4", "Ess #otimes G4 (FFT)", ctau,
                                *hexpoNP_cb, g4);

    hfft_cb2->setBufferFraction(fftBuffer);

    if (hybridNoEds) {
      // --- Hybrid without Eds: R = c1*G1 + c2*CB2 + (1-c1-c2)*G4 ---
      std::cout << "  Mode: hybridNoEds (no exponential double-sided)"
                << std::endl;
      std::cout << "  Numerical:  Ess x CB2 (FFT)" << std::endl;

      // --- Prompt component: G1 + CB2 [+ G4 if includeG4] ---
      // includeG4=true  → hybridReso = G1 + CB2 + G4  (G4 inside prompt)
      // includeG4=false → hybridReso = G1 + CB2        (G4 separate in
      // sngModel)
      if (includeG4) {
        hybridReso = new RooAddPdf(
            "hybridReso", "G1 + CB2 + G4 (hybrid resolution, no Eds)",
            RooArgList(g1, cb2, g4), RooArgList(frac1, frac2), true);
      } else {
        hybridReso =
            new RooAddPdf("hybridReso", "G1 + CB2 (hybrid resolution, no Eds)",
                          RooArgList(g1, cb2), RooArgList(frac1));
      }

      // --- NP: mirrors hybridReso fractions for linearity: Ess⊗R ---
      // includeG4=true  → hdec1 + hfft_cb2 + hdec4  (frac1, frac2 recursive)
      // includeG4=false → hdec1 + hfft_cb2           (frac1)
      if (includeG4) {
        hNPconv = new RooAddPdf(
            "hNPconv", "NP hybrid (c1*Decay1 + c2*FFT_CB + c3*Decay4, no Eds)",
            RooArgList(*hdec1, *hfft_cb2, *hdec4), RooArgList(frac1, frac2),
            true);
      } else {
        hNPconv = new RooAddPdf(
            "hNPconv", "NP hybrid (c1*Decay1 + c2*FFT_CB, no Eds)",
            RooArgList(*hdec1, *hfft_cb2), RooArgList(frac1));
      }

      // --- Final model (includeG4 from JSON) ---
      // includeG4=true  → G4 is inside hybridReso: f_B*NPconv + (1-f_B)*Prompt
      // includeG4=false → G4 is a separate 3rd component: a1*NP + a2*Prompt +
      // (1-a1-a2)*G4
      if (includeG4) {
        sngModel = new RooAddPdf(
            "sngModel", "f_B*NPconv + (1-f_B)*PromptReso[G1+CB2+G4] (no Eds)",
            RooArgList(*hNPconv, *hybridReso), RooArgList(f_B));
        std::cout << "  Model: f_B*(NP conv) + (1-f_B)*(Prompt[G1+CB2+G4])"
                  << std::endl;
      } else {
        sngModel = new RooAddPdf(
            "sngModel", "a1*NPconv + a2*PromptReso + (1-a1-a2)*G4 (no Eds)",
            RooArgList(*hNPconv, *hybridReso, g4), RooArgList(a1, a2), true);
        std::cout
            << "  Model: a1*(NP conv) + a2*(Prompt[G1+CB2]) + (1-a1-a2)*G4"
            << std::endl;
      }

    } else if (!hybridNoEds && !hybridAsymDecay) {
      // --- Full hybrid: R = c1*G1 + c2*CB2 + c3*Eds + c4*G4 ---
      hexpoNP_eds =
          new RooGenericPdf("hexpoNP_eds", "NP exp (hybrid, for Eds)",
                            "(@0>=0.0)*exp(-@0/@1)", RooArgList(ctau, tauB));

      // Double-sided exponential: same parameters as expDoubleSided.
      heds_conv = new RooGenericPdf("heds_conv", "Eds conv (hybrid)",
                                    "exp(-abs(@0-@2)*@1)",
                                    RooArgList(ctau, lambda_exp3, mean_exp3));

      hfft_eds = new RooFFTConvPdf("hfft_eds", "Ess #otimes Eds (FFT)", ctau,
                                   *hexpoNP_eds, *heds_conv);
      hfft_eds->setBufferFraction(fftBuffer);

      std::cout << "  Numerical:  Ess x CB2 (FFT), Ess x Eds (FFT)"
                << std::endl;

      // --- Prompt component: Mixed resolution model ---
      // includeG4=true  → R = G1 + CB2 + Eds + G4  (G4 inside prompt)
      // includeG4=false → R = G1 + CB2 + Eds        (G4 separate in sngModel)
      if (includeG4) {
        hybridReso = new RooAddPdf(
            "hybridReso", "G1 + CB2 + Eds + G4 (hybrid resolution, with Eds)",
            RooArgList(g1, cb2, expDoubleSided, g4),
            RooArgList(frac1, frac2, frac3), true);  // recursive fractions
      } else {
        hybridReso = new RooAddPdf(
            "hybridReso", "G1 + CB2 + Eds (hybrid resolution, with Eds)",
            RooArgList(g1, cb2, expDoubleSided), RooArgList(frac1, frac2),
            true);  // recursive fractions
      }

      // --- Combine all convolution products ---
      // NP = mirrors hybridReso fractions for linearity: Ess⊗R
      // includeG4=true  → c1*(Ess⊗G1) + c2*(Ess⊗CB2) + c3*(Ess⊗Eds) +
      // c4*(Ess⊗G4) includeG4=false → c1*(Ess⊗G1) + c2*(Ess⊗CB2) + c3*(Ess⊗Eds)
      if (includeG4) {
        hNPconv = new RooAddPdf(
            "hNPconv",
            "NP hybrid (c1*Decay1 + c2*FFT_CB + c3*FFT_Eds + c4*Decay4)",
            RooArgList(*hdec1, *hfft_cb2, *hfft_eds, *hdec4),
            RooArgList(frac1, frac2, frac3), true);  // recursive fractions
      } else {
        hNPconv = new RooAddPdf(
            "hNPconv", "NP hybrid (c1*Decay1 + c2*FFT_CB + c3*FFT_Eds)",
            RooArgList(*hdec1, *hfft_cb2, *hfft_eds), RooArgList(frac1, frac2),
            true);  // recursive fractions
      }

      // --- Final model (includeG4 from JSON) ---
      // includeG4=true  → G4 is inside hybridReso: f_B*NPconv + (1-f_B)*Prompt
      // includeG4=false → G4 is a separate 3rd component: a1*NP + a2*Prompt +
      // (1-a1-a2)*G4
      if (includeG4) {
        sngModel = new RooAddPdf(
            "sngModel",
            "f_B*NPconv + (1-f_B)*PromptReso[G1+CB2+Eds+G4] (with Eds)",
            RooArgList(*hNPconv, *hybridReso), RooArgList(f_B));
        std::cout << "  Model: f_B*(NP conv) + (1-f_B)*(Prompt[G1+CB2+Eds+G4])"
                  << std::endl;
      } else {
        sngModel = new RooAddPdf(
            "sngModel", "a1*NPconv + a2*PromptReso + (1-a1-a2)*G4 (with Eds)",
            RooArgList(*hNPconv, *hybridReso, g4), RooArgList(a1, a2), true);
        std::cout
            << "  Model: a1*(NP conv) + a2*(Prompt[G1+CB2+Eds]) + (1-a1-a2)*G4"
            << std::endl;
      }
    }

    // =================================================================
    // hybridAsymDecay: R ⊗ [f_pr·δ + (1-f_pr)·f_R·ExpRight +
    // (1-f_pr)·(1-f_R)·ExpLeft]
    // =================================================================
    // R = G1 + CB2 + Eds [+ G4]  — Eds IS inside the resolution
    // Distributing convolution linearity:
    //   f_pr·R  +  (1-f_pr)·f_R·(R⊗ExpRight)  +  (1-f_pr)·(1-f_R)·(R⊗ExpLeft)
    // where R⊗Exp = c1·(G1⊗Exp) + c2·(CB2⊗Exp) + c3·(Eds⊗Exp) [+ c4·(G4⊗Exp)]
    if (hybridAsymDecay) {
      std::cout << "\n===== HYBRID ASYMMETRIC DECAY CONVOLUTION ====="
                << std::endl;
      std::cout << "  Kernel: f_pr*delta + (1-f_pr)*f_R*ExpRight + "
                   "(1-f_pr)*(1-f_R)*ExpLeft"
                << std::endl;
      std::cout << "  R = G1 + CB2 + Eds [+ G4]  (Eds IS in resolution)"
                << std::endl;
      std::cout << "  ExpRight = Theta(tz>=0)*Exp(-tz/tauD),"
                   "  ExpLeft = Theta(tz<=0)*Exp(tz/tauLambda)"
                << std::endl;

      // --- Parameters from JSON ---
      {
        const auto& jA = cfg["asymDecay"];
        // tauLambda.active=false → 2-component model: f_pr*R +
        // (1-f_pr)*(R⊗ExpRight) tauLambda.active=true  → 3-component model:
        // adds (1-f_pr)*(1-f_R)*(R⊗ExpLeft)
        hybridAsymExpLeft = jA["tauLambda"]["active"].get<bool>();

        tauD_var = new RooRealVar(
            "tauD", "#tau_{D}", jA["tauD"]["value"].get<double>(),
            jA["tauD"]["min"].get<double>(), jA["tauD"]["max"].get<double>());
        f_pr_var = new RooRealVar(
            "f_pr", "f_{pr}", jA["f_pr"]["value"].get<double>(),
            jA["f_pr"]["min"].get<double>(), jA["f_pr"]["max"].get<double>());
        if (jA["tauD"]["fixed"].get<bool>()) tauD_var->setConstant(true);
        if (jA["f_pr"]["fixed"].get<bool>()) f_pr_var->setConstant(true);

        if (hybridAsymExpLeft) {
          tauLambda_var = new RooRealVar("tauLambda", "#tau_{#lambda}",
                                         jA["tauLambda"]["value"].get<double>(),
                                         jA["tauLambda"]["min"].get<double>(),
                                         jA["tauLambda"]["max"].get<double>());
          f_R_var = new RooRealVar(
              "f_R", "f_{R}", jA["f_R"]["value"].get<double>(),
              jA["f_R"]["min"].get<double>(), jA["f_R"]["max"].get<double>());
          if (jA["tauLambda"]["fixed"].get<bool>())
            tauLambda_var->setConstant(true);
          if (jA["f_R"]["fixed"].get<bool>()) f_R_var->setConstant(true);
        }
      }
      std::cout << "  ExpLeft component: "
                << (hybridAsymExpLeft ? "ACTIVE"
                                      : "DISABLED (delta + ExpRight only)")
                << std::endl;

      int fftBinsAsym = 50000;
      double fftBufferAsym = 0.65;
      ctau.setBins(fftBinsAsym, "cache");

      // --- GaussModels for analytical RooDecay ---
      if (!hgm1)
        hgm1 =
            new RooGaussModel("hgm1", "GM1 (hybrid)", ctau, reso_mean1, sigma1);
      if (includeG4 && !hgm4)
        hgm4 =
            new RooGaussModel("hgm4", "GM4 (hybrid)", ctau, reso_mean4, sigma4);

      // --- Prompt: R = G1 + CB2 + Eds [+ G4] (Eds inside resolution) ---
      if (includeG4) {
        hybridReso = new RooAddPdf("hybridReso",
                                   "G1 + CB2 + Eds + G4 (prompt, with Eds)",
                                   RooArgList(g1, cb2, expDoubleSided, g4),
                                   RooArgList(frac1, frac2, frac3), true);
      } else {
        hybridReso =
            new RooAddPdf("hybridReso", "G1 + CB2 + Eds (prompt, with Eds)",
                          RooArgList(g1, cb2, expDoubleSided),
                          RooArgList(frac1, frac2), true);
      }

      // CB2 copy for FFT (avoid cross-talk with standalone cb2)
      asymCB2_conv = new RooCrystalBall("asymCB2_conv", "CB2 conv (asym)", ctau,
                                        cb2MeanRef, cb2_sigma, cbAlpha1, cbN1,
                                        cbAlpha2Ref, cbN2Ref);
      // Eds copy for FFT (avoid cross-talk with expDoubleSided)
      asymEds_conv = new RooGenericPdf(
          "asymEds_conv", "Eds conv (asym)", "exp(-abs(@0-@2)*@1)",
          RooArgList(ctau, lambda_exp3, mean_exp3));

      // --- Right-side: R ⊗ Θ(τz≥0)·Exp(-τz/τD) ---
      asymDecR_G1 = new RooDecay("asymDecR_G1", "G1 #otimes ExpRight", ctau,
                                 *tauD_var, *hgm1, RooDecay::SingleSided);
      asymExpR_cb = new RooGenericPdf("asymExpR_cb", "ExpRight for CB2 FFT",
                                      "(@0>=0.0)*exp(-@0/@1)",
                                      RooArgList(ctau, *tauD_var));
      asymFFT_cb2_R =
          new RooFFTConvPdf("asymFFT_cb2_R", "CB2 #otimes ExpRight (FFT, asym)",
                            ctau, *asymExpR_cb, *asymCB2_conv);
      asymFFT_cb2_R->setBufferFraction(fftBufferAsym);
      // Separate ExpRight copy for Eds FFT (avoid normalization cross-talk)
      asymExpR_eds = new RooGenericPdf("asymExpR_eds", "ExpRight for Eds FFT",
                                       "(@0>=0.0)*exp(-@0/@1)",
                                       RooArgList(ctau, *tauD_var));
      asymFFT_eds_R =
          new RooFFTConvPdf("asymFFT_eds_R", "Eds #otimes ExpRight (FFT, asym)",
                            ctau, *asymExpR_eds, *asymEds_conv);
      asymFFT_eds_R->setBufferFraction(fftBufferAsym);

      if (includeG4) {
        asymDecR_G4 = new RooDecay("asymDecR_G4", "G4 #otimes ExpRight", ctau,
                                   *tauD_var, *hgm4, RooDecay::SingleSided);
        asymConvRight =
            new RooAddPdf("asymConvRight", "R #otimes ExpRight (with Eds)",
                          RooArgList(*asymDecR_G1, *asymFFT_cb2_R,
                                     *asymFFT_eds_R, *asymDecR_G4),
                          RooArgList(frac1, frac2, frac3), true);
      } else {
        asymConvRight = new RooAddPdf(
            "asymConvRight", "R #otimes ExpRight (with Eds)",
            RooArgList(*asymDecR_G1, *asymFFT_cb2_R, *asymFFT_eds_R),
            RooArgList(frac1, frac2), true);
      }

      if (hybridAsymExpLeft) {
        // --- Left-side: R ⊗ Θ(τz≤0)·Exp(τz/τλ) ---
        asymDecL_G1 = new RooDecay("asymDecL_G1", "G1 #otimes ExpLeft", ctau,
                                   *tauLambda_var, *hgm1, RooDecay::Flipped);
        asymExpL_cb = new RooGenericPdf("asymExpL_cb", "ExpLeft for CB2 FFT",
                                        "(@0<=0.0)*exp(@0/@1)",
                                        RooArgList(ctau, *tauLambda_var));
        asymFFT_cb2_L = new RooFFTConvPdf("asymFFT_cb2_L",
                                          "CB2 #otimes ExpLeft (FFT, asym)",
                                          ctau, *asymExpL_cb, *asymCB2_conv);
        asymFFT_cb2_L->setBufferFraction(fftBufferAsym);
        asymExpL_eds = new RooGenericPdf("asymExpL_eds", "ExpLeft for Eds FFT",
                                         "(@0<=0.0)*exp(@0/@1)",
                                         RooArgList(ctau, *tauLambda_var));
        asymFFT_eds_L = new RooFFTConvPdf("asymFFT_eds_L",
                                          "Eds #otimes ExpLeft (FFT, asym)",
                                          ctau, *asymExpL_eds, *asymEds_conv);
        asymFFT_eds_L->setBufferFraction(fftBufferAsym);

        if (includeG4) {
          asymDecL_G4 = new RooDecay("asymDecL_G4", "G4 #otimes ExpLeft", ctau,
                                     *tauLambda_var, *hgm4, RooDecay::Flipped);
          asymConvLeft =
              new RooAddPdf("asymConvLeft", "R #otimes ExpLeft (with Eds)",
                            RooArgList(*asymDecL_G1, *asymFFT_cb2_L,
                                       *asymFFT_eds_L, *asymDecL_G4),
                            RooArgList(frac1, frac2, frac3), true);
        } else {
          asymConvLeft = new RooAddPdf(
              "asymConvLeft", "R #otimes ExpLeft (with Eds)",
              RooArgList(*asymDecL_G1, *asymFFT_cb2_L, *asymFFT_eds_L),
              RooArgList(frac1, frac2), true);
        }

        // 3-component model: f_pr·R + (1-f_pr)·f_R·(R⊗ExpRight) +
        // (1-f_pr)·(1-f_R)·(R⊗ExpLeft)
        sngModel = new RooAddPdf(
            "sngModel",
            "f_{pr}R + f_{R}(R#otimesExpRight) + "
            "(1-f_{pr}-f_{R})(R#otimesExpLeft)"
            " (hybridAsymDecay)",
            RooArgList(*hybridReso, *asymConvRight, *asymConvLeft),
            RooArgList(*f_pr_var, *f_R_var), true);
      } else {
        // 2-component model: f_pr·R + (1-f_pr)·(R⊗ExpRight)
        sngModel = new RooAddPdf("sngModel",
                                 "f_{pr}R + (1-f_{pr})(R#otimesExpRight) "
                                 "(hybridAsymDecay, no ExpLeft)",
                                 RooArgList(*hybridReso, *asymConvRight),
                                 RooArgList(*f_pr_var));
      }
    }
  }
  // ========================================================================
  // Apply JSON config: override initial values, ranges, and fixed flags.
  // RooFit PDFs hold references to RooRealVars, so updating the vars here
  // (after PDF construction but before the fit) is fully correct.
  // ========================================================================
  {
    auto setFromCfg = [](RooRealVar& var, const nlohmann::json& node) {
      var.setRange(node["min"].get<double>(), node["max"].get<double>());
      var.setVal(node["value"].get<double>());
      if (node["fixed"].get<bool>()) var.setConstant(true);
    };

    const auto& jG1CB2 = cfg["reso_G1_CB2"];
    const auto& jEds = cfg["reso_Eds"];
    const auto& jG4 = cfg["reso_G4"];
    const auto& jNP = cfg["signal_NP"];

    setFromCfg(reso_mean1, jG1CB2["reso_mean1"]);
    if (useCb2OwnMean) setFromCfg(cb2_mean, jG1CB2["cb2_mean"]);
    setFromCfg(sigma1, jG1CB2["sigma_G1"]);
    setFromCfg(cb2_sigma, jG1CB2["sigma_CB2"]);
    setFromCfg(cbAlpha1, jG1CB2["cbAlpha1"]);
    setFromCfg(cbN1, jG1CB2["cbN1"]);
    if (!useCbRatio) {
      setFromCfg(cbAlpha2, jG1CB2["cbAlpha2"]);
      setFromCfg(cbN2, jG1CB2["cbN2"]);
    }
    setFromCfg(frac1, jG1CB2["frac1"]);

    setFromCfg(lambda_exp3, jEds["lambda_exp3"]);
    setFromCfg(mean_exp3, jEds["mean_exp3"]);
    setFromCfg(frac2, jEds["frac2"]);

    setFromCfg(reso_mean4, jG4["reso_mean4"]);
    setFromCfg(sigma4, jG4["sigma_G4"]);

    setFromCfg(tauB, jNP["tauB"]);
    setFromCfg(f_B, jNP["f_B"]);
    setFromCfg(a1, jNP["a1"]);
    setFromCfg(a2, jNP["a2"]);
  }

  if (useHybridConv && !sngModel) {
    std::cerr << "ERROR: Failed to build signal model (HybridConv mode)"
              << std::endl;
    return;
  }
  if (!useHybridConv && !useHistPdf && !useKeysPdf && !resoConv) {
    std::cerr << "ERROR: Failed to build resolution model (ResoConv mode)"
              << std::endl;
    return;
  }
  if (useHistPdf && !histPdf) {
    std::cerr << "ERROR: Failed to build RooHistPdf" << std::endl;
    return;
  }
  if (useKeysPdf && !keysPdf) {
    std::cerr << "ERROR: Failed to build RooKeysPdf" << std::endl;
    return;
  }

  // ========================================================================
  // 5a. Convert to binned data if requested
  // ========================================================================
  double binWidth = 0.1;  // mm (same as plotting)
  int nBins = static_cast<int>((ctauMax - ctauMin) / binWidth);

  RooDataHist* dataBinned = nullptr;
  RooAbsData* fitData = dataReduced;  // default: unbinned

  if (useBinnedData) {
    ctau.setBins(nBins);
    dataBinned = new RooDataHist("dataBinned", "Binned ctau data",
                                 RooArgSet(ctau), *dataReduced);
    fitData = dataBinned;
    std::cout << "Converted to binned data: " << nBins
              << " bins, binWidth = " << binWidth << " mm" << std::endl;
    std::cout << "  Binned entries: " << dataBinned->numEntries()
              << ", sum of weights: " << dataBinned->sumEntries() << std::endl;
  }

  // ========================================================================
  // 5b. Configure numerical integrator (like SUBARooFitter)
  // ========================================================================
  // Use adaptive Gauss-Kronrod integrator with tighter tolerances
  // This fixes the "integral did not converge" warnings for expo_NP
  RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-9);
  RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-9);
  RooAbsReal::defaultIntegratorConfig()
      ->getConfigSection("RooIntegrator1D")
      .setRealValue("maxSteps", 50);
  RooAbsReal::defaultIntegratorConfig()->method1D().setLabel(
      "RooAdaptiveGaussKronrodIntegrator1D");
  std::cout << "Integration config: AdaptiveGaussKronrod, eps=1e-9, maxSteps=50"
            << std::endl;

  RooFitResult* result = nullptr;
  RooAbsPdf* fittedPdf = nullptr;  // Track which PDF was actually fit

  std::cout << "\n===== SINGLE-SHOT FIT =====" << std::endl;
  std::cout << "  mainMode = " << mainMode << ",  subMode = " << subMode
            << std::endl;

  if (useHybridConv) {
    // --- HybridConv: fit the full sngModel (NP exp ⊗ R + prompt R) ---
    fittedPdf = sngModel;
    result = sngModel->fitTo(
        *dataReduced, RooFit::Extended(false), RooFit::SumW2Error(true),
        RooFit::Minimizer("Minuit2", "Migrad"), RooFit::Strategy(1),
        RooFit::Offset(true), RooFit::Optimize(true), RooFit::NumCPU(8),
        RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(1),
        RooFit::Save(true));
  } else if (useHistPdf) {
    // RooHistPdf has no free parameters — skip fitTo entirely to avoid the
    // "zero parameters" minimizer error and subsequent segfault on
    // result->Print.
    fittedPdf = histPdf;
    std::cout << "  RooHistPdf: non-parametric mode, skipping fitTo."
              << std::endl;
  } else if (useKeysPdf) {
    // RooKeysPdf has no free parameters — same reason as RooHistPdf above.
    fittedPdf = keysPdf;
    std::cout << "  RooKeysPdf: non-parametric mode, skipping fitTo."
              << std::endl;
  } else {
    // --- ResoConv: fit the resolution model directly (G1+CB2 [+Eds]) ---
    fittedPdf = resoConv;
    result = resoConv->fitTo(
        *dataReduced, RooFit::Extended(false), RooFit::SumW2Error(true),
        RooFit::Minimizer("Minuit2", "Migrad"), RooFit::Strategy(1),
        RooFit::Offset(true), RooFit::Optimize(true), RooFit::NumCPU(8),
        RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(1),
        RooFit::Save(true));
  }

  // fittedPdf = &resoModel;
  // result = resoModel.fitTo(
  //     *fitData, RooFit::Extended(false), RooFit::SumW2Error(true),
  //     RooFit::Minimizer("Minuit2", "Migrad"), RooFit::Strategy(1),
  //     RooFit::Offset(true), RooFit::Optimize(true), RooFit::NumCPU(8),
  //     RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(0),
  //     RooFit::Save(true));

  // fittedPdf = cbConv;
  // result = cbConv->fitTo(
  //     *dataReduced, RooFit::Extended(false), RooFit::SumW2Error(true),
  //     RooFit::Minimizer("Minuit2", "Migrad"), RooFit::Strategy(1),
  //     RooFit::Offset(true), RooFit::Optimize(true), RooFit::NumCPU(8),
  //     RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(0),
  //     RooFit::Save(true));

  if (result) {
    result->Print("v");
    result->correlationMatrix().Print();
  } else {
    std::cout << "  No fit result (non-parametric mode)." << std::endl;
  }

  // ========================================================================
  // 7. Plot
  // ========================================================================
  // binWidth and nBins already defined in section 5a

  TCanvas c("c", "Signal ctau fit", 900, 900);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  bool isHybridFit = useHybridConv;
  bool isNonParametric = (useHistPdf || useKeysPdf);

  // For non-parametric modes, pin ctau's default binning so that both
  // data->plotOn and pdf->plotOn share identical bin edges — required for
  // frame->pullHist("data","total") to match points bin-by-bin.
  if (isNonParametric) ctau.setBins(nBins);

  RooPlot* frame = ctau.frame(RooFit::Title("Signal c#tau fit (SPlot)"),
                              RooFit::Bins(nBins));

  // Data
  fitData->plotOn(frame, RooFit::Name("data"), RooFit::MarkerSize(0.5),
                  RooFit::DataError(RooAbsData::SumW2));

  // Plot the fitted PDF (whatever it was)
  if (fittedPdf) {
    // RooHistPdf: plain plotOn works — its normalization integral is
    // well-behaved. RooKeysPdf: the KDE can return NaN during the normalization
    // integral when ctau=0 (old default) was outside the data support; with the
    // default value now set to the range midpoint the integral is fine, so
    // plain plotOn works for both.  Keep the branch in case future debugging is
    // needed.
    fittedPdf->plotOn(frame, RooFit::Name("total"), RooFit::LineColor(kBlue),
                      RooFit::LineWidth(2));

    if (isNonParametric) {
      // Non-parametric modes have no sub-components to plot.
      // The single curve already drawn is the complete model.
      std::cout << "\nNon-parametric mode (" << mainModeStr
                << "): no sub-components to plot." << std::endl;
    } else if (isHybridFit) {
      // ============================================================
      // Hybrid-specific component plotting
      // ============================================================
      std::cout << "\nPlotting hybrid model components..." << std::endl;
      RooArgSet* allComps = fittedPdf->getComponents();

      if (hybridAsymDecay) {
        // ============================================================
        // hybridAsymDecay: 3-component model
        //   f_pr · R  +  (1-f_pr)·f_R · (R⊗ExpRight)
        //             +  (1-f_pr)·(1-f_R) · (R⊗ExpLeft)
        // Level 1 (solid, width 4): the 3 model components
        // Level 2 (dashed, width 2): sub-components of R⊗ExpRight  [red shades]
        //         (dash-dot, width 2): sub-components of R⊗ExpLeft  [blue
        //         shades] (dotted,   width 2): resolution sub-comps of R [green
        //         shades]
        // ============================================================
        std::cout << "  hybridAsymDecay plotting: 3 main + sub-components"
                  << std::endl;

        // --- Level 1: 3 main model components ---
        if (allComps->find("hybridReso"))
          fittedPdf->plotOn(frame, RooFit::Components("hybridReso"),
                            RooFit::LineColor(kGreen + 2),
                            RooFit::LineStyle(kSolid), RooFit::LineWidth(4),
                            RooFit::Name("asym_prompt"));
        if (allComps->find("asymConvRight"))
          fittedPdf->plotOn(frame, RooFit::Components("asymConvRight"),
                            RooFit::LineColor(kRed), RooFit::LineStyle(kSolid),
                            RooFit::LineWidth(4), RooFit::Name("asym_right"));
        if (hybridAsymExpLeft && allComps->find("asymConvLeft"))
          fittedPdf->plotOn(frame, RooFit::Components("asymConvLeft"),
                            RooFit::LineColor(kAzure + 1),
                            RooFit::LineStyle(kSolid), RooFit::LineWidth(4),
                            RooFit::Name("asym_left"));

        // --- Level 2: sub-components of R⊗ExpRight (dashed, red family) ---
        if (allComps->find("asymDecR_G1"))
          fittedPdf->plotOn(frame, RooFit::Components("asymDecR_G1"),
                            RooFit::LineColor(kRed - 7),
                            RooFit::LineStyle(kDashed), RooFit::LineWidth(2),
                            RooFit::Name("asym_R_G1"));
        if (allComps->find("asymFFT_cb2_R"))
          fittedPdf->plotOn(frame, RooFit::Components("asymFFT_cb2_R"),
                            RooFit::LineColor(kOrange + 1),
                            RooFit::LineStyle(kDashed), RooFit::LineWidth(2),
                            RooFit::Name("asym_R_CB2"));
        if (allComps->find("asymFFT_eds_R"))
          fittedPdf->plotOn(frame, RooFit::Components("asymFFT_eds_R"),
                            RooFit::LineColor(kPink + 1),
                            RooFit::LineStyle(kDashed), RooFit::LineWidth(2),
                            RooFit::Name("asym_R_Eds"));
        if (allComps->find("asymDecR_G4"))
          fittedPdf->plotOn(frame, RooFit::Components("asymDecR_G4"),
                            RooFit::LineColor(kMagenta + 1),
                            RooFit::LineStyle(kDashed), RooFit::LineWidth(2),
                            RooFit::Name("asym_R_G4"));

        // --- Level 2: sub-components of R⊗ExpLeft (dash-dot, blue family) ---
        if (hybridAsymExpLeft) {
          if (allComps->find("asymDecL_G1"))
            fittedPdf->plotOn(frame, RooFit::Components("asymDecL_G1"),
                              RooFit::LineColor(kAzure + 7),
                              RooFit::LineStyle(kDashDotted),
                              RooFit::LineWidth(2), RooFit::Name("asym_L_G1"));
          if (allComps->find("asymFFT_cb2_L"))
            fittedPdf->plotOn(frame, RooFit::Components("asymFFT_cb2_L"),
                              RooFit::LineColor(kCyan + 2),
                              RooFit::LineStyle(kDashDotted),
                              RooFit::LineWidth(2), RooFit::Name("asym_L_CB2"));
          if (allComps->find("asymFFT_eds_L"))
            fittedPdf->plotOn(frame, RooFit::Components("asymFFT_eds_L"),
                              RooFit::LineColor(kTeal + 2),
                              RooFit::LineStyle(kDashDotted),
                              RooFit::LineWidth(2), RooFit::Name("asym_L_Eds"));
          if (allComps->find("asymDecL_G4"))
            fittedPdf->plotOn(frame, RooFit::Components("asymDecL_G4"),
                              RooFit::LineColor(kBlue - 7),
                              RooFit::LineStyle(kDashDotted),
                              RooFit::LineWidth(2), RooFit::Name("asym_L_G4"));
        }

        // --- Level 2: resolution sub-components of R = hybridReso (dotted,
        // green family) ---
        if (allComps->find("reso_gauss1"))
          fittedPdf->plotOn(frame, RooFit::Components("reso_gauss1"),
                            RooFit::LineColor(kGreen - 6),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_g1"));
        if (allComps->find("cb2"))
          fittedPdf->plotOn(frame, RooFit::Components("cb2"),
                            RooFit::LineColor(kGreen + 3),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_cb2"));
        if (allComps->find("expDoubleSided"))
          fittedPdf->plotOn(frame, RooFit::Components("expDoubleSided"),
                            RooFit::LineColor(kSpring - 5),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_eds"));
        if (allComps->find("reso_gauss4"))
          fittedPdf->plotOn(frame, RooFit::Components("reso_gauss4"),
                            RooFit::LineColor(kGray + 2),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_g4"));

      } else {
        // ============================================================
        // hybridWOEds / hybridWEds: NP conv + prompt
        // ============================================================

        // --- Main components (solid, bold) ---
        if (allComps->find("hNPconv"))
          fittedPdf->plotOn(frame, RooFit::Components("hNPconv"),
                            RooFit::LineColor(kRed), RooFit::LineStyle(kSolid),
                            RooFit::LineWidth(3), RooFit::Name("np_total"));
        if (allComps->find("hybridReso"))
          fittedPdf->plotOn(frame, RooFit::Components("hybridReso"),
                            RooFit::LineColor(kGreen + 2),
                            RooFit::LineStyle(kSolid), RooFit::LineWidth(3),
                            RooFit::Name("p_total"));

        // --- NP sub-components (dashed) ---
        if (allComps->find("hdec1"))
          fittedPdf->plotOn(frame, RooFit::Components("hdec1"),
                            RooFit::LineColor(kRed - 7),
                            RooFit::LineStyle(kDashed), RooFit::LineWidth(2),
                            RooFit::Name("np_dec1"));
        if (allComps->find("hfft_cb2"))
          fittedPdf->plotOn(frame, RooFit::Components("hfft_cb2"),
                            RooFit::LineColor(kOrange + 1),
                            RooFit::LineStyle(kDashed), RooFit::LineWidth(2),
                            RooFit::Name("np_fft_cb2"));
        if (allComps->find("hfft_eds"))
          fittedPdf->plotOn(frame, RooFit::Components("hfft_eds"),
                            RooFit::LineColor(kMagenta),
                            RooFit::LineStyle(kDashed), RooFit::LineWidth(2),
                            RooFit::Name("np_fft_eds"));
        if (allComps->find("hdec4"))
          fittedPdf->plotOn(frame, RooFit::Components("hdec4"),
                            RooFit::LineColor(kViolet),
                            RooFit::LineStyle(kDashed), RooFit::LineWidth(2),
                            RooFit::Name("np_dec4"));

        // --- Prompt sub-components (dotted) ---
        if (allComps->find("reso_gauss1"))
          fittedPdf->plotOn(frame, RooFit::Components("reso_gauss1"),
                            RooFit::LineColor(kGreen - 6),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_g1"));
        if (allComps->find("cb2"))
          fittedPdf->plotOn(frame, RooFit::Components("cb2"),
                            RooFit::LineColor(kCyan + 2),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_cb2"));
        if (allComps->find("expDoubleSided"))
          fittedPdf->plotOn(frame, RooFit::Components("expDoubleSided"),
                            RooFit::LineColor(kPink + 1),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_eds"));
        if (allComps->find("reso_gauss4"))
          fittedPdf->plotOn(frame, RooFit::Components("reso_gauss4"),
                            RooFit::LineColor(kGray + 2),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_g4"));
        if (allComps->find("reso_gauss5"))
          fittedPdf->plotOn(frame, RooFit::Components("reso_gauss5"),
                            RooFit::LineColor(kOrange + 2),
                            RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
                            RooFit::Name("p_g5"));
      }

      delete allComps;

    } else {
      // Try to plot sub-components (only direct components, not nested ones)
      TString pdfName = fittedPdf->GetName();
      std::cout << "\nPlotting fitted PDF: " << pdfName << std::endl;

      // Define specific component names to plot based on model structure
      std::vector<TString> componentNames;

      // For RooAddPdf, we can get the direct component list
      if (fittedPdf->InheritsFrom(RooAddPdf::Class())) {
        RooAddPdf* addPdf = dynamic_cast<RooAddPdf*>(fittedPdf);
        if (addPdf) {
          const RooArgList& pdfList = addPdf->pdfList();
          std::cout << "Direct components in RooAddPdf: " << pdfList.getSize()
                    << std::endl;

          for (int i = 0; i < pdfList.getSize(); i++) {
            RooAbsArg* arg = pdfList.at(i);
            if (arg) {
              componentNames.push_back(arg->GetName());
            }
          }
        }
      }

      // Define colors for different component types
      std::vector<int> colors = {
          kRed,    kGreen + 2, kMagenta,  kOrange + 1, kCyan + 2,
          kViolet, kGray + 2,  kBlue + 2, kYellow + 2, kPink + 1};

      // Plot each direct component
      for (size_t i = 0; i < componentNames.size(); i++) {
        TString compName = componentNames[i];
        int color = colors[i % colors.size()];

        std::cout << "  Attempting to plot component: " << compName
                  << std::endl;

        // Create component label string
        TString compLabel = Form("%s_comp", compName.Data());

        // Try to plot this component
        fittedPdf->plotOn(frame, RooFit::Components(compName.Data()),
                          RooFit::LineColor(color), RooFit::LineStyle(kDashed),
                          RooFit::LineWidth(2), RooFit::Name(compLabel.Data()));

        // Verify it was plotted successfully
        if (frame->findObject(compLabel.Data())) {
          std::cout << "    Successfully plotted: " << compName << std::endl;
        } else {
          std::cout << "    Warning: Failed to plot: " << compName << std::endl;
        }
      }

      // Plot sub-components of resoModel (auto-detected)
      // Only in HybridConv mode, where resoModel is nested inside sngModel.
      // In Reso mode, resoConv shares the same components — already plotted
      // above.
      if (useHybridConv && resoModel.InheritsFrom(RooAddPdf::Class())) {
        const RooArgList& resoComps =
            dynamic_cast<RooAddPdf&>(resoModel).pdfList();
        int subColors[] = {kRed - 7,    kGreen - 6, kMagenta - 6, kCyan - 6,
                           kOrange - 3, kBlue - 7,  kYellow - 6,  kPink - 6};
        int nSubColors = sizeof(subColors) / sizeof(subColors[0]);
        std::cout << "  resoModel sub-components: " << resoComps.getSize()
                  << std::endl;
        for (int i = 0; i < resoComps.getSize(); i++) {
          RooAbsArg* sub = resoComps.at(i);
          if (!sub) continue;
          TString subName = sub->GetName();
          TString label = subName + "_sub";
          resoModel.plotOn(
              frame, RooFit::Components(subName.Data()),
              RooFit::Normalization(1.0 - f_B.getVal(), RooAbsReal::Relative),
              RooFit::LineColor(subColors[i % nSubColors]),
              RooFit::LineStyle(kDotted), RooFit::LineWidth(2),
              RooFit::Name(label.Data()));
          if (frame->findObject(label.Data())) {
            std::cout << "    Sub-component plotted: " << subName << std::endl;
          }
        }
      }
    }  // end else (non-hybrid plotting)
  }

  frame->SetMinimum(1.0);
  frame->SetMaximum(frame->GetMaximum() * 10.0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->Draw();

  // Legend (dynamic - only add entries that exist)
  TLegend leg(0.55, 0.40, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.023);
  leg.AddEntry(frame->findObject("data"), dataLabel, "ep");
  if (frame->findObject("total")) {
    leg.AddEntry(frame->findObject("total"), "Total signal", "l");
  }

  if (isNonParametric) {
    // Non-parametric: only the total curve entry (already added above).
    TString npLabel;
    if (useHistPdf)
      npLabel = Form("RooHistPdf (interp. %s)", subModeStr.c_str());
    else
      npLabel = Form("RooKeysPdf KDE (%s, #rho=%.2f)", keysMirrorStr.c_str(),
                     keysRhoFinal);
    if (frame->findObject("total"))
      leg.AddEntry(frame->findObject("total"), npLabel.Data(), "l");
  } else if (isHybridFit) {
    if (hybridAsymDecay) {
      // ---- hybridAsymDecay legend ----
      // Level 1: 3 main model components
      if (frame->findObject("asym_prompt"))
        leg.AddEntry(frame->findObject("asym_prompt"),
                     "Prompt: f_{pr} #times R", "l");
      if (frame->findObject("asym_right"))
        leg.AddEntry(frame->findObject("asym_right"),
                     "NP right: f_{R} #times R#otimesExp_{R}(#tau_{D})", "l");
      if (hybridAsymExpLeft && frame->findObject("asym_left"))
        leg.AddEntry(
            frame->findObject("asym_left"),
            "NP left: (1-f_{R}) #times R#otimesExp_{L}(#tau_{#lambda})", "l");
      // Level 2: R⊗ExpRight sub-components
      if (frame->findObject("asym_R_G1"))
        leg.AddEntry(frame->findObject("asym_R_G1"), "  G1 #otimes Exp_{R}",
                     "l");
      if (frame->findObject("asym_R_CB2"))
        leg.AddEntry(frame->findObject("asym_R_CB2"),
                     "  CB2 #otimes Exp_{R} (FFT)", "l");
      if (frame->findObject("asym_R_Eds"))
        leg.AddEntry(frame->findObject("asym_R_Eds"),
                     "  E_{ds} #otimes Exp_{R} (FFT)", "l");
      if (frame->findObject("asym_R_G4"))
        leg.AddEntry(frame->findObject("asym_R_G4"), "  G4 #otimes Exp_{R}",
                     "l");
      // Level 2: R⊗ExpLeft sub-components (only when active)
      if (hybridAsymExpLeft) {
        if (frame->findObject("asym_L_G1"))
          leg.AddEntry(frame->findObject("asym_L_G1"), "  G1 #otimes Exp_{L}",
                       "l");
        if (frame->findObject("asym_L_CB2"))
          leg.AddEntry(frame->findObject("asym_L_CB2"),
                       "  CB2 #otimes Exp_{L} (FFT)", "l");
        if (frame->findObject("asym_L_Eds"))
          leg.AddEntry(frame->findObject("asym_L_Eds"),
                       "  E_{ds} #otimes Exp_{L} (FFT)", "l");
        if (frame->findObject("asym_L_G4"))
          leg.AddEntry(frame->findObject("asym_L_G4"), "  G4 #otimes Exp_{L}",
                       "l");
      }
      // Level 2: resolution sub-components of R
      if (frame->findObject("p_g1"))
        leg.AddEntry(frame->findObject("p_g1"), "  R: G1", "l");
      if (frame->findObject("p_cb2"))
        leg.AddEntry(frame->findObject("p_cb2"), "  R: CB2", "l");
      if (frame->findObject("p_eds"))
        leg.AddEntry(frame->findObject("p_eds"), "  R: E_{ds}", "l");
      if (frame->findObject("p_g4"))
        leg.AddEntry(frame->findObject("p_g4"), "  R: G4 (bad assoc.)", "l");

    } else {
      // ---- hybridWOEds / hybridWEds legend ----
      if (frame->findObject("np_total"))
        leg.AddEntry(frame->findObject("np_total"), "Non-Prompt", "l");
      if (frame->findObject("p_total"))
        leg.AddEntry(frame->findObject("p_total"), "Prompt", "l");
      // NP sub-components
      if (frame->findObject("np_dec1"))
        leg.AddEntry(frame->findObject("np_dec1"),
                     "NP: E_{ss} #otimes G1 (Decay)", "l");
      if (frame->findObject("np_fft_cb2"))
        leg.AddEntry(frame->findObject("np_fft_cb2"),
                     "NP: E_{ss} #otimes CB2 (FFT)", "l");
      if (frame->findObject("np_fft_eds"))
        leg.AddEntry(frame->findObject("np_fft_eds"),
                     "NP: E_{ss} #otimes E_{ds} (FFT)", "l");
      if (frame->findObject("np_dec4"))
        leg.AddEntry(frame->findObject("np_dec4"),
                     "NP: E_{ss} #otimes G4 (Decay)", "l");
      // Prompt sub-components
      if (frame->findObject("p_g1"))
        leg.AddEntry(frame->findObject("p_g1"), "P: G1", "l");
      if (frame->findObject("p_cb2"))
        leg.AddEntry(frame->findObject("p_cb2"), "P: CB2", "l");
      if (frame->findObject("p_eds"))
        leg.AddEntry(frame->findObject("p_eds"), "P: E_{ds}", "l");
      if (frame->findObject("p_g4"))
        leg.AddEntry(frame->findObject("p_g4"), "G4: Bad Assoc.", "l");
    }

  } else {
    // Dynamically find and add all plotted components
    if (fittedPdf && fittedPdf->InheritsFrom(RooAddPdf::Class())) {
      RooAddPdf* addPdf = dynamic_cast<RooAddPdf*>(fittedPdf);
      if (addPdf) {
        const RooArgList& pdfList = addPdf->pdfList();

        // Only look for direct components that were plotted
        for (int i = 0; i < pdfList.getSize(); i++) {
          RooAbsArg* arg = pdfList.at(i);
          if (!arg) continue;

          TString compName = arg->GetName();
          TString compLabel = Form("%s_comp", compName.Data());

          // Check if this component was successfully plotted
          if (frame->findObject(compLabel.Data())) {
            // Clean up the name for display
            TString displayName = arg->GetTitle();
            if (displayName == "") displayName = compName;
            leg.AddEntry(frame->findObject(compLabel.Data()),
                         displayName.Data(), "l");
          }
        }
      }
    }

    // Legend entries for resoModel sub-components (auto-detected)
    // Only in HybridConv mode — in Reso mode, direct components are already
    // added above.
    if (fittedPdf && useHybridConv &&
        resoModel.InheritsFrom(RooAddPdf::Class())) {
      const RooArgList& resoComps =
          dynamic_cast<RooAddPdf&>(resoModel).pdfList();
      for (int i = 0; i < resoComps.getSize(); i++) {
        RooAbsArg* sub = resoComps.at(i);
        if (!sub) continue;
        TString label = TString(sub->GetName()) + "_sub";
        if (frame->findObject(label.Data())) {
          TString displayName = sub->GetTitle();
          if (displayName == "") displayName = sub->GetName();
          displayName += " (reso)";
          leg.AddEntry(frame->findObject(label.Data()), displayName.Data(),
                       "l");
        }
      }
    }
  }  // end else (non-hybrid legend)
  leg.Draw();

  // Parameter box (dynamic - show only fitted parameters)
  TLatex tex;
  tex.SetNDC();
  tex.SetTextSize(0.025);

  double yPos = 0.88;
  double yStep = 0.032;

  // Get all parameters (floating + fixed) from the fitted PDF
  if (result) {
    const RooArgList& floatPars = result->floatParsFinal();
    const RooArgList& constPars = result->constPars();

    // Display floating parameters (with errors)
    for (int i = 0; i < floatPars.getSize(); i++) {
      RooRealVar* var = dynamic_cast<RooRealVar*>(floatPars.at(i));
      if (var) {
        TString parName = var->GetTitle();
        if (parName == "") parName = var->GetName();
        tex.DrawLatex(0.15, yPos,
                      Form("%s = %.4g #pm %.4g", parName.Data(), var->getVal(),
                           var->getError()));
        yPos -= yStep;
      }
    }

    // Display fixed parameters (no errors)
    for (int i = 0; i < constPars.getSize(); i++) {
      RooRealVar* var = dynamic_cast<RooRealVar*>(constPars.at(i));
      if (var) {
        TString parName = var->GetTitle();
        if (parName == "") parName = var->GetName();
        tex.DrawLatex(0.15, yPos,
                      Form("%s = %.4g (fixed)", parName.Data(), var->getVal()));
        yPos -= yStep;
      }
    }

    // Fit info
    yPos -= yStep * 0.5;  // Extra space
    tex.DrawLatex(0.15, yPos, Form("Status: %d (0 = OK)", result->status()));
    yPos -= yStep;
  } else {
    tex.DrawLatex(0.15, yPos, "Non-parametric PDF (no fit parameters)");
    yPos -= yStep * 1.5;
  }
  // nFloatPars = 0 for non-parametric modes: chi2/ndf is computed at the
  // plotting resolution (binWidth grid) comparing data vs template — nBinsHist
  // is irrelevant since the dataset is unbinned and the template can exceed the
  // plotting resolution.
  int nFloatPars = result ? result->floatParsFinal().getSize() : 0;
  int ndfVal = frame->GetNbinsX() - nFloatPars;
  double chi2val = frame->chiSquare("total", "data", nFloatPars) * ndfVal;
  tex.DrawLatex(0.15, yPos,
                Form("#chi^{2}/ndf = %.2f / %d = %.2f", chi2val, ndfVal,
                     chi2val / ndfVal));

  // Pull histogram
  c.cd();
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();

  RooHist* pullHist = frame->pullHist("data", "total");
  RooPlot* pullFrame = ctau.frame(RooFit::Title(""));
  if (!pullHist) {
    std::cerr << "WARNING: pullHist is null — data/total curve binning mismatch"
              << std::endl;
  } else {
    pullFrame->addPlotable(pullHist, "P");
  }
  pullFrame->SetYTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.12);
  pullFrame->GetYaxis()->SetTitleOffset(0.35);
  pullFrame->GetYaxis()->SetLabelSize(0.1);
  pullFrame->GetYaxis()->SetNdivisions(505);
  pullFrame->GetXaxis()->SetTitleSize(0.12);
  pullFrame->GetXaxis()->SetTitleOffset(1.0);
  pullFrame->GetXaxis()->SetLabelSize(0.1);
  pullFrame->Draw();

  TLine pullLine(ctauMin, 0, ctauMax, 0);
  pullLine.SetLineColor(kRed);
  pullLine.SetLineStyle(kDashed);
  pullLine.Draw();

  // Save output files into the specified folder
  std::string outDir = std::string(gSystem->DirName(__FILE__)) +
                       "/../../plots/splot/" + inputFolder;
  gSystem->mkdir(outDir.c_str(), kTRUE);
  c.SaveAs(Form("%s/RedoSplotSignalFit.pdf", outDir.c_str()));

  // Save workspace with data, PDFs, and parameters for interactive exploration
  RooWorkspace ws("ws", "Signal fit workspace");
  ws.import(*dataReduced, RooFit::Rename("dataw_Jpsi"));
  if (fittedPdf) {
    ws.import(*fittedPdf, RooFit::RecycleConflictNodes());
  }
  ws.var("ctau")->setRange(ctauMin, ctauMax);

  TFile outFile(Form("%s/RedoSplotSignalFit.root", outDir.c_str()), "RECREATE");
  ws.Write();
  c.Write();
  outFile.Close();

  std::cout << "\nWorkspace saved to RedoSplotSignalFit.root" << std::endl;
  std::cout << "To explore interactively:" << std::endl;
  std::cout << "  TFile f(\"" << outDir << "/RedoSplotSignalFit.root\");"
            << std::endl;
  std::cout << "  RooWorkspace* ws = (RooWorkspace*)f.Get(\"ws\");"
            << std::endl;
  if (fittedPdf) {
    std::cout << "  ws->pdf(\"" << fittedPdf->GetName() << "\")->plotOn(frame);"
              << std::endl;
  }
  std::cout << "  auto* frame = ws->var(\"ctau\")->frame();" << std::endl;
  std::cout << "  ws->data(\"dataw_Jpsi\")->plotOn(frame);" << std::endl;
  std::cout << "  frame->Draw();" << std::endl;

  // ========================================================================
  // Save fit results as JSON (same structure as input config, with errors)
  // ========================================================================
  {
    // Helper: update a param node with post-fit value, range, fixed flag,
    // and — only for free parameters — an "error" field.
    auto updateParam = [](nlohmann::json& node, RooRealVar& var) {
      node["value"] = var.getVal();
      node["min"] = var.getMin();
      node["max"] = var.getMax();
      node["fixed"] = static_cast<bool>(var.isConstant());
      if (!var.isConstant())
        node["error"] = var.getError();
      else
        node.erase("error");
    };

    nlohmann::json outCfg = cfg;  // copy preserves all existing fields/comments

    auto& jG1CB2 = outCfg["reso_G1_CB2"];
    auto& jEds = outCfg["reso_Eds"];
    auto& jG4 = outCfg["reso_G4"];
    auto& jNP = outCfg["signal_NP"];

    updateParam(jG1CB2["reso_mean1"], reso_mean1);
    if (useCb2OwnMean) updateParam(jG1CB2["cb2_mean"], cb2_mean);
    updateParam(jG1CB2["sigma_G1"], sigma1);
    updateParam(jG1CB2["sigma_CB2"], cb2_sigma);
    updateParam(jG1CB2["cbAlpha1"], cbAlpha1);
    updateParam(jG1CB2["cbN1"], cbN1);
    if (!useCbRatio) {
      updateParam(jG1CB2["cbAlpha2"], cbAlpha2);
      updateParam(jG1CB2["cbN2"], cbN2);
    }
    updateParam(jG1CB2["frac1"], frac1);

    updateParam(jEds["lambda_exp3"], lambda_exp3);
    updateParam(jEds["mean_exp3"], mean_exp3);
    updateParam(jEds["frac2"], frac2);

    updateParam(jG4["reso_mean4"], reso_mean4);
    updateParam(jG4["sigma_G4"], sigma4);

    updateParam(jNP["tauB"], tauB);
    updateParam(jNP["f_B"], f_B);
    updateParam(jNP["a1"], a1);
    updateParam(jNP["a2"], a2);

    if (hybridAsymDecay && tauD_var) {
      auto& jA = outCfg["asymDecay"];
      updateParam(jA["tauD"], *tauD_var);
      updateParam(jA["f_pr"], *f_pr_var);
      if (hybridAsymExpLeft && tauLambda_var) {
        updateParam(jA["tauLambda"], *tauLambda_var);
        updateParam(jA["f_R"], *f_R_var);
      }
    }

    // Fit quality summary (extra block appended to the output JSON)
    outCfg["fit_result"]["status"] = result ? result->status() : -1;
    outCfg["fit_result"]["edm"] = result ? result->edm() : -1.0;
    outCfg["fit_result"]["chi2"] = chi2val;
    outCfg["fit_result"]["ndf"] = ndfVal;
    outCfg["fit_result"]["chi2_ndf"] = chi2val / ndfVal;

    std::string cfgBase = gSystem->BaseName(configFile);
    std::string outJsonPath = outDir + "/out_" + cfgBase;
    std::ofstream outJsonStream(outJsonPath);
    outJsonStream << outCfg.dump(2) << "\n";
    std::cout << "Fit result JSON saved to: " << outJsonPath << std::endl;
  }

  std::cout << "\n=== Fit Summary ===" << std::endl;
  std::cout << "Fitted PDF: " << (fittedPdf ? fittedPdf->GetName() : "Unknown")
            << std::endl;
  if (result) {
    std::cout << "Status: " << result->status() << " (0 = converged)"
              << std::endl;
    std::cout << "\nFitted Parameters:" << std::endl;
    const RooArgList& finalPars = result->floatParsFinal();
    for (int i = 0; i < finalPars.getSize(); i++) {
      RooRealVar* var = dynamic_cast<RooRealVar*>(finalPars.at(i));
      if (var) {
        std::cout << "  " << var->GetName() << " = " << var->getVal() << " +/- "
                  << var->getError() << std::endl;
      }
    }
  } else {
    std::cout << "Non-parametric mode: no fit result." << std::endl;
  }
  std::cout << "===================" << std::endl;

  // Cleanup
  f->Close();
  delete sngModel;

  if (useHybridConv) {
    delete hybridReso;
    delete hNPconv;
    delete hdec1;
    delete hdec4;
    delete hgm1;
    delete hgm4;
    delete hfft_cb2;
    delete hfft_eds;
    delete hexpoNP_cb;
    delete hexpoNP_eds;
    delete hcb2_conv;
    delete heds_conv;
    // Asymmetric decay cleanup
    delete tauD_var;
    delete tauLambda_var;
    delete f_pr_var;
    delete f_R_var;
    delete asymDecR_G1;
    delete asymDecL_G1;
    delete asymDecR_G4;
    delete asymDecL_G4;
    delete asymExpR_cb;
    delete asymExpL_cb;
    delete asymExpR_eds;
    delete asymExpL_eds;
    delete asymEds_conv;
    delete asymFFT_cb2_R;
    delete asymFFT_cb2_L;
    delete asymFFT_eds_R;
    delete asymFFT_eds_L;
    delete asymConvRight;
    delete asymConvLeft;
    delete asymCB2_conv;
  }
  delete dataBinned;
  // histPdf and histPdfData are owned by the workspace after ws.import — do not
  // delete
  delete keysPdf;
}

/// @}  // AddonRedoSplotSignalFit
