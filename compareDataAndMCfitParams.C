
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstdlib>

#include <TROOT.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TColor.h>

using namespace std;

// ============================================================
// CONFIGURATION
// ============================================================

bool use_pO = true;
bool use_OO = true;

bool plotData = true;
bool plotMC   = false;

bool skipIntegratedBin = true;

const char* outputDir = "plots";

// ============================================================
// OBSERVABLE DEFINITION
// ============================================================

struct Observable {
    string name;
    string label;
    string group;
    int color;
    double Ymin;
    double Ymax;
    bool logY;
};

vector<Observable> GetObservables(bool isOO)
{
    vector<Observable> obs = {
        {"fit chi2 (last performed)", "Fir #chi^{2}", "chi2", kBlack, 0, 100, false},
        {"fit chi2 (last performed)", "Fir #chi^{2} zoomed", "chi2zoom", kBlack, 0, 10, false},

        {"mean_tauzRes (ns)", "Mean #tau_{z} resolution", "mean", kBlack, -0.00001, 0.00001, false},

        {"sigma0_tauzRes (ns)", "#sigma_{0}", "sigma", kRed+1, 0.0001, 1.0, true},
        {"sigma1_tauzRes (ns)", "#sigma_{1}", "sigma", kBlue+1, 0.0001, 1.0, true},
        {"sigma2_tauzRes (ns)", "#sigma_{2}", "sigma", kGreen+2, 0.0001, 1.0, true},

        {"fGaus0_tauzRes", "f_{0}", "fraction", kRed+1, 0.01, 1.1, false},
        {"fGaus1_tauzRes", "f_{1}", "fraction", kBlue+1, 0.01, 1.1, false}
    };

    if (isOO) {
        obs.push_back(
            {"sigma3_tauzRes (ns)", "#sigma_{3}", "sigma", kMagenta+1, 0.0001, 1.0, true}
        );

        obs.push_back(
            {"fGaus2_tauzRes", "f_{2}", "fraction", kGreen+2, 0.01, 1, true}
        );
    }

    return obs;
}

// ============================================================
// DATA STRUCTURE
// ============================================================

struct DataPoint {
    double x;
    double ex;
    double y;
    double ey;
};

struct Table {
    vector<double> ptMin;
    vector<double> ptMax;
    map<string, vector<DataPoint>> values;
};

// ============================================================
// PARSE VALUE AND UNCERTAINTY
// Example: 0.002 (7.9892e-08)
// ============================================================

bool ParseValue(string input, double& value, double& error)
{
    size_t pos = input.find('(');

    string valStr = input.substr(0, pos);
    string errStr = "0";

    if (pos != string::npos) {
        size_t end = input.find(')', pos);
        if (end != string::npos)
            errStr = input.substr(pos + 1, end - pos - 1);
    }

    try {
        size_t n1, n2;

        value = stod(valStr, &n1);
        error = stod(errStr, &n2);

        if (!isfinite(value) || !isfinite(error))
            return false;

        return true;
    }
    catch (...) {
        return false;
    }
}

// ============================================================
// READ TSV
// ============================================================

Table ReadTable(const char* filename, const vector<Observable>& observables)
{
    Table table;

    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "[ERROR] Cannot open " << filename << endl;
        return table;
    }

    string line;

    if (!getline(file, line))
        return table;

    vector<string> headers;
    stringstream headerStream(line);
    string field;

    while (getline(headerStream, field, '\t'))
        headers.push_back(field);

    map<string, int> columns;

    for (int i = 0; i < (int)headers.size(); ++i) {
        string h = headers[i];

        if (!h.empty() && h.back() == '\r')
            h.pop_back();

        columns[h] = i;
    }

    if (!columns.count("pT_bin (GeV/c)")) {
        cerr << "[ERROR] Missing pT_bin column in "
             << filename << endl;
        return table;
    }

    for (const auto& obs : observables) {
        if (!columns.count(obs.name)) {
            cerr << "[WARNING] Missing column "
                 << obs.name << " in " << filename << endl;
        }
    }

    while (getline(file, line)) {

        if (line.empty())
            continue;

        vector<string> fields;
        stringstream ss(line);

        while (getline(ss, field, '\t'))
            fields.push_back(field);

        if (fields.size() < headers.size())
            continue;

        double low, high;

        string bin = fields[columns["pT_bin (GeV/c)"]];

        size_t comma = bin.find(',');
        size_t left  = bin.find('[');
        size_t right = bin.find(']');

        if (comma == string::npos ||
            left == string::npos ||
            right == string::npos)
            continue;

        try {
            low = stod(bin.substr(left + 1, comma - left - 1));
            high = stod(bin.substr(comma + 1, right - comma - 1));
        }
        catch (...) {
            continue;
        }

        if (skipIntegratedBin && low == 0 && high == 20)
            continue;

        double x = 0.5 * (low + high);
        double ex = 0.5 * (high - low);

        table.ptMin.push_back(low);
        table.ptMax.push_back(high);

        for (const auto& obs : observables) {

            if (!columns.count(obs.name))
                continue;

            double value, error;

            if (!ParseValue(fields[columns[obs.name]], value, error))
                continue;

            table.values[obs.name].push_back(
                {x, ex, value, error}
            );
        }
    }

    return table;
}

// ============================================================
// DRAW ONE OBSERVABLE GROUP
// ============================================================

void DrawGroup(
    const string& group,
    const vector<Observable>& observables,
    const Table& data,
    const Table& mc,
    bool hasData,
    bool hasMC,
    const string& system)
{
    TCanvas* canvas = new TCanvas(
        Form("c_%s_%s", system.c_str(), group.c_str()),
        Form("%s - %s", system.c_str(), group.c_str()),
        900, 650
    );

    TMultiGraph* mg = new TMultiGraph();

    TLegend* legend = new TLegend(0.68, 0.68, 0.89, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    for (const auto& obs : observables) {

        if (obs.group != group)
            continue;

        // DATA: solid line
        if (hasData && data.values.count(obs.name)) {

            const auto& points = data.values.at(obs.name);

            if (!points.empty()) {

                TGraphErrors* gr = new TGraphErrors(points.size());

                for (int i = 0; i < (int)points.size(); ++i) {
                    gr->SetPoint(i, points[i].x, points[i].y);
                    gr->SetPointError(i, points[i].ex, points[i].ey);
                }

                gr->SetLineColor(obs.color);
                gr->SetMarkerColor(obs.color);
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.9);
                gr->SetLineWidth(2);
                gr->SetLineStyle(1);

                mg->Add(gr, "P");

                legend->AddEntry(
                    gr,
                    Form("%s - Data", obs.label.c_str()),
                    "p"
                );
            }
        }

        // MC: dashed line
        if (hasMC && mc.values.count(obs.name)) {

            const auto& points = mc.values.at(obs.name);

            if (!points.empty()) {

                TGraphErrors* gr = new TGraphErrors(points.size());

                for (int i = 0; i < (int)points.size(); ++i) {
                    gr->SetPoint(i, points[i].x, points[i].y);
                    gr->SetPointError(i, points[i].ex, points[i].ey);
                }

                gr->SetLineColor(obs.color);
                gr->SetMarkerColor(obs.color);
                gr->SetMarkerStyle(24);
                gr->SetMarkerSize(0.9);
                gr->SetLineWidth(2);
                gr->SetLineStyle(2);

                mg->Add(gr, "P");

                legend->AddEntry(
                    gr,
                    Form("%s - MC", obs.label.c_str()),
                    "p"
                );
            }
        }
    }

    mg->SetTitle("");
    mg->Draw("A");

    mg->SetTitle(Form("%s", system.c_str()));
    mg->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    mg->GetYaxis()->SetTitle("Fit parameter");

    mg->GetXaxis()->SetLimits(0, 12);

    mg->GetXaxis()->SetTitleSize(0.045);
    mg->GetYaxis()->SetTitleSize(0.045);

    mg->GetXaxis()->SetLabelSize(0.04);
    mg->GetYaxis()->SetLabelSize(0.04);

    legend->Draw();

    canvas->SetTicks(1, 1);
    canvas->SetLeftMargin(0.13);
    canvas->SetBottomMargin(0.12);

    for (const auto& obs : observables) {
        if (obs.group == group) {
            mg->GetYaxis()->SetRangeUser(obs.Ymin, obs.Ymax);
            if (obs.logY)
                canvas->SetLogy();
            break;
        }
    }

    canvas->SaveAs(
        Form("%s/%s_%s.pdf",
             outputDir, system.c_str(), group.c_str())
    );

    canvas->SaveAs(
        Form("%s/%s_%s.png",
             outputDir, system.c_str(), group.c_str())
    );
}

// ============================================================
// MAIN
// ============================================================

void compareDataAndMCfitParams(
    const char* caseName_pO_data,
    const char* caseName_pO_MC,
    const char* caseName_OO_data,
    const char* caseName_OO_MC)
{
    gStyle->SetOptStat(0);

    // ========================================================
    // CONSTRUCT INPUT FILENAMES FROM ARGUMENTS
    // ========================================================

    TString file_pO_data = Form(
        "output/fitSummary_pO_%s.tsv",
        caseName_pO_data
    );

    TString file_pO_MC = Form(
        "output/fitSummary_pO_MC_%s.tsv",
        caseName_pO_MC
    );

    TString file_OO_data = Form(
        "output/fitSummary_OO_%s.tsv",
        caseName_OO_data
    );

    TString file_OO_MC = Form(
        "output/fitSummary_OO_MC_%s.tsv",
        caseName_OO_MC
    );

    cout << "\n[INFO] Input files:" << endl;
    cout << "pO data: " << file_pO_data << endl;
    cout << "pO MC:   " << file_pO_MC << endl;
    cout << "OO data: " << file_OO_data << endl;
    cout << "OO MC:   " << file_OO_MC << endl;

    // ========================================================
    // VALIDATION
    // ========================================================

    if (!use_pO && !use_OO) {
        cout << "[ERROR] Select at least one collision system." << endl;
        return;
    }

    if (!plotData && !plotMC) {
        cout << "[ERROR] Select data and/or MC." << endl;
        return;
    }

    gSystem->mkdir(outputDir, true);

    vector<string> groups = {
        "chi2", "chi2zoom", "mean", "sigma", "fraction"
    };

    // ========================================================
    // pO
    // ========================================================

    if (use_pO) {

        auto observables = GetObservables(false);

        Table data, mc;

        if (plotData)
            data = ReadTable(file_pO_data.Data(), observables);

        if (plotMC)
            mc = ReadTable(file_pO_MC.Data(), observables);

        for (const auto& group : groups) {
            DrawGroup(
                group, observables, data, mc,
                plotData, plotMC, "pO"
            );
        }
    }

    // ========================================================
    // OO
    // ========================================================

    if (use_OO) {

        auto observables = GetObservables(true);

        Table data, mc;

        if (plotData)
            data = ReadTable(file_OO_data.Data(), observables);

        if (plotMC)
            mc = ReadTable(file_OO_MC.Data(), observables);

        for (const auto& group : groups) {
            DrawGroup(
                group, observables, data, mc,
                plotData, plotMC, "OO"
            );
        }
    }

    cout << "\n[INFO] Plots saved in " << outputDir << endl;
}
