#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "TLatex.h"

void sethistostyle(TH1D *h1, int color = 1, int style = 20)
{
    h1->SetLineColor(color);
    h1->SetMarkerColor(color);
    h1->SetMarkerStyle(style);
    h1->SetMarkerSize(0.8);
    h1->SetLineWidth(2);
    h1->SetStats(0);
    h1->SetTitle("");
    h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h1->GetYaxis()->SetTitle("|#it{v}_{2}|");
    h1->GetYaxis()->SetTitleOffset(1.6);
}

/// Compares two files containing histograms and prints the differences.
/// The histograms to compare can be specified by the `histograms` parameter.
/// If no histograms are specified, it defaults to comparing "v2_vs_mass".
///
/// @param fileName1 The first file to compare.
/// @param fileName2 The second file to compare.
/// @param histograms A comma-separated list of histogram names to compare.
void compare_v2_mc( //
    const char *filename1 = "/Users/spolitan/alice/alice-correlations/oo/hh_all_mc/out_Wsubtraction_deltaeta1.2.root",
    const char *filename2 = "/Users/spolitan/alice/alice-correlations/oo/hh_all_mc/out_WOsubtraction_deltaeta1.2.root",
    const char *outputFileName = "/Users/spolitan/alice/alice-correlations/oo/hh_all_mc/comparison_v2scaled_hh_mc_deltaeta1.2_withsub", // output file name without .root
    const bool isPhi = false)
{

    std::vector<std::string> mult_bins = {"25-50", "50-150", "150-300"};
    std::vector<std::string> leg_labels = {
        "with LM subtraction",
        "without LM subtraction"};

    auto outFile = TFile::Open(Form("%s.root", outputFileName), "RECREATE");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.033);

    std::vector<int> colors = {
        TColor::GetColorTransparent(kRed + 1, 0.55), TColor::GetColorTransparent(kAzure + 4, 0.55), TColor::GetColorTransparent(kGreen + 2, 0.55), TColor::GetColorTransparent(kOrange + 1, 0.55)};

    // th1 vs mult
    std::vector<double> mult_values1 = {0., 10.0, 20.0};
    std::vector<double> mult_values2 = {25., 50.0, 150.0, 300.0};
    std::vector<double> pt_values = {0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0};
    int nptbins = 6; // Number of pT bins

    TCanvas *c = new TCanvas("c", "Comparison Canvas", 1800, 600);
    c->Divide(3, 1);
    for (int ihist = 1; ihist < 4; ++ihist)
    {
        const char *histogram1 = Form("v2scaled_vs_mass_mult%d", ihist);

        TH1D *hist1 = (TH1D *)TFile::Open(filename1)->Get(histogram1);
        TH1D *hist2 = (TH1D *)TFile::Open(filename2)->Get(histogram1);
        // TH1D *hist3 = (TH1D *)TFile::Open(filename3)->Get(histogram1);

        if (!hist1 || !hist2)
        {
            std::cerr << "Error: One of the histograms does not exist in the files." << std::endl;
            continue;
        }

        sethistostyle(hist1, colors[0], 20);
        sethistostyle(hist2, colors[1], 21);
        // sethistostyle(hist3, colors[2], 22);

        TLegend *legend = new TLegend(0.12, 0.6, 0.3, 0.8);
        legend->AddEntry(hist1, leg_labels[0].c_str(), "lp");
        // legend->AddEntry(fitFunc1, "[0] = " + TString::Format("%.4f #pm %.4f", fitFunc1->GetParameter(0), fitFunc1->GetParError(0)), "l");
        legend->AddEntry(hist2, leg_labels[1].c_str(), "lp");
        // legend->AddEntry(fitFunc2, "[0] = " + TString::Format("%.4f #pm %.4f", fitFunc2->GetParameter(0), fitFunc2->GetParError(0)), "l");
        //  legend->AddEntry(hist3, leg_labels[2].c_str(), "lp");
        legend->SetTextSize(0.04);
        legend->SetTextFont(42);
        legend->SetBorderSize(0);

        c->cd(ihist);
        hist1->GetYaxis()->SetRangeUser(-0.1, 0.3);
        hist1->GetYaxis()->SetDecimals(true);
        hist1->Draw("E");
        hist2->Draw("E same");

        TF1 *fitFunc1 = new TF1(Form("fpol0_1_%d", ihist), "pol0", 0.5, 8.0);
        fitFunc1->SetLineColor(colors[0]);
        fitFunc1->SetLineWidth(8);
        hist1->Fit(fitFunc1, "RQM0");
        fitFunc1->SetLineStyle(2);
        fitFunc1->Draw("same");

        TLatex *fitParam1 = new TLatex(0.15, 0.15, Form("[0] = %.4f #pm %.4f", fitFunc1->GetParameter(0), fitFunc1->GetParError(0)));
        fitParam1->SetNDC();
        fitParam1->SetTextSize(0.03);
        fitParam1->SetTextColor(colors[0]);
        fitParam1->Draw();

        if (ihist == 1)
        {
            legend->Draw();
            latex.DrawLatex(0.15, 0.85, "LHC25h3 MC (OO apass2)");
            latex.DrawLatex(0.55, 0.78, "|#it{#Delta#eta}| > 1.2");
        }
        latex.DrawLatex(0.55, 0.85, Form("Nch: [%s]", mult_bins[ihist - 1].c_str()));
    }
    c->Update();

    // Save the canvas to the output file
    outFile->cd();
    c->Write();
    c->SaveAs(Form("%s.pdf", outputFileName));
    // c2->Write();
    outFile->Close();
}
