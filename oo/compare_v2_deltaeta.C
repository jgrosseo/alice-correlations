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
    h1->GetYaxis()->SetTitle("#it{v}_{2}");
    h1->GetYaxis()->SetTitleOffset(1.6);
}

/// Compares two files containing histograms and prints the differences.
/// The histograms to compare can be specified by the `histograms` parameter.
/// If no histograms are specified, it defaults to comparing "v2_vs_mass".
///
/// @param fileName1 The first file to compare.
/// @param fileName2 The second file to compare.
/// @param histograms A comma-separated list of histogram names to compare.
void compare_v2_deltaeta(const char *fileName1 = "/Users/spolitan/alice/alice-correlations/oo/k0_hm_largeretacheck/v2ofK0sh_deltaeta1.3.root",
                         const char *fileName3 = "/Users/spolitan/alice/alice-correlations/oo/k0_hm_largeretacheck/v2ofK0sh_deltaeta1.5.root", // "/Users/spolitan/Downloads/stefano_phi_final/Invmassfit_nosub/v2ofphih_LMsubtraction.root",
                         const char *fileName2 = "/Users/spolitan/alice/alice-correlations/oo/k0_hm_largeretacheck/v2ofK0sh_deltaeta1.4.root",
                         const char *sclaing = "/Users/spolitan/alice/alice-correlations/oo/hh_lm_largeretacheck/",
                         const char *outputFileName = "/Users/spolitan/alice/alice-correlations/oo/k0_hm_largeretacheck/comparison_v2_k0s_deltaeta")
{

    std::vector<std::string> mult_bins = {"25-50", "50-100", "100-300"};
    std::vector<std::string> leg_labels = {"|#it{#Delta#eta}| > 1.2", "|#it{#Delta#eta}| > 1.4", "|#it{#Delta#eta}| > 1.5"};

    auto outFile = TFile::Open(Form("%s.root", outputFileName), "RECREATE");

    std::vector<int> colors = {
        TColor::GetColorTransparent(kRed + 1, 0.55), TColor::GetColorTransparent(kAzure + 4, 0.55), TColor::GetColorTransparent(kGreen + 2, 0.55), TColor::GetColorTransparent(kOrange + 1, 0.55)};

    // th1 vs mult
    std::vector<double> mult_values = {25., 50.0, 100.0, 300.0};
    std::vector<double> pt_values = {0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0};
    int nptbins = 6; // Number of pT bins
    std::vector<TH1D *> histograms1;
    std::vector<TH1D *> histograms2;
    std::vector<TH1D *> histograms3;

    // scaling
    TH1D *hscaling1 = nullptr;
    TH1D *hscaling2 = nullptr;
    TH1D *hscaling3 = nullptr;
    if (strlen(sclaing) != 0)
    {
        hscaling1 = (TH1D *)TFile::Open("/Users/spolitan/alice/alice-correlations/oo/hh_lm_largeretacheck/out_sub2hadron_deltaeta1.3.root")->Get("v2vsmult");
        hscaling2 = (TH1D *)TFile::Open("/Users/spolitan/alice/alice-correlations/oo/hh_lm_largeretacheck/out_sub2hadron_deltaeta1.4.root")->Get("v2vsmult");
        hscaling3 = (TH1D *)TFile::Open("/Users/spolitan/alice/alice-correlations/oo/hh_lm_largeretacheck/out_sub2hadron_deltaeta1.5.root")->Get("v2vsmult");
    }

    for (int i = 0; i < nptbins; ++i)
    {
        histograms1.push_back(new TH1D("h1", "Comparison of v2 vs mass", 3, mult_values.data()));
        histograms2.push_back(new TH1D("h2", "Comparison of v2 vs mass", 3, mult_values.data()));
        sethistostyle(histograms1.back(), kRed, 20);
        sethistostyle(histograms2.back(), kBlue, 21);
        histograms1.back()->GetXaxis()->SetTitle("Ntrk");
        histograms2.back()->GetXaxis()->SetTitle("Ntrk");
    }

    TCanvas *c = new TCanvas("c", "Comparison Canvas", 900, 600);
    c->Divide(3, 2);
    int padIndex = 1;
    for (const char *histogram : {"hv2_cent0", "hv2_cent1", "hv2_cent2"})
    {
        TH1D *hist1 = (TH1D *)TFile::Open(fileName1)->Get(histogram);
        TH1D *hist2 = (TH1D *)TFile::Open(fileName2)->Get(histogram);

        TH1D *hist3;
        if (fileName3 && strlen(fileName3) > 0)
        {
            hist3 = (TH1D *)TFile::Open(fileName3)->Get(histogram);
            hist3->SetDirectory(0); // Detach from file to avoid deletion when file is closed
            sethistostyle(hist3, colors[2], 22);
        }

        if (hscaling1 && hscaling2 && hscaling3)
        {
            hist1->Scale(1.0 / sqrt(hscaling1->GetBinContent(padIndex)));
            hist2->Scale(1.0 / sqrt(hscaling2->GetBinContent(padIndex)));
            hist3->Scale(1.0 / sqrt(hscaling3->GetBinContent(padIndex)));
        }

        if (!hist1 || !hist2)
        {
            std::cerr << "Error: One of the histograms does not exist in the files." << std::endl;
            continue;
        }

        sethistostyle(hist1, colors[0], 20);
        sethistostyle(hist2, colors[1], 21);

        for (int ptindex = 1; ptindex <= nptbins; ++ptindex)
        {
            histograms1[ptindex - 1]->SetBinContent(padIndex, hist1->GetBinContent(ptindex));
            histograms2[ptindex - 1]->SetBinContent(padIndex, hist2->GetBinContent(ptindex));
            histograms1[ptindex - 1]->SetBinError(padIndex, hist1->GetBinError(ptindex));
            histograms2[ptindex - 1]->SetBinError(padIndex, hist2->GetBinError(ptindex));
        }

        auto hist1_uncertainty = (TH1D *)hist1->Clone("hist1_uncertainty");
        auto hist2_uncertainty = (TH1D *)hist2->Clone("hist2_uncertainty");
        hist1_uncertainty->GetYaxis()->SetTitle("#it{v}_{2} Uncertainty");
        hist2_uncertainty->GetYaxis()->SetTitle("#it{v}_{2} Uncertainty");
        hist1_uncertainty->GetYaxis()->SetRangeUser(0, 0.02);
        hist2_uncertainty->GetYaxis()->SetRangeUser(0, 0.02);

        for (int i = 1; i <= hist1->GetNbinsX(); ++i)
        {
            double error1 = hist1->GetBinError(i);
            double error2 = hist2->GetBinError(i);
            hist1_uncertainty->SetBinContent(i, error1);
            hist2_uncertainty->SetBinContent(i, error2);
            hist1_uncertainty->SetBinError(i, 1.e-6); // Set a small error for visibility
            hist2_uncertainty->SetBinError(i, 1.e-6); // Set a small error for visibility
        }

        TLegend *legend = new TLegend(0.6, 0.6, 0.8, 0.8);
        legend->AddEntry(hist1, leg_labels[0].c_str(), "lp");
        legend->AddEntry(hist2, leg_labels[1].c_str(), "lp");
        if (fileName3 && strlen(fileName3) > 0 && padIndex != 1)
        {
            legend->AddEntry(hist3, leg_labels[2].c_str(), "lp");
        }
        legend->SetTextSize(0.04);
        legend->SetTextFont(42);
        legend->SetBorderSize(0);

        c->cd(padIndex);
        hist1->GetYaxis()->SetRangeUser(-0.01, 0.22);
        hist1->GetYaxis()->SetDecimals(true);
        hist1->Draw("E");
        hist2->Draw("E same");
        if (fileName3 && strlen(fileName3) > 0)
        {
            hist3->Draw("E same");
        }
        if (padIndex == 2)
            legend->Draw();

        auto ratio = (TH1D *)hist2->Clone("ratio");
        ratio->Divide(hist2, hist1, 1, 1, "E");
        ratio->SetTitle(" ");
        ratio->GetYaxis()->SetTitle("Ratio");
        ratio->GetXaxis()->SetTitle(hist1->GetXaxis()->GetTitle());
        ratio->GetYaxis()->SetRangeUser(.8, 1.2);
        ratio->SetMarkerStyle(20);

        auto ratio3 = (TH1D *)hist3->Clone("ratio3");
        if (fileName3 && strlen(fileName3) > 0)
        {
            ratio3->Divide(hist3, hist1, 1, 1, "E");
            sethistostyle(ratio3, colors[2], 22);
        }

        // fit the ratio
        // TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x", hist1->GetXaxis()->GetXmin(), hist1->GetXaxis()->GetXmax());
        // ratio->Fit(fitFunc, "QSE0");
        // fitFunc->SetLineColor(kRed);
        //
        c->cd(padIndex + 3);
        // set grid
        c->cd(padIndex + 3)->SetGridy();
        ratio->Draw("E");
        if (fileName3 && strlen(fileName3) > 0)
        {
            ratio3->Draw("E same");
        }
        // fitFunc->Draw("same");

        // Add the fit parameters to the legend
        // TLatex *fitParam0 = new TLatex(0.2, 0.8, Form("[0] = %.3f #pm %.3f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
        // TLatex *fitParam1 = new TLatex(0.2, 0.75, Form("[1] = %.3f #pm %.3f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
        // TLatex *fitChi2 = new TLatex(0.2, 0.7, Form("#chi^{2}/ndf = %.2f", fitFunc->GetChisquare() / fitFunc->GetNDF()));
        // fitParam0->SetNDC();
        // fitParam1->SetNDC();
        // fitChi2->SetNDC();
        // fitParam0->SetTextSize(0.04);
        // fitParam1->SetTextSize(0.04);
        // fitChi2->SetTextSize(0.04);
        // fitParam0->SetTextColor(kRed);
        // fitParam1->SetTextColor(kRed);
        // fitChi2->SetTextColor(kRed);
        //
        // fitParam0->Draw();
        // fitParam1->Draw();
        // fitChi2->Draw();

        // Add the uncertainty histograms
        // c->cd(padIndex + 6);
        // hist1_uncertainty->Draw("E");
        // hist2_uncertainty->Draw("E same");

        padIndex++;
    }
    c->Update();
    /*
    TCanvas *c2 = new TCanvas("c2", "Comparison Canvas 2", 900, 600);
    c2->Divide(3, 2);
    for (int i = 0; i < nptbins; ++i)
    {
        c2->cd(i + 1);
        if (i == 0)
        {
            TLegend *legend = new TLegend(0.6, 0.7, 0.8, 0.9);
            legend->AddEntry(histograms1[i], "w/ sub", "lp");
            legend->AddEntry(histograms2[i], "w/o sub", "lp");
            legend->SetTextSize(0.04);
            legend->SetTextFont(42);
            legend->SetBorderSize(0);
            legend->Draw();
        }
        histograms1[i]->SetTitle(Form("Comparison of v2 vs mass for pt bin %f - %f GeV/c", pt_values[i], pt_values[i + 1]));
        histograms1[i]->Draw("E");
        histograms2[i]->Draw("E same");
        histograms1[i]->GetYaxis()->SetRangeUser(-0.01, 0.02);
        histograms1[i]->GetYaxis()->SetTitle("#it{v}_{2}");
        histograms1[i]->GetYaxis()->SetTitleOffset(1.6);
    }

    */
    // Save the canvas to the output file
    outFile->cd();
    c->Write();
    c->SaveAs(Form("%s.pdf", outputFileName));
    // c2->Write();
    outFile->Close();
}
