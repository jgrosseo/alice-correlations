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
void compare_v2_ft0m( // const char *fileName1 = "/Users/spolitan/alice/alice-correlations/oo/k0_hm_ft0m/v2ofK0sh_ft0m.root",
                      // const char *fileName2 = "/Users/spolitan/Downloads/For_prottay_Stefano/HMextraction/v2ofK0sh.root",
    // const char *fileName1 = "/Users/spolitan/alice/alice-correlations/oo/phi_hm_ft0m/v2ofphih_LMsubtraction.root",
    const char *fileName1 = "/Users/spolitan/alice/alice-correlations/oo/lambda_hm_ft0m/combinedfit_fulldetKs_sourav5.root",
    const char *fileName2 = "/Users/spolitan/alice/alice-correlations/oo/combinedfit_fulldetKs_sourav5.root",
    // const char *fileName2 = "/Users/spolitan/Downloads/stefano_phi_final/Invmassfit_subtracted/v2ofphih_LMsubtraction.root",
    const char *sclaing1 = "/Users/spolitan/alice/alice-correlations/oo/hh_all_ft0m/out_sub2hadron_deltaeta1.2.root",
    const char *sclaing2 = "/Users/spolitan/Downloads/stefano_phi_final/hadron_extraction/out_sub2hadron.root",
    const char *outputFileName = "/Users/spolitan/alice/alice-correlations/oo/lambda_hm_ft0m/comparison_v2_lambda_ft0m",
    const bool isPhi = false)
{

    std::vector<std::string> mult_bins = {"25-50", "50-100", "100-300"};
    std::vector<std::string> leg_labels = {
        "10-20% (FT0M)",
        "0-10% (FT0M)",
        "Nch 50-150",
        "Nch 150-300",
    };

    auto outFile = TFile::Open(Form("%s.root", outputFileName), "RECREATE");

    std::vector<int> colors = {
        TColor::GetColorTransparent(kRed + 1, 0.55), TColor::GetColorTransparent(kAzure + 4, 0.55), TColor::GetColorTransparent(kGreen + 2, 0.55), TColor::GetColorTransparent(kOrange + 1, 0.55)};

    // th1 vs mult
    std::vector<double> mult_values1 = {0., 10.0, 20.0};
    std::vector<double> mult_values2 = {25., 50.0, 150.0, 300.0};
    std::vector<double> pt_values = {0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0};
    int nptbins = 6; // Number of pT bins
    std::vector<TH1D *> histograms1;
    std::vector<TH1D *> histograms2;
    std::vector<TH1D *> histograms3;

    // scaling
    TH1D *hscaling1 = nullptr;
    TH1D *hscaling2 = nullptr;
    TH1D *hscaling3 = nullptr;
    if (strlen(sclaing1) != 0 && strlen(sclaing2) != 0)
    {
        hscaling1 = (TH1D *)TFile::Open(sclaing1)->Get("v2vsmult");
        hscaling2 = (TH1D *)TFile::Open(sclaing2)->Get("v2vsmult");
    }

    const int centBin[] = {1, 0}; // 0-10%, 10-20% for ft0m
    TCanvas *c = new TCanvas("c", "Comparison Canvas", 900, 600);
    c->Divide(2, 2);
    int padIndex = 1;
    for (int ihist = 1; ihist < 3; ++ihist)
    {

        int centindex = centBin[ihist - 1];
        const char *histogram1 = Form("hv2_cent%d", centindex);
        const char *histogram2 = Form("hv2_cent%d", ihist);

        cout << "Comparing " << histogram1 << " and " << histogram2 << "\n";
        TH1D *hist1 = (TH1D *)TFile::Open(fileName1)->Get(histogram1);
        TH1D *hist2 = (TH1D *)TFile::Open(fileName2)->Get(histogram2);

        if (hscaling1 && hscaling2)
        {
            hist1->Scale(1.0 / sqrt(hscaling1->GetBinContent(centBin[ihist - 1] + 1)));
            if (!isPhi)
            {
                if (ihist != 2)
                    hist2->Scale(1.0 / sqrt(hscaling2->GetBinContent(ihist + 1)));
                else // for 150-300
                    hist2->Scale(1.0 / sqrt(hscaling2->GetBinContent(4)));
            }
        }

        if (!hist1 || !hist2)
        {
            std::cerr << "Error: One of the histograms does not exist in the files." << std::endl;
            continue;
        }

        sethistostyle(hist1, colors[ihist - 1], 20);
        sethistostyle(hist2, colors[ihist + 1], 21);

        TLegend *legend = new TLegend(0.6, 0.6, 0.8, 0.8);
        legend->AddEntry(hist1, leg_labels[ihist - 1].c_str(), "lp");
        legend->AddEntry(hist2, leg_labels[ihist + 1].c_str(), "lp");
        legend->SetTextSize(0.04);
        legend->SetTextFont(42);
        legend->SetBorderSize(0);

        c->cd(padIndex);
        hist1->GetYaxis()->SetRangeUser(-0.01, 0.22);
        hist1->GetYaxis()->SetDecimals(true);
        hist1->Draw("E");
        hist2->Draw("E same");
        legend->Draw();

        if (padIndex == 2)
        {
            TLatex *lat = new TLatex();
            lat->SetNDC();
            lat->SetTextSize(0.04);
            lat->DrawLatex(0.2, 0.85, "LM FT0M: 70-90%");
            lat->DrawLatex(0.2, 0.80, "LM Nch: 0-25");
        }

        auto ratio = (TH1D *)hist2->Clone("ratio");
        ratio->Divide(hist2, hist1, 1, 1, "E");
        ratio->SetTitle(" ");
        ratio->GetYaxis()->SetTitle("Ratio");
        ratio->GetXaxis()->SetTitle(hist1->GetXaxis()->GetTitle());
        ratio->GetYaxis()->SetRangeUser(.8, 1.2);
        ratio->SetMarkerStyle(20);

        // fit the ratio
        // TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x", hist1->GetXaxis()->GetXmin(), hist1->GetXaxis()->GetXmax());
        // ratio->Fit(fitFunc, "QSE0");
        // fitFunc->SetLineColor(kRed);
        //
        c->cd(padIndex + 2);
        // set grid
        c->cd(padIndex + 2)->SetGridy();
        ratio->Draw("E");
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
