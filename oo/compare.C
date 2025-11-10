#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"

void plot(const char *histogram = "v2_vs_mass", const char *fileName1 = "massandv2_sub.root", const char *fileName2 = "massandv2.root", TFile *outFile = nullptr)
{
    auto file1 = TFile::Open(fileName1);
    auto file2 = TFile::Open(fileName2);

    auto canvas = new TCanvas(Form("canvas_%s", histogram), Form("Canvas for %s", histogram), 800, 600);
    canvas->Divide(5, 3);
    auto canvas_ratio = new TCanvas(Form("canvas_ratio_%s", histogram), Form("Canvas for ratio of %s", histogram), 800, 600);
    canvas_ratio->Divide(5, 3);

    int currentPad = 1;

    int imult = 0;
    while (imult < 10)
    {
        int itrig = 0;
        while (true)
        {
            Printf("%d %d", imult, itrig);

            TString tmp;
            tmp.Form("%s_cent%d_pt%d", histogram, imult, itrig);

            auto graph1 = (TGraphErrors *)file1->Get(tmp);
            if (!graph1)
            {
                break;
            }
            else
            {
                graph1->SetTitle(Form("%s, mult: %d, trig: %d", histogram, imult, itrig));
            }

            canvas->cd(currentPad);
            graph1->Draw("AP");

            auto graph2 = (TGraphErrors *)file2->Get(tmp);
            if (graph2)
            {
                graph2->SetLineColor(2);
                graph2->SetMarkerColor(2);
                graph2->SetMarkerStyle(4);
                graph2->Draw("PSAME");
            }

            auto graphratio = (TGraphErrors *)graph1->Clone(Form("%s_ratio", tmp.Data()));
            graphratio->SetLineColor(4);
            graphratio->SetMarkerColor(4);
            graphratio->SetMarkerStyle(4);
            if (graphratio)
            {
                for (int i = 0; i < graphratio->GetN(); i++)
                {
                    double x, y, ex, ey;
                    graphratio->GetPoint(i, x, y);
                    ex = graphratio->GetErrorX(i);
                    ey = graphratio->GetErrorY(i);

                    double y2 = graph2->GetY()[i];
                    double ey2 = graph2->GetErrorY(i);

                    if (y2 != 0)
                    {
                        graphratio->SetPoint(i, x, y / y2);
                        graphratio->SetPointError(i, ex, TMath::Sqrt((ey / y) * (ey / y) + (ey2 / y2) * (ey2 / y2)));
                    }
                }
            }

            canvas_ratio->cd(currentPad);
            graphratio->Draw("AP");

            itrig++;
            currentPad++;
        }

        imult++;

        canvas->Update();
        canvas_ratio->Update();
    }

    // Save the canvases to files
    outFile->cd();
    canvas->Write();
    canvas_ratio->Write();
}

/// Compares two files containing histograms and prints the differences.
/// The histograms to compare can be specified by the `histograms` parameter.
/// If no histograms are specified, it defaults to comparing "v2_vs_mass".
///
/// @param fileName1 The first file to compare.
/// @param fileName2 The second file to compare.
/// @param histograms A comma-separated list of histogram names to compare.
void compare(const char *fileName1 = "/Users/spolitan/alice/alice-correlations/oo/massandv2lambdaSub.root", const char *fileName2 = "/Users/spolitan/alice/alice-correlations/oo/massandv2lambda.root", std::vector<const char *> histograms = {"v1_vs_mass", "v2_vs_mass", "v3_vs_mass", "v4_vs_mass", "v5_vs_mass", "Chi2/ndf_vs_mass", "c_vs_mass"},
             const char *outputFileName = "comparison_results_lambdas.root")
{
    auto outFile = TFile::Open(outputFileName, "RECREATE");
    if (!outFile || outFile->IsZombie())

        outFile->cd();
    for (const char *histogram : histograms)
    {
        plot(histogram, fileName1, fileName2, outFile);
    }
    outFile->Close();
}
