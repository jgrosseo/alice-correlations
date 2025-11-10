#include <iostream>
#include <cmath>

using namespace std;

TLatex *latexExp = new TLatex();
latexExp->SetTextSize(0.038);
latexExp->SetTextFont(42);
latexExp->SetNDC();
TLatex *latexDetail = new TLatex();
latexDetail->SetTextSize(0.032);
latexDetail->SetTextFont(42);
latexDetail->SetNDC();

gStyle->SetGridColor(kGray + 2);

TLegend *DrawLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{

  TLegend *legend = new TLegend(x1, y1, x2, y2);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);

  return legend;
}

TCanvas *DrawCanvas(TString opt = "c")
{
  TCanvas *c1 = new TCanvas(opt.Data(), opt.Data(), 20, 20, 600, 600);
  c1->cd(1);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.15);
  return c1;
}

void SetHistoStyle(TH1D *h1, Int_t MarkerColor = 1, Int_t MarkerStyle = 20, const char *yaxistitle = "")
{

  h1->SetTitle("");
  h1->SetMarkerStyle(MarkerStyle);
  h1->SetMarkerColor(MarkerColor);
  h1->SetLineColor(MarkerColor);
  h1->SetMarkerSize(1.4);
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h1->GetXaxis()->CenterTitle(false);
  h1->GetXaxis()->SetNdivisions(506);
  h1->GetYaxis()->SetNdivisions(505);
  h1->GetXaxis()->SetLabelOffset(0.015);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetTickLength(0.04);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTitleOffset(1.6);
  h1->GetYaxis()->SetDecimals(false);
  // h1->GetYaxis()->SetNdivisions(310);
  h1->GetYaxis()->SetLabelOffset(0.015);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.05);
  h1->GetYaxis()->SetTickLength(0.04);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitle(yaxistitle);
  // h1->GetYaxis()->SetTitle("Ratio^{#frac{#Lambda}{#bar{#Lambda}}}");
  h1->GetYaxis()->SetTitleFont(42);
}

void SetGraphStyle(TGraphErrors *g1, Int_t MarkerColor = 1, Int_t MarkerStyle = 20, const char *yaxistitle = "")
{

  g1->SetTitle("");
  g1->SetMarkerStyle(MarkerStyle);
  g1->SetMarkerColor(MarkerColor);
  g1->SetLineColor(MarkerColor);
  g1->SetMarkerSize(1.4);
  g1->SetLineWidth(2);
  g1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  g1->GetXaxis()->CenterTitle(false);
  g1->GetXaxis()->SetNdivisions(506);
  g1->GetYaxis()->SetNdivisions(505);
  g1->GetXaxis()->SetLabelOffset(0.015);
  g1->GetXaxis()->SetLabelFont(42);
  g1->GetXaxis()->SetTitleFont(42);
  g1->GetXaxis()->SetLabelSize(0.05);
  g1->GetXaxis()->SetTitleSize(0.05);
  g1->GetXaxis()->SetTickLength(0.04);
  g1->GetXaxis()->SetTitleOffset(1.2);
  g1->GetYaxis()->SetTitleOffset(1.6);
  g1->GetYaxis()->SetDecimals(false);
  // g1->GetYaxis()->SetNdivisions(310);
  g1->GetYaxis()->SetLabelOffset(0.015);
  g1->GetYaxis()->SetLabelFont(42);
  g1->GetYaxis()->SetLabelSize(0.05);
  g1->GetYaxis()->SetTickLength(0.04);
  g1->GetYaxis()->SetTitleSize(0.05);
  if (yaxistitle != nullptr)
  {
    g1->GetYaxis()->SetTitle(yaxistitle);
    // g1->GetYaxis()->SetTitle("Ratio^{#frac{#Lambda}{#bar{#Lambda}}}");
    g1->GetYaxis()->SetTitleFont(42);
  }
}

TH1D *Plothisto2(TH1D *h1, Int_t MarkerColor = 1, Int_t MarkerStyle = 20)
{

  h1->SetTitle("");
  h1->SetMarkerStyle(MarkerStyle);
  h1->SetMarkerColor(MarkerColor);
  h1->SetLineColor(MarkerColor);
  h1->SetMarkerSize(0.0);
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h1->GetXaxis()->CenterTitle(false);
  h1->GetXaxis()->SetNdivisions(506);
  h1->GetYaxis()->SetNdivisions(505);
  h1->GetXaxis()->SetLabelOffset(0.015);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetTickLength(0.04);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTitleOffset(1.85);
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetDecimals(false);
  // h1->GetYaxis()->SetNdivisions(310);
  h1->GetYaxis()->SetLabelOffset(0.015);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.05);
  h1->GetYaxis()->SetTickLength(0.04);
  h1->GetYaxis()->SetTitleSize(0.05);
  // h1->GetYaxis()->SetTitle("#chi^{2}/Ndf");
  h1->GetYaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
  h1->GetYaxis()->SetTitleFont(42);
  return h1;
}

void MyStyle()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFillColor(0);
  gStyle->SetLineColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleColor(1);
  gStyle->SetLineWidth(2);
  gStyle->SetErrorX(0.4);
}

TGraphErrors *DrawFrame(TGraphErrors *h, Int_t MCol, Int_t MSty, Bool_t mrk = 0)
{
  // h->GetXaxis()->SetTitle("d#it{N}_{ch}/d#eta");
  h->GetXaxis()->SetTitle("<#it{N}_{part}>");
  if (mrk)
  {
    h->SetMarkerColor(MCol);
    h->SetMarkerStyle(MSty);
    h->SetMarkerSize(1.4);
  }
  h->SetLineColor(MCol);
  h->SetLineWidth(2);
  h->GetXaxis()->CenterTitle(false);
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetLabelOffset(0.017);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTickLength(0.06);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->CenterTitle(false);
  h->GetYaxis()->SetDecimals(true);
  // h->GetYaxis()->SetNdivisions(310);
  h->GetYaxis()->SetLabelOffset(0.015);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTickLength(0.04);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleFont(42);

  return h;
}

void getPartV2vsPt(const char *inputFileName, std::vector<TH1D *> &hPartV2vsPt, const std::vector<int> &centBins, const std::vector<int> &colors, const std::vector<int> &markers, const char *particle, const char *v2ofparticle)
{
  TFile *inputFileRef = TFile::Open(v2ofparticle);
  if (!inputFileRef)
  {
    std::cerr << "Error opening reference file: " << v2ofparticle << std::endl;
    return;
  }

  TH1D *hRef = (TH1D *)inputFileRef->Get("v2vsmult");
  if (!hRef)
  {
    std::cerr << "Reference histogram hV2_vs_cent0 not found in file: " << v2ofparticle << std::endl;
    return;
  }
  hRef->SetDirectory(0); // Detach from file to avoid deletion when file is closed
  inputFileRef->Close();

  TFile *inputFile = TFile::Open(inputFileName);
  if (!inputFile)
  {
    std::cerr << "Error opening file: " << inputFileName << std::endl;
    return;
  }

  for (size_t i = 0; i < centBins.size() - 1; i++)
  {
    TString histName = Form("hv2_cent%i", int(i));
    TH1D *h = (TH1D *)inputFile->Get(histName.Data());
    if (!h)
    {
      std::cerr << "Histogram " << histName.Data() << " not found in file." << std::endl;
      continue;
    }
    hPartV2vsPt.push_back(h);
    hPartV2vsPt[i]->SetDirectory(0); // Detach from file to avoid deletion when file is closed
    if (particle != std::string("Phi"))
    {
      std::cout << "Scaling histogram " << histName.Data() << " by 1/sqrt(" << hRef->GetBinContent(3) << ")" << std::endl;
      hPartV2vsPt[i]->Scale(1.0 / sqrt(hRef->GetBinContent(4))); // Normalize to reference histogram
    }
    SetHistoStyle(hPartV2vsPt[i], colors[i], markers[i], Form("#it{v}_{2}(%s)", particle));
  }
}

void v2plot(const char *infileLambda = "/Users/spolitan/alice/alice-correlations/oo/combinedfit_fulldetKs_sourav5.root",
            const char *infileAntiLambda = "", //"/Users/spolitan/alice/alice-correlations/oo/antilambda_full_pass2_021025_lm_050/combinedfit_fulldetAntiLambda_sourav.root",
            const char *infileK0s = "/Users/spolitan/Downloads/v2ofK0sh.root",
            const char *infilePhi = "/Users/spolitan/Downloads/stefano_phi_final/Invmassfit_subtracted/v2ofphih_LMsubtraction.root",
            const char *v2ofparticle = "/Users/spolitan/Downloads/stefano_phi_final/hadron_extraction/out_sub2hadron.root",
            const char *outdir = "/Users/spolitan/alice/alice-correlations/oo/lambda_full_pass2_021025_lm_050/",
            const string &suffix = "full_all")
{

  MyStyle();

  // Style and binning
  std::vector<int> centBins = {25, 50, 150, 300};
  // std::vector<double> colors = {
  //     kRed + 1,
  //     kAzure + 4,
  //     kGreen + 2,
  //     kOrange + 1,
  // };
  // color with trasparency

  // Define colors with transparency (alpha = 0.6)
  std::vector<int> colors = {
      TColor::GetColorTransparent(kRed + 1, 0.55), TColor::GetColorTransparent(kAzure + 4, 0.55), TColor::GetColorTransparent(kGreen + 2, 0.55), TColor::GetColorTransparent(kOrange + 1, 0.55)};

  std::vector<int> markers = {
      kFullSquare,
      kFullCircle,
      kFullDiamond,
      kOpenSquare};
  double xOffset = 0.20;

  // Load files if provided
  std::vector<TH1D *> hPartV2vsPtLambda;
  std::vector<TH1D *> hPartV2vsPtAntiLambda;
  std::vector<TH1D *> hPartV2vsPtK0s;
  std::vector<TH1D *> hPartV2vsPtPhi;
  if (strlen(infileLambda) > 0)
  {
    getPartV2vsPt(infileLambda, hPartV2vsPtLambda, centBins, colors, markers, "Lambda", v2ofparticle);
  }
  if (strlen(infileAntiLambda) > 0)
  {
    getPartV2vsPt(infileAntiLambda, hPartV2vsPtAntiLambda, centBins, colors, markers, "AntiLambda", v2ofparticle);
  }
  if (strlen(infileK0s) > 0)
  {
    getPartV2vsPt(infileK0s, hPartV2vsPtK0s, centBins, colors, markers, "K0s", v2ofparticle);
  }
  if (strlen(infilePhi) > 0)
  {
    getPartV2vsPt(infilePhi, hPartV2vsPtPhi, centBins, colors, markers, "Phi", v2ofparticle);
  }

  TFile *outputFile = TFile::Open(Form("%s/v2_plots_%s.root", outdir, suffix.c_str()), "RECREATE");

  // Draw hPartV2vsPt for each particle separately
  if (!hPartV2vsPtLambda.empty())
  {
    TLegend *lp = DrawLegend(0.2, 0.55, 0.3, 0.75);
    lp->SetHeader("Nch");
    TCanvas *cLambda = DrawCanvas("cLambda");
    for (size_t i = 0; i < hPartV2vsPtLambda.size(); i++)
    {
      if (i == 0)
      {
        hPartV2vsPtLambda[i]->GetYaxis()->SetRangeUser(0, 0.3);
        hPartV2vsPtLambda[i]->GetYaxis()->SetTitle("#it{v}_{2} (#Lambda+#bar{#Lambda})");
        hPartV2vsPtLambda[i]->GetXaxis()->SetRangeUser(0.2, 20);
        hPartV2vsPtLambda[i]->Draw("p");
      }
      else
      {
        hPartV2vsPtLambda[i]->Draw("p same");
      }
      lp->AddEntry(hPartV2vsPtLambda[i], Form("%d-%d", centBins[i], centBins[i + 1]), "p");
      hPartV2vsPtLambda[i]->Write(Form("v2_lambda_cent%i", int(i)));
    }
    lp->Draw();
    latexExp->DrawLatex(xOffset, 0.88, "ALICE, O#minusO, #sqrt{#it{s}}_{NN}=5.36 TeV");
    latexDetail->DrawLatex(xOffset, 0.83, "LHC25ae_cpass2");
    latexDetail->DrawLatex(xOffset, 0.78, "|#Delta#it{#eta}| > 1.2");
    cLambda->Write();
    cLambda->SaveAs(Form("%s/v2_lambda_centrality_%s.pdf", outdir, suffix.c_str()));
  }

  if (!hPartV2vsPtAntiLambda.empty())
  {
    TLegend *lp = DrawLegend(0.2, 0.55, 0.3, 0.75);
    lp->SetHeader("Nch");
    TCanvas *cAntiLambda = DrawCanvas("cAntiLambda");
    for (size_t i = 0; i < hPartV2vsPtAntiLambda.size(); i++)
    {
      if (i == 0)
      {
        hPartV2vsPtAntiLambda[i]->GetYaxis()->SetRangeUser(0, 0.3);
        hPartV2vsPtAntiLambda[i]->GetYaxis()->SetTitle("#it{v}_{2} (#Lambda+#bar{#Lambda})");
        hPartV2vsPtAntiLambda[i]->GetXaxis()->SetRangeUser(0.2, 20);
        hPartV2vsPtAntiLambda[i]->Draw("p");
      }
      else
      {
        hPartV2vsPtAntiLambda[i]->Draw("p same");
      }
      lp->AddEntry(hPartV2vsPtAntiLambda[i], Form("%d-%d", centBins[i], centBins[i + 1]), "p");
      hPartV2vsPtAntiLambda[i]->Write(Form("v2_antilambda_cent%i", int(i)));
    }
    lp->Draw();
    latexExp->DrawLatex(xOffset, 0.88, "ALICE, O#minusO, #sqrt{#it{s}}_{NN}=5.36 TeV");
    latexDetail->DrawLatex(xOffset, 0.83, "LHC25ae_cpass2");
    latexDetail->DrawLatex(xOffset, 0.78, "|#Delta#it{#eta}| > 1.2");
    cAntiLambda->Write();
    cAntiLambda->SaveAs(Form("%s/v2_antilambda_centrality_%s.pdf", outdir, suffix.c_str()));
  }

  if (!hPartV2vsPtK0s.empty())
  {
    TLegend *lp = DrawLegend(0.2, 0.55, 0.3, 0.75);
    lp->SetHeader("Nch");
    TCanvas *cK0s = DrawCanvas("cK0s");
    for (size_t i = 0; i < hPartV2vsPtK0s.size(); i++)
    {
      if (i == 0)
      {
        hPartV2vsPtK0s[i]->GetYaxis()->SetRangeUser(0, 0.3);
        hPartV2vsPtK0s[i]->GetXaxis()->SetRangeUser(0.2, 20);
        hPartV2vsPtK0s[i]->Draw("p");
      }
      else
      {
        hPartV2vsPtK0s[i]->Draw("p same");
      }
      lp->AddEntry(hPartV2vsPtK0s[i], Form("%d-%d", centBins[i], centBins[i + 1]), "p");
      hPartV2vsPtK0s[i]->Write(Form("v2_k0s_cent%i", int(i)));
    }
    lp->Draw();
    latexExp->DrawLatex(xOffset, 0.88, "ALICE, O#minusO, #sqrt{#it{s}}_{NN}=5.36 TeV");
    latexDetail->DrawLatex(xOffset, 0.83, "LHC25ae_cpass2");
    latexDetail->DrawLatex(xOffset, 0.78, "|#Delta#it{#eta}| > 1.2");
    cK0s->Write();
    cK0s->SaveAs(Form("%s/v2_k0s_centrality_%s.pdf", outdir, suffix.c_str()));
  }

  if (!hPartV2vsPtPhi.empty())
  {
    TLegend *lp = DrawLegend(0.2, 0.55, 0.3, 0.75);
    lp->SetHeader("Nch");
    TCanvas *cPhi = DrawCanvas("cPhi");
    for (size_t i = 0; i < hPartV2vsPtPhi.size(); i++)
    {
      if (i == 0)
      {
        hPartV2vsPtPhi[i]->GetYaxis()->SetRangeUser(0, 0.3);
        hPartV2vsPtPhi[i]->GetYaxis()->SetTitle("#it{v}_{2} (#phi)");
        hPartV2vsPtPhi[i]->GetXaxis()->SetRangeUser(0.2, 20);
        hPartV2vsPtPhi[i]->Draw("p");
      }
      else
      {
        hPartV2vsPtPhi[i]->Draw("p same");
      }
      lp->AddEntry(hPartV2vsPtPhi[i], Form("Nch: %i-%i%%", centBins[i], centBins[i + 1]), "p");
      hPartV2vsPtPhi[i]->Write(Form("v2_phi_cent%i", int(i)));
    }
    lp->Draw();
    latexExp->DrawLatex(xOffset, 0.88, "ALICE, O#minusO, #sqrt{#it{s}}_{NN}=5.36 TeV");
    latexDetail->DrawLatex(xOffset, 0.83, "LHC25ae_cpass2");
    latexDetail->DrawLatex(xOffset, 0.78, "|#Delta#it{#eta}| > 1.2");
    cPhi->Write();
    cPhi->SaveAs(Form("%s/v2_phi_centrality_%s.pdf", outdir, suffix.c_str()));
  }

  // dummy tgraph to plot prottay's points from x, y map (k0s)
  // define map
  std::map<double, double> xytmp = {
      {0.7553956834532372, 0.060886699507389175}, {1.2589928057553954, 0.10463054187192117}, {1.762589928057554, 0.12650246305418716}, {2.5077081192189103, 0.13566502463054186}, {3.509763617677287, 0.13448275862068965}, {4.5015416238437815, 0.12413793103448276}};

  TGraph *gtmp = new TGraph();
  for (const auto &xy : xytmp)
  {
    gtmp->SetPoint(gtmp->GetN(), xy.first, xy.second);
  }
  gtmp->SetMarkerSize(1.0);
  gtmp->SetMarkerStyle(21);
  gtmp->SetMarkerColor(kBlack);
  gtmp->SetLineColor(kBlack);
  gtmp->SetLineWidth(2);

  // dummy tgraph to plot prottay's points from x, y map (lambda)
  // define map
  std::map<double, double> xytmp2 = {{0.7502569373072969, 0.03871921182266008}, {1.2589928057553954, 0.0771428571428571}, {1.762589928057554, 0.1223645320197044}, {2.5077081192189103, 0.15783251231527093}, {3.509763617677287, 0.1850246305418719}, {4.5015416238437815, 0.1921182266009852}};
  TGraph *gtmp2 = new TGraph();
  for (const auto &xy : xytmp2)
  {
    gtmp2->SetPoint(gtmp2->GetN(), xy.first, xy.second);
  }
  gtmp2->SetMarkerSize(1.0);
  gtmp2->SetMarkerStyle(22);
  // gtmp2->SetMarkerColor(kGray + 2);
  // gtmp2->SetLineColor(kGray + 2);

  // Draw all particles together for one centrality range
  int centralityIndex = 2; // Change this to select different centrality ranges
  TCanvas *cAll = DrawCanvas("cAll");
  // cAll->SetGrid();

  TLegend *lpAllCent = DrawLegend(0.20, 0.55, 0.45, 0.75);
  if (!hPartV2vsPtPhi.empty())
  {
    hPartV2vsPtPhi[centralityIndex]->SetLineColor(colors[2]);
    hPartV2vsPtPhi[centralityIndex]->SetMarkerColor(colors[2]);
    hPartV2vsPtPhi[centralityIndex]->SetMarkerStyle(kFullCross);
    hPartV2vsPtPhi[centralityIndex]->GetYaxis()->SetTitle("#it{v}_{2}");
    hPartV2vsPtPhi[centralityIndex]->GetYaxis()->SetNdivisions(510);
    hPartV2vsPtPhi[centralityIndex]->GetYaxis()->SetDecimals(true);
    hPartV2vsPtPhi[centralityIndex]->GetYaxis()->SetRangeUser(0, 0.3);
    hPartV2vsPtPhi[centralityIndex]->GetXaxis()->SetRangeUser(0.2, 8);
    hPartV2vsPtPhi[centralityIndex]->Draw("p");
    lpAllCent->AddEntry(hPartV2vsPtPhi[2], "#phi", "p");
  }
  if (!hPartV2vsPtK0s.empty())
  {
    hPartV2vsPtK0s[centralityIndex]->SetLineColor(colors[1]);
    hPartV2vsPtK0s[centralityIndex]->SetMarkerColor(colors[1]);
    hPartV2vsPtK0s[centralityIndex]->SetMarkerStyle(kFullCircle);
    hPartV2vsPtK0s[centralityIndex]->GetYaxis()->SetTitle("#it{v}_{2}");
    hPartV2vsPtK0s[centralityIndex]->GetYaxis()->SetNdivisions(510);
    hPartV2vsPtK0s[centralityIndex]->GetYaxis()->SetDecimals(true);
    hPartV2vsPtK0s[centralityIndex]->GetYaxis()->SetRangeUser(0, 0.3);
    hPartV2vsPtK0s[centralityIndex]->GetXaxis()->SetRangeUser(0.2, 20);
    hPartV2vsPtK0s[centralityIndex]->Draw("p same");
    lpAllCent->AddEntry(hPartV2vsPtK0s[centralityIndex], "K0s", "p");
  }
  if (!hPartV2vsPtAntiLambda.empty())
  {
    hPartV2vsPtAntiLambda[centralityIndex]->SetLineColor(colors[3]);
    hPartV2vsPtAntiLambda[centralityIndex]->SetMarkerColor(colors[3]);
    hPartV2vsPtAntiLambda[centralityIndex]->SetMarkerStyle(kFullTriangleUp);
    hPartV2vsPtAntiLambda[centralityIndex]->GetXaxis()->SetRangeUser(0.2, 20);
    hPartV2vsPtAntiLambda[centralityIndex]
        ->Draw("p same");
    lpAllCent->AddEntry(hPartV2vsPtAntiLambda[centralityIndex], "#bar{#Lambda}", "p");
  }
  if (!hPartV2vsPtLambda.empty())
  {
    hPartV2vsPtLambda[centralityIndex]->SetLineColor(colors[0]);
    hPartV2vsPtLambda[centralityIndex]->SetMarkerColor(colors[0]);
    hPartV2vsPtLambda[centralityIndex]->SetMarkerStyle(kFullSquare);
    hPartV2vsPtLambda[centralityIndex]->GetXaxis()->SetRangeUser(0.2, 20);
    hPartV2vsPtLambda[centralityIndex]->Draw("p same");
    lpAllCent->AddEntry(hPartV2vsPtLambda[centralityIndex], "#Lambda", "p");
  }
  // gtmp->Draw("p same");
  // gtmp2->Draw("p same");
  lpAllCent->Draw();
  latexExp->DrawLatex(xOffset, 0.88, "ALICE, O#minusO #sqrt{#it{s}}_{NN}=5.36 TeV");
  latexDetail->DrawLatex(xOffset, 0.83, "LHC25ae_cpass2");
  latexDetail->DrawLatex(xOffset, 0.78, Form("#it{N}_{ch}: %d-%d, |#Delta#it{#eta}| > 1.2", centBins[centralityIndex], centBins[centralityIndex + 1]));
  cAll->Write();
  cAll->SaveAs(Form("%s/v2_all_particles_centrality_%s.pdf", outdir, suffix.c_str()));

  // Draw ratios to K0s
  TCanvas *cAllratio = new TCanvas("cAllratio", "cAllratio", 800, 800);
  cAllratio->SetLeftMargin(0.15);
  cAllratio->SetRightMargin(0.05);
  cAllratio->SetTopMargin(0.05);
  cAllratio->SetBottomMargin(0.15);
  cAllratio->cd();
  cAllratio->SetGrid();
  lpAllCent->Draw();
  latexExp->DrawLatex(xOffset, 0.88, "ALICE, O#minusO, #sqrt{#it{s}}_{NN}=5.36 TeV");
  latexDetail->DrawLatex(xOffset, 0.83, "LHC25ae_cpass2");
  latexDetail->DrawLatex(xOffset, 0.78, Form("Nch: %d-%d, |#Delta#it{#eta}| > 1.2", centBins[centralityIndex], centBins[centralityIndex + 1]));
  TLine *line = new TLine(0, 1, 6, 1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  // add ratio to gtmp if K0s is available
  if (!hPartV2vsPtK0s.empty() && gtmp->GetN() > 0)
  {
    TGraph *gRatio = (TGraph *)gtmp->Clone("gRatio");
    for (int i = 0; i < gRatio->GetN(); i++)
    {
      double x, y;
      gRatio->GetPoint(i, x, y);
      // find corresponding bin in K0s histogram
      int bin = hPartV2vsPtK0s[centralityIndex]->GetXaxis()->FindBin(x);
      double k0sV2 = hPartV2vsPtK0s[centralityIndex]->GetBinContent(bin);
      if (k0sV2 != 0)
      {
        gRatio->SetPoint(i, x, y / k0sV2);
      }
      else
      {
        gRatio->SetPoint(i, x, 0);
      }
    }
    gRatio->SetMarkerSize(1.0);
    gRatio->SetMarkerStyle(21);
    gRatio->SetMarkerColor(kBlack);
    gRatio->SetLineColor(kBlack);
    gRatio->SetLineWidth(2);
    gRatio->Draw("p");
  }
  // add ratio to gtmp2 if Lambda is available
  if (!hPartV2vsPtLambda.empty() && gtmp2->GetN() > 0)
  {
    TGraph *gRatio2 = (TGraph *)gtmp2->Clone("gRatio2");
    for (int i = 0; i < 1; i++)
    {
      double x, y;
      gRatio2->GetPoint(i, x, y);
      // find corresponding bin in K0s histogram
      int bin = hPartV2vsPtK0s[centralityIndex]->GetXaxis()->FindBin(x);
      double k0sV2 = hPartV2vsPtK0s[centralityIndex]->GetBinContent(bin);
      if (k0sV2 != 0)
      {
        gRatio2->SetPoint(i, x, y / k0sV2);
      }
      else
      {
        gRatio2->SetPoint(i, x, 0);
      }
    }
    gRatio2->SetMarkerSize(1.0);
    gRatio2->SetMarkerStyle(22);
    gRatio2->SetMarkerColor(kGray + 2);
    gRatio2->SetLineColor(kGray + 2);
    gRatio2->SetLineWidth(2);
    gRatio2->Draw("p same");
  }
  line->Draw("same");
  cAllratio->Write();
  cAllratio->SaveAs(Form("%s/v2_all_particles_centrality_%s_ratio.pdf", outdir, suffix.c_str()));

  // Draw all particles together for one centrality range vs CMS
  // split ALICE and CMS legends
  TCanvas *cAllCMS = DrawCanvas("cAllCMS");
  TLegend *lpAllAlice = DrawLegend(0.20, 0.76, 0.45, 0.82);
  lpAllAlice->SetTextSize(0.035);
  TLegend *lpAllCMS = DrawLegend(0.20, 0.56, 0.45, 0.62);
  lpAllCMS->SetTextSize(0.035);

  TFile *inputFileCMSK0 = TFile::Open("/Users/spolitan/alice/alice-correlations/oo/cms_kzero_pPb_nch185_250_nobbcorr.root");
  TGraphErrors *hCMSK0 = (TGraphErrors *)inputFileCMSK0->Get("Table 1/Graph1D_y1");
  TFile *inputFileCMSLambda = TFile::Open("/Users/spolitan/alice/alice-correlations/oo/cms_lambda_pPb_nch185_250_nobbcorr.root");
  TGraphErrors *hCMSLambda = (TGraphErrors *)inputFileCMSLambda->Get("Table 2/Graph1D_y1");

  outputFile->cd();

  centralityIndex = 2; // Change this to select different centrality ranges

  SetGraphStyle(hCMSK0, kAzure + 4, 24, "#it{v}_{2}(K^{0}_{S})");
  SetGraphStyle(hCMSLambda, kRed + 1, 25, "#it{v}_{2}(#Lambda)");
  if (!hPartV2vsPtK0s.empty())
  {
    hPartV2vsPtK0s[centralityIndex]->SetLineColor(colors[1]);
    hPartV2vsPtK0s[centralityIndex]->SetMarkerColor(colors[1]);
    hPartV2vsPtK0s[centralityIndex]->SetMarkerStyle(kFullCircle);
    hPartV2vsPtK0s[centralityIndex]->GetYaxis()->SetRangeUser(0, 0.3);
    hPartV2vsPtK0s[centralityIndex]->GetXaxis()->SetRangeUser(0, 20);
    hPartV2vsPtK0s[centralityIndex]->GetYaxis()->SetTitle("#it{v}_{2}");
    hPartV2vsPtK0s[centralityIndex]->Draw("p");
    lpAllAlice->AddEntry(hPartV2vsPtK0s[centralityIndex], "K^{0}_{S}", "p");
  }
  if (!hPartV2vsPtLambda.empty())
  {
    hPartV2vsPtLambda[centralityIndex]->SetLineColor(colors[0]);
    hPartV2vsPtLambda[centralityIndex]->SetMarkerColor(colors[0]);
    hPartV2vsPtLambda[centralityIndex]->SetMarkerStyle(kFullSquare);
    hPartV2vsPtLambda[centralityIndex]->Draw("p same");
    lpAllAlice->AddEntry(hPartV2vsPtLambda[centralityIndex], "#Lambda", "p");
  }
  if (!hPartV2vsPtPhi.empty())
  {
    hPartV2vsPtPhi[centralityIndex]->SetLineColor(colors[2]);
    hPartV2vsPtPhi[centralityIndex]->SetMarkerColor(colors[2]);
    hPartV2vsPtPhi[centralityIndex]->SetMarkerStyle(kFullCross);
    hPartV2vsPtPhi[centralityIndex]->Draw("p same");
    lpAllAlice->AddEntry(hPartV2vsPtPhi[centralityIndex], "#phi", "p");
  }
  hCMSK0->Draw("epz same");
  hCMSLambda->Draw("epz same");

  xOffset = 0.20;
  latexExp->DrawLatex(xOffset, 0.88, "ALICE, O#minusO #sqrt{#it{s}}_{NN} = 5.36 TeV");
  latexDetail->DrawLatex(xOffset, 0.84, Form("Nch: %d-%d, |#Delta#it{#eta}| > 1.2", centBins[centralityIndex], centBins[centralityIndex + 1]));
  latexExp->DrawLatex(xOffset, 0.68, "CMS, pPb #sqrt{#it{s}}_{NN} = 5.02 TeV");
  latexDetail->DrawLatex(xOffset, 0.64, "Nch: 180-250, |#Delta#it{#eta}| > 1");

  lpAllCMS->AddEntry(hCMSK0, "K^{0}_{S}", "p");
  lpAllCMS->AddEntry(hCMSLambda, "#Lambda", "p");
  lpAllAlice->Draw();
  lpAllCMS->Draw();
  cAllCMS->SaveAs(Form("%s/v2_all_particles_centrality_%s_CMScomp.pdf", outdir, suffix.c_str()));
  cAllCMS->Write();

  // Close the output file
  outputFile->Close();
}
