static const double PT_EDGES[] = {0.5, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 5, 8};
static const int N_PT_BINS = int(sizeof(PT_EDGES) / sizeof(double)) - 1;

/*static const double MASS_EDGES[] = {0.996,1.004,1.01,1.012,1.014,1.016,1.018,1.02,1.022,1.024,1.026,1.028,1.03,1.036,1.042,1.048,1.054,1.06};
 */

static const double MASS_EDGES[] = {0.466, 0.47, 0.474, 0.478, 0.482, 0.486, 0.488, 0.49, 0.492, 0.494, 0.496, 0.498, 0.5, 0.502, 0.504, 0.506, 0.51, 0.514, 0.518, 0.524};

static const int N_MASS_BINS = int(sizeof(MASS_EDGES) / sizeof(double)) - 1;

static const int MULT_LIST[] = {25, 50, 150, 300}; // labels only
static const int N_MULT = int(sizeof(MULT_LIST) / sizeof(int)) - 1;
static const int LM_INDEX = 0; // low-mult template index

static const char *FILE_LIST_LM[] = {
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading0.root",
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading1.root",
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading2.root",
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading3.root",
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading4.root",
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading5.root",
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading6.root",
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading7.root",
    "./k0_lm_largeretacheck/dphi_corr_k0s_ptleading8.root"};
static const char *FILE_LIST_HM[] = {
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading0.root",
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading1.root",
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading2.root",
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading3.root",
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading4.root",
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading5.root",
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading6.root",
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading7.root",
    "./k0_hm_largeretacheck/dphi_corr_k0s_ptleading8.root"};

static_assert(N_PT_BINS == int(sizeof(FILE_LIST_HM) / sizeof(FILE_LIST_HM[0])),
              "Number of pT bins must match number of files in FILE_LIST");

static const double ABS_ETA_MAX = 0.8;                                                                    // project only |Δη| > ABS_ETA_MIN
static const char *OUT_FILE = Form("./out_SRjetpeak_WOsubtraction_deltaeta%.1f_ratio.root", ABS_ETA_MAX); // output root file

// ===================== Helpers =====================
static TString RangeLabel(double a, double b, int prec = 3)
{
    TString fmt = TString::Format("%%.%df–%%.%df", prec, prec);
    return TString::Format(fmt.Data(), a, b);
}

void setHistoStyle(TH1 *h, int color, int marker)
{
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(1.);
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->SetStats(kFALSE);
}

static TH1 *FindInvMass(TFile *fin, int ptIndex, int multIndex, TString &matchedName)
{
    if (!fin)
        return nullptr;
    const char *patterns[] = {"hmass_%d_%d", "hmass_pt%d_mult%d", "invMass_%d_%d", "invMass_pt%d_mult%d"};
    for (auto pat : patterns)
    {
        matchedName.Form(pat, ptIndex, multIndex);
        if (auto h = dynamic_cast<TH1 *>(fin->Get(matchedName)))
            return h;
    }
    matchedName = "";
    return nullptr;
}

// ---- Δη -> Δφ projection using inverse-variance weights (w = 1/σ^2) ----
// y(phi) = sum_i w_i y_i / sum_i w_i,   σ(phi) = sqrt(1 / sum_i w_i)
static TH1D *ProjectPhiWeighted(const TH2 *h2, const char *name, double absEtaMin)
{
    if (!h2)
        return nullptr;
    const TAxis *ax = h2->GetXaxis(); // Δφ
    const TAxis *ay = h2->GetYaxis(); // Δη
    TH1D *out = new TH1D(name, h2->GetTitle(), ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
    out->Sumw2();

    for (int j = 1; j <= ax->GetNbins(); ++j)
    {
        const double dphi = ax->GetBinCenter(j);
        if (!(dphi > -1.0 && dphi < 1.0))
        {
            std::cout << "Skip dphi: " << dphi << " for near-side only" << std::endl;
            continue; // near-side only
        }
        double sw = 0.0;  // sum of weights
        double swy = 0.0; // sum of w*y
        for (int i = 1; i <= ay->GetNbins(); ++i)
        {
            const double deta = ay->GetBinCenter(i);
            if (std::abs(deta) > absEtaMin)
            {
                std::cout << "Skip deta: " << deta << " for absEtaMin: " << absEtaMin << std::endl;
                continue; // short-range only
            }
            const double y = h2->GetBinContent(j, i);
            const double e = h2->GetBinError(j, i);
            if (!(e > 0.0) || !std::isfinite(e))
                continue;
            const double w = 1.0 / (e * e); // stat. error^2 as weight (inverse-variance)
            sw += w;
            swy += w * y;
        }
        if (sw > 0.0)
        {
            out->SetBinContent(j, swy / sw);
            out->SetBinError(j, std::sqrt(1.0 / sw));
        }
        else
        {
            out->SetBinContent(j, 0.0);
            out->SetBinError(j, 0.0);
        }
    }
    return out;
}

void compareSRjetpeak()
{

    std::vector<int> colors = {
        TColor::GetColorTransparent(kRed + 1, 0.55), TColor::GetColorTransparent(kAzure + 4, 0.55), TColor::GetColorTransparent(kGreen + 2, 0.55), TColor::GetColorTransparent(kOrange + 1, 0.55)};

    TFile *outFile = TFile::Open(OUT_FILE, "RECREATE");

    // th1 vs inv mass in different pt bins in different mult bins for the sigma
    TH1D *hsigma_LM[N_PT_BINS];
    TH1D *hsigma_HM[N_PT_BINS];

    TH1D *halpha_LM[N_PT_BINS];
    TH1D *halpha_HM[N_PT_BINS];
    TH1D *hbeta_LM[N_PT_BINS];
    TH1D *hbeta_HM[N_PT_BINS];
    TH1D *hchi2_LM[N_PT_BINS];
    TH1D *hchi2_HM[N_PT_BINS];

    TCanvas *c3 = new TCanvas(Form("c_sigma_pt%d", N_PT_BINS), "c_sigma", 900, 600);
    c3->Divide(3, 3);
    TCanvas *calpha_vsmass2 = new TCanvas("calpha_SRjetpeak", "calpha_SRjetpeak", 900, 600);
    calpha_vsmass2->Divide(3, 3);
    TCanvas *cbeta_vsmass2 = new TCanvas("cbeta_SRjetpeak", "cbeta_SRjetpeak", 900, 600);
    cbeta_vsmass2->Divide(3, 3);
    TCanvas *cchi2_vsmass2 = new TCanvas("cchi2_SRjetpeak", "cchi2_SRjetpeak", 900, 600);
    cchi2_vsmass2->Divide(3, 3);

    TCanvas *c_vsmass2[N_PT_BINS];
    TCanvas *cratio_vsmass[N_PT_BINS];
    for (int i = 0; i < N_PT_BINS; ++i)
    {
        c_vsmass2[i] = new TCanvas(Form("c_SRjetpeak_pt%d", i), Form("c_SRjetpeak_pt%d", i), 1200, 800);
        c_vsmass2[i]->Divide(5, 4);
        cratio_vsmass[i] = new TCanvas(Form("cratio_SRjetpeak_pt%d", i), Form("cratio_SRjetpeak_pt%d", i), 1200, 800);
        cratio_vsmass[i]->Divide(5, 4);
    }

    for (int i = 0; i < N_PT_BINS; ++i)
    {

        // TCanvas *c_vsmass = new TCanvas(Form("c_SRjetpeak_pt%d", N_PT_BINS), "c_SRjetpeak", 1200, 800);
        // c_vsmass->Divide(5, 5);
        //  load LM and HM file
        TFile *fLM = TFile::Open(FILE_LIST_LM[i], "READ");
        TFile *fHM = TFile::Open(FILE_LIST_HM[i], "READ");
        if (!fLM || fLM->IsZombie())
        {
            std::cerr << "Skip " << FILE_LIST_LM[i] << "\n";
            continue;
        }
        if (!fHM || fHM->IsZombie())
        {
            std::cerr << "Skip " << FILE_LIST_HM[i] << "\n";
            continue;
        }

        // Loop over mass bins
        hsigma_LM[i] = new TH1D(Form("hsigma_LM_pt%d", i), "LM sigma vs mass bin", N_MASS_BINS, MASS_EDGES);
        hsigma_HM[i] = new TH1D(Form("hsigma_HM_pt%d", i), "HM sigma vs mass bin", N_MASS_BINS, MASS_EDGES);

        halpha_LM[i] = new TH1D(Form("halpha_LM_pt%d", i), "LM alpha vs mass bin", N_MASS_BINS, MASS_EDGES);
        halpha_HM[i] = new TH1D(Form("halpha_HM_pt%d", i), "HM alpha vs mass bin", N_MASS_BINS, MASS_EDGES);
        hbeta_LM[i] = new TH1D(Form("hbeta_LM_pt%d", i), "LM beta vs mass bin", N_MASS_BINS, MASS_EDGES);
        hbeta_HM[i] = new TH1D(Form("hbeta_HM_pt%d", i), "HM beta vs mass bin", N_MASS_BINS, MASS_EDGES);
        hchi2_LM[i] = new TH1D(Form("hchi2_LM_pt%d", i), "LM chi2 vs mass bin", N_MASS_BINS, MASS_EDGES);
        hchi2_HM[i] = new TH1D(Form("hchi2_HM_pt%d", i), "HM chi2 vs mass bin", N_MASS_BINS, MASS_EDGES);
        for (int jm = 0; jm < N_MASS_BINS; ++jm)
        {
            // load LM map
            TString inName_LM = TString::Format("dphi_%d_%d_%d_%d", i, 0, LM_INDEX, jm); // assoc index assumed 0
            TH2 *h2_LM = (TH2 *)fLM->Get(inName_LM);
            if (!h2_LM)
                continue;

            // load HM map
            TString inName_HM = TString::Format("dphi_%d_%d_%d_%d", i, 0, 2, jm); // assoc index assumed 0
            TH2 *h2_HM = (TH2 *)fHM->Get(inName_HM);
            if (!h2_HM)
                continue;

            auto h1_HM_SR = h2_HM->ProjectionY(Form("h1_HM_SR_mass%d_pt%d", jm, i), h2_HM->GetXaxis()->FindBin(-1.3), h2_HM->GetXaxis()->FindBin(1.3));
            auto h1_LM_SR = h2_LM->ProjectionY(Form("h1_LM_SR_mass%d_pt%d", jm, i), h2_LM->GetXaxis()->FindBin(-1.3), h2_LM->GetXaxis()->FindBin(1.3));
            h1_HM_SR->SetDirectory(0);
            h1_LM_SR->SetDirectory(0);

            auto h1_HM_LM_ratio = (TH1D *)h1_HM_SR->Clone(Form("h1_HM_LM_ratio_mass%d_pt%d", jm, i));
            h1_HM_LM_ratio->Scale(h1_LM_SR->Integral() / h1_HM_SR->Integral());
            h1_HM_LM_ratio->Divide(h1_LM_SR);
            outFile->cd();
            cratio_vsmass[i]->cd(jm + 1);
            h1_HM_LM_ratio->GetYaxis()->SetTitle("HM / LM");
            h1_HM_LM_ratio->GetXaxis()->SetTitle("#Delta#eta");
            h1_HM_LM_ratio->SetTitle(Form("pT %3.2f-%3.2f GeV/c, mass %3.3f-%3.3f GeV/c^{2}",
                                          PT_EDGES[i], PT_EDGES[i + 1],
                                          MASS_EDGES[jm], MASS_EDGES[jm + 1]));
            // h1_HM_LM_ratio->GetYaxis()->SetRangeUser(h1_HM_LM_ratio->GetMinimum() * 0.8, h1_HM_LM_ratio->GetMaximum() * 1.2);
            h1_HM_LM_ratio->GetXaxis()->SetRangeUser(-1.3, 1.3);
            setHistoStyle(h1_HM_LM_ratio, colors[0], 21);
            setHistoStyle(h1_HM_SR, colors[0], 20);
            setHistoStyle(h1_LM_SR, colors[1], 24);
            // h1_HM_LM_ratio->Draw("E Same");
            h1_LM_SR->DrawNormalized("E Same");
            h1_HM_SR->DrawNormalized("E SAME");
            h1_HM_LM_ratio->Write();

            // fit generalised gaussian to get sigma
            TF1 *generalised_gaus_lm = new TF1(Form("fgaus_LM_mass%d_pt%d", jm, i), "[0] + ([1]/(2*[2]*TMath::Gamma(1/[3])))*TMath::Exp(-TMath::Power(TMath::Abs(x/[2]),[3]))", -1.0, 1.0);
            generalised_gaus_lm->SetParameters(0.0, h1_LM_SR->GetMaximum(), 0.2, 2.0);
            generalised_gaus_lm->SetParLimits(2, 0.05, 1.0);
            generalised_gaus_lm->SetParLimits(3, 1.0, 10.0);
            generalised_gaus_lm->SetLineColor(kAzure + 4);
            TF1 *generalised_gaus_hm = new TF1(Form("fgaus_HM_mass%d_pt%d", jm, i), "[0] + ([1]/(2*[2]*TMath::Gamma(1/[3])))*TMath::Exp(-TMath::Power(TMath::Abs(x/[2]),[3]))", -1.0, 1.0);
            generalised_gaus_hm->SetParameters(0.0, h1_HM_SR->GetMaximum(), 0.2, 2.0);
            generalised_gaus_hm->SetParLimits(2, 0.05, 1.0);
            generalised_gaus_hm->SetParLimits(3, 1.0, 10.0);

            // normalisation (optional)
            //  h1_LM_SR->Scale(1.0 / h1_LM_SR->Integral());
            //  h1_HM_SR->Scale(1.0 / h1_HM_SR->Integral());
            h1_LM_SR->Fit(generalised_gaus_lm, "RQM0");
            double sigma_LM = TMath::Sqrt(generalised_gaus_lm->GetParameter(2) * generalised_gaus_lm->GetParameter(2) * TMath::Gamma(3.0 / generalised_gaus_lm->GetParameter(3)) / TMath::Gamma(1.0 / generalised_gaus_lm->GetParameter(3)));
            h1_HM_SR->Fit(generalised_gaus_hm, "RQM0");
            double sigma_HM = TMath::Sqrt(generalised_gaus_hm->GetParameter(2) * generalised_gaus_hm->GetParameter(2) * TMath::Gamma(3.0 / generalised_gaus_hm->GetParameter(3)) / TMath::Gamma(1.0 / generalised_gaus_hm->GetParameter(3)));
            hsigma_LM[i]->SetBinContent(jm + 1, sigma_LM);
            hsigma_LM[i]->SetBinError(jm + 1, 1.e-9); // h1_LM_SR->GetRMS());
            hsigma_HM[i]->SetBinContent(jm + 1, sigma_HM);
            hsigma_HM[i]->SetBinError(jm + 1, 1.e-9); // h1_HM_SR->GetRMS());
            halpha_LM[i]->SetBinContent(jm + 1, generalised_gaus_lm->GetParameter(3));
            halpha_LM[i]->SetBinError(jm + 1, generalised_gaus_lm->GetParError(3));
            halpha_HM[i]->SetBinContent(jm + 1, generalised_gaus_hm->GetParameter(3));
            halpha_HM[i]->SetBinError(jm + 1, generalised_gaus_hm->GetParError(3));
            hbeta_LM[i]->SetBinContent(jm + 1, generalised_gaus_lm->GetParameter(2));
            hbeta_LM[i]->SetBinError(jm + 1, generalised_gaus_lm->GetParError(2));
            hbeta_HM[i]->SetBinContent(jm + 1, generalised_gaus_hm->GetParameter(2));
            hbeta_HM[i]->SetBinError(jm + 1, generalised_gaus_hm->GetParError(2));
            hchi2_LM[i]->SetBinContent(jm + 1, generalised_gaus_lm->GetChisquare() / generalised_gaus_lm->GetNDF());
            hchi2_HM[i]->SetBinContent(jm + 1, generalised_gaus_hm->GetChisquare() / generalised_gaus_hm->GetNDF());

            h1_HM_SR->GetYaxis()->SetTitle("1/N_{trig} dN/d#Delta#phi (|#Delta#phi|<1.3)");
            h1_HM_SR->GetXaxis()->SetTitle("#Delta#eta");
            h1_HM_SR->SetTitle(Form("pT %3.2f-%3.2f GeV/c, mass %3.3f-%3.3f GeV/c^{2}",
                                    PT_EDGES[i], PT_EDGES[i + 1],
                                    MASS_EDGES[jm], MASS_EDGES[jm + 1]));
            h1_LM_SR->GetYaxis()->SetTitle("1/N_{trig} dN/d#Delta#phi (|#Delta#phi|<1.3)");
            h1_LM_SR->GetXaxis()->SetTitle("#Delta#eta");
            h1_LM_SR->SetTitle(Form("pT %3.2f-%3.2f GeV/c, mass %3.3f-%3.3f GeV/c^{2}",
                                    PT_EDGES[i], PT_EDGES[i + 1],
                                    MASS_EDGES[jm], MASS_EDGES[jm + 1]));
            h1_HM_SR->GetYaxis()->SetRangeUser(0.0, TMath::Max(h1_HM_SR->GetMaximum(), h1_LM_SR->GetMaximum()) * 1.5);
            h1_HM_SR->GetXaxis()->SetRangeUser(-1.3, 1.3);
            setHistoStyle(h1_HM_SR, colors[0], 20);
            setHistoStyle(h1_LM_SR, colors[1], 24);

            c_vsmass2[i]->cd(jm + 1);
            h1_HM_SR->Draw("E Same");
            h1_LM_SR->Draw("E SAME");
            generalised_gaus_hm->Draw("SAME");
            generalised_gaus_lm->Draw("SAME");
            if (jm == 0)
            {
                auto legend = new TLegend(0.6, 0.6, 0.88, 0.88);
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                legend->AddEntry(h1_HM_SR, "HM (150-300)", "lep");
                legend->AddEntry(h1_LM_SR, "LM (0-25)", "lep");
                legend->AddEntry(generalised_gaus_hm, Form("HM #sigma=%.3f", sigma_HM), "l");
                legend->AddEntry(generalised_gaus_lm, Form("HM Chi2/NDF=%.2f", generalised_gaus_hm->GetChisquare() / generalised_gaus_hm->GetNDF()), "");
                legend->AddEntry(generalised_gaus_lm, Form("LM #sigma=%.3f", sigma_LM), "l");
                legend->AddEntry(generalised_gaus_lm, Form("LM Chi2/NDF=%.2f", generalised_gaus_lm->GetChisquare() / generalised_gaus_lm->GetNDF()), "");
                legend->Draw();
            }
            else
            {
                auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                legend->AddEntry(generalised_gaus_hm, Form("HM #sigma=%.3f", sigma_HM), "l");
                legend->AddEntry(generalised_gaus_lm, Form("HM Chi2/NDF=%.2f", generalised_gaus_hm->GetChisquare() / generalised_gaus_hm->GetNDF()), "");
                legend->AddEntry(generalised_gaus_lm, Form("LM #sigma=%.3f", sigma_LM), "l");
                legend->AddEntry(generalised_gaus_lm, Form("LM Chi2/NDF=%.2f", generalised_gaus_lm->GetChisquare() / generalised_gaus_lm->GetNDF()), "");
                legend->Draw();
            }
            c_vsmass2[i]->cd(jm + 1)->Update();

            outFile->cd();
            h2_LM->Write();
            h2_HM->Write();
            h1_LM_SR->Write();
            h1_HM_SR->Write();
            generalised_gaus_lm->Write();
            generalised_gaus_hm->Write();
        } // mass bin

        c_vsmass2[i]->SaveAs(Form("k0_srjetpeak_deltaeta%.1f_pt%d.pdf", ABS_ETA_MAX, i));
        cratio_vsmass[i]->SaveAs(Form("k0_srjetpeak_ratio_deltaeta%.1f_pt%d.pdf", ABS_ETA_MAX, i));

        // sigma vs mass
        outFile->cd();
        hsigma_LM[i]->Write();
        hsigma_HM[i]->Write();
        c3->cd();
        c3->cd(i + 1);
        c3->SetGridy();
        hsigma_LM[i]->GetYaxis()->SetTitle("#sigma_{#Delta#eta}");
        hsigma_LM[i]->GetXaxis()->SetTitle("#it{M}(#pi#pi) (GeV/c^{2})");
        hsigma_LM[i]->GetYaxis()->SetRangeUser(0, max(hsigma_LM[i]->GetMaximum() * 1.2, hsigma_HM[i]->GetMaximum()) * 1.2);
        hsigma_LM[i]->SetTitle(Form("pT %3.2f-%3.2f GeV/c",
                                    PT_EDGES[i], PT_EDGES[i + 1]));
        setHistoStyle(hsigma_LM[i], colors[1], 24);
        setHistoStyle(hsigma_HM[i], colors[0], 20);
        hsigma_LM[i]->Draw("E");
        hsigma_HM[i]->Draw("E SAME");
        // fit with pol0 outside mass peak region
        TF1 *fpol0_LM = new TF1(Form("fpol0_LM_pt%d_sb_left", i), "pol0", MASS_EDGES[0], 0.486);
        TF1 *fpol0_LM_right = new TF1(Form("fpol0_LM_pt%d_sb_right", i), "pol0", 0.506, MASS_EDGES[N_MASS_BINS]);
        hsigma_LM[i]->Fit(fpol0_LM, "RQM0");
        hsigma_LM[i]->Fit(fpol0_LM_right, "RQM0+");
        TF1 *fpol0_HM = new TF1(Form("fpol0_HM_pt%d_sb", i), "pol0", MASS_EDGES[0], 0.486);
        TF1 *fpol0_HM_right = new TF1(Form("fpol0_HM_pt%d_sb_right", i), "pol0", 0.506, MASS_EDGES[N_MASS_BINS]);
        hsigma_HM[i]->Fit(fpol0_HM, "RQM0");
        hsigma_HM[i]->Fit(fpol0_HM_right, "RQM0+");
        // fit peak region only
        TF1 *fpol0_LM_peak = new TF1(Form("fpol0_LM_pt%d_peak", i), "pol0", 0.486, 0.506);
        TF1 *fpol0_HM_peak = new TF1(Form("fpol0_HM_pt%d_peak", i), "pol0", 0.486, 0.506);
    
        hsigma_LM[i]->Fit(fpol0_LM_peak, "RQM0");
        hsigma_HM[i]->Fit(fpol0_HM_peak, "RQM0");

        fpol0_LM->SetLineColor(colors[1]);
        fpol0_LM->SetLineStyle(2);
        fpol0_HM->SetLineColor(colors[0]);
        fpol0_HM->SetLineStyle(2);
        fpol0_LM_peak->SetLineColor(colors[1]);
        fpol0_LM_peak->SetLineStyle(1);
        fpol0_HM_peak->SetLineColor(colors[0]);
        fpol0_HM_peak->SetLineStyle(1);
        fpol0_LM_right->SetLineColor(colors[1]);
        fpol0_LM_right->SetLineStyle(2);
        fpol0_HM_right->SetLineColor(colors[0]);
        fpol0_HM_right->SetLineStyle(2);
        fpol0_LM_peak->Draw("SAME");
        fpol0_HM_peak->Draw("SAME");
        fpol0_HM->Draw("SAME");
        fpol0_LM->Draw("SAME");
        fpol0_HM_right->Draw("SAME");
        fpol0_LM_right->Draw("SAME");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.15, 0.85, Form("LM #sigma=%.3f", fpol0_LM->GetParameter(0)));
        latex.DrawLatex(0.15, 0.80, Form("HM #sigma=%.3f", fpol0_HM->GetParameter(0)));
        latex.DrawLatex(0.4, 0.85, Form("LM peak #sigma=%.3f", fpol0_LM_peak->GetParameter(0)));
        latex.DrawLatex(0.4, 0.80, Form("HM peak #sigma=%.3f", fpol0_HM_peak->GetParameter(0)));
        latex.DrawLatex(0.65, 0.85, Form("LM peak #sigma=%.3f", fpol0_LM_right->GetParameter(0)));
        latex.DrawLatex(0.65, 0.80, Form("HM peak #sigma=%.3f", fpol0_HM_right->GetParameter(0)));

        if (i == 0)
        {
            auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(hsigma_HM[i], "HM (150-300)", "lep");
            legend->AddEntry(hsigma_LM[i], "LM (0-25)", "lep");
            legend->Draw();
        }

        // alpha vs mass
        outFile->cd();
        halpha_LM[i]->Write();
        halpha_HM[i]->Write();
        calpha_vsmass2->cd();
        calpha_vsmass2->cd(i + 1);
        calpha_vsmass2->SetGridy();
        halpha_LM[i]->GetYaxis()->SetTitle("#alpha");
        halpha_LM[i]->GetXaxis()->SetTitle("#it{M}(#pi#pi) (GeV/c^{2})");
        halpha_LM[i]->GetYaxis()->SetRangeUser(0, 10.0);
        halpha_LM[i]->SetTitle(Form("pT %3.2f-%3.2f GeV/c",
                                    PT_EDGES[i], PT_EDGES[i + 1]));
        setHistoStyle(halpha_LM[i], colors[1], 24);
        setHistoStyle(halpha_HM[i], colors[0], 20);
        halpha_LM[i]->Draw("E");
        halpha_HM[i]->Draw("E SAME");
        if (i == 0)
        {
            auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(halpha_HM[i], "HM (150-300)", "lep");
            legend->AddEntry(halpha_LM[i], "LM (0-25)", "lep");
            legend->Draw();
        }

        // beta vs mass
        outFile->cd();
        hbeta_LM[i]->Write();
        hbeta_HM[i]->Write();
        cbeta_vsmass2->cd();
        cbeta_vsmass2->cd(i + 1);
        cbeta_vsmass2->SetGridy();
        hbeta_LM[i]->GetYaxis()->SetTitle("#beta");
        hbeta_LM[i]->GetXaxis()->SetTitle("#it{M}(#pi#pi) (GeV/c^{2})");
        hbeta_LM[i]->GetYaxis()->SetRangeUser(0, 1.0);
        hbeta_LM[i]->SetTitle(Form("pT %3.2f-%3.2f GeV/c",
                                   PT_EDGES[i], PT_EDGES[i + 1]));
        setHistoStyle(hbeta_LM[i], colors[1], 24);
        setHistoStyle(hbeta_HM[i], colors[0], 20);
        hbeta_LM[i]->Draw("E");
        hbeta_HM[i]->Draw("E SAME");
        if (i == 0)
        {
            auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(hbeta_HM[i], "HM (150-300)", "lep");
            legend->AddEntry(hbeta_LM[i], "LM (0-25)", "lep");
            legend->Draw();
        }

        // chi2 vs mass
        outFile->cd();
        hchi2_LM[i]->Write();
        hchi2_HM[i]->Write();
        cchi2_vsmass2->cd();
        cchi2_vsmass2->cd(i + 1);
        cchi2_vsmass2->SetGridy();
        hchi2_LM[i]->GetYaxis()->SetTitle("#chi^{2}/NDF");
        hchi2_LM[i]->GetXaxis()->SetTitle("#it{M}(#pi#pi) (GeV/c^{2})");
        hchi2_LM[i]->GetYaxis()->SetRangeUser(0, max(hchi2_LM[i]->GetMaximum(), hchi2_HM[i]->GetMaximum()) * 1.2);
        hchi2_LM[i]->SetTitle(Form("pT %3.2f-%3.2f GeV/c",
                                   PT_EDGES[i], PT_EDGES[i + 1]));
        setHistoStyle(hchi2_LM[i], colors[1], 24);
        setHistoStyle(hchi2_HM[i], colors[0], 20);
        hchi2_LM[i]->Draw("E");
        hchi2_HM[i]->Draw("E SAME");
        if (i == 0)
        {
            auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(hchi2_HM[i], "HM (150-300)", "lep");
            legend->AddEntry(hchi2_LM[i], "LM (0-25)", "lep");
            legend->Draw();
        }

        // c->Write();
        // c2->Write();
        // fLM->Close();
        // fHM->Close();
    } // pt bin
    c3->SaveAs(Form("k0_sigma_vs_mass_deltaeta%.1f.pdf", ABS_ETA_MAX));
    calpha_vsmass2->SaveAs(Form("k0_alpha_vs_mass_deltaeta%.1f.pdf", ABS_ETA_MAX));
    cbeta_vsmass2->SaveAs(Form("k0_beta_vs_mass_deltaeta%.1f.pdf", ABS_ETA_MAX));
    cchi2_vsmass2->SaveAs(Form("k0_chi2_vs_mass_deltaeta%.1f.pdf", ABS_ETA_MAX));
    outFile->cd();
    c3->Write();

    outFile->Close();
}