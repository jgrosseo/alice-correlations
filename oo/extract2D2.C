#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <map>

Float_t gpTMin = 0.5;
Float_t gpTMax = 8;
Float_t gZVtxRange = -7;

gStyle->SetPalette(kRainBow);
gROOT->SetBatch(kTRUE);

void SetupRanges(CorrelationContainer *obj)
{
    obj->setEtaRange(0, 0);
    obj->setPtRange(gpTMin, gpTMax);
    if (gZVtxRange > 0)
        obj->setZVtxRange(-gZVtxRange + 0.01, gZVtxRange - 0.01);
}

struct BinningConfig
{
    std::vector<float> leadingPtArr;
    std::vector<float> assocPtArr;
    std::vector<float> centralityArr;
    std::vector<double> massArr;
    std::string outPrefix;
};

std::map<std::string, BinningConfig>
getBinningConfigs()
{
    std::map<std::string, BinningConfig> binmap;

    // Lambda
    binmap["Lambda"] = {
        {0.5, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 5, 8},
        {0.5, 5},
        {0, 10, 20}, // 50, 150, 300},
        {1.1, 1.102, 1.104, 1.106, 1.108, 1.11, 1.112, 1.114, 1.116, 1.118, 1.12, 1.122, 1.124, 1.126, 1.13},
        "lambda"};
    // Antiambda
    binmap["Antilambda"] = {
        {0.5, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 5, 8},
        {0.5, 5},
        {0, 25, 50, 150, 300},
        {1.1, 1.102, 1.104, 1.106, 1.108, 1.11, 1.112, 1.114, 1.116, 1.118, 1.12, 1.122, 1.124, 1.126, 1.13},
        "antilambda"};
    // K0s
    binmap["K0s"] = {
        {0.5, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 5, 8},
        {0.5, 8},
        //{25, 50, 150, 30}, // nch
        {0, 10, 20}, // v0m
        {0.466, 0.47, 0.474, 0.478, 0.482, 0.486, 0.488, 0.49, 0.492, 0.494, 0.496, 0.498, 0.5, 0.502, 0.504, 0.506, 0.51, 0.514, 0.518, 0.524},
        "k0s"};
    // phi
    binmap["Phi"] = {
        {0.5, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 5, 8},
        {0.5, 8},
        {70, 90},
        {0.996, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06},
        "phi"};

    binmap["Hadrons"] = {
        {0.5, 8.0},
        {0.5, 8.0},
        {0, 25}, // 50, 150, 300},
        {1.07, 1.072, 1.074, 1.076, 1.078, 1.08, 1.082, 1.084, 1.086, 1.088, 1.09, 1.092, 1.094, 1.096, 1.098, 1.1, 1.102, 1.104, 1.106, 1.108, 1.11, 1.112, 1.114, 1.116, 1.118, 1.12, 1.122, 1.124, 1.126, 1.128, 1.13, 1.132, 1.134, 1.136, 1.138, 1.14, 1.142, 1.144, 1.146, 1.148, 1.15, 1.152, 1.154, 1.156, 1.158, 1.16, 1.162, 1.164, 1.166, 1.168, 1.17},
        "hadrons"};
    return binmap;
}

// Plot a TH2 as a 3D surface on a given canvas pad, with options
void PlotTH2OnPad(
    TH1 *h,
    TPad *pad,
    double x1 = -999, double x2 = -999, // X axis range, -999 disables
    double y1 = -999, double y2 = -999, // Y axis range
    double z1 = -999, double z2 = -999, // Z axis range
    double phi = 30, double theta = 30  // Rotation: phi (azimuth), theta (elevation)
)
{
    if (!h || !pad)
        return;
    pad->cd();
    h->SetStats(0);
    h->SetLineWidth(1);
    h->SetLineColor(kBlack);
    h->SetLineStyle(8);

    // Set the grid empty
    pad->SetGridx(0);
    pad->SetGridy(0);

    // Set axis ranges if requested
    if (x1 != -999 || x2 != -999)
    {
        h->GetXaxis()->SetRangeUser(x1 == -999 ? h->GetXaxis()->GetXmin() : x1,
                                    x2 == -999 ? h->GetXaxis()->GetXmax() : x2);
    }
    if (y1 != -999 || y2 != -999)
    {
        h->GetYaxis()->SetRangeUser(y1 == -999 ? h->GetYaxis()->GetXmin() : y1,
                                    y2 == -999 ? h->GetYaxis()->GetXmax() : y2);
    }
    if (z1 != -999 || z2 != -999)
    {
        h->SetMinimum(z1 == -999 ? h->GetMinimum() : z1);
        h->SetMaximum(z2 == -999 ? h->GetMaximum() : z2);
    }

    h->GetXaxis()->SetTitle("#Delta#varphi (rad)");
    h->GetYaxis()->SetTitle("#Delta#eta");
    h->GetZaxis()->SetTitle("Counts");

    h->GetXaxis()->SetTitleOffset(1.4);
    h->GetYaxis()->SetTitleOffset(1.4);
    h->GetZaxis()->SetTitleOffset(1.4);
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetZaxis()->SetMaxDigits(3);
    h->GetXaxis()->SetDecimals();
    h->GetYaxis()->SetDecimals();
    h->GetZaxis()->SetDecimals();

    pad->SetTheta(theta);
    pad->SetPhi(phi);

    h->DrawCopy("SURF1");
}

// Unified GetSumOfRatios (with optional mass selection)
void GetSumOfRatiosUnified(
    CorrelationContainer *h, CorrelationContainer *hMixed, TH1 **hist,
    CorrelationContainer::CFStep step,
    Float_t centralityBegin, Float_t centralityEnd,
    Float_t ptBegin, Float_t ptEnd,
    Bool_t normalizePerTrigger = kTRUE,
    Bool_t useMass = kFALSE, double massBegin = 0.0, double massEnd = 0.0)
{
    h->setCentralityRange(centralityBegin, centralityEnd);
    hMixed->setCentralityRange(centralityBegin, centralityEnd);

    if (useMass)
    {
        std::cout << "********************* Using mass range: " << massBegin << " - " << massEnd << " *********************" << std::endl;
        h->getPairHist()->getTHn(6)->GetAxis(6)->SetRangeUser(massBegin, massEnd);
        hMixed->getPairHist()->getTHn(6)->GetAxis(6)->SetRangeUser(massBegin, massEnd);
        h->getTriggerHist()->getTHn(6)->GetAxis(3)->SetRangeUser(massBegin, massEnd);
        hMixed->getTriggerHist()->getTHn(6)->GetAxis(3)->SetRangeUser(massBegin, massEnd);
    }
    else
    {
        h->getPairHist()->getTHn(6)->GetAxis(6)->SetRangeUser(0.0, 0.0);
        hMixed->getPairHist()->getTHn(6)->GetAxis(6)->SetRangeUser(0.0, 0.0);
        h->getTriggerHist()->getTHn(6)->GetAxis(3)->SetRangeUser(0.0, 0.0);
        hMixed->getTriggerHist()->getTHn(6)->GetAxis(3)->SetRangeUser(0.0, 0.0);
    }

    *hist = h->getSumOfRatios(hMixed, step, ptBegin, ptEnd, normalizePerTrigger);

    TString str;
    str.Form("%.1f < p_{T,trig} < %.1f", ptBegin, ptEnd);
    TString str2;
    str2.Form("%.2f < p_{T,assoc} < %.2f", gpTMin - 0.01, gpTMax + 0.01);
    TString str3;
    if (useMass)
        str3.Form("%.4f < M < %.4f", massBegin, massEnd);
    else
        str3.Form("%.4f < M < %.4f", 0.0, 0.0);

    TString newTitle;
    newTitle.Form("%s - %s - %s - %.0f-%.0f", str.Data(), str2.Data(), str3.Data(), centralityBegin, centralityEnd);
    if ((*hist))
        (*hist)->SetTitle(newTitle);

    // Printout for debugging
    std::cout << "[GetSumOfRatiosUnified] Centrality: " << centralityBegin << " - " << centralityEnd
              << ", Pt: " << ptBegin << " - " << ptEnd
              << ", Assoc (global): " << gpTMin << " - " << gpTMax;
    if (useMass)
        std::cout << ", Mass: " << massBegin << " - " << massEnd;
    std::cout << std::endl;
}

// Main function - output file is per leading pT bin, user provides particle name
void extract2D2(const char *fileNamePbPb = "/Users/spolitan/cernbox/AnalysisResults_full_pass2_525561_lambda_ft0m.root", const char *outdir = "./lambda_hm_ft0m", const char *folder = "correlation-task", const char *particleName = "Lambda", bool qaplot = false)
{
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();

    auto binmap = getBinningConfigs();
    std::string particle(particleName);
    if (binmap.find(particle) == binmap.end())
    {
        std::cerr << "ERROR: Unknown particle type '" << particle << "'. Available: ";
        for (auto const &p : binmap)
            std::cerr << p.first << " ";
        std::cerr << std::endl;
        return;
    }
    const BinningConfig &cfg = binmap[particle];

    float leadingPtArrC[32], assocPtArrC[32], centralityArrC[32];
    double massArrC[64];
    int maxLeadingPt = cfg.leadingPtArr.size() - 1;
    int maxAssocPt = cfg.assocPtArr.size() - 1;
    int maxCentrality = cfg.centralityArr.size() - 1;
    int maxMass = cfg.massArr.size() - 1;
    for (int i = 0; i < maxLeadingPt + 1; ++i)
        leadingPtArrC[i] = cfg.leadingPtArr[i];
    for (int i = 0; i < maxAssocPt + 1; ++i)
        assocPtArrC[i] = cfg.assocPtArr[i];
    for (int i = 0; i < maxCentrality + 1; ++i)
        centralityArrC[i] = cfg.centralityArr[i];
    for (int i = 0; i < maxMass + 1; ++i)
        massArrC[i] = cfg.massArr[i];

    // Check the binning configuration
    auto *inputFile = TFile::Open(fileNamePbPb);
    // Axis 0: axis0 - #Delta#eta
    // Axis 1 : axis1 - p_{T}(GeV / c)
    // Axis 2 : axis2 - p_{T}(GeV / c)
    // Axis 3 : axis3 - multiplicity / centrality
    // Axis 4 : axis4 - #Delta #varphi(rad)
    // Axis 5 : axis5 - z - vtx(cm)    // Axis 6 : axis6 - m(GeV / c ^ 2)
    auto h = (CorrelationContainer *)inputFile->Get(Form("%s/sameEvent", folder));
    auto thn = h->getPairHist()->getTHn(6);
    if (!thn)
    {
        std::cerr << "FATAL: THnBase is null!\n";
        return;
    }
    else
    {
        int nAxes = thn->GetNdimensions();

        if (thn->GetAxis(2)->GetXmin() > cfg.leadingPtArr[0] || thn->GetAxis(2)->GetXmax() < cfg.leadingPtArr[maxLeadingPt])
        {
            std::cerr << "WARNING: Axis 1 range [" << thn->GetAxis(1)->GetXmin() << ", " << thn->GetAxis(1)->GetXmax()
                      << "] is outside the configured leading pT range [" << cfg.leadingPtArr[0] << ", "
                      << cfg.leadingPtArr[maxLeadingPt] << "]\n";
        }
        if (thn->GetAxis(1)->GetXmin() > cfg.assocPtArr[0] || thn->GetAxis(1)->GetXmax() < cfg.assocPtArr[maxAssocPt])
        {
            std::cerr << "WARNING: Axis 2 range [" << thn->GetAxis(2)->GetXmin() << ", " << thn->GetAxis(2)->GetXmax()
                      << "] is outside the configured assoc pT range [" << cfg.assocPtArr[0] << ", "
                      << cfg.assocPtArr[maxAssocPt] << "]\n";
        }
        if (thn->GetAxis(3)->GetXmin() > cfg.centralityArr[0] || thn->GetAxis(3)->GetXmax() < cfg.centralityArr[maxCentrality])
        {
            std::cerr << "WARNING: Axis 3 range [" << thn->GetAxis(3)->GetXmin() << ", " << thn->GetAxis(3)->GetXmax()
                      << "] is outside the configured centrality range [" << cfg.centralityArr[0] << ", "
                      << cfg.centralityArr[maxCentrality] << "]\n";
        }
        if (thn->GetAxis(6)->GetXmin() > cfg.massArr[0] || thn->GetAxis(6)->GetXmax() < cfg.massArr[maxMass])
        {
            std::cerr << "WARNING: Axis 6 range [" << thn->GetAxis(6)->GetXmin() << ", " << thn->GetAxis(6)->GetXmax()
                      << "] is outside the configured mass range [" << cfg.massArr[0] << ", "
                      << cfg.massArr[maxMass] << "]\n";
        }

        ///*
        if (qaplot)
        {
            TFile outfileQa("qa_k0_check.root", "RECREATE");
            auto hMapSame = thn->Projection(0, 4);
            // auto hMapMixed = thnMixed->Projection(0, 4);
            hMapSame->Write("same_deltaEtaDeltaPhi");
            // hMapMixed->Write("mixed_deltaEtaDeltaPhi");

            thn->GetAxis(1)->SetRangeUser(0., 1.0);
            thn->GetAxis(2)->SetRangeUser(0., 1.0);
            thn->GetAxis(3)->SetRangeUser(0, 100);
            // thnMixed->GetAxis(1)->SetRangeUser(0, 1.0);
            // thnMixed->GetAxis(2)->SetRangeUser(0, 1.0);
            // thnMixed->GetAxis(3)->SetRangeUser(0, 50);
            for (int i = 0; i < nAxes; ++i)
            {
                auto hAx = thn->Projection(i);
                hAx->SetLineColor(kBlue);
                hAx->Write(Form("same_axis%d", i));
                // auto hAxMixed = thnMixed->Projection(i);
                // hAxMixed->SetLineColor(kRed);
                // hAxMixed->Write(Form("mixed_axis%d", i));
            }

            outfileQa.Close();
            return;
        }
        // */
    }

    gStyle->SetOptStat(1111111);
    gSystem->Load("libO2PhysicsAnalysisCore.so");
    gSystem->Load("libO2PhysicsPWGCFCore.so");
    gSystem->Load("libO2PhysicsTwoPartCorrCore.so");

    CorrelationContainer::CFStep step = CorrelationContainer::kCFStepReconstructed;
    Bool_t normalizePerTrigger = kTRUE;

    Float_t *meanCent = new Float_t[maxCentrality];
    TH1F *proj = inputFile->Get<TH1F>(Form("%s/multiplicity", folder));
    for (int i = 0; i < maxCentrality; ++i)
    {
        proj->GetXaxis()->SetRangeUser(cfg.centralityArr[i] + 0.001, cfg.centralityArr[i + 1] - 0.001);
        meanCent[i] = proj->GetMean();
    }

    int total_iterations = maxLeadingPt * maxAssocPt * maxCentrality * maxMass;
    int iter_count = 0;

    for (int iLeadingPt = 0; iLeadingPt < maxLeadingPt; ++iLeadingPt)
    {
        auto tbin_start = clock::now();

        TString outFileName;
        outFileName.Form("%s/dphi_corr_%s_ptleading%d.root", outdir, cfg.outPrefix.c_str(), iLeadingPt);
        gSystem->Exec(Form("mkdir -p %s", outdir)); // Ensure output directory exists
        TFile file(outFileName, "RECREATE");

        std::cout << "\n\n\n####################### Starting leading pT bin " << iLeadingPt << "/" << maxLeadingPt << " for particle " << particle << " #######################\n";

        TTree axes("axes", "Axes final binning");
        UInt_t NleadingPt = maxLeadingPt;
        UInt_t NassocPt = maxAssocPt;
        UInt_t Ncentrality = maxCentrality;
        UInt_t Nmass = maxMass;
        std::cout << "NleadingPt: " << NleadingPt << ", NassocPt: " << NassocPt
                  << ", Ncentrality: " << Ncentrality << ", Nmass: " << Nmass << std::endl;
        std::cout << "Mass bins: ";
        for (int i = 0; i < maxMass + 1; ++i)
        {
            std::cout << cfg.massArr[i] << " ";
        }
        std::cout << std::endl;
        axes.Branch("NleadingPt", &NleadingPt, "NleadingPt/i");
        axes.Branch("leadingPt", leadingPtArrC, Form("leadingPt[%d]/F", maxLeadingPt + 1));
        axes.Branch("NassocPt", &NassocPt, "NassocPt/i");
        axes.Branch("assocPt", assocPtArrC, Form("assocPt[%d]/F", maxAssocPt + 1));
        axes.Branch("Ncentrality", &Ncentrality, "Ncentrality/i");
        axes.Branch("centrality", centralityArrC, Form("centrality[%d]/F", maxCentrality + 1));
        axes.Branch("Nmass", &Nmass, "Nmass/i");
        axes.Branch("mass", massArrC, Form("mass[%d]/D", maxMass + 1));
        axes.Branch("meanCent", meanCent, Form("meanCent[%u]/F", maxCentrality));
        axes.Fill();
        axes.Write();

        for (int iAssocPt = 0; iAssocPt < maxAssocPt; ++iAssocPt) // Loop over assoc pT bins
        {
            if (cfg.assocPtArr[iAssocPt] >= cfg.leadingPtArr[iLeadingPt + 1]) // Skip if assoc pT is not valid for leading pT
                continue;

            gpTMin = cfg.assocPtArr[iAssocPt] + 0.01;
            gpTMax = cfg.assocPtArr[iAssocPt + 1] - 0.01;

            for (int mult = 0; mult < maxCentrality; ++mult) // Loop over centrality
            {
                TCanvas *cMass = new TCanvas(Form("cmass_%d_%d_%d", iLeadingPt, iAssocPt, mult),
                                             "All mass cuts", 1200, 1000);
                if (qaplot)
                {
                    int nPads = maxMass;
                    int nRows = std::ceil(std::sqrt(nPads));
                    int nCols = std::ceil(double(nPads) / nRows);
                    cMass->Divide(nCols, nRows, 0.001, 0.001);
                }

                // Set up the correlation containers
                auto hcopy = (CorrelationContainer *)inputFile->Get(Form("%s/sameEvent", folder));
                auto hMixedcopy = (CorrelationContainer *)inputFile->Get(Form("%s/mixedEvent", folder));
                SetupRanges(hcopy);
                SetupRanges(hMixedcopy);

                if (iAssocPt == 0) // First assoc pT bin, get integrated inv. mass
                {
                    auto thnOriginal = hcopy->getPairHist()->getTHn(6);
                    thnOriginal->GetAxis(3)->SetRangeUser(cfg.centralityArr[mult] + 0.01, cfg.centralityArr[mult + 1] - 0.01);
                    thnOriginal->GetAxis(2)->SetRangeUser(cfg.leadingPtArr[iLeadingPt] + 0.01, cfg.leadingPtArr[iLeadingPt + 1] - 0.01);
                    auto hMixedcopy2 = hMixedcopy->getPairHist()->getTHn(6);
                    hMixedcopy2->GetAxis(3)->SetRangeUser(cfg.centralityArr[mult] + 0.01, cfg.centralityArr[mult + 1] - 0.01);
                    hMixedcopy2->GetAxis(2)->SetRangeUser(cfg.leadingPtArr[iLeadingPt] + 0.01, cfg.leadingPtArr[iLeadingPt + 1] - 0.01);
                    auto hmass2 = thnOriginal->Projection(6);
                    auto hMixedMass = hMixedcopy2->Projection(6);
                    if (qaplot)
                    {
                        auto hMapSame = thnOriginal->Projection(0, 4);
                        auto thnMixed = hMixedcopy->getPairHist()->getTHn(6);
                        thnMixed->GetAxis(3)->SetRangeUser(cfg.centralityArr[mult] + 0.01, cfg.centralityArr[mult + 1] - 0.01);
                        thnMixed->GetAxis(2)->SetRangeUser(cfg.leadingPtArr[iLeadingPt] + 0.01, cfg.leadingPtArr[iLeadingPt + 1] - 0.01);
                        auto hMapMixed = thnMixed->Projection(0, 4);
                        hMapSame->SetName(Form("hDeltaEtaDeltaPhiMapSameMass_%d_%d", iLeadingPt, mult));
                        hMapMixed->SetName(Form("hDeltaEtaDeltaPhiMapMixedMass_%d_%d", iLeadingPt, mult));
                        hMapSame->SetTitle("DeltaEta vs DeltaPhi Map Same Event");
                        hMapMixed->SetTitle("DeltaEta vs DeltaPhi Map Mixed Event");
                        hMapSame->SetStats(0);
                        hMapMixed->SetStats(0);
                        file.cd();
                        hMapSame->Write();
                        hMapMixed->Write();
                        TCanvas cMap(Form("cMap_%d_%d", iLeadingPt, mult), "DeltaEta vs DeltaPhi Map", 800, 600);
                        cMap.Divide(2, 1);
                        cMap.cd(1);
                        hMapSame->Draw("SURF1");
                        cMap.cd(2);
                        hMapMixed->Draw("SURF1SAME");
                        cMap.SaveAs(Form("%s/deltaEtaDeltaPhiMap_%d_%d.pdf", outdir, iLeadingPt, mult));
                    }
                    if (hmass2)
                    {
                        file.cd();
                        hmass2->SetName(Form("hmass_%d_%d", iLeadingPt, mult));
                        hmass2->Write();
                    }
                    if (hMixedMass)
                    {
                        file.cd();
                        hMixedMass->SetName(Form("hmassMixed_%d_%d", iLeadingPt, mult));
                        hMixedMass->Write();
                    }
                }

                // Get the sum of ratios for the current assoc pT and centrality
                /*
                TH1 *hist0 = 0;
                GetSumOfRatiosUnified(
                    hcopy, hMixedcopy, &hist0, step,
                    cfg.centralityArr[mult] + 0.01, cfg.centralityArr[mult + 1] - 0.01,
                    cfg.leadingPtArr[iLeadingPt] + 0.01, cfg.leadingPtArr[iLeadingPt + 1] - 0.01,
                    normalizePerTrigger, kFALSE);
                    if (hist0)
                    {
                        file.cd();
                        hist0->SetName(Form("dphi_%d_%d_%d", iLeadingPt, iAssocPt, mult));
                        hist0->Write();
                    }
                    delete hist0;
                */
                delete hcopy;
                delete hMixedcopy;

                for (int mass = 0; mass < maxMass; ++mass) // Loop over mass bins
                {
                    auto hcopy2 = (CorrelationContainer *)inputFile->Get(Form("%s/sameEvent", folder));
                    auto hMixedcopy2 = (CorrelationContainer *)inputFile->Get(Form("%s/mixedEvent", folder));
                    SetupRanges(hcopy2);
                    SetupRanges(hMixedcopy2);

                    TH1 *hist1 = 0;
                    GetSumOfRatiosUnified(
                        hcopy2, hMixedcopy2, &hist1, step,
                        cfg.centralityArr[mult] + 0.01, cfg.centralityArr[mult + 1] - 0.01,
                        cfg.leadingPtArr[iLeadingPt] + 0.01, cfg.leadingPtArr[iLeadingPt + 1] - 0.01,
                        normalizePerTrigger, kTRUE, cfg.massArr[mass] + 0.0001, cfg.massArr[mass + 1] - 0.0001);
                    if (hist1)
                    {
                        file.cd();
                        hist1->SetName(Form("dphi_%d_%d_%d_%d", iLeadingPt, iAssocPt, mult, mass));
                        hist1->Write();

                        // Draw on canvas pad
                        if (qaplot)
                        {
                            cMass->cd(mass + 1);
                            PlotTH2OnPad(hist1, (TPad *)cMass->cd(mass + 1), -1.8, -4.4, -1.4, 1.4, hist1->GetMaximum() * 0.4, hist1->GetMaximum() * 1.01,
                                         140, 40);
                            // Save the canvas with all mass plots
                            TString pdfname = Form("%s/dphi_%s_ptleading%d_assoc%d_cent%d_allmass.pdf",
                                                   outdir, cfg.outPrefix.c_str(), iLeadingPt, iAssocPt, mult);
                            cMass->SaveAs(pdfname);
                        }
                    }
                    delete hist1;
                    delete hcopy2;
                    delete hMixedcopy2;
                }
                file.cd();
                if (qaplot)
                {
                    cMass->SetName(Form("cmass_%d_%d_%d", iLeadingPt, iAssocPt, mult));
                    cMass->Write();
                    delete cMass;
                }
            } // mult
        } // iAssocPt

        auto tbin_end = clock::now();
        auto elapsed_bin = std::chrono::duration_cast<std::chrono::minutes>(tbin_end - tbin_start).count();
        std::cout << "File " << outFileName << " written! Leading pT bin time: " << elapsed_bin << " min" << std::endl;
        file.Close();
    } // iLeadingPt

    auto t1 = clock::now();
    auto total_elapsed = std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count();
    std::cout << "All done for particle " << particle << "! Total elapsed time: " << total_elapsed << " min" << std::endl;
    delete[] meanCent;
}