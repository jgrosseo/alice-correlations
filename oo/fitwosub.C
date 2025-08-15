#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
#include <cstdio>

static TH1D *gpsimMass;
static TGraphErrors *gpvsb;

static const double massLimits[] = {0.47, 0.52};

void setTH1Style(TH1 *h, const char *xtitle, const char *ytitle, const int color = kBlack)
{
	h->SetMarkerStyle(kFullCircle);
	h->SetMarkerColor(color);
	h->SetLineColor(color);
	h->SetMarkerSize(1.0);
	h->SetLineWidth(2);
	h->GetXaxis()->SetTitle(xtitle);
	h->GetYaxis()->SetTitle(ytitle);
}

void setTGraphStyle(TGraph *g, const char *xtitle, const char *ytitle)
{
	g->SetMarkerStyle(kFullCircle);
	g->SetMarkerColor(kBlack);
	g->SetLineColor(kBlack);
	g->SetMarkerSize(1.0);
	g->SetLineWidth(2);
	g->GetXaxis()->SetTitle(xtitle);
	g->GetYaxis()->SetTitle(ytitle);
}

std::vector<std::string> GetAllRootFiles(const std::string &dirpath, const bool debug = false)
{
	std::vector<std::string> rootfiles;
	DIR *dir = opendir(dirpath.c_str());
	if (!dir)
	{
		std::cerr << "Cannot open directory: " << dirpath << std::endl;
		return rootfiles;
	}
	struct dirent *entry;
	while ((entry = readdir(dir)) != nullptr)
	{
		std::string fname = entry->d_name;
		if (fname.size() > 5 &&
			fname.substr(0, 4) == "dphi" &&
			fname.substr(fname.size() - 5) == ".root")
		{
			if (debug)
				std::cout << "Found file: " << fname << std::endl;
			rootfiles.push_back(dirpath + "/" + fname);
		}
	}
	closedir(dir);
	return rootfiles;
}

template <bool massFit, bool vnFit>
double SimFit(const double *par)
{
	// invariant mass fit
	double chi2_1 = 0.0, chi2_2 = 0.0;
	if (massFit)
	{
		/*
		int a = gpsimMass->GetXaxis()->FindBin(massLimits[0]);
		int b = gpsimMass->GetXaxis()->FindBin(massLimits[1]);
		for(int i = a; i < b; ++i){
		  double x = gpsimMass->GetBinCenter(i);
		  //double delta = gpsimMass->GetBinContent(i)
		  //-(par[0]+par[1]*x+par[2]*x*x+par[3]*TMath::Gaus(x,par[4],par[5]));
		  double delta = gpsimMass->GetBinContent(i)
		-(par[0]+par[1]*x+par[2]*x*x+par[3]*x*x*x+par[4]*TMath::Gaus(x,par[5],par[6])+par[7]*TMath::Gaus(x,par[8],par[9]));
		  double sigErr = gpsimMass->GetBinError(i);
		  chi2_1 += delta*delta/(sigErr*sigErr);
		  }*/
	}

	// vsb(m) = vs * (Ns/(Nb+Ns))(m) + vb(m) * (Nb/(Nb+Ns))(m), vb(m) = a*m+b
	if (vnFit)
	{
		for (uint i = 0; i < gpvsb->GetN(); ++i)
		{
			double x = gpvsb->GetX()[i];
			double Nbkg = par[0] + par[1] * x + par[2] * x * x + par[3] * x * x * x;
			// double Nsig = par[3]*TMath::Gaus(x,par[4],par[5]);
			double Nsig = par[4] * TMath::Gaus(x, par[5], par[6]) + par[7] * TMath::Gaus(x, par[8], par[9]);
			// double delta = gpvsb->GetY()[i]
			//-(par[6]*Nsig/(Nsig+Nbkg)+(par[7]*x+par[8])*Nbkg/(Nsig+Nbkg)); //vn^B modeled with linear dependency, according to CMS
			double delta = gpvsb->GetY()[i] - (par[10] * Nsig / (Nsig + Nbkg) + (par[11] * x + par[12]) * Nbkg / (Nsig + Nbkg)); // vn^B modeled with linear dependency, according to CMS
			double vsbErr = gpvsb->GetEY()[i];
			chi2_2 += delta * delta / (vsbErr * vsbErr);
		}
	}
	// chi2_2 /= double(gpvsb->GetN());
	return chi2_1 + chi2_2;
}

TF1 *extract(TH1D *phf1, const char *pnamepf, bool write = false)
{
	const TString fitName = Form("fit_%s", pnamepf);

	auto fit = new TF1(fitName, "[0]+[1]*(1.0 + 2*[2]*TMath::Cos(x) + 2*[3]*TMath::Cos(2*x) + 2*[4]*TMath::Cos(3*x) + 2*[5]*TMath::Cos(4*x) + 2*[6]*TMath::Cos(5*x))");
	fit->SetParNames("czyam", "c", "v1", "v2", "v3", "v4", "v5");

	// fit->FixParameter(5,0.0); // v4 := 0
	// fit->FixParameter(6,0.0); // v5 := 0
	fit->FixParameter(0, 0.0); // offset := 0
	
	//fit->SetParLimits(1, 0.0,1.0);
	fit->SetParLimits(2, -1.0,1.0);
	fit->SetParLimits(3, 0.0,1.0);
	fit->SetParLimits(4, 0.0,1.0);
	//fit->SetParLimits(5, 0.0,1.0);
	//fit->SetParLimits(6, 0.0,1.0);
	

	// new TCanvas;
	const char *option = "0SE"; //"0QSE";
	TFitResultPtr r = phf1->Fit(fit, option, "", -TMath::Pi() / 2.0, 3.0 / 2.0 * TMath::Pi());

	// fit->Print("V");

	if (write)
	{
		fit->SetRange(-TMath::Pi() / 2.0, 3.0 / 2.0 * TMath::Pi());
		fit->Write();
	}

	return fit;
}

TF1 *extractSubtraction(TH1D *phf1, const char *pnamepf, TF1 *subtract, bool write = false)
{
	const TString fitName = Form("fit%s", pnamepf);

	auto fit = new TF1(fitName, "[0]+[1]*(1 + 2*[2]*TMath::Cos(x)  + 2*[3] *TMath::Cos(2*x) + 2*[4] *TMath::Cos(3*x) + 2*[5] *TMath::Cos(4*x) + 2*[6]*TMath::Cos(5*x)) + "
								"[7]*([8]+[9]*(1 + 2*[10]*TMath::Cos(x) + 2*[11]*TMath::Cos(2*x) + 2*[12]*TMath::Cos(3*x) + 2*[13]*TMath::Cos(4*x) + 2*[14]*TMath::Cos(5*x)))");
	const char *parNames[] = {"czyam", "c", "v1", "v2", "v3", "v4", "v5", "scaling",
							  "czyam_lm", "c_lm", "v1_lm", "v2_lm", "v3_lm", "v4_lm", "v5_lm"};
	for (int i = 0; i < 15; i++)
		fit->SetParName(i, parNames[i]);
	/*
	fit->FixParameter(2,0.0); // v1 := 0
	fit->FixParameter(5,0.0); // v4 := 0
	fit->FixParameter(6,0.0); // v5 := 0*/
	fit->FixParameter(0, 0.0); // offset := 0
	fit->SetParameter(7, 1.6);

	// fit->FixParameter(7, 1);
	fit->SetParameter(7, 1);
	fit->SetParLimits(7, 0.9, 1.1);

	for (int i = 0; i < 7; i++)
		fit->FixParameter(i + 8, subtract->GetParameter(i));

	// new TCanvas;
	TString option = "0SE"; //"0SE";
	printf("Pass 1 (ignore the warnings):\n");
	phf1->Fit(fit, option + "Q", "", -TMath::Pi() / 2.0, 3.0 / 2.0 * TMath::Pi());
	printf("Pass 2 (check the output):\n");
	TFitResultPtr r = phf1->Fit(fit, option, "", -TMath::Pi() / 2.0, 3.0 / 2.0 * TMath::Pi());

	// fit->Print("V");

	if (write)
	{
		fit->SetRange(-TMath::Pi() / 2.0, 3.0 / 2.0 * TMath::Pi());
		fit->Write();
	}

	return fit;
}

void fit(const char *pinFileName = "infile.root",
		 int itrig = 0,
		 const char *outFileName = "out.root",
		 bool subtraction = false, 
		 double absDeltaEtaMin = 1.2, double absDeltaEtaMax = 1.8,
		 double absDeltaPhiMax = 1.3)
{
	TFile *pf = new TFile(pinFileName, "read");
	if (!pf || pf->IsZombie())
	{
		std::cerr << "Error opening input file: " << pinFileName << std::endl;
		return;
	}

	TTree *paxes = (TTree *)pf->Get("axes");
	UInt_t Ncent, Nmass, Ntrig, Nassoc;
	paxes->SetBranchAddress("Ncentrality", &Ncent);
	paxes->SetBranchAddress("Nmass", &Nmass);
	paxes->SetBranchAddress("NleadingPt", &Ntrig);
	paxes->SetBranchAddress("NassocPt", &Nassoc);
	paxes->GetEntry(0);
	float *pcentBins = new float[Ncent + 1];
	//double *pmassBins = new double[Nmass + 1];
	float *pmeanCents = new float[Ncent + 1];
	float *passocPt = new float[Nassoc + 1];
	float *ptrigPt = new float[Ntrig + 1];
	paxes->SetBranchAddress("centrality", pcentBins);
	//paxes->SetBranchAddress("mass", pmassBins);
	paxes->SetBranchAddress("meanCent", pmeanCents);
	paxes->SetBranchAddress("assocPt", passocPt);
	paxes->SetBranchAddress("leadingPt", ptrigPt);
	paxes->GetEntry(0);

	Double_t pmassBins[51] = {0.45,0.452,0.454,0.456,0.458,0.46,0.462,0.464,0.466,0.468,0.47,0.472,0.474,0.476,0.478,0.48,0.482,0.484,0.486,0.488,0.49,0.492,0.494,0.496,0.498,0.5,0.502,0.504,0.506,0.508,0.51,0.512,0.514,0.516,0.518,0.52,0.522,0.524,0.526,0.528,0.53,0.532,0.534,0.536,0.538,0.54,0.542,0.544,0.546,0.548,0.55};

	auto absDeltaEtaMinOrig = absDeltaEtaMin;
	auto absDeltaEtaMaxOrig = absDeltaEtaMax;
	absDeltaEtaMin += 0.001;
	absDeltaEtaMax -= 0.001;

	ROOT::Math::Minimizer *pmin = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
	pmin->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
	pmin->SetMaxIterations(10000);		 // for GSL
	pmin->SetTolerance(0.001);
	pmin->SetPrintLevel(1);

	auto v22NoMassTemplate = new TH3F("v22NoMassTemplate", ";pT,trig;pT,assoc;multiplicity,v22", int(Ntrig), ptrigPt, int(Nassoc), passocPt, int(Ncent), pcentBins);
	auto v22NoMassLMSubtraction = (TH3 *)v22NoMassTemplate->Clone("v22NoMassLMSubtraction");

	// Cherck if the output file already exists
	TFile *outFile = TFile::Open(outFileName, "UPDATE");
	if (!outFile || outFile->IsZombie())
	{
		std::cerr << "Error opening output file: " << outFileName << std::endl;
		return;
	}
	if (!outFile->cd())
	{
		std::cerr << "Error changing directory to output file: " << outFileName << std::endl;
		return;
	}
	outFile->cd();

	std::cout << "Output file opened successfully: " << outFileName << std::endl;
	
	uint iassoc = 0; // fixed

	double TempSub, TempSubErr;

	std::vector<TF1*> lowestFit;
	std::vector<TH2*> lowestHist;

	for (uint imult = 0; imult < Ncent; ++imult) // Loop over centrality
	{
		std::vector<TGraphErrors*> gParams;

		// Get invariant mass histogram
		gpsimMass = (TH1D *)pf->Get(Form("hmass_%u_%u", itrig, imult))->Clone(Form("hmass_%u_%u_%d", itrig, iassoc, imult));
		gpsimMass->SetDirectory(0);
		setTH1Style(gpsimMass, "M_{inv} (GeV/c^{2})", "Counts");
		gpsimMass->Write();
		TGraphErrors *pvsbMass = new TGraphErrors(Nmass); // mass dependent V2 Fourier coefficient

		for (int imass = 0; imass < int(Nmass); ++imass) // Loop over mass bins
		{
			const TString masspf = imass < 0 ? "_" : Form("_%d", imass);
			const TString namepf = Form("_%u_%u_%u%s", itrig, iassoc, imult, masspf.Data());
			const TString namepf2 = Form("%u%u%d%u", itrig, iassoc, imass, imult);
			TH2D *ph = (TH2D *)pf->Get(Form("dphi%s", namepf.Data())); // 2d histogram to be normalized to DeltaEta
			if (!ph)
			{
				printf("No histograms corresponding mult bin %u. (itrig=%u, iassoc=%u, imass=%d)\n", imult, itrig, iassoc, imass);
				continue;
			}
			std::cout << "itrig = " << itrig
						<< ", iassoc = " << iassoc
						<< ", imult = " << imult
						<< ", imass = " << imass
						<< ", masspf = '" << masspf << "'"
						<< ", namepf = '" << namepf << "'" << std::endl;

			int a = ph->GetYaxis()->FindBin(absDeltaEtaMin);
			int b = ph->GetYaxis()->FindBin(absDeltaEtaMax);
			TH1D *pp = ph->ProjectionX(Form("proj_dphi_Pos%s", namepf.Data()), a, b, "e"); // positive side DeltaPhi long-range
			setTH1Style(pp, "#Delta#varphi (rad)", Form("dN/d#Delta#varphi (%d < #Delta#eta < %d)", int(absDeltaEtaMin), int(absDeltaEtaMax)), kRed + 1);

			a = ph->GetYaxis()->FindBin(-absDeltaEtaMax);
			b = ph->GetYaxis()->FindBin(-absDeltaEtaMin);
			TH1D *pn = ph->ProjectionX(Form("proj_dphi_Neg%s", namepf.Data()), a, b, "e"); // negative side DeltaPhi long-range
			setTH1Style(pn, "#Delta#varphi (rad)", Form("dN/d#Delta#varphi (%d < #Delta#eta < %d)", int(-absDeltaEtaMax), int(-absDeltaEtaMin)), kAzure + 4);
			// pp->Write();
			// pn->Write();

			TH1D *phf1 = (TH1D *)pp->Clone(Form("projdphi_%u_%u_%d_%u", itrig, iassoc, imass, imult)); // combined DeltaPhi distribution (no CZYAM)
			setTH1Style(phf1, "#Delta#varphi (rad)", Form("dN/d#Delta#varphi (%d < #Delta#eta < %d)", int(absDeltaEtaMin), int(absDeltaEtaMax)));
			phf1->Add(pp, pn, 0.5, 0.5);
			phf1->Rebin(2);
			phf1->Scale(0.5);

			double normPhi = 2.0 * (absDeltaEtaMaxOrig - absDeltaEtaMinOrig);
			phf1->Scale(1.0 / normPhi); // normalize to long-range width (small delta-factor)
			// phf1->Write();

			a = ph->GetXaxis()->FindBin(-absDeltaPhiMax);
			b = ph->GetXaxis()->FindBin(absDeltaPhiMax);
			TH1D *pd = ph->ProjectionY(Form("projdeta_%u_%u_%d_%u", itrig, iassoc, imass, imult)); // near side DeltaEta projection
			setTH1Style(pd, "#Delta#eta", Form("dN/d#Delta#eta (%d < #Delta#varphi < %d)", int(-absDeltaPhiMax), int(absDeltaPhiMax)));
			double normEta = 2.0 * (absDeltaPhiMax);
			pd->Scale(1.0 / normEta);
			// pd->Write();

			double masscenter = (pmassBins[imass] + pmassBins[imass + 1]) / 2.0;

			TF1* fit = nullptr;

			if (imult == 0 || !subtraction) {
				fit = extract(phf1, namepf2.Data());

				if (subtraction) {
					lowestFit.push_back(fit);
					lowestHist.push_back(ph);
				}
			} else {
				fit = extractSubtraction(phf1, namepf.Data(), lowestFit[imass]);
				float scaling = fit->GetParameter(7);

				// draw a subtracted 2d histogram with the scale factor
				if (false) {
					auto clone = (TH2*) ph->Clone(Form("%s_sub", ph->GetName()));
					clone->Add(lowestHist[imass], -1.0 * scaling);
					clone->GetYaxis()->SetRangeUser(-1.59, 1.59);
					clone->Write();
				}					
			}

			// add all fit params to a graph
			for (int i=0; i<fit->GetNpar()+1+subtraction; i++) {
				if (gParams.size() <= i) {
					auto name = "NONAME";
					if (i < fit->GetNpar())
						name = fit->GetParName(i);
					else if (i == fit->GetNpar())
						name = "Chi2/ndf";
					else if (i == fit->GetNpar()+1)
						name = "v2_lmsub";
					
					auto graph = new TGraphErrors();
					graph->SetName(Form("%s_vs_mass_cent%d_pt%d", name, imult, itrig));
					setTGraphStyle(graph, "M_{inv} (GeV/c^{2})", name);
					gParams.push_back(graph);
				}
				auto graph = gParams.at(i);

				float value = 0;
				float error = 0;
				if (i < fit->GetNpar()) { // all parameters
					value = fit->GetParameter(i);
					error = fit->GetParError(i);
				} else if (i == fit->GetNpar()) { // Chi2/ndf
					value = fit->GetChisquare() / fit->GetNDF();
				} else if (i == fit->GetNpar()+1 && imult > 0) { // LM subtraction
					// convert to LM subtraction by rescaling
					float factor = fit->GetParameter(1) / (fit->GetParameter(1) + fit->GetParameter(7) * fit->GetParameter(9));
					value = fit->GetParameter(3) * factor;
					error = fit->GetParError(3) * factor; // TODO uncertainty propagation of param 1 and 9					
				}

				graph->SetPoint(imass, masscenter, value);
				graph->SetPointError(imass, 0.0, error); // no error in mass axis
			}
		}

		std::cout << "Writing v2 vs mass histograms to output file." << std::endl;
		//const int Nbins = Nmass + 1; // number of bins in the mass histogram
		const int Nbins = 49 + 1; // number of bins in the mass histogram

		auto hV2Mass = new TH1D(Form("hV2_vs_mass_cent%d_pt%d", imult,itrig),
								Form("v2 vs mass (centrality %d and pt %d)", imult, itrig),
										Nbins, pmassBins);
		TGraphErrors *g = gParams.at(3);
		
		if (!g) {
			std::cerr << "Warning: gv2Mass[" << imult << "][" << itrig << "] is null. Skipping." << std::endl;
			continue;
		}
		for (int i = 0; i < g->GetN(); ++i)
		{
			double x = g->GetX()[i];
			double y = g->GetY()[i];
			double ey = g->GetEY()[i];
			
			int bin = hV2Mass->FindBin(x);
			hV2Mass->SetBinContent(bin, y);
			hV2Mass->SetBinError(bin, ey);
		}
		// gv2Mass->Write();
		hV2Mass->Write();

		for (auto graph : gParams) {
			graph->Write();
		}		
	}

	outFile->Close();

	pf->Close();
	delete pf;
}

void fitwosub(bool subtraction = false, const char *inFilePath = "./k0s", const char *outFileName = "./massandv2.root", double absDeltaEtaMin = 1.2, double absDeltaEtaMax = 1.8, double absDeltaPhiMax = 1.3)
{
	std::cout << "Starting fitwosubt with input file path: " << inFilePath << std::endl;
	// Loading all the files called dphi*.root from the input directory
	std::vector<std::string> files = GetAllRootFiles(inFilePath, true);
	// Ordering files by name
	std::sort(files.begin(), files.end());
	std::cout << "Found " << files.size() << " files to process." << std::endl;
	// Check if output file already exists
	if (TFile::Open(outFileName))
	{
		std::cout << "WARNING: Output file already exists: " << outFileName << " - it will be overwritten!" << std::endl;
		TFile::Open(outFileName)->Close(); // Close the file to avoid issues
		remove(outFileName);			   // Remove the file
		std::cout << "Removed existing output file: " << outFileName << std::endl;
	}
	for (int iFile = 0; iFile < files.size(); ++iFile) // Uncomment the condition to process all files: for (int iFile = 0; iFile < (int)
	{
		const std::string &f = files[iFile];
		std::cout << "Processing file " << (iFile + 1) << " of " << files.size() << ": " << f << std::endl;
		fit(f.c_str(), iFile, outFileName, subtraction, absDeltaEtaMin, absDeltaEtaMax, absDeltaPhiMax);
	}
	std::cout << "Finished processing all files." << std::endl;
}
