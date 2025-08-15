void drawSubtracted(const char* fileName = "massandv2_sub.root")
{
    auto file = TFile::Open(fileName);

    auto canvas = new TCanvas("c1", "c1", 1500, 1000);
    canvas->Divide(5, 3);
    int currentPad = 1;

    int imult = 1;
    while (imult != -1) {
        int itrig = 0;
        while (true) {
            Printf("%d %d", imult, itrig);

            TString tmp;
            tmp.Form("dphi_%d_0_%d_25_sub", itrig, imult);

            auto hist = (TH2*) file->Get(tmp);
            if (!hist) {
                if (itrig == 0)
                    imult = -2; // to entirely break
                break;
            }

            canvas->cd(currentPad);
            hist->Draw("SURF1");
            hist->GetYaxis()->SetRangeUser(-1.39, 1.39);

            currentPad++;
            itrig++;
        }
        imult++;
    }
}