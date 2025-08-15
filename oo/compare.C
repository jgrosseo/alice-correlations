void compare(const char* histogram = "v2_vs_mass", const char* fileName1 = "massandv2.root", const char* fileName2 = "massandv2_sub.root")
{
    auto file1 = TFile::Open(fileName1);
    auto file2 = TFile::Open(fileName2);

    auto canvas = new TCanvas("c1", "c1", 1500, 1000);
    canvas->Divide(5, 3);
    int currentPad = 1;

    int imult = 0;
    while (imult != -1) {
        int itrig = 0;
        while (true) {
            Printf("%d %d", imult, itrig);

            TString tmp;
            tmp.Form("%s_cent%d_pt%d", histogram, imult, itrig);

            auto graph1 = (TGraphErrors*) file1->Get(tmp);
            if (!graph1) {
                if (itrig == 0)
                    imult = -2; // to entirely break
                break;
            }

            canvas->cd(currentPad);
            graph1->Draw("AP");

            auto graph2 = (TGraphErrors*) file2->Get(tmp);
            graph2->SetLineColor(2);
            graph2->SetMarkerColor(2);
            graph2->SetMarkerStyle(4);
            graph2->Draw("PSAME");

            currentPad++;
            itrig++;
        }
        imult++;
    }
}