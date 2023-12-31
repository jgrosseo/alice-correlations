Float_t gpTMin = 0.51;
Float_t gpTMax = 49.99;
Float_t gZVtxRange = -1;

void SetupRanges(CorrelationContainer* obj)
{
  obj->setEtaRange(0, 0);
  obj->setPtRange(gpTMin, gpTMax);
  if (gZVtxRange > 0)
    obj->setZVtxRange(-gZVtxRange+0.01, gZVtxRange-0.01);
}

void GetSumOfRatios(CorrelationContainer* h, CorrelationContainer* hMixed, TH1** hist, CorrelationContainer::CFStep step, Float_t centralityBegin, Float_t centralityEnd, Float_t ptBegin, Float_t ptEnd, Bool_t normalizePerTrigger = kTRUE)
{
  Printf("GetSumOfRatios | step %d | %.1f-%.1f%% | %.1f - %.1f GeV/c | %.1f - %.1f GeV/c", step, centralityBegin, centralityEnd, gpTMin, gpTMax, ptBegin, ptEnd);

  h->setCentralityRange(centralityBegin, centralityEnd);
  hMixed->setCentralityRange(centralityBegin, centralityEnd);
  *hist  = h->getSumOfRatios(hMixed, step, ptBegin, ptEnd, normalizePerTrigger);
  
  TString str;
  str.Form("%.1f < p_{T,trig} < %.1f", ptBegin - 0.01, ptEnd + 0.01);
  
  TString str2;
  str2.Form("%.2f < p_{T,assoc} < %.2f", gpTMin - 0.01, gpTMax + 0.01);
    
  TString newTitle;
  newTitle.Form("%s - %s - %.0f-%.0f", str.Data(), str2.Data(), centralityBegin, centralityEnd);
  if ((*hist))
    (*hist)->SetTitle(newTitle);
}

void extract2D(const char* fileNamePbPb = "AnalysisResults.root", const char* outputFile = "dphi_corr.root", const char* folder = "correlation-task")
{
  gStyle->SetOptStat(1111111);

  CorrelationContainer::CFStep step = CorrelationContainer::kCFStepReconstructed;
  Bool_t normalizePerTrigger = kTRUE; // don't do this if histograms are to be merged with other periods later -> Use MergeDPhiFiles below

  auto file = TFile::Open(outputFile, "RECREATE");
  file->Close();
  
  Int_t maxLeadingPt = 2;
  Float_t leadingPtArr[] = { 1.0, 2.0, 3.0 };

  Int_t maxAssocPt = 2;
  Float_t assocPtArr[] =  { 1.0, 2.0, 3.0 };
  
  Int_t maxCentrality = 4;
  Float_t centralityArr[] = { 0, 4, 8, 16, 80 };
  
  auto inputFile = TFile::Open(fileNamePbPb);
  auto h      = (CorrelationContainer*) inputFile->Get(Form("%s/sameEvent", folder));
  auto hMixed = (CorrelationContainer*) inputFile->Get(Form("%s/mixedEvent", folder));
  
//   hMixed->getTriggerHist()->FillParent();
//   hMixed->getTriggerHist()->getTHn(6)->Projection(2)->Draw();
//   return;
  
  auto eventHist = h->getEventCount();
  Printf("Events with centrality: %f", eventHist->Integral(eventHist->GetXaxis()->FindBin(6.), eventHist->GetXaxis()->FindBin(6.), eventHist->GetYaxis()->FindBin(0.), eventHist->GetYaxis()->FindBin(99.9)));
  Printf("Events: %f", eventHist->ProjectionX()->Integral(eventHist->GetXaxis()->FindBin(6.), eventHist->GetXaxis()->FindBin(6.)));
//   eventHist->ProjectionX()->Draw();
//   return;

  /*h->setPtRange(0.51, 1.99);
  h->setCentralityRange(0, 49.9);

  TH2* sameTwoD = (TH2*) h->getPerTriggerYield(step, 2.01, 3.99);
  new TCanvas;
  sameTwoD->Draw("SURF1");
  
//   return;

  TH2* mixedTwoD = (TH2*) hMixed->getPerTriggerYield(step, 2.01, 3.99);
  new TCanvas;
  mixedTwoD->Draw("SURF1");

  auto ratio = (TH2*) sameTwoD->Clone("ratio");
  ratio->Divide(mixedTwoD);

  auto c3 = new TCanvas;
  ratio->GetYaxis()->SetRangeUser(-1.59, 1.59);
  ratio->Draw("SURF1");
  c3->SaveAs("pertriggeryield.png");*/


  // new TCanvas;
  // ((TH3*) inputFile->Get("correlation-task/yields"))->GetXaxis()->SetRangeUser(0.01, 49.99);
  // ((TH3*) inputFile->Get("correlation-task/yields"))->Project3D("zy")->Draw("COLZ");

  //return;
  
  for (Int_t i=0; i<maxLeadingPt; i++)
  {
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      if (assocPtArr[j] >= leadingPtArr[i+1])
        continue;

      gpTMin = assocPtArr[j] + 0.01;
      gpTMax = assocPtArr[j+1] - 0.01;
        
      SetupRanges(h);
      SetupRanges(hMixed);
  
      for (Int_t mult=0; mult<maxCentrality; mult++)
      {
        TH1* hist1 = 0;
        GetSumOfRatios(h, hMixed, &hist1, step, centralityArr[mult], centralityArr[mult+1] - 1.0, leadingPtArr[i] + 0.01, leadingPtArr[i+1] - 0.01, normalizePerTrigger); 

        file = TFile::Open(outputFile, "UPDATE");
      
        if (hist1)
        {
          hist1->SetName(Form("dphi_%d_%d_%d", i, j, mult));
          hist1->Write();
        }
      
        file->Close();

        delete hist1;
      }
    }
    
    TH1* triggers = h->getTriggersAsFunctionOfMultiplicity(step, leadingPtArr[i] + 0.01, leadingPtArr[i+1] - 0.01);
    triggers->SetName(Form("triggers_%d", i));
    TString str;
    str.Form("%.1f < p_{T,trig} < %.1f", leadingPtArr[i], leadingPtArr[i+1]);
    triggers->SetTitle(str);

    file = TFile::Open(outputFile, "UPDATE");
    triggers->Write();
    file->Close();    
  }
    
  delete h;
  delete hMixed;
}

