// =============================================================================
// analyze_cherenkov.C (ULTIMATE EDITION - PUBLICATION QUALITY)
// =============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <iostream>
#include <cmath>

// ── Kích thước ranh giới chuẩn từ Geant4 ──
const double B_SD1  = 50.0;
const double B_SD2  = 53.0; 
const double B_SIPM = 53.5; 

void StyleHist(TH1* h, int color, int width=2, int style=1) { 
    h->SetLineColor(color); h->SetLineWidth(width); h->SetLineStyle(style); h->SetStats(0); 
}

// Hàm trải phẳng 3D thành 2D Origami
void FillOrigami(TH2D* H[6], double x, double y, double z, double bound) {
    double ax = std::abs(x), ay = std::abs(y), az = std::abs(z);
    double dY = std::abs(ay - bound), dZ = std::abs(az - bound), dX = std::abs(ax - bound);
    double dMin = std::min({dY, dZ, dX});
    if (dMin > 0.5) return;
    if (dMin == dY) { if(y>0) H[0]->Fill(x,z); else H[1]->Fill(x,z); }
    else if (dMin == dZ) { if(z>0) H[2]->Fill(x,y); else H[3]->Fill(x,-y); }
    else { if(x>0) H[5]->Fill(y,z); else H[4]->Fill(-y,z); }
}

void DrawOrigami(TPad* pad, TH2D* H[6], const char* title) {
    pad->cd(); 
    TPaveText *pt = new TPaveText(0.1, 0.9, 0.9, 0.98, "NDC");
    pt->AddText(title); pt->SetFillColor(0); pt->SetBorderSize(0); pt->Draw();
    
    TPad* grid = new TPad("grid", "grid", 0, 0, 1, 0.9); grid->Draw(); grid->cd();
    grid->Divide(3, 4, 0.001, 0.001);
    int kPads[] = {2, 8, 5, 11, 7, 9}; // Top, Bot, Frt, Bck, Lft, Rgt
    for(int i=0; i<6; ++i) { grid->cd(kPads[i]); H[i]->Draw("COLZ"); }
}

void analyze_cherenkov() {
    TString filepath = "/home/huyhoang/software/Geant4_hoang/Cerenkov/build/cherenkov_water_data.root";
    TFile* fIn = TFile::Open(filepath, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cout << "Error: Không tìm thấy file " << filepath << std::endl; return;
    }

    TTree* tTruth = (TTree*)fIn->Get("Truth");
    TTree* tIn    = (TTree*)fIn->Get("InnerData");
    TTree* tOut   = (TTree*)fIn->Get("OuterData");
    TTree* tSum   = (TTree*)fIn->Get("EventSummary");
    TTree* tDet   = (TTree*)fIn->Get("DetectedSignal");

    gStyle->SetOptStat(0); gStyle->SetPalette(kRainBow); 

    // =========================================================================
    // CANVAS 1: PRIMARY PHYSICS & CALORIMETRY
    // =========================================================================
    TCanvas* c1 = new TCanvas("c1", "1. Primary Physics", 1500, 500); c1->Divide(3,1);
    
    c1->cd(1);
    TH1D* hPriE = new TH1D("hPriE", "Muon Initial Energy; Energy (GeV); Events", 50, 0, 10);
    StyleHist(hPriE, kBlue); hPriE->SetFillColorAlpha(kBlue, 0.3);
    tTruth->Project("hPriE", "Initial_E_MeV/1000.0", "PDG_Code==13 || PDG_Code==-13");
    hPriE->Draw("HIST");

    c1->cd(2);
    TH1D* hEdep = new TH1D("hEdep", "Energy Deposit in Water; dE/dx (MeV); Events", 50, 0, 100);
    StyleHist(hEdep, kRed); hEdep->SetFillColorAlpha(kRed, 0.3);
    tTruth->Project("hEdep", "Edep_Water_MeV", "Edep_Water_MeV > 0");
    hEdep->Draw("HIST");

    c1->cd(3);
    TH2D* hCalo = new TH2D("hCalo", "Generated Photons vs Energy Deposit; Edep in Water (MeV); Generated Photons", 50, 0, 50, 50, 0, 20000);
    tTruth->Project("hCalo", "OptPhotonsGen:Edep_Water_MeV", "Edep_Water_MeV > 0");
    hCalo->Draw("COLZ");

    // =========================================================================
    // CANVAS 2: OPTICAL PIPELINE (EFFICIENCY)
    // =========================================================================
    TCanvas* c2 = new TCanvas("c2", "2. Photon Pipeline", 800, 600); c2->SetGridy();
    
    double sumGen = 0, sumSD1 = 0, sumSD2 = 0, sumDet = 0;
    int gw, ga, w2a, a2c, det;
    tSum->SetBranchAddress("GenWater", &gw); tSum->SetBranchAddress("GenAcrylic", &ga);
    tSum->SetBranchAddress("W2A_Unique", &w2a); tSum->SetBranchAddress("A2C_Unique", &a2c);
    tSum->SetBranchAddress("Detected_Total", &det);
    for(int i=0; i<tSum->GetEntries(); ++i) {
        tSum->GetEntry(i);
        sumGen += (gw + ga); sumSD1 += w2a; sumSD2 += a2c; sumDet += det;
    }

    TH1D* hPipe = new TH1D("hPipe", "Optical Efficiency Pipeline; Stage; Unique Photons", 4, 0, 4);
    hPipe->GetXaxis()->SetBinLabel(1, "1. Generated"); hPipe->SetBinContent(1, sumGen);
    hPipe->GetXaxis()->SetBinLabel(2, "2. Pass SD1 (Water->Acrylic)"); hPipe->SetBinContent(2, sumSD1);
    hPipe->GetXaxis()->SetBinLabel(3, "3. Pass SD2 (Acrylic->Gel)"); hPipe->SetBinContent(3, sumSD2);
    hPipe->GetXaxis()->SetBinLabel(4, "4. Detected (SiPM)"); hPipe->SetBinContent(4, sumDet);
    hPipe->SetFillColor(kAzure+2); hPipe->SetBarWidth(0.6); hPipe->SetBarOffset(0.2);
    hPipe->SetMarkerSize(1.5); hPipe->Draw("b text");

    // =========================================================================
    // CANVAS 3: SPECTROSCOPY (WAVELENGTH EVOLUTION)
    // =========================================================================
    TCanvas* c3 = new TCanvas("c3", "3. Wavelength Spectroscopy", 1000, 600);
    
    TH1D* hwSD1 = new TH1D("hwSD1", "Spectral Evolution; Wavelength (nm); Photons", 80, 200, 800);
    TH1D* hwSD2 = new TH1D("hwSD2", "Spectral Evolution; Wavelength (nm); Photons", 80, 200, 800);
    TH1D* hwDet = new TH1D("hwDet", "Spectral Evolution; Wavelength (nm); Photons", 80, 200, 800);
    StyleHist(hwSD1, kBlue, 2); StyleHist(hwSD2, kRed, 2); StyleHist(hwDet, kGreen+2, 2);
    hwSD1->SetFillColorAlpha(kBlue, 0.1); hwSD2->SetFillColorAlpha(kRed, 0.2); hwDet->SetFillColorAlpha(kGreen+2, 0.3);

    tIn->Project("hwSD1", "Wavelength_nm");
    tOut->Project("hwSD2", "Wavelength_nm");
    tDet->Project("hwDet", "Wavelength_nm");

    hwSD1->Draw("HIST"); hwSD2->Draw("HIST SAME"); hwDet->Draw("HIST SAME");
    TLegend* lSpec = new TLegend(0.65, 0.65, 0.88, 0.88);
    lSpec->AddEntry(hwSD1, "At SD1 (Before Acrylic)", "f");
    lSpec->AddEntry(hwSD2, "At SD2 (Survived UV Cut-off)", "f");
    lSpec->AddEntry(hwDet, "Detected (Filtered by SiPM PDE)", "f");
    lSpec->Draw();

    // =========================================================================
    // CANVAS 4: ANGULAR DISTRIBUTIONS (TIR ANALYSIS)
    // =========================================================================
    TCanvas* c4 = new TCanvas("c4", "4. Angular Distributions", 1000, 500); c4->Divide(2,1);
    
    c4->cd(1);
    TH1D* hAng1 = new TH1D("hAng1", "Incident Angle at SD1 (Water->Acrylic); Angle (deg); Photons", 90, 0, 90);
    StyleHist(hAng1, kBlue); hAng1->SetFillColorAlpha(kBlue, 0.3);
    tIn->Project("hAng1", "AngleToPlane_deg");
    hAng1->Draw("HIST");

    c4->cd(2);
    TH1D* hAng2 = new TH1D("hAng2", "Incident Angle at SD2 (Acrylic->Gel); Angle (deg); Photons", 90, 0, 90);
    StyleHist(hAng2, kRed); hAng2->SetFillColorAlpha(kRed, 0.3);
    tOut->Project("hAng2", "AngleToPlane_deg");
    hAng2->Draw("HIST");

    // =========================================================================
    // CANVAS 5: TIMING & TRAPPING (BOUNCE)
    // =========================================================================
    TCanvas* c5 = new TCanvas("c5", "5. Timing & Bounces", 1500, 500); c5->Divide(3,1);

    c5->cd(1);
    TH1D* hTime = new TH1D("hTime", "Photon Arrival Time at SiPM; Time (ns); Photons", 100, 0, 10);
    StyleHist(hTime, kOrange+2); hTime->SetFillColorAlpha(kOrange+1, 0.4);
    tDet->Project("hTime", "Time_ns");
    hTime->Draw("HIST");

    c5->cd(2);
    TH1D* hBounce = new TH1D("hBounce", "Trapping: Number of Bounces before Detection; Bounces; Photons", 20, 0, 20);
    StyleHist(hBounce, kMagenta+2); hBounce->SetFillColorAlpha(kMagenta, 0.4);
    tDet->Project("hBounce", "BounceCount");
    hBounce->Draw("HIST");

    c5->cd(3);
    TH2D* hTimeBounce = new TH2D("hTimeBounce", "Correlation: Arrival Time vs Bounces; Arrival Time (ns); Bounce Count", 50, 0, 10, 20, 0, 20);
    tDet->Project("hTimeBounce", "BounceCount:Time_ns");
    hTimeBounce->Draw("COLZ");

    // =========================================================================
    // CANVAS 6: 2D ORIGAMI HIT MAPS
    // =========================================================================
    TCanvas* c6 = new TCanvas("c6", "6. Spatial Hit Maps", 1800, 600); c6->Divide(3,1);
    
    TH2D *hMapSD1[6], *hMapSD2[6], *hMapDet[6];
    const char* fN[] = {"Top","Bot","Frt","Bck","Lft","Rgt"};
    for(int i=0; i<6; ++i) {
        hMapSD1[i] = new TH2D(Form("M1_%s",fN[i]), Form("SD1 %s",fN[i]), 40,-B_SD1,B_SD1, 40,-B_SD1,B_SD1); hMapSD1[i]->SetStats(0);
        hMapSD2[i] = new TH2D(Form("M2_%s",fN[i]), Form("SD2 %s",fN[i]), 40,-B_SD2,B_SD2, 40,-B_SD2,B_SD2); hMapSD2[i]->SetStats(0);
        hMapDet[i] = new TH2D(Form("M3_%s",fN[i]), Form("SiPM %s",fN[i]), 40,-B_SIPM,B_SIPM, 40,-B_SIPM,B_SIPM); hMapDet[i]->SetStats(0);
    }

    double hx, hy, hz;
    tIn->SetBranchAddress("HitX_mm", &hx); tIn->SetBranchAddress("HitY_mm", &hy); tIn->SetBranchAddress("HitZ_mm", &hz);
    for(int i=0; i<tIn->GetEntries(); ++i) { tIn->GetEntry(i); FillOrigami(hMapSD1, hx, hy, hz, B_SD1); }
    
    tOut->SetBranchAddress("HitX_mm", &hx); tOut->SetBranchAddress("HitY_mm", &hy); tOut->SetBranchAddress("HitZ_mm", &hz);
    for(int i=0; i<tOut->GetEntries(); ++i) { tOut->GetEntry(i); FillOrigami(hMapSD2, hx, hy, hz, B_SD2); }
    
    tDet->SetBranchAddress("HitX_mm", &hx); tDet->SetBranchAddress("HitY_mm", &hy); tDet->SetBranchAddress("HitZ_mm", &hz);
    for(int i=0; i<tDet->GetEntries(); ++i) { tDet->GetEntry(i); FillOrigami(hMapDet, hx, hy, hz, B_SIPM); }

    DrawOrigami((TPad*)c6->cd(1), hMapSD1, "Hits at SD1 (Water->Acrylic)");
    DrawOrigami((TPad*)c6->cd(2), hMapSD2, "Hits at SD2 (Acrylic->Gel)");
    DrawOrigami((TPad*)c6->cd(3), hMapDet, "Absorbed by SiPM");

    std::cout << "\n[SUCCESS] Đã tạo thành công 6 Canvas phân tích toàn diện!\n";
}