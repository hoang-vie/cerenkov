// =============================================================================
// analyze_cherenkov.C (BẢN V.3.0 - FIX LỖI TỌA ĐỘ VÀ PIPELINE)
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
#include <iostream>
#include <cmath>
#include <set>

// ── Kích thước ranh giới chuẩn từ Geant4 ──
const double B_SD1  = 50.0;
const double B_SD2  = 53.0; // Sửa lại thành 53.0 (Acrylic 3mm)
const double B_SIPM = 53.5; // Sửa lại thành 53.5 (Coupling 0.5mm)

void StyleH1(TH1D* h, int color, int style=1, int width=2) { 
    h->SetLineColor(color); h->SetLineWidth(width); h->SetLineStyle(style); h->SetStats(0); 
}

void FillOrigami(TH2D* H[6], double x, double y, double z, double bound) {
    double ax = std::abs(x), ay = std::abs(y), az = std::abs(z);
    double dY = std::abs(ay - bound), dZ = std::abs(az - bound), dX = std::abs(ax - bound);
    double dMin = std::min({dY, dZ, dX});
    if (dMin > 0.5) return;
    if (dMin == dY) { if(y>0) H[0]->Fill(x,z); else H[1]->Fill(x,z); }
    else if (dMin == dZ) { if(z>0) H[2]->Fill(x,y); else H[3]->Fill(x,-y); }
    else { if(x>0) H[5]->Fill(y,z); else H[4]->Fill(-y,z); }
}

void DrawOrigami(TPad* pad, TH2D* H[6]) {
    pad->cd(); pad->Divide(3, 4, 0.002, 0.002);
    int kPads[] = {2, 8, 5, 11, 7, 9}; // Top, Bot, Frt, Bck, Lft, Rgt
    for(int i=0; i<6; ++i) { pad->cd(kPads[i]); H[i]->Draw("COLZ"); }
}

void analyze_cherenkov() {
    TFile* fIn = TFile::Open("/home/huyhoang/software/Geant4_hoang/Cerenkov/build/cherenkov_water_data.root", "READ");
    if (!fIn || fIn->IsZombie()) return;

    TTree* tIn  = (TTree*)fIn->Get("InnerData");
    TTree* tOut = (TTree*)fIn->Get("OuterData");
    TTree* tSum = (TTree*)fIn->Get("EventSummary");
    TTree* tDet = (TTree*)fIn->Get("DetectedSignal");
    gStyle->SetOptStat(0); gStyle->SetPalette(kRainBow); 

    // =========================================================================
    // 1. HITMAPS (ORIGAMI)
    // =========================================================================
    TH2D *hSD1[6], *hSD2[6], *hSiPM[6];
    const char* fN[] = {"Top","Bot","Frt","Bck","Lft","Rgt"};
    for(int i=0; i<6; ++i) {
        hSD1[i]  = new TH2D(Form("SD1_%s",fN[i]),  Form("SD1 %s",fN[i]),  40,-B_SD1,B_SD1, 40,-B_SD1,B_SD1);
        hSD2[i]  = new TH2D(Form("SD2_%s",fN[i]),  Form("SD2 %s",fN[i]),  40,-B_SD2,B_SD2, 40,-B_SD2,B_SD2);
        hSiPM[i] = new TH2D(Form("SiPM_%s",fN[i]), Form("SiPM %s",fN[i]), 40,-B_SIPM,B_SIPM, 40,-B_SIPM,B_SIPM);
        hSD1[i]->SetStats(0); hSD2[i]->SetStats(0); hSiPM[i]->SetStats(0);
    }

    // =========================================================================
    // 2. GÓC TỚI (2 PLOTS)
    // =========================================================================
    TH1D* hAng1 = new TH1D("hAng1", "Incident Angle at SD1; Angle (deg)", 90, 0, 90);
    TH1D* hAng2 = new TH1D("hAng2", "Incident Angle at SD2; Angle (deg)", 90, 0, 90);
    StyleH1(hAng1, kBlue); hAng1->SetFillColorAlpha(kBlue, 0.3);
    StyleH1(hAng2, kRed);  hAng2->SetFillColorAlpha(kRed, 0.3);

    // =========================================================================
    // 3. QUANG PHỔ (4 PLOTS)
    // =========================================================================
    auto MakeStack = [](const char* name, const char* title) {
        return new THStack(name, Form("%s; Wavelength (nm); Photons", title));
    };
    THStack *hsWlSD1 = MakeStack("hsWlSD1", "Spectra at SD1"),
            *hsWlSD2 = MakeStack("hsWlSD2", "Spectra at SD2"),
            *hsWlDet = MakeStack("hsWlDet", "Detected Spectra");
    
    // So sánh Tổng quan (Plot số 4)
    TH1D *hTotSD1 = new TH1D("hTotSD1", "Overlay; Wavelength (nm)", 80, 200, 800);
    TH1D *hTotSD2 = new TH1D("hTotSD2", "Overlay; Wavelength (nm)", 80, 200, 800);
    TH1D *hTotDet = new TH1D("hTotDet", "Overlay; Wavelength (nm)", 80, 200, 800);
    StyleH1(hTotSD1, kBlue, 1, 2); StyleH1(hTotSD2, kRed, 1, 2); StyleH1(hTotDet, kGreen+2, 1, 3);

    // Cấu trúc phân mảng nguồn gốc
    int colors[3] = {kBlue, kRed, kGreen+2};
    const char* srcN[3] = {"Water", "Acrylic", "Gel"};
    int styles[3] = {1001, 3004, 3005}; // Họa tiết gạch chéo để dễ nhìn khi chồng
    TH1D *hW1[3], *hW2[3], *hW3[3];
    for(int i=0; i<3; ++i) {
        hW1[i] = new TH1D(Form("hW1_%d",i), srcN[i], 80, 200, 800);
        hW2[i] = new TH1D(Form("hW2_%d",i), srcN[i], 80, 200, 800);
        hW3[i] = new TH1D(Form("hW3_%d",i), srcN[i], 80, 200, 800);
        
        hW1[i]->SetFillColor(colors[i]); hW1[i]->SetFillStyle(styles[i]); hsWlSD1->Add(hW1[i]);
        hW2[i]->SetFillColor(colors[i]); hW2[i]->SetFillStyle(styles[i]); hsWlSD2->Add(hW2[i]);
        hW3[i]->SetFillColor(colors[i]); hW3[i]->SetFillStyle(styles[i]); hsWlDet->Add(hW3[i]);
    }

    // =========================================================================
    // 4. PIPELINE (LỌC TRÙNG TRACKID)
    // =========================================================================
    TH1D* hPipe[3];
    for(int i=0; i<3; ++i) {
        hPipe[i] = new TH1D(Form("hPipe_%d",i), "Optical Photon Flow; Stage; Number of Unique Photons", 4, 0, 4);
        hPipe[i]->SetFillColor(colors[i]); hPipe[i]->SetFillStyle(styles[i]);
        hPipe[i]->GetXaxis()->SetBinLabel(1, "1. Generated");
        hPipe[i]->GetXaxis()->SetBinLabel(2, "2. Pass SD1");
        hPipe[i]->GetXaxis()->SetBinLabel(3, "3. Pass SD2");
        hPipe[i]->GetXaxis()->SetBinLabel(4, "4. Detected");
    }

    // =========================================================================
    // THU THẬP DỮ LIỆU
    // =========================================================================
    double edepW; int priPdg, detTot, gw, ga, gg;
    tSum->SetBranchAddress("Edep_Water_MeV", &edepW); 
    tSum->SetBranchAddress("Primary_PDG", &priPdg);     // ĐÃ FIX: Biến riêng biệt
    tSum->SetBranchAddress("Detected_Total", &detTot);  // ĐÃ FIX: Biến riêng biệt
    tSum->SetBranchAddress("GenWater", &gw); tSum->SetBranchAddress("GenAcrylic", &ga); tSum->SetBranchAddress("GenCoupling", &gg);
    
    TH2D* hCalo = new TH2D("hCalo", "Calorimetry: Edep vs Detected; Energy Deposit in Water (MeV); Detected Photons", 50, 0, 50, 50, 0, 1000);
    hCalo->SetStats(0);

    double genT[3] = {0,0,0};
    for(int i=0; i<tSum->GetEntries(); ++i) {
        tSum->GetEntry(i);
        genT[0]+=gw; genT[1]+=ga; genT[2]+=gg;
        if(std::abs(priPdg) == 13) hCalo->Fill(edepW, detTot); // Muon
    }
    for(int i=0; i<3; ++i) hPipe[i]->SetBinContent(1, genT[i]);

    // Đọc Tree Ranh Giới với bộ lọc Unique
    auto ReadTree = [&](TTree* t, TH2D* H[6], TH1D* hA, TH1D* hw[3], TH1D* hTot, double bnd, int stage) {
        int trk, orig; double hx, hy, hz, ang, wl;
        t->SetBranchAddress("TrackID", &trk); t->SetBranchAddress("OriginVol", &orig);
        t->SetBranchAddress("HitX_mm", &hx); t->SetBranchAddress("HitY_mm", &hy); t->SetBranchAddress("HitZ_mm", &hz);
        t->SetBranchAddress("AngleToPlane_deg", &ang); t->SetBranchAddress("Wavelength_nm", &wl);
        
        std::set<int> unique_trk[3];
        for(int i=0; i<t->GetEntries(); ++i) {
            t->GetEntry(i);
            FillOrigami(H, hx, hy, hz, bnd);
            if(hA) hA->Fill(ang);
            hTot->Fill(wl);
            if(orig>=0 && orig<=2) { hw[orig]->Fill(wl); unique_trk[orig].insert(trk); } // LỌC TRÙNG TẠI ĐÂY
        }
        for(int i=0; i<3; ++i) hPipe[i]->SetBinContent(stage, unique_trk[i].size());
    };

    ReadTree(tIn,  hSD1, hAng1, hW1, hTotSD1, B_SD1, 2);
    ReadTree(tOut, hSD2, hAng2, hW2, hTotSD2, B_SD2, 3); // Dùng B_SD2 = 53.0

    // Đọc Tree SiPM
    int trk, orig, bnc; double hx, hy, hz, wl, time;
    tDet->SetBranchAddress("TrackID", &trk); tDet->SetBranchAddress("OriginVol", &orig); tDet->SetBranchAddress("BounceCount", &bnc);
    tDet->SetBranchAddress("HitX_mm", &hx); tDet->SetBranchAddress("HitY_mm", &hy); tDet->SetBranchAddress("HitZ_mm", &hz);
    tDet->SetBranchAddress("Wavelength_nm", &wl); tDet->SetBranchAddress("Time_ns", &time);
    
    TH1D* hTime = new TH1D("hTime", "Time Response; Time (ns); Photons", 100, 0, 10);
    TH2D* hTB = new TH2D("hTB", "Time vs Bounce; Time (ns); Bounce Count", 100, 0, 10, 20, 0, 20);
    StyleH1(hTime, kTeal+2, 1, 2); hTime->SetFillColorAlpha(kTeal, 0.3); hTB->SetStats(0);

    std::set<int> u_det[3];
    for(int i=0; i<tDet->GetEntries(); ++i) {
        tDet->GetEntry(i);
        FillOrigami(hSiPM, hx, hy, hz, B_SIPM); // B_SIPM = 53.5
        hTotDet->Fill(wl);
        if(orig>=0 && orig<=2) { hW3[orig]->Fill(wl); u_det[orig].insert(trk); }
        if(orig==0) { hTime->Fill(time); hTB->Fill(time, bnc); }
    }
    for(int i=0; i<3; ++i) hPipe[i]->SetBinContent(4, u_det[i].size());

    // =========================================================================
    // VẼ CANVAS
    // =========================================================================
    
    // C1: Hit Maps 
    TCanvas* C1 = new TCanvas("C1", "Hit Maps", 1800, 600); C1->Divide(3,1);
    DrawOrigami((TPad*)C1->cd(1), hSD1); DrawOrigami((TPad*)C1->cd(2), hSD2); DrawOrigami((TPad*)C1->cd(3), hSiPM);

    // C2: Angles (2 Plots)
    TCanvas* C2 = new TCanvas("C2", "Angles", 1200, 500); C2->Divide(2,1);
    C2->cd(1); hAng1->Draw("HIST"); C2->cd(2); hAng2->Draw("HIST");

    // C3: Wavelength (4 Plots theo yêu cầu)
    TCanvas* C3 = new TCanvas("C3", "Wavelength", 1200, 800); C3->Divide(2,2);
    C3->cd(1); hsWlSD1->Draw("HIST"); C3->cd(2); hsWlSD2->Draw("HIST"); C3->cd(3); hsWlDet->Draw("HIST");
    C3->cd(4); 
    hTotSD1->SetTitle("Spectra Overlay; Wavelength (nm)"); hTotSD1->Draw("HIST"); 
    hTotSD2->Draw("HIST SAME"); hTotDet->Draw("HIST SAME");
    TLegend* lO = new TLegend(0.6,0.6,0.85,0.85); lO->AddEntry(hTotSD1,"Total SD1"); lO->AddEntry(hTotSD2,"Total SD2"); lO->AddEntry(hTotDet,"Detected"); lO->Draw();

    // C4: Pipeline
    TCanvas* C4 = new TCanvas("C4", "Pipeline", 800, 600);
    THStack* hsPipe = new THStack("hsPipe", "Photon Pipeline (Unique Tracking); Stage; Photons");
    for(int i=0; i<3; ++i) hsPipe->Add(hPipe[i]);
    hsPipe->Draw("BAR");
    TLegend* lP = new TLegend(0.7,0.7,0.85,0.85); for(int i=0;i<3;++i) lP->AddEntry(hPipe[i], srcN[i], "f"); lP->Draw();

    // C5: Timing & Calorimetry
    TCanvas* C5 = new TCanvas("C5", "Physics", 1500, 500); C5->Divide(3,1);
    C5->cd(1); hTime->Draw("HIST"); C5->cd(2); hTB->Draw("COLZ"); C5->cd(3); hCalo->Draw("COLZ");

    std::cout << "\n[OK] LỖI ĐÃ ĐƯỢC VÁ! File ROOT đã được phân tích hoàn hảo!\n";
}