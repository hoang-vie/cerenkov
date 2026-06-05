// =============================================================================
// analyze_cherenkov.C
// Geant4 Cherenkov simulation analysis  -  ROOT/C++ macro
//
// Data Flow:
//   Tree 4 "OpticalTrace" : Hit-maps, Incident Angles, Reflections, Pipeline
//   Tree 2 "InnerData"    : Total Wavelength Spectrum at SD1
//   Tree 3 "OuterData"    : Total Wavelength Spectrum at SD2
//   Tree 5 "EventSummary" : N_GenPhotons baseline
//
// Surface convention (6 Faces - Separated Top/Bottom):
//   1 = Top (+Y) | 2 = Bottom (-Y) | 3 = Front (+Z) 
//   4 = Back (-Z)| 5 = Left (-X)   | 6 = Right (+X)
// =============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TMath.h>
#include <iostream>
#include <cmath>
#include <set>
#include <utility>
#include <algorithm>
#include <vector>

// =============================================================================
// Geometry constants
// =============================================================================
static const double BOUND_SD1  = 50.0;  // mm  (Water half-size)
static const double BOUND_SD2  = 52.0;  // mm  (Acrylic half-size)
static const double TOLER      = 0.5;   // mm  (Tolerance)

int SurfaceID(double hX, double hY, double hZ, int bID) {
    const double bound = (bID == 1) ? BOUND_SD1 : BOUND_SD2;
    const double dY = std::abs(std::abs(hY) - bound);
    const double dZ = std::abs(std::abs(hZ) - bound);
    const double dX = std::abs(std::abs(hX) - bound);
    const double dMin = std::min({dY, dZ, dX});
    
    if (dMin > TOLER) return -1;
    if (dMin == dY)   return (hY > 0) ? 1 : 2; 
    if (dMin == dZ)   return (hZ > 0) ? 3 : 4; 
    return (hX > 0) ? 6 : 5;                   
}

TVector3 FaceNormal(int sid) {
    switch (sid) {
        case 1: return TVector3( 0, 1, 0);
        case 2: return TVector3( 0,-1, 0);
        case 3: return TVector3( 0, 0, 1);
        case 4: return TVector3( 0, 0,-1);
        case 5: return TVector3(-1, 0, 0);
        case 6: return TVector3( 1, 0, 0);
        default: return TVector3(0, 0, 0);
    }
}

void StyleH1(TH1D* h, int color) { h->SetLineColor(color); h->SetLineWidth(2); h->SetStats(0); }

void StyleH2(TH2D* h) {
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.06); h->GetXaxis()->SetTitleOffset(0.85);
    h->GetYaxis()->SetTitleSize(0.06); h->GetYaxis()->SetTitleOffset(0.85);
    h->GetZaxis()->SetLabelSize(0.04);
}


void SyncZScale(const std::vector<TH2D*>& histos) {
    double maxZ = 1.0;
    for (auto* h : histos) {
        if (h) maxZ = std::max(maxZ, h->GetMaximum());
    }

    for (auto* h : histos) {
        if (h) { h->SetMaximum(maxZ); h->SetMinimum(0); }
    }
}

TLegend* FaceLegend(double x1, double y1, double x2, double y2, std::vector<TH1D*>& h) {
    TLegend* leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(1); leg->SetFillStyle(0); leg->SetTextSize(0.036);
    leg->AddEntry(h[0], "All faces",  "l");
    leg->AddEntry(h[1], "Top (+Y)",   "l");
    leg->AddEntry(h[2], "Bot (-Y)",   "l");
    leg->AddEntry(h[3], "Front (+Z)", "l");
    leg->AddEntry(h[4], "Back (-Z)",  "l");
    leg->AddEntry(h[5], "Left (-X)",  "l");
    leg->AddEntry(h[6], "Right (+X)", "l");
    return leg;
}

void DrawOrigamiHalf(TPad* pad, TH2D* hTop, TH2D* hBot, TH2D* hFront, TH2D* hBack, TH2D* hLeft, TH2D* hRight) {
    pad->cd(); pad->Divide(3, 4, 0.002, 0.002);
    const int kCOLZ_pad[] = {2, 8, 5, 11, 7, 9}; 
    TH2D* histos[]        = {hTop, hBot, hFront, hBack, hLeft, hRight};
    for (int i = 0; i < 6; ++i) {
        pad->cd(kCOLZ_pad[i]); gPad->SetRightMargin(0.18); histos[i]->Draw("COLZ");
    }
}

// =============================================================================
// MAIN ANALYSIS MACRO
// =============================================================================
void analyze_cherenkov() {

    // CHÚ Ý: Bạn có thể sửa tên file ở đây để test ("data_air.root" hoặc "data_gel.root")
    const TString kInFile  = "/home/huyhoang/software/Geant4_hoang/Cerenkov/build/cherenkov_water_data.root";
    const TString kOutFile = "cherenkov_water_data.root";

    TFile* fIn = TFile::Open(kInFile, "READ");
    if (!fIn || fIn->IsZombie()) { std::cerr << "[ERROR] Cannot open: " << kInFile << std::endl; return; }
    std::cout << "[INFO] Opened: " << kInFile << std::endl;

    TTree* tBnd = (TTree*)fIn->Get("OpticalTrace");
    if (!tBnd) tBnd = (TTree*)fIn->Get("Optical_Boundaries");
    TTree* tSum = (TTree*)fIn->Get("EventSummary");
    TTree* tIn  = (TTree*)fIn->Get("InnerData");
    TTree* tOut = (TTree*)fIn->Get("OuterData");

    if (!tBnd) { std::cerr << "[ERROR] Cannot find Tree OpticalTrace!\n"; return; }
    const bool hasAng = (tBnd->GetBranch("AngleToNormal_deg") != nullptr);

    // ==========================================================================
    // STEP 1  -  EventSummary
    // ==========================================================================
    long long sumGen = 0, sumHitSD1 = 0, sumEscSD2 = 0;
    if (tSum) {
        int evID = 0, nGen = 0, nHit = 0, nEsc = 0;
        tSum->SetBranchAddress("EventID", &evID);
        if (tSum->GetBranch("N_GenPhotons")) {
            tSum->SetBranchAddress("N_GenPhotons", &nGen); tSum->SetBranchAddress("N_HitSD1", &nHit); tSum->SetBranchAddress("N_EscapedSD2", &nEsc);
        } else {
            tSum->SetBranchAddress("N_Gen", &nGen); tSum->SetBranchAddress("N_W2A", &nHit); tSum->SetBranchAddress("N_A2Air", &nEsc);
        }
        for (Long64_t i = 0; i < tSum->GetEntries(); ++i) { tSum->GetEntry(i); sumGen += nGen; sumHitSD1 += nHit; sumEscSD2 += nEsc; }
        tSum->ResetBranchAddresses();
    }

    // ==========================================================================
    // STEP 2  -  Book Histograms
    // ==========================================================================
    const double BMAX = std::max(BOUND_SD1, BOUND_SD2) + 5.0; 
    const int    NBXY = 60, NBW = 100, NBA = 90;
    const int kFaceColor[7] = {kBlack, kCyan+2, kBlue, kRed, kGreen+2, kMagenta, kOrange+1};

    TH1D* hWaveSD1 = new TH1D("hWaveSD1", "Wavelength at SD1 (Total); #lambda (nm); Unique photons", NBW, 200, 800);
    TH1D* hWaveSD2 = new TH1D("hWaveSD2", "Wavelength at SD2 (Total); #lambda (nm); Unique photons", NBW, 200, 800);
    StyleH1(hWaveSD1, kBlue); hWaveSD1->SetFillColorAlpha(kBlue, 0.3);
    StyleH1(hWaveSD2, kRed);  hWaveSD2->SetFillColorAlpha(kRed, 0.3);

    auto MakeH1Set = [&](const char* tag, const char* title) -> std::vector<TH1D*> {
        const char* sfx[7] = {"_Tot","_Top","_Bot","_Frnt","_Back","_Left","_Rgt"};
        std::vector<TH1D*> v(7);
        for (int i = 0; i < 7; ++i) {
            v[i] = new TH1D(Form("%s%s", tag, sfx[i]), Form("%s; Angle to normal (deg); Unique photons", title), NBA, 0, 90);
            StyleH1(v[i], kFaceColor[i]);
        }
        return v;
    };
    auto aSD1 = MakeH1Set("aSD1","Incident angle at SD1");
    auto aSD2 = MakeH1Set("aSD2","Incident angle at SD2");

    TH2D* h1Top   = new TH2D("h1Top",  "SD1 Top (+Y); X (mm); Z (mm)",   NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h1Bot   = new TH2D("h1Bot",  "SD1 Bot (-Y); X (mm); Z (mm)",   NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h1Front = new TH2D("h1Front","SD1 Front (+Z); X (mm); Y (mm)", NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h1Back  = new TH2D("h1Back", "SD1 Back (-Z); X (mm); -Y (mm)", NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h1Left  = new TH2D("h1Left", "SD1 Left (-X); -Y (mm); Z (mm)", NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h1Right = new TH2D("h1Right","SD1 Right (+X); Y (mm); Z (mm)", NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);

    TH2D* h2Top   = new TH2D("h2Top",  "SD2 Top (+Y); X (mm); Z (mm)",   NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h2Bot   = new TH2D("h2Bot",  "SD2 Bot (-Y); X (mm); Z (mm)",   NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h2Front = new TH2D("h2Front","SD2 Front (+Z); X (mm); Y (mm)", NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h2Back  = new TH2D("h2Back", "SD2 Back (-Z); X (mm); -Y (mm)", NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h2Left  = new TH2D("h2Left", "SD2 Left (-X); -Y (mm); Z (mm)", NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);
    TH2D* h2Right = new TH2D("h2Right","SD2 Right (+X); Y (mm); Z (mm)", NBXY,-BMAX,BMAX, NBXY,-BMAX,BMAX);

    for (auto* h : {h1Top,h1Bot,h1Front,h1Back,h1Left,h1Right, h2Top,h2Bot,h2Front,h2Back,h2Left,h2Right}) StyleH2(h);

    long long rfSD1[6][2] = {}, rfSD2[6][2] = {};
    auto FaceIdx = [](int sid) -> int { return (sid >= 1 && sid <= 6) ? sid - 1 : -1; };
    using Key = std::pair<int,int>;

    // ==========================================================================
    // STEP 3  -  Wavelength Extraction
    // ==========================================================================
    if (tIn) {
        double wl = 0; int tid = 0, eid = 0; tIn->SetBranchAddress("Wavelength_nm", &wl); tIn->SetBranchAddress("TrackID", &tid); tIn->SetBranchAddress("EventID", &eid);
        std::set<Key> seenWL1;
        for(Long64_t i=0; i<tIn->GetEntries(); ++i) { tIn->GetEntry(i); if(seenWL1.insert({eid, tid}).second) hWaveSD1->Fill(wl); }
        tIn->ResetBranchAddresses();
    }
    if (tOut) {
        double wl = 0; int tid = 0, eid = 0; tOut->SetBranchAddress("Wavelength_nm", &wl); tOut->SetBranchAddress("TrackID", &tid); tOut->SetBranchAddress("EventID", &eid);
        std::set<Key> seenWL2;
        for(Long64_t i=0; i<tOut->GetEntries(); ++i) { tOut->GetEntry(i); if(seenWL2.insert({eid, tid}).second) hWaveSD2->Fill(wl); }
        tOut->ResetBranchAddresses();
    }

    // ==========================================================================
    // STEP 4  -  Trace boundaries
    // ==========================================================================
    int evID=0, trkID=0, bndID=0, isRef=0; double hX=0, hY=0, hZ=0, inX=0, inY=0, inZ=0, angDirect=0;
    tBnd->SetBranchAddress("EventID", &evID); tBnd->SetBranchAddress("TrackID", &trkID); tBnd->SetBranchAddress("BoundaryID", &bndID);
    tBnd->SetBranchAddress("IsReflected", &isRef); tBnd->SetBranchAddress("HitX_mm", &hX); tBnd->SetBranchAddress("HitY_mm", &hY); tBnd->SetBranchAddress("HitZ_mm", &hZ);
    if (!hasAng) { tBnd->SetBranchAddress("InDirX", &inX); tBnd->SetBranchAddress("InDirY", &inY); tBnd->SetBranchAddress("InDirZ", &inZ); } 
    else { tBnd->SetBranchAddress("AngleToNormal_deg", &angDirect); }

    std::set<Key> seenSD1, seenSD2, transSD1, transSD2;
    for (Long64_t i = 0; i < tBnd->GetEntries(); ++i) {
        tBnd->GetEntry(i);
        if (bndID != 1 && bndID != 2) continue;  
        const int sid = SurfaceID(hX, hY, hZ, bndID); if (sid < 0) continue; 
        const Key uid = {evID, trkID}; const int fidx = FaceIdx(sid);

        double angle = hasAng ? angDirect : std::acos(std::min(std::abs(TVector3(inX, inY, inZ).Unit().Dot(FaceNormal(sid))), 1.0)) * TMath::RadToDeg();

        if (bndID == 1 && seenSD1.insert(uid).second) {
            aSD1[0]->Fill(angle); aSD1[fidx+1]->Fill(angle);
            if(sid==1) h1Top->Fill(hX,hZ); else if(sid==2) h1Bot->Fill(hX,hZ); else if(sid==3) h1Front->Fill(hX,hY); else if(sid==4) h1Back->Fill(hX,-hY); else if(sid==5) h1Left->Fill(-hY,hZ); else if(sid==6) h1Right->Fill(hY,hZ);
            rfSD1[fidx][1]++; if (isRef == 1) rfSD1[fidx][0]++; else transSD1.insert(uid); 
        }
        else if (bndID == 2 && seenSD2.insert(uid).second) {
            aSD2[0]->Fill(angle); aSD2[fidx+1]->Fill(angle);
            if(sid==1) h2Top->Fill(hX,hZ); else if(sid==2) h2Bot->Fill(hX,hZ); else if(sid==3) h2Front->Fill(hX,hY); else if(sid==4) h2Back->Fill(hX,-hY); else if(sid==5) h2Left->Fill(-hY,hZ); else if(sid==6) h2Right->Fill(hY,hZ);
            rfSD2[fidx][1]++; if (isRef == 1) rfSD2[fidx][0]++; else transSD2.insert(uid);
        }
    }
    tBnd->ResetBranchAddresses();

    // ==========================================================================
    // STEP 5  -  Efficiency Calculation
    // ==========================================================================
    const long long N_Gen      = (sumGen > 0) ? sumGen : (long long)seenSD1.size();
    const long long N_HitSD1   = (long long)seenSD1.size();
    const long long N_TransSD1 = (long long)transSD1.size();
    const long long N_HitSD2   = (long long)seenSD2.size();
    const long long N_TransSD2 = (long long)transSD2.size();

    auto Pct = [](long long num, long long den) -> TString { return (den<=0) ? "N/A" : TString(Form("%.3f %%", 100.0*num/den)); };

    std::cout << "\n=========== CHI TIET TRUYEN QUA VA PHAN XA ===========\n";
    std::cout << "  1. Bắt đầu (EventSummary): " << N_Gen      << " photons\n";
    std::cout << "  2. Đập vào Acrylic (SD1) : " << N_HitSD1   << " (" << Pct(N_HitSD1,  N_Gen)      << ")\n";
    std::cout << "  3. Lọt qua SD1           : " << N_TransSD1 << " (" << Pct(N_TransSD1,N_HitSD1)   << " truyền qua SD1)\n";
    std::cout << "  4. Đập vào lớp vỏ (SD2)  : " << N_HitSD2   << "\n";
    std::cout << "  5. LỌT QUA VỎ (TransSD2) : " << N_TransSD2 << " (" << Pct(N_TransSD2,N_HitSD2)   << " TRUYỀN QUA VỎ MỚI!)\n";
    std::cout << "======================================================\n\n";

    // ==========================================================================
    // STEP 6  -  Draw Canvases
    // ==========================================================================
    gStyle->SetOptStat(0); gStyle->SetPalette(kRainBow); gStyle->SetTitleSize(0.048, "t");

    std::vector<TH2D*> allMaps = {
        h1Top, h1Bot, h1Front, h1Back, h1Left, h1Right,
        h2Top, h2Bot, h2Front, h2Back, h2Left, h2Right
    };
    SyncZScale(allMaps);

    TCanvas* cOrig = new TCanvas("cOrig","Origami Hit Maps - SD1 (left) | SD2 (right)", 1800, 900);
    TPad* pLeft  = new TPad("pLeft",  "SD1", 0.00, 0.0, 0.50, 1.0); pLeft->Draw();
    TPad* pRight = new TPad("pRight", "SD2", 0.50, 0.0, 1.00, 1.0); pRight->Draw();
    DrawOrigamiHalf(pLeft, h1Top, h1Bot, h1Front, h1Back, h1Left, h1Right);
    DrawOrigamiHalf(pRight, h2Top, h2Bot, h2Front, h2Back, h2Left, h2Right);

    TCanvas* cWave = new TCanvas("cWave","Wavelength Distributions", 1400, 500);
    cWave->Divide(2, 1);
    cWave->cd(1); gPad->SetGridy(); hWaveSD1->Draw("HIST");
    cWave->cd(2); gPad->SetGridy(); hWaveSD2->Draw("HIST");

    auto DrawAngleHalf = [&](TCanvas* cv, int padTot, int padSub, std::vector<TH1D*>& aSet, const char* sdLabel) {
        cv->cd(padTot); gPad->SetGridy();
        aSet[0]->SetTitle(Form("%s  -  All faces; Angle to normal (deg); Unique photons", sdLabel));
        aSet[0]->Draw("HIST"); for (int f=1; f<7; ++f) aSet[f]->Draw("HIST SAME");
        FaceLegend(0.64,0.42,0.92,0.92, aSet)->Draw();

        cv->cd(padSub); ((TPad*)gPad)->Divide(3, 2, 0.004, 0.004);
        const char* fnames[6] = {"Top (+Y)","Bottom (-Y)","Front (+Z)","Back (-Z)","Left (-X)","Right (+X)"};
        for (int f=0; f<6; ++f) {
            ((TPad*)gPad)->cd(f+1); gPad->SetGridy();
            aSet[f+1]->SetTitle(Form("%s  -  %s; Angle (deg); Counts", sdLabel, fnames[f]));
            aSet[f+1]->Draw("HIST");
        }
    };
    TCanvas* cAng = new TCanvas("cAng","Incident Angle Distributions",1600,700);
    cAng->Divide(2, 2); 
    DrawAngleHalf(cAng, 1, 2, aSD1, "SD1 (Water#rightarrowAcrylic)");
    DrawAngleHalf(cAng, 3, 4, aSD2, "SD2 (Acrylic#rightarrowCoupling)");

    TCanvas* cPipe = new TCanvas("cPipe","Photon Survival Pipeline",1200,650);
    cPipe->SetBottomMargin(0.22); cPipe->SetRightMargin(0.06); cPipe->SetGridy();
    TH1D* hPipe = new TH1D("hPipe","Optical Photon Survival Pipeline; ; Unique photons", 5, 0, 5);
    hPipe->SetFillColor(kAzure+1); hPipe->SetBarWidth(0.55); hPipe->SetBarOffset(0.225);
    const char* stageLabels[5] = {"Generated","Hit SD1","Trans SD1","Hit SD2","Trans SD2 (Escaped)"};
    const long long stageVals[5] = {N_Gen, N_HitSD1, N_TransSD1, N_HitSD2, N_TransSD2};
    for (int b = 1; b <= 5; ++b) { hPipe->SetBinContent(b, stageVals[b-1]); hPipe->GetXaxis()->SetBinLabel(b, stageLabels[b-1]); }
    hPipe->GetXaxis()->SetLabelSize(0.050); hPipe->Draw("B TEXT90");

    // TỐI ƯU C5: HIỂN THỊ RÕ RÀNG PHẢN XẠ VÀ TRUYỀN QUA
    TCanvas* cRef = new TCanvas("cRef","Reflection vs Transmission per Face",1600,600);
    cRef->Divide(2, 1);
    const char* faceName[6] = {"Top","Bottom","Front","Back","Left","Right"};
    auto DrawRefTransPanel = [&](long long rf[][2], const char* title) {
        TH1D* h = new TH1D(title, Form("%s; Face; Ratio (%%)", title), 6, 0, 6);
        h->SetStats(0); h->GetYaxis()->SetRangeUser(0.0, 130.0); h->SetFillColor(kAzure+2);
        for (int f=0; f<6; ++f) {
            h->GetXaxis()->SetBinLabel(f+1, faceName[f]);
            if(rf[f][1]>0) h->SetBinContent(f+1, 100.0*rf[f][0]/rf[f][1]); // Vẽ cột Phản xạ
        }
        gPad->SetGridy(); gPad->SetBottomMargin(0.15); h->Draw("BAR");
        TLatex tval; tval.SetTextSize(0.038); tval.SetTextAlign(22);
        for (int f=0; f<6; ++f) {
            double fracRef = (rf[f][1]>0) ? 100.0*rf[f][0]/rf[f][1] : 0.0;
            double fracTrans = (rf[f][1]>0) ? 100.0 - fracRef : 0.0;
            // In thẳng lên đồ thị: Ref = ...%, Trans = ...%
            tval.DrawLatex(h->GetXaxis()->GetBinCenter(f+1), h->GetBinContent(f+1) + 8.0, 
                           Form("#color[2]{Ref: %.1f%%}\n#color[4]{Trans: %.1f%%}", fracRef, fracTrans));
        }
        return h;
    };
    cRef->cd(1); DrawRefTransPanel(rfSD1, "SD1 (Water#rightarrowAcrylic)");
    cRef->cd(2); DrawRefTransPanel(rfSD2, "SD2 (Acrylic#rightarrowCoupling)");

    // ==========================================================================
    // STEP 7  -  Save
    // ==========================================================================
    TFile* fOut = TFile::Open(kOutFile, "RECREATE");
    cOrig->Write("C1_OrigamiHitMaps"); cWave->Write("C2_Wavelength");
    cAng->Write("C3_IncidentAngle"); cPipe->Write("C4_PhotonPipeline"); cRef->Write("C5_ReflectionFraction");
    
    fOut->mkdir("HitMaps/SD1")->cd(); for (auto* h : {h1Top,h1Bot,h1Front,h1Back,h1Left,h1Right}) h->Write();
    fOut->mkdir("HitMaps/SD2")->cd(); for (auto* h : {h2Top,h2Bot,h2Front,h2Back,h2Left,h2Right}) h->Write();
    fOut->cd(); hPipe->Write("Hist_PhotonPipeline");
    fOut->Close(); fIn->Close();
    std::cout << "[DONE] Output saved to: " << kOutFile << "\n";
}

