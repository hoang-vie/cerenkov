// =============================================================================
// analyze_v6_integrated.C — Cherenkov water-tank analysis (schema v5: 8 Trees)
// Includes: True Particle PDG Tracking (Muon vs Electron), Absolute Time-Cut,
//           Secondary/Primary Ratio Discriminator, Waveform Tail QA,
//           and ALL original analysis plots with refined Physics Titles/Labels.
// =============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

// =============================================================================
// STYLE HELPERS
// =============================================================================
void ApplyStyles1D(TH1* h, int color, int fillStyle = 1001, int lineWidth = 2) {
    h->SetLineColor(color); h->SetLineWidth(lineWidth); h->SetFillColorAlpha(color, 0.25);
    h->SetFillStyle(fillStyle); h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.04);  h->GetYaxis()->SetLabelSize(0.04);
}
void ApplyStyles2D(TH2* h) {
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.04);  h->GetYaxis()->SetLabelSize(0.04);
}

// =============================================================================
// Origami: unfold the 6 cube faces onto one canvas
// =============================================================================
void FillOrigami(TH2D* hFaces[6], double x, double y, double z) {
    double aX = std::abs(x), aY = std::abs(y), aZ = std::abs(z);
    if      (aY >= aX && aY >= aZ && y <  0) hFaces[0]->Fill(x, z);    // BOTTOM (-Y)
    else if (aZ >= aX && aZ >= aY && z >= 0) hFaces[1]->Fill(x, y);    // FRONT  (+Z)
    else if (aX >= aY && aX >= aZ && x >= 0) hFaces[2]->Fill(y, z);    // RIGHT  (+X)
    else if (aZ >= aX && aZ >= aY && z <  0) hFaces[3]->Fill(x, -y);   // BACK   (-Z)
    else if (aX >= aY && aX >= aZ && x <  0) hFaces[4]->Fill(-y, z);   // LEFT   (-X)
    else                                      hFaces[5]->Fill(x, -z);  // TOP    (+Y)
}

void DrawOrigamiCanvas(TCanvas* cOut, const TString& pdfOut, TH2D* hFaces[6], const TString& title) {
    cOut->Clear();
    double h = 0.23; double w = h * (650.0 / 1200.0);
    double cx = 0.5; double cy = 0.5 - h / 2.0;
    TPad* pads[6] = {
        new TPad("pBot",   "B",  cx - w/2, cy - h/2,   cx + w/2, cy + h/2),
        new TPad("pFront", "F",  cx - w/2, cy + h/2,   cx + w/2, cy + 3*h/2),
        new TPad("pRight", "R",  cx + w/2, cy - h/2,   cx + 3*w/2, cy + h/2),
        new TPad("pBack",  "Ba", cx - w/2, cy - 3*h/2, cx + w/2, cy - h/2),
        new TPad("pLeft",  "L",  cx - 3*w/2, cy - h/2, cx - w/2, cy + h/2),
        new TPad("pTop",   "T",  cx - w/2, cy + 3*h/2, cx + w/2, cy + 5*h/2)
    };

    double gMax = 1;
    for (int i = 0; i < 6; i++) if (hFaces[i] && hFaces[i]->GetMaximum() > gMax) gMax = hFaces[i]->GetMaximum();
    for (int i = 0; i < 6; i++) {
        cOut->cd(); pads[i]->Draw(); pads[i]->cd();
        gPad->SetFixedAspectRatio(kTRUE);
        gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.18); gPad->SetTopMargin(0.15);
        if (hFaces[i]) { hFaces[i]->SetStats(0); hFaces[i]->SetMinimum(0); hFaces[i]->SetMaximum(gMax); hFaces[i]->Draw("COLZ"); }
    }
    cOut->cd();
    TLatex tex; tex.SetNDC(); tex.SetTextFont(62); tex.SetTextSize(0.04);
    tex.DrawLatex(0.05, 0.90, title);
    cOut->Print(pdfOut);
}

void DrawFunnelWithPct(TCanvas* cOut, const TString& pdfOut, TH1D* hF) {
    cOut->Clear();
    hF->SetMaximum(hF->GetMaximum() * 1.3);
    hF->Draw("HIST BAR");
    double base = hF->GetBinContent(1);
    for (int i = 1; i <= hF->GetNbinsX(); i++) {
        double v = hF->GetBinContent(i);
        double p = (base > 0) ? (v / base * 100.0) : 0;
        TLatex* tV = new TLatex(i - 0.5, v + hF->GetMaximum() * 0.02, Form("%.2f", v));
        tV->SetTextAlign(21); tV->SetTextSize(0.032); tV->Draw();
        TLatex* tP = new TLatex(i - 0.5, v + hF->GetMaximum() * 0.07, Form("(%.1f%%)", p));
        tP->SetTextAlign(21); tP->SetTextSize(0.032); tP->SetTextColor(kRed + 2); tP->Draw();
    }
    cOut->Print(pdfOut);
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================
void analyze(const char* inputName = "visual_data.root",
                 double timeCutWall_ns = 0.50,
                 long long maxEvents = -1,
                 double waterHalfXZ = 50.0) {

    TFile* fIn = TFile::Open(inputName, "READ");
    if (!fIn || fIn->IsZombie()) { std::cerr << "Error: cannot open " << inputName << "\n"; return; }

    TString base(inputName);
    { Ssiz_t slash = base.Last('/'); if (slash != kNPOS) base.Remove(0, slash + 1); }
    TString pdfOut  = "report_" + base; pdfOut.ReplaceAll(".root", ".pdf");
    TString rootOut = "hists_"  + base;

    std::cout << ">>> Processing file: " << inputName << "\n";
    std::cout << ">>> Time Cut Window for Early Wall Hitmap: " << timeCutWall_ns << " ns (relative to trigger).\n";
    if (maxEvents > 0) std::cout << ">>> TEST MODE: limited to " << maxEvents << " events.\n";

    TFile* fOut = new TFile(rootOut, "RECREATE");

    TTree* t0 = (TTree*)fIn->Get("Secondary");
    TTree* t1 = (TTree*)fIn->Get("PhotonGen");
    TTree* t2 = (TTree*)fIn->Get("Boundary_A2C");
    TTree* t3 = (TTree*)fIn->Get("EventSummary");
    TTree* t4 = (TTree*)fIn->Get("DetectedSignal");
    TTree* t5 = (TTree*)fIn->Get("WallHitmap");
    TTree* t6 = (TTree*)fIn->Get("PrimarySummary");
    TTree* t7 = (TTree*)fIn->Get("PrimaryEdepPerVolume");
    if (!t0 || !t1 || !t2 || !t3 || !t4 || !t5 || !t6 || !t7) {
        std::cerr << "Error: missing core tree(s).\n";
        return;
    }

    TTree* t8 = (TTree*)fIn->Get("RunConfig");

    // [v6] "Water*" -> "Core*": khoi loi (radiator) gio co the la Water HOAC Acrylic,
    // xem DetectorConstruction::SetCoreMaterial / Tree8.CoreMaterial trong TreeInfo.md.
    double coreFullX = -1, coreFullY = -1, coreFullZ = -1;
    double scintX=-1, scintY=-1, scintZ=-1, sipmZ=-1, acrThick=-1, tyvThick=-1;
    TString tyvCfgStr = "N/A";
    TString coreMaterialStr = "Unknown";  // [MOI] doc tu Tree8.CoreMaterial - PHAI kiem tra truoc khi
                                           // dien giai bat ky so lieu Cherenkov nao (xem TreeInfo.md)
    // [MOI] doc tu Tree8.UseAcrylicShellInt. Mac dinh TRUE de tuong thich file ROOT CU
    // (sinh ra truoc khi them cot nay) - cac file do luon co vo Acrylic, khong co lua chon khac.
    bool hasAcrylicShell = true;
    bool hasRunConfig = false;

    if (t8 && t8->GetEntries() > 0) {
        // Moi bien co gia tri mac dinh AN TOAN ngay tu khai bao: neu 1 branch bi doi ten/thieu,
        // SetBranchAddress se bao loi va bien se GIU nguyen gia tri mac dinh nay thay vi rac stack.
        double wx = -1, wy = -1, wz = -1;
        double scX = -1, scY = -1, scZ = -1, spZ = -1, acT = -1, tyT = -1;
        int tyvInt = -1;
        int shellInt = 1; // mac dinh 1 (co vo) neu branch khong ton tai (file ROOT cu)
        std::string* coreMatPtr = nullptr;  // G4Analysis SColumn -> branch kieu std::string trong ROOT

        t8->SetBranchAddress("CoreFullX_mm", &wx);
        t8->SetBranchAddress("CoreFullY_mm", &wy);
        t8->SetBranchAddress("CoreFullZ_mm", &wz);
        t8->SetBranchAddress("ScintX_mm", &scX);
        t8->SetBranchAddress("ScintY_mm", &scY);
        t8->SetBranchAddress("ScintZ_mm", &scZ);
        t8->SetBranchAddress("SipmZ_mm", &spZ);
        t8->SetBranchAddress("AcrylicThickness_mm", &acT);
        t8->SetBranchAddress("TyvekThickness_mm", &tyT);
        t8->SetBranchAddress("TyvekConfigInt", &tyvInt);
        if (t8->GetBranch("CoreMaterial")) t8->SetBranchAddress("CoreMaterial", &coreMatPtr);
        // [MOI] Branch nay CHI ton tai o file ROOT sinh sau khi them fUseAcrylicShell vao code -
        // kiem tra GetBranch() truoc, tranh loi "missing branch" tren file ROOT cu.
        if (t8->GetBranch("UseAcrylicShellInt")) t8->SetBranchAddress("UseAcrylicShellInt", &shellInt);
        t8->GetEntry(0);

        coreFullX = wx; coreFullY = wy; coreFullZ = wz;
        scintX = scX; scintY = scY; scintZ = scZ;
        sipmZ = spZ; acrThick = acT; tyvThick = tyT;
        tyvCfgStr = (tyvInt == 1) ? "Inside" : "Outside";
        if (coreMatPtr) coreMaterialStr = *coreMatPtr;
        hasAcrylicShell = (shellInt != 0);
        hasRunConfig = true;

        std::cout << ">>> Core material for this run: " << coreMaterialStr << "\n";
        std::cout << ">>> Acrylic shell used for this run: " << (hasAcrylicShell ? "YES" : "NO (Tyvek/Gel wrap Core directly)") << "\n";
    }

    if (hasRunConfig && coreFullX > 0) {
        waterHalfXZ = coreFullX / 2.0;
    }

    gStyle->SetPalette(kTemperatureMap); gStyle->SetOptTitle(1);
    TCanvas* cOut = new TCanvas("cOut", "Analysis", 1200, 650);
    cOut->SetLeftMargin(0.12); cOut->SetBottomMargin(0.13); cOut->SetRightMargin(0.10); cOut->SetTopMargin(0.10);
    cOut->Print(pdfOut + "[");

    // =========================================================================
    // PARTICLE QA & DISCRIMINATION HISTOGRAMS
    // =========================================================================
    TH1D* hSec_PDG = new TH1D("hSec_PDG", "Secondary Charged Particles Produced (Tree0);Particle Species;Track Count", 5, 0, 5);
    ApplyStyles1D(hSec_PDG, kOrange + 2);
    hSec_PDG->GetXaxis()->SetBinLabel(1, "e^{-} / e^{+}");
    hSec_PDG->GetXaxis()->SetBinLabel(2, "#gamma (Photons)");
    hSec_PDG->GetXaxis()->SetBinLabel(3, "#mu^{-} / #mu^{+}");
    hSec_PDG->GetXaxis()->SetBinLabel(4, "#pi^{-} / #pi^{+}");
    hSec_PDG->GetXaxis()->SetBinLabel(5, "Others");

    TH1D* hdEdx_TrueMuon = new TH1D("hdEdx_TrueMuon", "Energy Loss Rate (dE/dx) in Core: True Muons (|PDG| == 13);dE/dx (MeV/cm);Primary Tracks", 150, 0, 5);
    ApplyStyles1D(hdEdx_TrueMuon, kBlue);

    TH1D* hdEdx_PrimaryEl = new TH1D("hdEdx_PrimaryEl", "Energy Loss Rate (dE/dx) in Core: Primary Electrons (|PDG| == 11);dE/dx (MeV/cm);Primary Tracks", 150, 0, 5);
    ApplyStyles1D(hdEdx_PrimaryEl, kRed);

    TH1D* hSecToPrimRatio = new TH1D("hSecToPrimRatio", "Particle Identification Discriminator: Secondary / Primary Photon Ratio;Ratio N_{secondary} / N_{primary};Events", 100, 0, 0.5);
    ApplyStyles1D(hSecToPrimRatio, kViolet + 2);

    // =========================================================================
    // GROUP A: TRUTH -- dE/dx and Cherenkov emission angle
    // =========================================================================
    const double E_MIN_MIP = 300.0; // MeV

    TH2D* hBetheBloch2D = new TH2D("hBetheBloch2D", "Bethe-Bloch Kinematics Validation: Initial Kinetic Energy vs. dE/dx in Core;Primary Particle Energy E_{k} (MeV);dE/dx (MeV/cm)", 100, 0, 2000, 100, 0, 10);
    ApplyStyles2D(hBetheBloch2D);

    TH2D* hFT = new TH2D("hFT", "Frank-Tamm Yield Validation: Cherenkov Photons vs. Path Length in Core;Path Length L_{core} (mm);Generated Photon Count", 100, 0, 150, 100, 0, 12000);
    ApplyStyles2D(hFT);

    // [MOI] Nhung histogram gan voi VolumeID/OriginVol==1 (Acrylic) se RONG (0 entries) khi
    // hasAcrylicShell=false - KHONG PHAI loi, chi la Tree7/Tree1 khong con dong nao voi
    // VolumeID/OriginVol==1 nua vi Acrylic khong ton tai trong hinh hoc run nay. Gan hau to
    // "[No Acrylic shell]" vao title de nguoi xem PDF khong nham la thieu du lieu do bug.
    TString acrSuffix = hasAcrylicShell ? "" : " [No Acrylic shell this run]";

    TH1D* hdEdx_MuAcrylic = new TH1D("hdEdx_MuAcrylic", "Primary Energy Loss Rate (dE/dx) in Acrylic (VolumeID == 1)" + acrSuffix + ";dE/dx (MeV/cm);Primary Tracks", 150, 0, 8);
    ApplyStyles1D(hdEdx_MuAcrylic, kOrange + 1);

    TH1D* hdEdx_MuWater = new TH1D("hdEdx_MuWater", "Primary Energy Loss Rate (dE/dx) in Core (VolumeID == 0);dE/dx (MeV/cm);Primary Tracks", 150, 0, 5);
    ApplyStyles1D(hdEdx_MuWater, kBlue);

    TH1D* hdEdx_MuWaterTotal = new TH1D("hdEdx_MuWaterTotal", "Total Energy Loss Rate in Core (Primary + Secondaries, N_{primary} = 1);dE/dx (MeV/cm);Events", 150, 0, 5);
    ApplyStyles1D(hdEdx_MuWaterTotal, kRed);

    TH1D* hNumPrimaryPerEvt = new TH1D("hNumPrimaryPerEvt", "Primary Multiplicity Distribution per Event;Number of Primary Particles;Events", 10, -0.5, 9.5);
    ApplyStyles1D(hNumPrimaryPerEvt, kGray + 2);

    TH1D* hFT_Occupancy = new TH1D("hFT_Occupancy", "Primary Path Length Distribution in Core;Path Length L_{core} (mm);Primary Tracks", 50, 0, 150);
    ApplyStyles1D(hFT_Occupancy, kGray + 1);

    TH1D* hAng_MuWater = new TH1D("hAng_MuWater", "Cherenkov Emission Angle -- Primary Source, Core Origin;Cherenkov Angle #theta_{C} (deg);Photon Count", 90, 0, 90); ApplyStyles1D(hAng_MuWater, kBlue);
    TH1D* hAng_SeWater = new TH1D("hAng_SeWater", "Cherenkov Emission Angle -- Secondary Source, Core Origin;Cherenkov Angle #theta_{C} (deg);Photon Count", 90, 0, 90); ApplyStyles1D(hAng_SeWater, kCyan+2);
    TH1D* hAng_MuAcr   = new TH1D("hAng_MuAcr",   "Cherenkov Emission Angle -- Primary Source, Acrylic Origin" + acrSuffix + ";Cherenkov Angle #theta_{C} (deg);Photon Count", 90, 0, 90); ApplyStyles1D(hAng_MuAcr, kRed);
    TH1D* hAng_SeAcr   = new TH1D("hAng_SeAcr",   "Cherenkov Emission Angle -- Secondary Source, Acrylic Origin" + acrSuffix + ";Cherenkov Angle #theta_{C} (deg);Photon Count", 90, 0, 90); ApplyStyles1D(hAng_SeAcr, kOrange);
    TH1D* hAng_MuGel   = new TH1D("hAng_MuGel",   "Cherenkov Emission Angle -- Primary Source, Gel Origin;Cherenkov Angle #theta_{C} (deg);Photon Count", 90, 0, 90); ApplyStyles1D(hAng_MuGel, kGreen+2);
    TH1D* hAng_SeGel   = new TH1D("hAng_SeGel",   "Cherenkov Emission Angle -- Secondary Source, Gel Origin;Cherenkov Angle #theta_{C} (deg);Photon Count", 90, 0, 90); ApplyStyles1D(hAng_SeGel, kMagenta);

    TH1D* hTTS = new TH1D("hTTS", "First-Photon Arrival Delay Relative to Trigger;t_{first} - t_{trigger} (ns);Events", 150, 0, 15);
    ApplyStyles1D(hTTS, kRed - 4);

    // =========================================================================
    // GROUP B: EFFICIENCY -- Funnel chain & Yield
    // =========================================================================
    TH1D* hFunnelMu = new TH1D("hFunnelMu", "Photon Transport Efficiency Chain -- Primary Source;Transport Stage;Mean Photons / Event", 4, 0, 4); ApplyStyles1D(hFunnelMu, kBlue);
    TH1D* hFunnelSe = new TH1D("hFunnelSe", "Photon Transport Efficiency Chain -- Secondary Source;Transport Stage;Mean Photons / Event", 4, 0, 4); ApplyStyles1D(hFunnelSe, kRed);

    const double WVL_MIN = 300.0; const double WVL_MAX = 750.0;
    TH1D* hFunnelMu_F = new TH1D("hFunnelMu_F", Form("Filtered Transport Efficiency Chain (%.0f-%.0f nm) -- Primary Source;Transport Stage;Mean Photons / Event", WVL_MIN, WVL_MAX), 4, 0, 4); ApplyStyles1D(hFunnelMu_F, kAzure + 2);
    TH1D* hFunnelSe_F = new TH1D("hFunnelSe_F", Form("Filtered Transport Efficiency Chain (%.0f-%.0f nm) -- Secondary Source;Transport Stage;Mean Photons / Event", WVL_MIN, WVL_MAX), 4, 0, 4); ApplyStyles1D(hFunnelSe_F, kOrange + 7);

    for (auto* h : {hFunnelMu, hFunnelSe, hFunnelMu_F, hFunnelSe_F}) {
        h->GetXaxis()->SetBinLabel(1, "1. Generated");
        h->GetXaxis()->SetBinLabel(2, "2. Reach Wall");
        h->GetXaxis()->SetBinLabel(3, "3. Extracted (A#rightarrowC)");
        h->GetXaxis()->SetBinLabel(4, "4. Detected");
    }

    TH1D* hYield_PEMeV = new TH1D("hYield_PEMeV", "Photoelectron Yield per Energy Deposit;N_{pe} / E_{dep} (p.e./MeV);Events", 100, 0, 300); ApplyStyles1D(hYield_PEMeV, kCyan + 2);
    TH1D* hEdepPerSecondary = new TH1D("hEdepPerSecondary", "Mean Energy Deposit per Secondary Particle (Core Volume);E_{dep} / Particle (MeV);Events", 100, 0, 5); ApplyStyles1D(hEdepPerSecondary, kPink + 2);

    // =========================================================================
    // GROUP C: DETECTED SIGNAL & SPECTRA
    // =========================================================================
    TH1D* hWvl_MuWater = new TH1D("hWvl_MuWater", "Detected Spectrum -- Primary Source, Core Origin;Wavelength #lambda (nm);Photon Count", 100, 200, 800); ApplyStyles1D(hWvl_MuWater, kBlue);
    TH1D* hWvl_SeWater = new TH1D("hWvl_SeWater", "Detected Spectrum -- Secondary Source, Core Origin;Wavelength #lambda (nm);Photon Count", 100, 200, 800); ApplyStyles1D(hWvl_SeWater, kCyan+2);
    TH1D* hWvl_MuAcr   = new TH1D("hWvl_MuAcr",   "Detected Spectrum -- Primary Source, Acrylic Origin" + acrSuffix + ";Wavelength #lambda (nm);Photon Count", 100, 200, 800); ApplyStyles1D(hWvl_MuAcr, kRed);
    TH1D* hWvl_SeAcr   = new TH1D("hWvl_SeAcr",   "Detected Spectrum -- Secondary Source, Acrylic Origin" + acrSuffix + ";Wavelength #lambda (nm);Photon Count", 100, 200, 800); ApplyStyles1D(hWvl_SeAcr, kOrange);
    TH1D* hWvl_Gel     = new TH1D("hWvl_Gel",     "Detected Spectrum -- Origin: Coupling Gel;Wavelength #lambda (nm);Photon Count", 100, 200, 800); ApplyStyles1D(hWvl_Gel, kGreen+2);

    TH1D* hWvlGen_MuWater = new TH1D("hWvlGen_MuWater", "Generated Spectrum -- Primary Source, Core Origin;Wavelength #lambda (nm);Photon Count", 100, 200, 800); ApplyStyles1D(hWvlGen_MuWater, kGray + 2);
    TH1D* hAIn_Mu = new TH1D("hAIn_Mu", "SiPM Surface Incident Angle -- Primary Source;Incident Angle #alpha_{in} (deg);Photon Count", 90, 0, 90); ApplyStyles1D(hAIn_Mu, kBlue);
    TH1D* hAIn_Se = new TH1D("hAIn_Se", "SiPM Surface Incident Angle -- Secondary Source;Incident Angle #alpha_{in} (deg);Photon Count", 90, 0, 90); ApplyStyles1D(hAIn_Se, kRed);

    TH1D* hAIn_Pct = new TH1D("hAIn_Pct", "SiPM Incident Angle Distribution (All Detected Photons);Incident Angle #alpha_{in} (deg);Fraction of Detected Photons (%)", 90, 0, 90); ApplyStyles1D(hAIn_Pct, kViolet + 1);

    TH1D* hTime_Mu_Dir = new TH1D("hTM_D", "Primary Source -- Direct Path;Arrival Time t (ns);Photon Count", 150, 0, 10); ApplyStyles1D(hTime_Mu_Dir, kBlue);
    TH1D* hTime_Mu_Bnc = new TH1D("hTM_B", "Primary Source -- Reflected Path;Arrival Time t (ns);Photon Count", 150, 0, 10); ApplyStyles1D(hTime_Mu_Bnc, kTeal+2, 0, 2); hTime_Mu_Bnc->SetLineStyle(2);
    TH1D* hTime_Se_Dir = new TH1D("hTD_D", "Secondary Source -- Direct Path;Arrival Time t (ns);Photon Count", 150, 0, 10); ApplyStyles1D(hTime_Se_Dir, kRed);
    TH1D* hTime_Se_Bnc = new TH1D("hTD_B", "Secondary Source -- Reflected Path;Arrival Time t (ns);Photon Count", 150, 0, 10); ApplyStyles1D(hTime_Se_Bnc, kOrange+1, 0, 2); hTime_Se_Bnc->SetLineStyle(2);
    TH1D* hTime_Tot     = new TH1D("hTime_Tot", "Detected Signal Timing Breakdown;Arrival Time t (ns);Photon Count", 150, 0, 10); ApplyStyles1D(hTime_Tot, kBlack, 3004, 3);

    TH2D* hWlVsAngle = new TH2D("hWlVsAngle", "Wavelength vs. SiPM Incident Angle;Incident Angle #alpha_{in} (deg);Wavelength #lambda (nm)", 90, 0, 90, 100, 200, 800); ApplyStyles2D(hWlVsAngle);

    // =========================================================================
    // GROUP D: HARDWARE METRICS
    // =========================================================================
    TH1D* hOccupancy = new TH1D("hOccupancy", "Mean Channel Occupancy (24 SiPMs);SiPM Channel Index;Mean Hits / Event", 24, -0.5, 23.5); ApplyStyles1D(hOccupancy, kSpring + 2);
    TH1D* hOccSpread = new TH1D("hOccSpread", "Single-Channel Hit Distribution Spread;Hits / Event / Channel;Occurrences", 40, 0, 40); ApplyStyles1D(hOccSpread, kOrange + 7);
    TH1D* hMultiplicity = new TH1D("hMultiplicity", "SiPM Channel Multiplicity per Event;Fired SiPM Channel Count;Events", 25, -0.5, 24.5); ApplyStyles1D(hMultiplicity, kMagenta + 1);

    TProfile* hAsymX = new TProfile("hAsymX", "Light Collection Asymmetry vs. Primary Impact X;Impact Position X (mm);Asymmetry A_{X} = (N_{X+}-N_{X-})/(N_{X+}+N_{X-})", 25, -waterHalfXZ, waterHalfXZ); hAsymX->SetMarkerStyle(20); hAsymX->SetMarkerColor(kBlue); hAsymX->SetLineColor(kBlue);
    TProfile* hAsymZ = new TProfile("hAsymZ", "Light Collection Asymmetry vs. Primary Impact Z;Impact Position Z (mm);Asymmetry A_{Z} = (N_{Z+}-N_{Z-})/(N_{Z+}+N_{Z-})", 25, -waterHalfXZ, waterHalfXZ); hAsymZ->SetMarkerStyle(21); hAsymZ->SetMarkerColor(kGreen + 2); hAsymZ->SetLineColor(kGreen + 2);

    TH1D* hTIRonly = new TH1D("hTIRonly", "Total Internal Reflections (TIR) per Photon;Number of TIR Bounces;Detected Photon Count", 60, 0, 60); ApplyStyles1D(hTIRonly, kBlue);
    TH1D* hTYVonly = new TH1D("hTYVonly", "Tyvek Diffuse Reflections per Photon;Number of Tyvek Reflections;Detected Photon Count", 60, 0, 60); ApplyStyles1D(hTYVonly, kOrange + 1);

    TH1D* hTIRLayer = new TH1D("hTIRLayer", "Reflection Occurrences by Geometric Interface;Interface Layer;Detected Photons", 4, 0, 4); ApplyStyles1D(hTIRLayer, kPink + 1);
    // [MOI] bit0 mang 2 y nghia tuy fUseAcrylicShell cua run nay (xem TreeInfo.md Tree4):
    // co vo Acrylic -> "Core-Acrylic" (TIR that); khong co vo -> "Core-Gel" (TIR that, truc tiep).
    hTIRLayer->GetXaxis()->SetBinLabel(1, hasAcrylicShell ? "Core-Acrylic (TIR)" : "Core-Gel direct (TIR)");
    hTIRLayer->GetXaxis()->SetBinLabel(2, "Core-Tyvek (diffuse)");
    hTIRLayer->GetXaxis()->SetBinLabel(3, hasAcrylicShell ? "Acrylic-Gel" : "N/A (no shell)");
    // [SUA NHAN] bit3 KHONG phai TIR tai Gel-Sensor (bien dielectric_metal khong the
    // sinh TotalInternalReflection that su) - ban chat la TIR NGUOC tu Gel ve lai
    // Acrylic/Core. Ten cot ROOT "TIR_Layer_Mask" bit3 giu nguyen vi tuong thich nguoc,
    // nhung nhan hien thi o day sua lai cho dung vat ly. Xem TreeInfo.md muc Tree 4.
    hTIRLayer->GetXaxis()->SetBinLabel(4, "Gel#rightarrowAcrylic/Core (TIR back)");

    TH1D* hAngleAnomaly = new TH1D("hAngleAnomaly", "SiPM Incident Angle Consistency QA;Incident Angle #alpha_{in} (deg);Photon Count", 100, 0, 180); ApplyStyles1D(hAngleAnomaly, kGray + 1);

    // =========================================================================
    // GROUP E & F: 2D SPECTRA & ORIGAMI
    // =========================================================================
    TH2D* hWlVsTime = new TH2D("hWlVsTime", "Optical Dispersion: Wavelength vs. Arrival Time;Arrival Time t (ns);Wavelength #lambda (nm)", 150, 0, 25, 100, 200, 800); ApplyStyles2D(hWlVsTime);
    TH2D* hTrap = new TH2D("hTrap", "Waveform Tail QA: Tyvek Reflections vs. Arrival Time;Arrival Time t (ns);Tyvek Reflection Count", 150, 0, 25, 60, 0, 60); ApplyStyles2D(hTrap);

    double panW = 2.0 * waterHalfXZ;
    TH2D* hPanorama = new TH2D("hPanorama", "Unfolded Side-Wall Hit Map (Primary Source);Perimeter Position U (mm);Vertical Position V (mm)", 400, 0, 4 * panW, 100, -panW/2, panW/2); ApplyStyles2D(hPanorama);

    TH2D* hTruthOrigami[6]; TH2D* hDetOrigami[6]; TH2D* hEarlyMuOrigami[6]; TH2D* hEarlyElOrigami[6];
    TString fN[6] = {"1. Bottom (-Y)", "2. Front (+Z)", "3. Right (+X)", "4. Back (-Z)", "5. Left (-X)", "6. Top (+Y)"};
    for (int i = 0; i < 6; i++) {
        hTruthOrigami[i]   = new TH2D(Form("hT%d", i), fN[i] + ";Local U (mm);Local V (mm)", 100, -waterHalfXZ, waterHalfXZ, 100, -waterHalfXZ, waterHalfXZ);
        hDetOrigami[i]     = new TH2D(Form("hD%d", i), fN[i] + ";Local U (mm);Local V (mm)", 100, -waterHalfXZ, waterHalfXZ, 100, -waterHalfXZ, waterHalfXZ);
        hEarlyMuOrigami[i] = new TH2D(Form("hEMu%d", i), fN[i] + ";Local U (mm);Local V (mm)", 100, -waterHalfXZ, waterHalfXZ, 100, -waterHalfXZ, waterHalfXZ);
        hEarlyElOrigami[i] = new TH2D(Form("hEEl%d", i), fN[i] + ";Local U (mm);Local V (mm)", 100, -waterHalfXZ, waterHalfXZ, 100, -waterHalfXZ, waterHalfXZ);
    }

    // =========================================================================
    // PASS 0 / Tree6 & Tree0 (Particle & Track PDG Mapping)
    // =========================================================================
    auto PrimaryKey = [](int e, int t) -> long long { return (long long)e * 100000LL + t; };
    std::unordered_set<long long> isPrimaryTrack;
    std::unordered_map<int, int> numPrimaryPerEvt;
    std::unordered_map<long long, double> PrimaryEnergyMap;
    std::unordered_map<long long, int> trackPDGMap;

    long long nT6 = t6->GetEntries();
    int eid6, trk6, pdg6; double vE6, enX6, enY6, enZ6, exX6, exY6, exZ6;
    t6->SetBranchAddress("EventID", &eid6); t6->SetBranchAddress("TrackID", &trk6);
    t6->SetBranchAddress("PDG_Code", &pdg6); t6->SetBranchAddress("VertexE_MeV", &vE6);
    t6->SetBranchAddress("EntryX_mm", &enX6); t6->SetBranchAddress("EntryY_mm", &enY6); t6->SetBranchAddress("EntryZ_mm", &enZ6);
    t6->SetBranchAddress("ExitX_mm", &exX6); t6->SetBranchAddress("ExitY_mm", &exY6); t6->SetBranchAddress("ExitZ_mm", &exZ6);

    for (long long i = 0; i < nT6; i++) {
        t6->GetEntry(i);
        long long key6 = PrimaryKey(eid6, trk6);
        isPrimaryTrack.insert(key6);
        PrimaryEnergyMap[key6] = vE6;
        trackPDGMap[key6] = pdg6;
        numPrimaryPerEvt[eid6]++;
    }
    for (auto& kv : numPrimaryPerEvt) hNumPrimaryPerEvt->Fill(kv.second);

    long long nT0 = t0->GetEntries();
    int eid0, trk0, pdg0;
    t0->SetBranchAddress("EventID", &eid0); t0->SetBranchAddress("TrackID", &trk0); t0->SetBranchAddress("PDG_Code", &pdg0);
    for (long long i = 0; i < nT0; i++) {
        t0->GetEntry(i);
        trackPDGMap[PrimaryKey(eid0, trk0)] = pdg0;
        int absP = std::abs(pdg0);
        if      (absP == 11) hSec_PDG->Fill(0.5);
        else if (absP == 22) hSec_PDG->Fill(1.5);
        else if (absP == 13) hSec_PDG->Fill(2.5);
        else if (absP == 211) hSec_PDG->Fill(3.5);
        else                 hSec_PDG->Fill(4.5);
    }

    // =========================================================================
    // PASS 2 / Tree3 (EventSummary)
    // =========================================================================
    long long nT3raw = t3->GetEntries();
    long long nT3 = (maxEvents > 0 && maxEvents < nT3raw) ? maxEvents : nT3raw;
    std::unordered_map<int, double> topTimeMap;

    int eid3 = -1, reachW = 0, eA = 0, dTot = 0, hitSiPM[24] = {0}, numMu3 = 0, secCountCore3 = 0;
    double eM = 0, lM = 0, eSec = 0, tTop = 0, firstHitT = -1, hitX3 = 0, hitZ3 = 0;
    t3->SetBranchAddress("EventID", &eid3);
    // [v6] "Edep_Primary_Water_MeV"/"Len_Primary_Water_mm"/"SecondaryCountWater" -> doi ten thanh
    // "..._Core_..."/"SecondaryCountCore" (xem TreeInfo.md muc "Nguyen tac thiet ke").
    t3->SetBranchAddress("Edep_Primary_Core_MeV", &eM); t3->SetBranchAddress("Len_Primary_Core_mm", &lM);
    t3->SetBranchAddress("Edep_Secondary_MeV", &eSec);
    t3->SetBranchAddress("ReachWall", &reachW); t3->SetBranchAddress("ExtractA2C", &eA);
    t3->SetBranchAddress("Detected_Total", &dTot);
    t3->SetBranchAddress("Time_TopScint_ns", &tTop); t3->SetBranchAddress("FirstHitTime_ns", &firstHitT);
    t3->SetBranchAddress("Primary_HitX_mm", &hitX3); t3->SetBranchAddress("Primary_HitZ_mm", &hitZ3);
    t3->SetBranchAddress("NumPrimaryTracks", &numMu3);
    t3->SetBranchAddress("SecondaryCountCore", &secCountCore3);
    for (int c = 0; c < 24; c++) t3->SetBranchAddress(Form("Hit_SiPM_%d", c), &hitSiPM[c]);

    double occSum[24] = {0};

    for (long long i = 0; i < nT3; i++) {
        t3->GetEntry(i);
        topTimeMap[eid3] = tTop;

        // Luu y: gia dinh hat so cap DUY NHAT co TrackID==1 khi numMu3==1 (dung trong da so
        // truong hop vi Geant4 danh so track sinh dau tien = 1, nhung khong duoc TreeInfo.md
        // dam bao tuong minh) - neu can chac chan 100%, nen lay TrackID truc tiep tu Tree6
        // thay vi hardcode "1" o day.
        long long keyMu = PrimaryKey(eid3, 1);
        bool isMIP = (PrimaryEnergyMap.count(keyMu) && PrimaryEnergyMap[keyMu] >= E_MIN_MIP);

        if (numMu3 == 1 && lM > 0 && isMIP) {
            hdEdx_MuWaterTotal->Fill((eM + eSec) / (lM / 10.0));
        }

        if ((eM + eSec) > 0) hYield_PEMeV->Fill((double)dTot / (eM + eSec));
        if (secCountCore3 > 0) hEdepPerSecondary->Fill(eSec / secCountCore3);
        if (dTot > 0) hTTS->Fill(firstHitT - tTop);

        int nFiring = 0;
        int sumPX=0, sumNX=0, sumPZ=0, sumNZ=0;
        for (int c = 0; c < 24; c++) {
            occSum[c] += hitSiPM[c];
            hOccSpread->Fill(hitSiPM[c]);
            if (hitSiPM[c] > 0) nFiring++;
            if (c>=0 && c<=3) sumPX += hitSiPM[c];
            if (c>=4 && c<=7) sumNX += hitSiPM[c];
            if (c>=8 && c<=11) sumPZ += hitSiPM[c];
            if (c>=12 && c<=15) sumNZ += hitSiPM[c];
        }
        hMultiplicity->Fill(nFiring);
        if (sumPX+sumNX>0) hAsymX->Fill(hitX3, (double)(sumPX-sumNX)/(sumPX+sumNX));
        if (sumPZ+sumNZ>0) hAsymZ->Fill(hitZ3, (double)(sumPZ-sumNZ)/(sumPZ+sumNZ));
    }
    double nEvents = (double)nT3;
    if (nEvents == 0) { std::cerr << "No triggered events found!\n"; return; }
    for (int c = 0; c < 24; c++) hOccupancy->SetBinContent(c+1, occSum[c]/nEvents);

    // =========================================================================
    // PASS 1 / Tree1 (PhotonGen)
    // =========================================================================
    long long nT1 = t1->GetEntries();
    int eid1, trk1, par1, cPdg1, orig1;
    double tAng1, wvl1;
    t1->SetBranchAddress("EventID", &eid1); t1->SetBranchAddress("TrackID", &trk1);
    t1->SetBranchAddress("ParentID", &par1); t1->SetBranchAddress("CreatorPDG", &cPdg1);
    t1->SetBranchAddress("OriginVol", &orig1); t1->SetBranchAddress("TrueAngle_deg", &tAng1);
    t1->SetBranchAddress("Wavelength_nm", &wvl1);

    std::unordered_map<int, long long> genMuMap, genSeMap;
    std::unordered_map<int, long long> genMuMap_F, genSeMap_F;
    std::map<std::pair<int,int>, bool> photonIsFromPrimary;
    std::map<std::pair<int,int>, int>  photonOriginVol;
    std::map<std::pair<int,int>, double> photonWvlMap;
    std::unordered_map<long long, long long> phoGenWaterByPrimaryTrack;

    for (long long i = 0; i < nT1; i++) {
        t1->GetEntry(i);
        if (topTimeMap.find(eid1) == topTimeMap.end()) continue;

        bool fromPrimary = isPrimaryTrack.count(PrimaryKey(eid1, par1)) > 0;
        photonIsFromPrimary[{eid1, trk1}] = fromPrimary;
        photonOriginVol[{eid1, trk1}]  = orig1;
        photonWvlMap[{eid1, trk1}]     = wvl1;

        if (fromPrimary) genMuMap[eid1]++; else genSeMap[eid1]++;
        if (wvl1 >= WVL_MIN && wvl1 <= WVL_MAX) {
            if (fromPrimary) genMuMap_F[eid1]++; else genSeMap_F[eid1]++;
        }

        if (fromPrimary && orig1 == 0) phoGenWaterByPrimaryTrack[PrimaryKey(eid1, par1)]++;

        if      (orig1 == 0 &&  fromPrimary) hAng_MuWater->Fill(tAng1);
        else if (orig1 == 0 && !fromPrimary) hAng_SeWater->Fill(tAng1);
        else if (orig1 == 1 &&  fromPrimary) hAng_MuAcr->Fill(tAng1);
        else if (orig1 == 1 && !fromPrimary) hAng_SeAcr->Fill(tAng1);
        else if (orig1 == 3 &&  fromPrimary) hAng_MuGel->Fill(tAng1);
        else if (orig1 == 3 && !fromPrimary) hAng_SeGel->Fill(tAng1);

        if (orig1 == 0 && fromPrimary) hWvlGen_MuWater->Fill(wvl1);
    }

    long long sGenMu = 0, sGenSe = 0;
    for (auto& kv : genMuMap) sGenMu += kv.second;
    for (auto& kv : genSeMap) sGenSe += kv.second;
    hFunnelMu->SetBinContent(1, (double)sGenMu / nEvents);
    hFunnelSe->SetBinContent(1, (double)sGenSe / nEvents);

    long long sGenMu_F = 0, sGenSe_F = 0;
    for (auto& kv : genMuMap_F) sGenMu_F += kv.second;
    for (auto& kv : genSeMap_F) sGenSe_F += kv.second;
    hFunnelMu_F->SetBinContent(1, (double)sGenMu_F / nEvents);
    hFunnelSe_F->SetBinContent(1, (double)sGenSe_F / nEvents);

    // =========================================================================
    // PASS 0.5 / Tree7 (PrimaryEdepPerVolume - Specific PDG Analysis)
    // =========================================================================
    long long nT7 = t7->GetEntries();
    int eid7, trk7, vol7;
    double edep7, len7;
    t7->SetBranchAddress("EventID", &eid7); t7->SetBranchAddress("TrackID", &trk7);
    t7->SetBranchAddress("VolumeID", &vol7);
    t7->SetBranchAddress("Edep_MeV", &edep7); t7->SetBranchAddress("Len_mm", &len7);

    std::unordered_map<int, double> edepPrimaryWaterByEvt, lenPrimaryWaterByEvt;

    for (long long i = 0; i < nT7; i++) {
        t7->GetEntry(i);
        if (topTimeMap.find(eid7) == topTimeMap.end()) continue;

        long long key7 = PrimaryKey(eid7, trk7);
        double vE = PrimaryEnergyMap.count(key7) ? PrimaryEnergyMap[key7] : 0.0;
        int pdg = trackPDGMap.count(key7) ? trackPDGMap[key7] : 0;

        if (vol7 == 0) {
            if (len7 > 0) {
                double dedx = edep7 / (len7 / 10.0);
                hdEdx_MuWater->Fill(dedx);
                hBetheBloch2D->Fill(vE, dedx);

                if (std::abs(pdg) == 13)      hdEdx_TrueMuon->Fill(dedx);
                else if (std::abs(pdg) == 11) hdEdx_PrimaryEl->Fill(dedx);
            }
            edepPrimaryWaterByEvt[eid7] += edep7;
            lenPrimaryWaterByEvt[eid7]  += len7;
            hFT_Occupancy->Fill(len7);

            long long phoGen = phoGenWaterByPrimaryTrack.count(key7) ? phoGenWaterByPrimaryTrack[key7] : 0;
            if (vE >= E_MIN_MIP) {
                hFT->Fill(len7, phoGen);
            }
        }
        else if (vol7 == 1) {
            if (len7 > 0) hdEdx_MuAcrylic->Fill(edep7 / (len7 / 10.0));
        }
    }

    long long PrimaryCacheFail = 0;
    for (long long i = 0; i < nT3; i++) {
        t3->GetEntry(i);
        double edepFromT7 = edepPrimaryWaterByEvt.count(eid3) ? edepPrimaryWaterByEvt[eid3] : 0.0;
        double lenFromT7  = lenPrimaryWaterByEvt.count(eid3)  ? lenPrimaryWaterByEvt[eid3]  : 0.0;
        if (std::abs(eM - edepFromT7) > 1e-6 || std::abs(lM - lenFromT7) > 1e-6) PrimaryCacheFail++;
    }

    // =========================================================================
    // PASS 3 / Tree5 (WallHitmap - Fixed Absolute Time Cut)
    // =========================================================================
    long long reachMu = 0, reachSe = 0;
    long long reachMu_F = 0, reachSe_F = 0;
    long long reachMuDirect = 0, reachSeDirect = 0;
    int eid5, trk5, par5; double x5, y5, z5, time5;
    t5->SetBranchAddress("EventID", &eid5); t5->SetBranchAddress("TrackID", &trk5);
    t5->SetBranchAddress("ParentID", &par5);
    t5->SetBranchAddress("HitX_mm", &x5); t5->SetBranchAddress("HitY_mm", &y5); t5->SetBranchAddress("HitZ_mm", &z5);
    t5->SetBranchAddress("Time_ns", &time5);

    for (long long i = 0; i < t5->GetEntries(); i++) {
        t5->GetEntry(i);
        if (topTimeMap.find(eid5) == topTimeMap.end()) continue;

        FillOrigami(hTruthOrigami, x5, y5, z5);

        // Fixed Time Cut relative to Scintillator Top Trigger
        double dt_absolute = time5 - topTimeMap[eid5];
        if (dt_absolute <= timeCutWall_ns) {
            long long parentKey = PrimaryKey(eid5, par5);
            int parentPDG = trackPDGMap.count(parentKey) ? trackPDGMap[parentKey] : 0;
            if (std::abs(parentPDG) == 13)      FillOrigami(hEarlyMuOrigami, x5, y5, z5);
            else if (std::abs(parentPDG) == 11) FillOrigami(hEarlyElOrigami, x5, y5, z5);
        }

        auto it = photonIsFromPrimary.find({eid5, trk5});
        bool fromPrimary_viaMap = (it != photonIsFromPrimary.end()) ? it->second : false;
        if (fromPrimary_viaMap) reachMu++; else reachSe++;

        auto itW = photonWvlMap.find({eid5, trk5});
        double wvl5 = (itW != photonWvlMap.end()) ? itW->second : -1.0;
        if (wvl5 >= WVL_MIN && wvl5 <= WVL_MAX) {
            if (fromPrimary_viaMap) reachMu_F++; else reachSe_F++;
        }

        bool fromPrimary_direct = isPrimaryTrack.count(PrimaryKey(eid5, par5)) > 0;
        if (fromPrimary_direct) reachMuDirect++; else reachSeDirect++;
    }
    hFunnelMu->SetBinContent(2, (double)reachMu / nEvents);
    hFunnelSe->SetBinContent(2, (double)reachSe / nEvents);

    hFunnelMu_F->SetBinContent(2, (double)reachMu_F / nEvents);
    hFunnelSe_F->SetBinContent(2, (double)reachSe_F / nEvents);

    // =========================================================================
    // PASS 4 / Tree2 (Boundary_A2C)
    // =========================================================================
    long long extrMu = 0, extrSe = 0;
    long long extrMu_F = 0, extrSe_F = 0;
    long long extrMuDirect = 0, extrSeDirect = 0;
    long long pathDirect = 0, pathBounced = 0;
    int eid2, trk2, par2, pathFlag2; double wvl2;
    t2->SetBranchAddress("EventID", &eid2); t2->SetBranchAddress("TrackID", &trk2);
    t2->SetBranchAddress("PathFlag", &pathFlag2); t2->SetBranchAddress("ParentID", &par2);
    t2->SetBranchAddress("Wavelength_nm", &wvl2);

    for (long long i = 0; i < t2->GetEntries(); i++) {
        t2->GetEntry(i);
        if (topTimeMap.find(eid2) == topTimeMap.end()) continue;

        auto it = photonIsFromPrimary.find({eid2, trk2});
        bool fromPrimary_viaMap = (it != photonIsFromPrimary.end()) ? it->second : false;

        if (fromPrimary_viaMap) extrMu++; else extrSe++;

        if (wvl2 >= WVL_MIN && wvl2 <= WVL_MAX) {
            if (fromPrimary_viaMap) extrMu_F++; else extrSe_F++;
        }

        bool fromPrimary_direct = isPrimaryTrack.count(PrimaryKey(eid2, par2)) > 0;
        if (fromPrimary_direct) extrMuDirect++; else extrSeDirect++;

        if (pathFlag2 == 1) pathDirect++; else pathBounced++;
    }
    hFunnelMu->SetBinContent(3, (double)extrMu / nEvents);
    hFunnelSe->SetBinContent(3, (double)extrSe / nEvents);

    hFunnelMu_F->SetBinContent(3, (double)extrMu_F / nEvents);
    hFunnelSe_F->SetBinContent(3, (double)extrSe_F / nEvents);

    // =========================================================================
    // PASS 5 / Tree4 (DetectedSignal)
    // =========================================================================
    long long detMu = 0, detSe = 0;
    long long detMu_F = 0, detSe_F = 0;
    int eid4, trk4, par4, cPdg4, tir4, tyv4, sipmId4, layerMask4;
    double t4ns, x4, y4, z4, wvl4, an4;
    t4->SetBranchAddress("EventID", &eid4); t4->SetBranchAddress("TrackID", &trk4);
    t4->SetBranchAddress("ParentID", &par4); t4->SetBranchAddress("CreatorPDG", &cPdg4);
    t4->SetBranchAddress("GlobalTime_ns", &t4ns);
    t4->SetBranchAddress("Wavelength_nm", &wvl4); t4->SetBranchAddress("Angle_In_deg", &an4);
    t4->SetBranchAddress("HitX_mm", &x4); t4->SetBranchAddress("HitY_mm", &y4); t4->SetBranchAddress("HitZ_mm", &z4);
    t4->SetBranchAddress("SiPM_ID", &sipmId4);
    t4->SetBranchAddress("TIR_Bounces", &tir4); t4->SetBranchAddress("Tyvek_Bounces", &tyv4);
    t4->SetBranchAddress("TIR_Layer_Mask", &layerMask4);

    long long angleAnomalyCount = 0;
    std::unordered_map<int, int> detPrimaryPerEvt, detSecondaryPerEvt;

    for (long long i = 0; i < t4->GetEntries(); i++) {
        t4->GetEntry(i);

        if (topTimeMap.find(eid4) == topTimeMap.end()) continue;
        double trueTime = t4ns - topTimeMap[eid4];

        FillOrigami(hDetOrigami, x4, y4, z4);

        bool fromPrimary = isPrimaryTrack.count(PrimaryKey(eid4, par4)) > 0;
        auto itOv = photonOriginVol.find({eid4, trk4});
        int  ov   = (itOv != photonOriginVol.end()) ? itOv->second : -1;

        if (fromPrimary) { detMu++; detPrimaryPerEvt[eid4]++; }
        else             { detSe++; detSecondaryPerEvt[eid4]++; }

        bool isDir = (tir4 + tyv4 == 0);

        if (wvl4 >= WVL_MIN && wvl4 <= WVL_MAX) {
            if (fromPrimary) detMu_F++; else detSe_F++;
        }

        if      (ov == 0 &&  fromPrimary) hWvl_MuWater->Fill(wvl4);
        else if (ov == 0 && !fromPrimary) hWvl_SeWater->Fill(wvl4);
        else if (ov == 1 &&  fromPrimary) hWvl_MuAcr->Fill(wvl4);
        else if (ov == 1 && !fromPrimary) hWvl_SeAcr->Fill(wvl4);
        else if (ov == 3)               hWvl_Gel->Fill(wvl4);

        if (fromPrimary) hAIn_Mu->Fill(an4); else hAIn_Se->Fill(an4);

        hAIn_Pct->Fill(an4);
        hAngleAnomaly->Fill(an4);
        if (an4 > 90.0) angleAnomalyCount++;
        hWlVsAngle->Fill(an4, wvl4);

        hTime_Tot->Fill(trueTime);
        if (fromPrimary && isDir)        hTime_Mu_Dir->Fill(trueTime);
        else if (fromPrimary && !isDir)  hTime_Mu_Bnc->Fill(trueTime);
        else if (!fromPrimary && isDir)  hTime_Se_Dir->Fill(trueTime);
        else                           hTime_Se_Bnc->Fill(trueTime);

        hWlVsTime->Fill(trueTime, wvl4);
        hTrap->Fill(trueTime, tyv4);
        hTIRonly->Fill(tir4);
        hTYVonly->Fill(tyv4);

        if (layerMask4 & 0x1) hTIRLayer->Fill(0.5);
        if (layerMask4 & 0x2) hTIRLayer->Fill(1.5);
        if (layerMask4 & 0x4) hTIRLayer->Fill(2.5);
        if (layerMask4 & 0x8) hTIRLayer->Fill(3.5);

        if (fromPrimary) {
            int face = sipmId4 / 4;
            double U = -999, V = y4;
            if      (face == 2) U = 0*panW + (x4 + waterHalfXZ);
            else if (face == 0) U = 1*panW + (-z4 + waterHalfXZ);
            else if (face == 3) U = 2*panW + (-x4 + waterHalfXZ);
            else if (face == 1) U = 3*panW + (z4 + waterHalfXZ);
            if (U > -999) hPanorama->Fill(U, V);
        }
    }
    hFunnelMu->SetBinContent(4, (double)detMu / nEvents);
    hFunnelSe->SetBinContent(4, (double)detSe / nEvents);

    hFunnelMu_F->SetBinContent(4, (double)detMu_F / nEvents);
    hFunnelSe_F->SetBinContent(4, (double)detSe_F / nEvents);

    // Particle PID Ratio Calculation (Secondary / Primary Photons per Event)
    for (auto& kv : topTimeMap) {
        int e = kv.first;
        double nPrim = detPrimaryPerEvt.count(e) ? detPrimaryPerEvt[e] : 0.0;
        double nSec  = detSecondaryPerEvt.count(e) ? detSecondaryPerEvt[e] : 0.0;
        if (nPrim > 0) hSecToPrimRatio->Fill(nSec / nPrim);
    }

    // =========================================================================
    // DRAW & EXPORT PDF REPORT
    // =========================================================================

    // --- 0. EXECUTIVE SUMMARY ---
    cOut->Clear();
    TLatex dash; dash.SetNDC(); dash.SetTextFont(62); dash.SetTextSize(0.045); dash.SetTextColor(kTeal+3);
    dash.DrawLatex(0.06, 0.94, "EXECUTIVE SUMMARY - Detector Performance & Particle QA");
    dash.SetTextFont(42); dash.SetTextSize(0.028); dash.SetTextColor(kBlack);
    double meanYield = hYield_PEMeV->GetMean();
    double ttsRMS = hTTS->GetRMS(), ttsMean = hTTS->GetMean();
    double avgMulti = hMultiplicity->GetMean();
    double avgNumPrimary = hNumPrimaryPerEvt->GetMean();
    double meanDEdxWater = hdEdx_MuWater->GetMean();
    double meanDEdxTotal = hdEdx_MuWaterTotal->GetMean();
    double meanDEdxAcrylic = hdEdx_MuAcrylic->GetMean();

    // [MOI] Hien thi vat lieu loi NGAY TREN DAU dashboard - TreeInfo.md yeu cau doc cot nay
    // TRUOC KHI dien giai bat ky so lieu Cherenkov nao (nguong phat/goc Cherenkov phu thuoc vat lieu).
    dash.SetTextFont(62); dash.SetTextColor(kRed+1);
    dash.DrawLatex(0.06, 0.965, Form("Core Material (Tree8.CoreMaterial): %s", coreMaterialStr.Data()));
    dash.SetTextFont(42); dash.SetTextColor(kBlack);

    if (hasRunConfig) {
        dash.DrawLatex(0.08, 0.89, Form("#bullet Geometry: Core(%s) %.0fx%.0fx%.0f mm | Scint %.0fx%.0fx%.0f mm | "
                        "Acrylic Shell: %s (%.2fmm) | Tyvek %.2fmm | SiPM thick %.2fmm",
                        coreMaterialStr.Data(), coreFullX, coreFullY, coreFullZ, scintX, scintY, scintZ,
                        hasAcrylicShell ? "ON" : "OFF", acrThick, tyvThick, sipmZ));
    }
    dash.DrawLatex(0.08, 0.84, Form("#bullet Triggered events (coincidence): %.0f", nEvents));
    dash.DrawLatex(0.08, 0.79, Form("#bullet Mean primary particles/event: %.3f", avgNumPrimary));
    // [MOI] Khi khong co vo Acrylic, hdEdx_MuAcrylic RONG (0 entries) -> GetMean()=0 gay hieu nham
    // la "co Acrylic nhung dE/dx = 0". Hien thi "N/A (no shell)" thay vi so 0.000 trong truong hop nay.
    TString acrDEdxStr = hasAcrylicShell ? Form("%.3f MeV/cm", meanDEdxAcrylic) : TString("N/A (no shell)");
    dash.DrawLatex(0.08, 0.74, Form("#bullet Mean dE/dx -- Core (%s): %.3f MeV/cm | Acrylic: %s",
                    coreMaterialStr.Data(), meanDEdxWater, acrDEdxStr.Data()));
    dash.DrawLatex(0.08, 0.69, Form("#bullet Mean dE/dx Core, Primary + Secondary: %.3f MeV/cm", meanDEdxTotal));
    dash.DrawLatex(0.08, 0.64, Form("#bullet Internal consistency check (Edep/Len): %s", (PrimaryCacheFail==0 ? "MATCH" : "MISMATCH")));
    dash.DrawLatex(0.08, 0.59, Form("#bullet Mean photon yield: %.2f p.e./MeV", meanYield));
    dash.DrawLatex(0.08, 0.54, Form("#bullet First-photon arrival time -- RMS: %.3f ns | Mean: %.3f ns", ttsRMS, ttsMean));
    dash.DrawLatex(0.08, 0.49, Form("#bullet Mean SiPM channel multiplicity: %.2f / 24", avgMulti));
    dash.DrawLatex(0.08, 0.44, Form("#bullet Generation-to-detection efficiency -- Primary: %.3f %% | Secondary: %.3f %%",
                    (sGenMu>0)?(100.0*detMu/sGenMu):0, (sGenSe>0)?(100.0*detSe/sGenSe):0));
    dash.DrawLatex(0.08, 0.39, Form("#bullet Fraction of detected photons from Primary source: %.1f %%", (detMu+detSe>0)?(100.0*detMu/(detMu+detSe)):0));
    dash.DrawLatex(0.08, 0.34, Form("#bullet Early Wall Hitmap Time Window: dt <= %.2f ns (rel. to trigger)", timeCutWall_ns));
    cOut->Print(pdfOut);

    // --- NEW PLOTS: Particle QA & Discriminators ---
    cOut->Clear(); hSec_PDG->Draw("HIST BAR"); cOut->Print(pdfOut);

    cOut->Clear();
    hdEdx_TrueMuon->Draw("HIST");
    hdEdx_PrimaryEl->Draw("HIST SAME");
    TLegend* lPDG = new TLegend(0.55, 0.70, 0.88, 0.85); lPDG->SetBorderSize(0);
    lPDG->AddEntry(hdEdx_TrueMuon, "True Muons (|PDG| == 13)", "f");
    lPDG->AddEntry(hdEdx_PrimaryEl, "Primary Electrons (|PDG| == 11)", "f");
    lPDG->Draw();
    cOut->Print(pdfOut);

    cOut->Clear(); hSecToPrimRatio->Draw("HIST"); cOut->Print(pdfOut);

    // --- A. TRUTH ---
    cOut->Clear(); hdEdx_MuWater->Draw("HIST"); cOut->Print(pdfOut);
    cOut->Clear(); hdEdx_MuWaterTotal->Draw("HIST"); cOut->Print(pdfOut);
    cOut->Clear(); hdEdx_MuAcrylic->Draw("HIST"); cOut->Print(pdfOut);

    cOut->Clear(); cOut->SetRightMargin(0.15); hBetheBloch2D->Draw("COLZ"); cOut->Print(pdfOut); cOut->SetRightMargin(0.10);

    cOut->Clear(); hNumPrimaryPerEvt->Draw("HIST"); cOut->Print(pdfOut);

    cOut->Clear(); cOut->Divide(3, 2);
    cOut->cd(1); hAng_MuWater->Draw("HIST"); cOut->cd(2); hAng_SeWater->Draw("HIST");
    cOut->cd(3); hAng_MuAcr->Draw("HIST");   cOut->cd(4); hAng_SeAcr->Draw("HIST");
    cOut->cd(5); hAng_MuGel->Draw("HIST");   cOut->cd(6); hAng_SeGel->Draw("HIST");
    cOut->Print(pdfOut); cOut->Clear();

    cOut->Clear(); cOut->SetRightMargin(0.15);
    hFT->Draw("COLZ");
    TF1* fitFT = new TF1("fitFT", "[0]*x", 0, 150);
    fitFT->SetLineColor(kRed); fitFT->SetLineWidth(3);
    hFT->Fit(fitFT, "Q", "", 100, 150);
    fitFT->Draw("SAME");

    TLatex ftTex; ftTex.SetNDC(); ftTex.SetTextSize(0.035); ftTex.SetTextColor(kRed+1);
    ftTex.DrawLatex(0.15, 0.82, Form("Linear Fit (MIP Regime): N_{#gamma} = %.2f #times L_{water} [mm^{-1}]", fitFT->GetParameter(0)));
    cOut->Print(pdfOut); cOut->SetRightMargin(0.10);

    // --- B. Funnel per source + yield ---
    DrawFunnelWithPct(cOut, pdfOut, hFunnelMu);
    DrawFunnelWithPct(cOut, pdfOut, hFunnelSe);

    DrawFunnelWithPct(cOut, pdfOut, hFunnelMu_F);
    DrawFunnelWithPct(cOut, pdfOut, hFunnelSe_F);

    cOut->Clear(); hYield_PEMeV->Draw("HIST"); cOut->Print(pdfOut);
    cOut->Clear(); hEdepPerSecondary->Draw("HIST"); cOut->Print(pdfOut);

    // --- C. Detected spectra ---
    cOut->Clear(); cOut->Divide(2, 2);
    cOut->cd(1); hWvl_MuWater->Draw("HIST"); cOut->cd(2); hWvl_SeWater->Draw("HIST");
    cOut->cd(3); hWvl_MuAcr->Draw("HIST");   cOut->cd(4); hWvl_SeAcr->Draw("HIST");
    cOut->Print(pdfOut); cOut->Clear();
    cOut->Clear(); hWvl_Gel->Draw("HIST"); cOut->Print(pdfOut);

    cOut->Clear(); hAIn_Mu->Draw("HIST"); hAIn_Se->Draw("HIST SAME");
    TLegend* lAI = new TLegend(0.6,0.75,0.88,0.88); lAI->SetBorderSize(0);
    lAI->AddEntry(hAIn_Mu,"Primary Source","f"); lAI->AddEntry(hAIn_Se,"Secondary Source","f"); lAI->Draw();
    cOut->Print(pdfOut);

    {
        double totalDetected = hAIn_Pct->Integral();
        if (totalDetected > 0) {
            for (int b = 1; b <= hAIn_Pct->GetNbinsX(); b++) {
                double raw = hAIn_Pct->GetBinContent(b);
                hAIn_Pct->SetBinContent(b, 100.0 * raw / totalDetected);
            }
        }
        hAIn_Pct->SetMinimum(0.0);
        cOut->Clear(); hAIn_Pct->Draw("HIST");
        TLatex totLbl; totLbl.SetNDC(); totLbl.SetTextSize(0.03);
        totLbl.DrawLatex(0.55, 0.85, Form("Total Detected Photons (100%%) = %.0f", totalDetected));
        cOut->Print(pdfOut);

        TH1D* hAIn_PctCum = (TH1D*)hAIn_Pct->Clone("hAIn_PctCum");
        hAIn_PctCum->SetTitle("Cumulative Distribution of SiPM Incident Angle;Incident Angle #alpha_{in} (deg);Cumulative Fraction of Detected Photons (%)");
        double running = 0.0;
        for (int b = 1; b <= hAIn_PctCum->GetNbinsX(); b++) {
            running += hAIn_Pct->GetBinContent(b);
            hAIn_PctCum->SetBinContent(b, running);
        }
        hAIn_PctCum->SetMinimum(0.0);
        hAIn_PctCum->SetMaximum(105.0);
        ApplyStyles1D(hAIn_PctCum, kGreen + 2);
        cOut->Clear(); hAIn_PctCum->Draw("HIST");
        { TLine* ln100 = new TLine(0,100,90,100); ln100->SetLineColor(kGray+1); ln100->SetLineStyle(2); ln100->Draw(); }
        cOut->Print(pdfOut);
    }

    cOut->Clear(); hTime_Tot->SetMaximum(hTime_Tot->GetMaximum()*1.4);
    hTime_Tot->Draw("HIST"); hTime_Mu_Dir->Draw("HIST SAME"); hTime_Mu_Bnc->Draw("HIST SAME");
    hTime_Se_Dir->Draw("HIST SAME"); hTime_Se_Bnc->Draw("HIST SAME");
    TLegend* lT = new TLegend(0.55,0.6,0.88,0.88); lT->SetBorderSize(0);
    lT->AddEntry(hTime_Tot,"Total Signal","l"); lT->AddEntry(hTime_Mu_Dir,"Primary Photon -- Direct","f");
    lT->AddEntry(hTime_Mu_Bnc,"Primary Photon -- Reflected","l"); lT->AddEntry(hTime_Se_Dir,"Secondary Photon -- Direct","f");
    lT->AddEntry(hTime_Se_Bnc,"Secondary Photon -- Reflected","l"); lT->Draw(); cOut->Print(pdfOut);

    cOut->Clear(); cOut->SetRightMargin(0.15); hWlVsAngle->Draw("COLZ"); cOut->Print(pdfOut);
    cOut->SetRightMargin(0.10);

    // --- D. Hardware ---
    cOut->Clear(); hOccupancy->SetMinimum(0); hOccupancy->Draw("HIST"); cOut->Print(pdfOut);
    cOut->Clear(); gPad->SetLogy(1); hOccSpread->Draw("HIST"); cOut->Print(pdfOut); gPad->SetLogy(0);
    cOut->Clear(); hMultiplicity->Draw("HIST"); cOut->Print(pdfOut);
    cOut->Clear(); hAsymX->Draw("PE"); cOut->Print(pdfOut);
    cOut->Clear(); hAsymZ->Draw("PE"); cOut->Print(pdfOut);

    cOut->Clear(); hTIRonly->Draw("HIST"); hTYVonly->Draw("HIST SAME");
    TLegend* lB = new TLegend(0.5,0.75,0.88,0.88); lB->SetBorderSize(0);
    lB->AddEntry(hTIRonly,"Total Internal Reflection (TIR)","f"); lB->AddEntry(hTYVonly,"Tyvek Diffuse Reflection","f");
    lB->Draw(); cOut->Print(pdfOut);

    cOut->Clear(); hTIRLayer->Draw("HIST BAR"); cOut->Print(pdfOut);

    cOut->Clear(); hAngleAnomaly->Draw("HIST");
    { TLine* ln90 = new TLine(90,0,90,hAngleAnomaly->GetMaximum()); ln90->SetLineColor(kRed); ln90->SetLineStyle(2); ln90->SetLineWidth(2); ln90->Draw(); }
    cOut->Print(pdfOut);

    // --- E. 2D & QA ---
    cOut->Clear(); cOut->SetRightMargin(0.15); hWlVsTime->Draw("COLZ"); cOut->Print(pdfOut);
    cOut->Clear(); hTrap->Draw("COLZ"); cOut->Print(pdfOut);
    cOut->SetRightMargin(0.10);

    cOut->Clear(); cOut->SetRightMargin(0.13); hPanorama->Draw("COLZ");
    { for (int k=1;k<=3;k++){ TLine* ln=new TLine(k*panW,-panW/2,k*panW,panW/2); ln->SetLineStyle(2); ln->SetLineColor(kWhite); ln->SetLineWidth(2); ln->Draw(); }
      const char* faceNames[4]={"+Z","+X","-Z","-X"};
      for (int k=0;k<4;k++){ TLatex* tl=new TLatex((k+0.5)*panW,-panW*0.55,Form("Face %s",faceNames[k])); tl->SetTextAlign(23); tl->SetTextSize(0.035); tl->Draw(); } }
    cOut->Print(pdfOut); cOut->SetRightMargin(0.10);

    TH1D* hTransfer = (TH1D*)hWvl_MuWater->Clone("hTransfer");
    hTransfer->Reset();
    hTransfer->SetTitle("Photon Detection Efficiency vs. Wavelength (Primary Source, Core Origin);Wavelength #lambda (nm);Detection Efficiency");
    const double kMinGenForRatio = 10.0;
    for (int b = 1; b <= hWvl_MuWater->GetNbinsX(); b++) {
        double num = hWvl_MuWater->GetBinContent(b);
        double den = hWvlGen_MuWater->GetBinContent(b);
        hTransfer->SetBinContent(b, (den >= kMinGenForRatio) ? (num / den) : 0.0);
    }
    hTransfer->SetMinimum(0.0);
    hTransfer->SetMaximum(1.0);
    ApplyStyles1D(hTransfer, kViolet+1);
    cOut->Clear(); hTransfer->Draw("HIST"); cOut->Print(pdfOut);

    TH1D* hGenOccupancy = (TH1D*)hWvlGen_MuWater->Clone("hGenOccupancy");
    hGenOccupancy->SetTitle("Wavelength of generated photon;Wavelength #lambda (nm);Generated Photon ");
    ApplyStyles1D(hGenOccupancy, kGray + 2);
    cOut->Clear(); gPad->SetLogy(1); hGenOccupancy->Draw("HIST"); cOut->Print(pdfOut); gPad->SetLogy(0);

    // --- F. Origami ---
    DrawOrigamiCanvas(cOut, pdfOut, hTruthOrigami, "Truth HitMap");

    // EARLY TIME-CUT HITMAPS (FIXED ABSOLUTE TIME REFERENCE)
    DrawOrigamiCanvas(cOut, pdfOut, hEarlyMuOrigami, Form("Early Wall Hit Map Muon (dt <= %.2f ns)", timeCutWall_ns));
    DrawOrigamiCanvas(cOut, pdfOut, hEarlyElOrigami, Form("Early Wall Hit Map Electron (dt <= %.2f ns)", timeCutWall_ns));

    DrawOrigamiCanvas(cOut, pdfOut, hDetOrigami,   "Detected Photon Hit Map ");

    cOut->Print(pdfOut + "]");
    fOut->Write(); fOut->Close(); fIn->Close();

    std::cout << "\n>>> PDF report written: " << pdfOut << "\n";
    std::cout << ">>> Histogram file written: " << rootOut << "\n\n";
}