#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TSystem.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>

struct StoppedMuInfo {
    double vertexKE_MeV = 0.0;
    double vertexTime_ns = 0.0;
    double properTime_ns = 0.0;
    double endTime_ns = 0.0;
    double endY_mm = 0.0;
    double totalDeltaKE_MeV = 0.0;

    bool coreDecay = false;
    bool stopped = false;
};

void validate_detector_muon_lifetime_stopped(
    const char* inputName="validation_detector_muon_lifetime_run0.root")
{
    constexpr double tauRef_us = 2.1969811;

    // VertexKE is kinetic energy, therefore the Michel positron
    // kinetic-energy endpoint is ~52.32 MeV, not 52.8 MeV.
    constexpr double michelKEMax_MeV = 52.32;

    TFile f(inputName,"READ");

    auto* tr = dynamic_cast<TTree*>(f.Get("Track"));
    auto* tv = dynamic_cast<TTree*>(f.Get("TrackVolume"));
    auto* vm = dynamic_cast<TTree*>(f.Get("VolumeMap"));

    if (!tr || !tv || !vm) {
        std::cerr << "Missing Track, TrackVolume, or VolumeMap\n";
        return;
    }

    // ------------------------------------------------------------
    // Resolve Core VolumeID.
    // ------------------------------------------------------------
    int volumeID = -1;
    char volumeName[256] = {};

    vm->SetBranchAddress("VolumeID",&volumeID);
    vm->SetBranchAddress("VolumeName",volumeName);

    int coreID = -1;

    for (Long64_t i=0;i<vm->GetEntries();++i) {
        vm->GetEntry(i);

        if (std::string(volumeName) == "Core") {
            coreID = volumeID;
            break;
        }
    }

    if (coreID < 0) {
        std::cerr << "Could not resolve Core VolumeID\n";
        return;
    }

    // ------------------------------------------------------------
    // Track branches.
    // ------------------------------------------------------------
    int event=0, track=0, parent=0, pdg=0;
    int endVolumeID=-1;

    float vertexKE=0;
    float vertexTime=0;
    float properNs=0;
    float endTime=0;
    float endY=0;

    char creatorProc[256] = {};
    char endProc[256] = {};

    tr->SetBranchAddress("EventID",&event);
    tr->SetBranchAddress("TrackID",&track);
    tr->SetBranchAddress("ParentID",&parent);
    tr->SetBranchAddress("PDG",&pdg);

    tr->SetBranchAddress("VertexKE_MeV",&vertexKE);
    tr->SetBranchAddress("VertexTime_ns",&vertexTime);

    tr->SetBranchAddress("EndProperTime_ns",&properNs);
    tr->SetBranchAddress("EndTime_ns",&endTime);
    tr->SetBranchAddress("EndY_mm",&endY);

    tr->SetBranchAddress("EndVolumeID",&endVolumeID);
    tr->SetBranchAddress("CreatorProcessName",creatorProc);
    tr->SetBranchAddress("EndProcessName",endProc);

    std::map<std::pair<int,int>,StoppedMuInfo> muons;

    // Pass 1: primary mu+ truth.
    for (Long64_t i=0;i<tr->GetEntries();++i) {
        tr->GetEntry(i);

        if (parent != 0 || pdg != -13)
            continue;

        StoppedMuInfo m;

        m.vertexKE_MeV = vertexKE;
        m.vertexTime_ns = vertexTime;
        m.properTime_ns = properNs;
        m.endTime_ns = endTime;
        m.endY_mm = endY;

        m.coreDecay =
            (endVolumeID == coreID &&
             std::string(endProc).find("Decay") != std::string::npos);

        muons[{event,track}] = m;
    }

    // ------------------------------------------------------------
    // Sum primary-muon DeltaKE over ALL volumes.
    // ------------------------------------------------------------
    int tvEvent=0, tvTrack=0;
    float tvDKE=0;

    tv->SetBranchAddress("EventID",&tvEvent);
    tv->SetBranchAddress("TrackID",&tvTrack);
    tv->SetBranchAddress("DeltaKE_MeV",&tvDKE);

    for (Long64_t i=0;i<tv->GetEntries();++i) {
        tv->GetEntry(i);

        auto it = muons.find({tvEvent,tvTrack});
        if (it == muons.end())
            continue;

        it->second.totalDeltaKE_MeV += tvDKE;
    }

    // Stopped definition:
    // primary lost >=99.5% of its initial kinetic energy before decay.
    long long nPrimary = muons.size();
    long long nCoreDecay = 0;
    long long nStoppedCoreDecay = 0;

    for (auto& kv : muons) {
        auto& m = kv.second;

        m.stopped =
            m.vertexKE_MeV > 0.0 &&
            m.totalDeltaKE_MeV >= 0.995*m.vertexKE_MeV;

        if (m.coreDecay)
            ++nCoreDecay;

        if (m.coreDecay && m.stopped)
            ++nStoppedCoreDecay;
    }

    TH1D hProper(
        "hDetectorMuProperStopped",
        "Stopped #mu^{+}: proper lifetime;"
        "Proper decay time (#mus);Decays",
        160,0,12);

    TH1D hDelayed(
        "hDetectorMuDelayedStopped",
        "Stopped #mu^{+} Michel delay;"
        "t_{e^{+} birth}-t_{#mu birth} (#mus);Events",
        160,0,12);

    TH1D hMichel(
        "hDetectorMuMichelStopped",
        "Stopped-#mu^{+} Michel e^{+} kinetic energy;"
        "Kinetic energy (MeV);e^{+}",
        160,0,53);

    TH1D hStopY(
        "hDetectorMuStopYStopped",
        "Stopped #mu^{+} decay position in Water Core;"
        "End y (mm);Decays",
        120,-60,60);

    // Fill stopped-muon lifetime truth.
    for (const auto& kv : muons) {
        const auto& m = kv.second;

        if (!m.coreDecay || !m.stopped)
            continue;

        if (std::isfinite(m.properTime_ns) &&
            m.properTime_ns >= 0)
            hProper.Fill(m.properTime_ns/1000.0);

        if (std::isfinite(m.endY_mm))
            hStopY.Fill(m.endY_mm);
    }

    // ------------------------------------------------------------
    // Direct Michel positrons from STOPPED primary mu+ only.
    // ------------------------------------------------------------
    long long nMichelMatched = 0;
    double maxMichelKE = 0.0;

    for (Long64_t i=0;i<tr->GetEntries();++i) {
        tr->GetEntry(i);

        if (pdg != -11)
            continue;

        auto it = muons.find({event,parent});
        if (it == muons.end())
            continue;

        const auto& mu = it->second;

        if (!mu.coreDecay || !mu.stopped)
            continue;

        if (std::string(creatorProc).find("Decay") == std::string::npos)
            continue;

        ++nMichelMatched;

        hMichel.Fill(vertexKE);

        if (vertexKE > maxMichelKE)
            maxMichelKE = vertexKE;

        const double dt_us =
            (double(vertexTime) - mu.vertexTime_ns)/1000.0;

        if (std::isfinite(dt_us) && dt_us >= 0)
            hDelayed.Fill(dt_us);
    }

    if (hProper.GetEntries() < 1000) {
        std::cerr
            << "Too few stopped Core decays: "
            << hProper.GetEntries() << "\n";
        return;
    }

    // ------------------------------------------------------------
    // Lifetime fits.
    // ------------------------------------------------------------
    TF1 fProper("fProperStopped",
                "[0]*exp(-x/[1])",
                0.1,10.0);

    fProper.SetParameters(
        hProper.GetMaximum(),
        tauRef_us);

    hProper.Fit(
        &fProper,"Q0","",
        0.1,10.0);

    TF1 fDelayed("fDelayedStopped",
                 "[0]*exp(-x/[1])",
                 0.1,10.0);

    fDelayed.SetParameters(
        hDelayed.GetMaximum(),
        tauRef_us);

    hDelayed.Fit(
        &fDelayed,"Q0","",
        0.1,10.0);

    const double tauProper =
        fProper.GetParameter(1);

    const double tauProperErr =
        fProper.GetParError(1);

    const double pullProper =
        tauProperErr > 0
        ? (tauProper-tauRef_us)/tauProperErr
        : std::numeric_limits<double>::quiet_NaN();

    const double tauDelayed =
        fDelayed.GetParameter(1);

    const double tauDelayedErr =
        fDelayed.GetParError(1);

    const double pullDelayed =
        tauDelayedErr > 0
        ? (tauDelayed-tauRef_us)/tauDelayedErr
        : std::numeric_limits<double>::quiet_NaN();

    // ------------------------------------------------------------
    // Validation.
    // ------------------------------------------------------------
    const bool passCore =
        nPrimary > 0 &&
        double(nCoreDecay)/nPrimary > 0.98;

    const bool passStopped =
        nStoppedCoreDecay > 0 &&
        double(nStoppedCoreDecay)/nCoreDecay > 0.99;

    const bool passProper =
        std::fabs(pullProper) < 3.0;

    // Since this branch is kinetic energy, expected endpoint ~52.32 MeV.
    const bool passMichelEndpoint =
        maxMichelKE > 52.0 &&
        maxMichelKE < 52.5;

    const bool passDelayed =
        std::fabs(pullDelayed) < 3.0;

    const bool passMatching =
        nMichelMatched == nStoppedCoreDecay;

    const bool allPass =
        passCore &&
        passStopped &&
        passProper &&
        passMichelEndpoint &&
        passDelayed &&
        passMatching;

    std::cout
        << "\nSTOPPED-MUON DETECTOR LIFETIME VALIDATION\n";

    std::cout
        << "Input file                  : "
        << inputName << "\n";

    std::cout
        << "Primary mu+                 : "
        << nPrimary << "\n";

    std::cout
        << "Core decays                 : "
        << nCoreDecay << "\n";

    std::cout
        << "Stopped Core decays         : "
        << nStoppedCoreDecay << "\n";

    std::cout
        << "Decay-in-flight Core decays : "
        << (nCoreDecay-nStoppedCoreDecay)
        << "\n";

    std::cout
        << "Matched stopped Michel e+   : "
        << nMichelMatched << "\n";

    std::cout
        << "Reference tau               : "
        << tauRef_us << " us\n";

    std::cout << std::setprecision(8)
        << "Stopped proper-time tau     : "
        << tauProper
        << " +/- "
        << tauProperErr
        << " us\n"

        << "Stopped proper-time pull    : "
        << pullProper << "\n"

        << "Stopped Michel-delay tau    : "
        << tauDelayed
        << " +/- "
        << tauDelayedErr
        << " us\n"

        << "Stopped Michel-delay pull   : "
        << pullDelayed << "\n"

        << "Max stopped Michel e+ KE    : "
        << maxMichelKE
        << " MeV\n";

    std::cout
        << "Core decay fraction         : "
        << (passCore ? "PASS":"FAIL") << "\n";

    std::cout
        << "Stopped-muon selection      : "
        << (passStopped ? "PASS":"FAIL") << "\n";

    std::cout
        << "Proper lifetime             : "
        << (passProper ? "PASS":"FAIL") << "\n";

    std::cout
        << "Michel endpoint (~52.32)    : "
        << (passMichelEndpoint ? "PASS":"FAIL")
        << "\n";

    std::cout
        << "Michel delayed timing       : "
        << (passDelayed ? "PASS":"FAIL") << "\n";

    std::cout
        << "Michel ancestry/matching    : "
        << (passMatching ? "PASS":"FAIL")
        << "\n";

    std::cout
        << "STOPPED-MUON LIFETIME VALIDATION: "
        << (allPass ? "PASS":"FAIL")
        << "\n";

    gSystem->mkdir("validation_plots",kTRUE);

    TCanvas c(
        "cStoppedMuonLifetime",
        "Stopped muon lifetime validation",
        1300,900);

    c.Divide(2,2);

    c.cd(1);
    hProper.Draw("HIST");
    fProper.Draw("SAME");

    c.cd(2);
    hDelayed.Draw("HIST");
    fDelayed.Draw("SAME");

    c.cd(3);
    hMichel.Draw("HIST");

    TLine endpoint(
        michelKEMax_MeV,
        0,
        michelKEMax_MeV,
        hMichel.GetMaximum()*1.05);

    endpoint.SetLineStyle(2);
    endpoint.Draw("SAME");

    c.cd(4);
    hStopY.Draw("HIST");

    c.SaveAs(
        "validation_plots/"
        "detector_muon_lifetime_stopped.png");

    c.SaveAs(
        "validation_plots/"
        "detector_muon_lifetime_stopped.pdf");
}
