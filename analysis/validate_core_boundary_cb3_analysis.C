#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {
using Key = std::pair<int,int>; // (EventID, PhotonID)

const char* FaceName(int face) {
    switch (face) {
        case 0: return "+X";
        case 1: return "-X";
        case 2: return "+Z";
        case 3: return "-Z";
        case 4: return "+Y";
        case 5: return "-Y";
        default: return "?";
    }
}

struct PhotonInfo {
    int parent = -1;
    int source = -1;
    int creatorType = -1;
    int creatorSubType = -1;
    int originVolume = -1;
    float wavelength = NAN;
};

struct CrossingInfo {
    int face = -1;
    int parent = -1;
    int source = -1;
    int root = -1;
    int originVolume = -1;
    float u = NAN;
    float v = NAN;
    float cosTheta = NAN;
    float wavelength = NAN;
    float time = NAN;
};

double pct(long long a, long long b) {
    return b ? 100.0 * double(a) / double(b) : 0.0;
}
}

void validate_core_boundary_cb3_analysis(
    const char* inputName="validation_core_boundary_cb2_phase_100_run0.root")
{
    TFile f(inputName, "READ");
    if (f.IsZombie()) {
        std::cerr << "Cannot open " << inputName << "\n";
        return;
    }

    auto* cb = dynamic_cast<TTree*>(f.Get("CoreBoundaryCrossing"));
    auto* ph = dynamic_cast<TTree*>(f.Get("Photon"));
    auto* si = dynamic_cast<TTree*>(f.Get("SensorIncident"));
    auto* ds = dynamic_cast<TTree*>(f.Get("DetectedSignal"));
    auto* vm = dynamic_cast<TTree*>(f.Get("VolumeMap"));

    if (!cb || !ph || !si || !ds || !vm) {
        std::cerr << "Missing required tree(s)\n";
        return;
    }

    // ---------------- Volume map ----------------
    int vmID = -1;
    char vmName[256] = {};
    vm->SetBranchAddress("VolumeID", &vmID);
    vm->SetBranchAddress("VolumeName", vmName);

    std::map<int,std::string> volumeNames;
    int coreVolumeID = -1;

    for (Long64_t i=0; i<vm->GetEntries(); ++i) {
        vm->GetEntry(i);
        volumeNames[vmID] = vmName;
        if (std::string(vmName) == "Core")
            coreVolumeID = vmID;
    }

    if (coreVolumeID < 0) {
        std::cerr << "Could not resolve Core VolumeID\n";
        return;
    }

    // ---------------- Photon lookup ----------------
    int pEvent=0, pPhoton=0, pParent=0, pSource=0;
    int pCreatorType=0, pCreatorSubType=0, pOrigin=0;
    float pWavelength=0;

    ph->SetBranchAddress("EventID", &pEvent);
    ph->SetBranchAddress("PhotonID", &pPhoton);
    ph->SetBranchAddress("ParentTrackID", &pParent);
    ph->SetBranchAddress("SourceTrackID", &pSource);
    ph->SetBranchAddress("CreatorProcessType", &pCreatorType);
    ph->SetBranchAddress("CreatorProcessSubType", &pCreatorSubType);
    ph->SetBranchAddress("OriginVolumeID", &pOrigin);
    ph->SetBranchAddress("Wavelength_nm", &pWavelength);

    std::map<Key,PhotonInfo> photonInfo;
    std::set<Key> allPhotons;
    std::set<Key> coreOriginPhotons;

    for (Long64_t i=0; i<ph->GetEntries(); ++i) {
        ph->GetEntry(i);
        Key k{pEvent,pPhoton};
        photonInfo[k] = {pParent,pSource,pCreatorType,pCreatorSubType,pOrigin,pWavelength};
        allPhotons.insert(k);
        if (pOrigin == coreVolumeID)
            coreOriginPhotons.insert(k);
    }

    // ---------------- CoreBoundaryCrossing ----------------
    int cEvent=0, cPhoton=0, cParent=0, cSource=0, cRoot=0, cOrigin=0, cFace=0;
    float cU=0, cV=0, cCos=0, cWave=0, cTime=0;

    cb->SetBranchAddress("EventID", &cEvent);
    cb->SetBranchAddress("PhotonID", &cPhoton);
    cb->SetBranchAddress("ParentTrackID", &cParent);
    cb->SetBranchAddress("SourceTrackID", &cSource);
    cb->SetBranchAddress("RootPrimaryTrackID", &cRoot);
    cb->SetBranchAddress("OriginVolumeID", &cOrigin);
    cb->SetBranchAddress("FaceID", &cFace);
    cb->SetBranchAddress("LocalU_mm", &cU);
    cb->SetBranchAddress("LocalV_mm", &cV);
    cb->SetBranchAddress("CosThetaNormal", &cCos);
    cb->SetBranchAddress("Wavelength_nm", &cWave);
    cb->SetBranchAddress("Time_ns", &cTime);

    std::set<Key> crossingPhotons;
    std::set<Key> coreOriginCrossingPhotons;
    std::map<Key,CrossingInfo> crossingInfo;
    long long faceCounts[6] = {};
    long long duplicateCrossings = 0;
    long long badFace = 0;
    long long badCos = 0;
    long long ancestryMismatch = 0;
    std::map<int,long long> crossingOriginCounts;
    std::map<std::pair<int,int>,long long> creatorCounts;

    TH2D* hUV[6];
    for (int face=0; face<6; ++face) {
        hUV[face] = new TH2D(
            Form("hUV_%d",face),
            Form("%s first-exit hitmap;Local U (mm);Local V (mm)",FaceName(face)),
            100,-50,50,100,-50,50);
    }

    TH1D hCos("hCoreExitCos",
              "Core first-exit angle;cos(#theta_{normal});Crossings",
              100,0,1);
    TH1D hWave("hCoreExitWave",
               "Core first-exit wavelength;Wavelength (nm);Crossings",
               140,200,900);
    TH1D hTime("hCoreExitTime",
               "Core first-exit time;Time (ns);Crossings",
               160,0,20);

    for (Long64_t i=0; i<cb->GetEntries(); ++i) {
        cb->GetEntry(i);
        Key k{cEvent,cPhoton};

        if (!crossingPhotons.insert(k).second)
            ++duplicateCrossings;

        crossingInfo[k] = {
            cFace,cParent,cSource,cRoot,cOrigin,cU,cV,cCos,cWave,cTime
        };

        if (cFace >=0 && cFace < 6) {
            ++faceCounts[cFace];
            hUV[cFace]->Fill(cU,cV);
        } else {
            ++badFace;
        }

        if (!std::isfinite(cCos) || cCos <= 0.0 || cCos > 1.00001)
            ++badCos;

        hCos.Fill(cCos);
        hWave.Fill(cWave);
        hTime.Fill(cTime);

        ++crossingOriginCounts[cOrigin];

        if (cOrigin == coreVolumeID)
            coreOriginCrossingPhotons.insert(k);

        auto pit = photonInfo.find(k);
        if (pit != photonInfo.end()) {
            const auto& p = pit->second;
            if (p.parent != cParent || p.source != cSource || p.originVolume != cOrigin)
                ++ancestryMismatch;
            ++creatorCounts[{p.creatorType,p.creatorSubType}];
        }
    }

    // ---------------- SensorIncident unique photons ----------------
    int sEvent=0, sPhoton=0;
    si->SetBranchAddress("EventID",&sEvent);
    si->SetBranchAddress("PhotonID",&sPhoton);

    std::set<Key> sensorPhotons;
    std::set<Key> sensorFromCoreExit;
    std::set<Key> sensorFromCoreOriginExit;

    for (Long64_t i=0; i<si->GetEntries(); ++i) {
        si->GetEntry(i);
        Key k{sEvent,sPhoton};
        sensorPhotons.insert(k);
        if (crossingPhotons.count(k))
            sensorFromCoreExit.insert(k);
        if (coreOriginCrossingPhotons.count(k))
            sensorFromCoreOriginExit.insert(k);
    }

    // ---------------- DetectedSignal unique photons ----------------
    int dEvent=0, dPhoton=0;
    ds->SetBranchAddress("EventID",&dEvent);
    ds->SetBranchAddress("PhotonID",&dPhoton);

    std::set<Key> detectedPhotons;
    std::set<Key> detectedFromCoreExit;
    std::set<Key> detectedFromCoreOriginExit;

    for (Long64_t i=0; i<ds->GetEntries(); ++i) {
        ds->GetEntry(i);
        Key k{dEvent,dPhoton};
        detectedPhotons.insert(k);
        if (crossingPhotons.count(k))
            detectedFromCoreExit.insert(k);
        if (coreOriginCrossingPhotons.count(k))
            detectedFromCoreOriginExit.insert(k);
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nCB3 CORE-BOUNDARY ANALYSIS SANITY\n";
    std::cout << "Input                         : " << inputName << "\n";
    std::cout << "Core VolumeID                 : " << coreVolumeID << "\n\n";

    std::cout << "=== UNIQUE-PHOTON PIPELINE ===\n";
    std::cout << "All Photon                    : " << allPhotons.size() << "\n";
    std::cout << "Core-origin Photon            : " << coreOriginPhotons.size() << "\n";
    std::cout << "CoreBoundary firstExit        : " << crossingPhotons.size() << "\n";
    std::cout << "Core-origin firstExit         : " << coreOriginCrossingPhotons.size() << "\n";
    std::cout << "Unique SensorIncident photons: " << sensorPhotons.size() << "\n";
    std::cout << "  matched to any firstExit    : " << sensorFromCoreExit.size() << "\n";
    std::cout << "  matched Core-origin exit    : " << sensorFromCoreOriginExit.size() << "\n";
    std::cout << "Unique Detected photons      : " << detectedPhotons.size() << "\n";
    std::cout << "  matched to any firstExit    : " << detectedFromCoreExit.size() << "\n";
    std::cout << "  matched Core-origin exit    : " << detectedFromCoreOriginExit.size() << "\n\n";

    std::cout << "=== CONDITIONAL EFFICIENCIES ===\n";
    std::cout << "Core-origin escape            : "
              << pct(coreOriginCrossingPhotons.size(),coreOriginPhotons.size())
              << " %\n";
    std::cout << "Collection after Core exit    : "
              << pct(sensorFromCoreOriginExit.size(),coreOriginCrossingPhotons.size())
              << " %\n";
    std::cout << "Detection after sensor reach  : "
              << pct(detectedFromCoreOriginExit.size(),sensorFromCoreOriginExit.size())
              << " %\n";
    std::cout << "End-to-end Core-origin detect : "
              << pct(detectedFromCoreOriginExit.size(),coreOriginPhotons.size())
              << " %\n\n";

    std::cout << "=== FACE OCCUPANCY ===\n";
    for (int face=0; face<6; ++face) {
        std::cout << std::setw(2) << FaceName(face)
                  << " : " << std::setw(8) << faceCounts[face]
                  << "  (" << pct(faceCounts[face],crossingPhotons.size()) << " %)\n";
    }

    std::cout << "\n=== DISTRIBUTION SUMMARY ===\n";
    std::cout << "mean cos(theta_normal)        : " << hCos.GetMean() << "\n";
    std::cout << "RMS  cos(theta_normal)        : " << hCos.GetRMS() << "\n";
    std::cout << "mean wavelength               : " << hWave.GetMean() << " nm\n";
    std::cout << "RMS  wavelength               : " << hWave.GetRMS() << " nm\n";
    std::cout << "mean first-exit time          : " << hTime.GetMean() << " ns\n";
    std::cout << "RMS  first-exit time          : " << hTime.GetRMS() << " ns\n";

    std::cout << "\n=== ORIGIN VOLUME BREAKDOWN ===\n";
    for (const auto& kv : crossingOriginCounts) {
        auto nit = volumeNames.find(kv.first);
        const std::string name = (nit == volumeNames.end()) ? "UNKNOWN" : nit->second;
        std::cout << "VolumeID=" << kv.first
                  << " (" << name << ") : "
                  << kv.second
                  << " (" << pct(kv.second,crossingPhotons.size()) << " %)\n";
    }

    std::vector<std::pair<long long,std::pair<int,int>>> creatorSorted;
    for (const auto& kv : creatorCounts)
        creatorSorted.push_back({kv.second,kv.first});
    std::sort(creatorSorted.rbegin(),creatorSorted.rend());

    std::cout << "\n=== TOP CREATOR PROCESS CODES ===\n";
    for (size_t i=0; i<std::min<size_t>(10,creatorSorted.size()); ++i) {
        const auto& x = creatorSorted[i];
        std::cout << "type=" << x.second.first
                  << " subtype=" << x.second.second
                  << " : " << x.first
                  << " (" << pct(x.first,crossingPhotons.size()) << " %)\n";
    }

    const bool structuralPass =
        duplicateCrossings == 0 &&
        badFace == 0 &&
        badCos == 0 &&
        ancestryMismatch == 0 &&
        crossingPhotons.size() == static_cast<size_t>(cb->GetEntries());

    std::cout << "\n=== STRUCTURAL SANITY ===\n";
    std::cout << "duplicate firstExit keys      : " << duplicateCrossings << "\n";
    std::cout << "bad FaceID                    : " << badFace << "\n";
    std::cout << "bad cos(theta)                : " << badCos << "\n";
    std::cout << "Photon/Crossing ancestry diff : " << ancestryMismatch << "\n";
    std::cout << "CB3 STRUCTURAL SANITY         : "
              << (structuralPass ? "PASS" : "FAIL") << "\n";

    // ---------------- Plots ----------------
    gSystem->mkdir("validation_plots", kTRUE);

    TCanvas cFaces("cCoreBoundaryFaces","Core boundary first-exit hitmaps",1500,900);
    cFaces.Divide(3,2);
    for (int face=0; face<6; ++face) {
        cFaces.cd(face+1);
        hUV[face]->Draw("COLZ");
    }
    cFaces.SaveAs("validation_plots/core_boundary_cb3_faces.png");
    cFaces.SaveAs("validation_plots/core_boundary_cb3_faces.pdf");

    TCanvas cDist("cCoreBoundaryDist","Core boundary first-exit distributions",1500,450);
    cDist.Divide(3,1);
    cDist.cd(1); hCos.Draw("HIST");
    cDist.cd(2); hWave.Draw("HIST");
    cDist.cd(3); hTime.Draw("HIST");
    cDist.SaveAs("validation_plots/core_boundary_cb3_distributions.png");
    cDist.SaveAs("validation_plots/core_boundary_cb3_distributions.pdf");
}
