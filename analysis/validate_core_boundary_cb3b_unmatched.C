#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {
using Key = std::pair<int,int>; // EventID, PhotonID

struct PInfo {
    int origin = -999;
    int creatorType = -999;
    int creatorSubType = -999;
    int reachedWall = -1;
    int nTIR = -1;
};

double pct(long long a, long long b) {
    return b ? 100.0 * double(a) / double(b) : 0.0;
}

std::string volName(int id, const std::map<int,std::string>& names) {
    auto it = names.find(id);
    return it == names.end() ? "UNKNOWN" : it->second;
}
}

void validate_core_boundary_cb3b_unmatched(
    const char* inputName="validation_core_boundary_cb2_phase_100_run0.root")
{
    TFile f(inputName, "READ");
    if (f.IsZombie()) {
        std::cerr << "Cannot open " << inputName << "\n";
        return;
    }

    auto* ph = dynamic_cast<TTree*>(f.Get("Photon"));
    auto* cb = dynamic_cast<TTree*>(f.Get("CoreBoundaryCrossing"));
    auto* si = dynamic_cast<TTree*>(f.Get("SensorIncident"));
    auto* ds = dynamic_cast<TTree*>(f.Get("DetectedSignal"));
    auto* vm = dynamic_cast<TTree*>(f.Get("VolumeMap"));

    if (!ph || !cb || !si || !ds || !vm) {
        std::cerr << "Missing required tree(s)\n";
        return;
    }

    // ---------- Volume names ----------
    int vID=-1;
    char vName[256] = {};
    vm->SetBranchAddress("VolumeID",&vID);
    vm->SetBranchAddress("VolumeName",vName);

    std::map<int,std::string> volumeNames;
    for (Long64_t i=0;i<vm->GetEntries();++i) {
        vm->GetEntry(i);
        volumeNames[vID] = vName;
    }

    // ---------- Photon lookup ----------
    int pEvent=0, pPhoton=0, pOrigin=0;
    int pCreatorType=0, pCreatorSubType=0;
    int pReachedWall=-1, pNTIR=-1;

    ph->SetBranchAddress("EventID",&pEvent);
    ph->SetBranchAddress("PhotonID",&pPhoton);
    ph->SetBranchAddress("OriginVolumeID",&pOrigin);
    ph->SetBranchAddress("CreatorProcessType",&pCreatorType);
    ph->SetBranchAddress("CreatorProcessSubType",&pCreatorSubType);

    const bool haveReachedWall = ph->GetBranch("ReachedCoreWall");
    const bool haveNTIR = ph->GetBranch("NTIR");

    if (haveReachedWall)
        ph->SetBranchAddress("ReachedCoreWall",&pReachedWall);
    if (haveNTIR)
        ph->SetBranchAddress("NTIR",&pNTIR);

    std::map<Key,PInfo> pinfo;

    for (Long64_t i=0;i<ph->GetEntries();++i) {
        ph->GetEntry(i);
        PInfo x;
        x.origin = pOrigin;
        x.creatorType = pCreatorType;
        x.creatorSubType = pCreatorSubType;
        x.reachedWall = haveReachedWall ? pReachedWall : -1;
        x.nTIR = haveNTIR ? pNTIR : -1;
        pinfo[{pEvent,pPhoton}] = x;
    }

    // ---------- CoreBoundary unique keys ----------
    int cEvent=0, cPhoton=0;
    cb->SetBranchAddress("EventID",&cEvent);
    cb->SetBranchAddress("PhotonID",&cPhoton);

    std::set<Key> exited;
    long long duplicateExit = 0;
    for (Long64_t i=0;i<cb->GetEntries();++i) {
        cb->GetEntry(i);
        if (!exited.insert({cEvent,cPhoton}).second)
            ++duplicateExit;
    }

    // ---------- Sensor unique keys ----------
    int sEvent=0, sPhoton=0;
    si->SetBranchAddress("EventID",&sEvent);
    si->SetBranchAddress("PhotonID",&sPhoton);

    std::set<Key> sensor;
    for (Long64_t i=0;i<si->GetEntries();++i) {
        si->GetEntry(i);
        sensor.insert({sEvent,sPhoton});
    }

    // ---------- Detected unique keys ----------
    int dEvent=0, dPhoton=0;
    ds->SetBranchAddress("EventID",&dEvent);
    ds->SetBranchAddress("PhotonID",&dPhoton);

    std::set<Key> detected;
    for (Long64_t i=0;i<ds->GetEntries();++i) {
        ds->GetEntry(i);
        detected.insert({dEvent,dPhoton});
    }

    std::map<int,long long> sensorByOrigin, sensorMatchedByOrigin, sensorUnmatchedByOrigin;
    std::map<int,long long> detectedByOrigin, detectedMatchedByOrigin, detectedUnmatchedByOrigin;

    std::map<std::pair<int,int>,long long> unmatchedSensorCreator;
    long long unmatchedSensorReachedWall = 0;
    long long unmatchedSensorHasTIR = 0;
    long long missingPhotonLookupSensor = 0;

    std::vector<Key> unmatchedSensorExamples;

    long long sensorMatched = 0;
    for (const auto& k : sensor) {
        auto pit = pinfo.find(k);
        if (pit == pinfo.end()) {
            ++missingPhotonLookupSensor;
            continue;
        }

        const auto& p = pit->second;
        ++sensorByOrigin[p.origin];

        if (exited.count(k)) {
            ++sensorMatched;
            ++sensorMatchedByOrigin[p.origin];
        } else {
            ++sensorUnmatchedByOrigin[p.origin];
            ++unmatchedSensorCreator[{p.creatorType,p.creatorSubType}];
            if (p.reachedWall > 0) ++unmatchedSensorReachedWall;
            if (p.nTIR > 0) ++unmatchedSensorHasTIR;
            if (unmatchedSensorExamples.size() < 20)
                unmatchedSensorExamples.push_back(k);
        }
    }

    long long detectedMatched = 0;
    long long missingPhotonLookupDetected = 0;

    for (const auto& k : detected) {
        auto pit = pinfo.find(k);
        if (pit == pinfo.end()) {
            ++missingPhotonLookupDetected;
            continue;
        }

        const auto& p = pit->second;
        ++detectedByOrigin[p.origin];

        if (exited.count(k)) {
            ++detectedMatched;
            ++detectedMatchedByOrigin[p.origin];
        } else {
            ++detectedUnmatchedByOrigin[p.origin];
        }
    }

    const long long sensorUnmatched = sensor.size() - sensorMatched;
    const long long detectedUnmatched = detected.size() - detectedMatched;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nCB3b SENSOR/DETECTED FIRST-EXIT DECOMPOSITION\n";
    std::cout << "Input                         : " << inputName << "\n";
    std::cout << "CoreBoundary unique keys      : " << exited.size() << "\n";
    std::cout << "CoreBoundary duplicate keys   : " << duplicateExit << "\n\n";

    std::cout << "=== SENSOR INCIDENT ===\n";
    std::cout << "Unique sensor photons         : " << sensor.size() << "\n";
    std::cout << "Matched firstExit             : " << sensorMatched
              << " (" << pct(sensorMatched,sensor.size()) << " %)\n";
    std::cout << "Unmatched firstExit           : " << sensorUnmatched
              << " (" << pct(sensorUnmatched,sensor.size()) << " %)\n";
    std::cout << "Missing Photon lookup         : " << missingPhotonLookupSensor << "\n\n";

    std::cout << "Sensor by OriginVolumeID:\n";
    for (const auto& kv : sensorByOrigin) {
        const auto id = kv.first;
        std::cout << "  " << id << " (" << volName(id,volumeNames) << ")"
                  << " total=" << kv.second
                  << " matched=" << sensorMatchedByOrigin[id]
                  << " unmatched=" << sensorUnmatchedByOrigin[id]
                  << "\n";
    }

    std::cout << "\n=== UNMATCHED SENSOR DIAGNOSTICS ===\n";
    std::cout << "ReachedCoreWall > 0           : " << unmatchedSensorReachedWall
              << " / " << sensorUnmatched << "\n";
    std::cout << "NTIR > 0                      : " << unmatchedSensorHasTIR
              << " / " << sensorUnmatched << "\n";

    std::vector<std::pair<long long,std::pair<int,int>>> creatorSorted;
    for (const auto& kv : unmatchedSensorCreator)
        creatorSorted.push_back({kv.second,kv.first});
    std::sort(creatorSorted.rbegin(),creatorSorted.rend());

    std::cout << "Top unmatched creator codes:\n";
    for (size_t i=0;i<std::min<size_t>(10,creatorSorted.size());++i) {
        const auto& x = creatorSorted[i];
        std::cout << "  type=" << x.second.first
                  << " subtype=" << x.second.second
                  << " count=" << x.first << "\n";
    }

    if (!unmatchedSensorExamples.empty()) {
        std::cout << "Example unmatched keys (EventID,PhotonID):\n";
        for (const auto& k : unmatchedSensorExamples)
            std::cout << "  (" << k.first << "," << k.second << ")\n";
    }

    std::cout << "\n=== DETECTED ===\n";
    std::cout << "Unique detected photons       : " << detected.size() << "\n";
    std::cout << "Matched firstExit             : " << detectedMatched
              << " (" << pct(detectedMatched,detected.size()) << " %)\n";
    std::cout << "Unmatched firstExit           : " << detectedUnmatched
              << " (" << pct(detectedUnmatched,detected.size()) << " %)\n";
    std::cout << "Missing Photon lookup         : " << missingPhotonLookupDetected << "\n\n";

    std::cout << "Detected by OriginVolumeID:\n";
    for (const auto& kv : detectedByOrigin) {
        const auto id = kv.first;
        std::cout << "  " << id << " (" << volName(id,volumeNames) << ")"
                  << " total=" << kv.second
                  << " matched=" << detectedMatchedByOrigin[id]
                  << " unmatched=" << detectedUnmatchedByOrigin[id]
                  << "\n";
    }

    const bool pass =
        duplicateExit == 0 &&
        missingPhotonLookupSensor == 0 &&
        missingPhotonLookupDetected == 0 &&
        sensorMatched + sensorUnmatched == static_cast<long long>(sensor.size()) &&
        detectedMatched + detectedUnmatched == static_cast<long long>(detected.size());

    std::cout << "\nCB3b DECOMPOSITION SANITY      : "
              << (pass ? "PASS" : "FAIL") << "\n";
}
