#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <map>
#include <string>

void validate_electron_brems(
    const char* inputName="validation_electron_brems_run0.root")
{
    TFile f(inputName, "READ");

    auto* tr = dynamic_cast<TTree*>(f.Get("Track"));
    auto* tv = dynamic_cast<TTree*>(f.Get("TrackVolume"));

    if (!tr || !tv) {
        std::cerr << "Missing Track or TrackVolume tree\n";
        return;
    }

    int trackID=0, parentID=0, pdg=0;
    std::string* creator=nullptr;

    tr->SetBranchAddress("TrackID",&trackID);
    tr->SetBranchAddress("ParentID",&parentID);
    tr->SetBranchAddress("PDG",&pdg);
    tr->SetBranchAddress("CreatorProcessName",&creator);

    long long nPrimaryE = 0;
    long long nBremsGamma = 0;
    long long nIonElectrons = 0;
    std::map<std::string,long long> creators;

    for (Long64_t i=0;i<tr->GetEntries();++i) {
        tr->GetEntry(i);

        if (parentID==0 && pdg==11)
            ++nPrimaryE;

        if (parentID>0 && creator) {
            creators[*creator]++;

            if (pdg==22 &&
                (creator->find("eBrem") != std::string::npos ||
                 creator->find("Brems") != std::string::npos))
                ++nBremsGamma;

            if (pdg==11 &&
                creator->find("eIoni") != std::string::npos)
                ++nIonElectrons;
        }
    }

    int tvTrack=0, vol=0;
    float edep=0, dke=0, len=0;

    tv->SetBranchAddress("TrackID",&tvTrack);
    tv->SetBranchAddress("VolumeID",&vol);
    tv->SetBranchAddress("Edep_MeV",&edep);
    tv->SetBranchAddress("DeltaKE_MeV",&dke);
    tv->SetBranchAddress("Length_mm",&len);

    double sumLen=0;
    double sumDKE=0;
    long long nCorePrimaryRows=0;

    for (Long64_t i=0;i<tv->GetEntries();++i) {
        tv->GetEntry(i);

        if (tvTrack != 1) continue;
        if (vol != 1) continue;

        sumLen += len;
        sumDKE += dke;
        ++nCorePrimaryRows;
    }

    std::cout << "Primary electrons=" << nPrimaryE << "\n";
    std::cout << "Brems gamma secondaries=" << nBremsGamma << "\n";
    std::cout << "Ionization electron secondaries=" << nIonElectrons << "\n";
    std::cout << "Primary-core TrackVolume rows=" << nCorePrimaryRows << "\n";

    if (nCorePrimaryRows>0) {
        std::cout << "Mean primary path="
                  << (sumLen/nCorePrimaryRows)/10.0
                  << " cm\n";
        std::cout << "Mean primary DeltaKE="
                  << sumDKE/nCorePrimaryRows
                  << " MeV\n";
    }

    std::cout << "\nCreator processes:\n";
    for (const auto& kv : creators)
        std::cout << "  "
                  << kv.first << " : "
                  << kv.second << "\n";

    const bool havePrimaries = nPrimaryE > 0;
    const bool haveBrems     = nBremsGamma > 0;
    const bool haveIoni      = nIonElectrons > 0;

    std::cout
        << "Electron ionization active: "
        << (haveIoni ? "PASS" : "FAIL") << "\n";

    std::cout
        << "Electron bremsstrahlung active: "
        << (haveBrems ? "PASS" : "FAIL") << "\n";

    std::cout
        << "ELECTRON TRANSPORT/BREMS VALIDATION: "
        << ((havePrimaries && haveBrems && haveIoni)
                ? "PASS"
                : "FAIL")
        << "\n";
}
