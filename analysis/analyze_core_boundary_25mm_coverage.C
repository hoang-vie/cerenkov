#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {
using Key = std::pair<int,int>; // (EventID, RootPrimaryTrackID)

struct Crossing {
    int face;
    float u;
    float v;
};

struct Scenario {
    std::string name;
    std::set<int> faces;
};

const char* FaceName(int f){
    static const char* n[]={"+X","-X","+Z","-Z","+Y","-Y"};
    return (f>=0 && f<6) ? n[f] : "?";
}

double pct(long long a,long long b){
    return b ? 100.0*double(a)/double(b) : 0.0;
}
}

void analyze_core_boundary_25mm_coverage(
    const char* inputName="study_cry_core_boundary_fullall_1000_run0.root",
    double sensorSize_mm=25.0)
{
    TFile f(inputName,"READ");
    auto* cb=(TTree*)f.Get("CoreBoundaryCrossing");
    if(!cb){
        std::cerr<<"Missing CoreBoundaryCrossing\n";
        return;
    }

    int event=0,root=0,face=0;
    float u=0,v=0;
    cb->SetBranchAddress("EventID",&event);
    cb->SetBranchAddress("RootPrimaryTrackID",&root);
    cb->SetBranchAddress("FaceID",&face);
    cb->SetBranchAddress("LocalU_mm",&u);
    cb->SetBranchAddress("LocalV_mm",&v);

    const double half = 0.5*sensorSize_mm;

    std::map<Key,std::vector<Crossing>> perPrimary;
    std::array<long long,6> fullFace{};
    std::array<long long,6> centralPatch{};
    long long total=0,badFace=0;

    for(Long64_t i=0;i<cb->GetEntries();++i){
        cb->GetEntry(i);
        if(face<0 || face>=6){
            ++badFace;
            continue;
        }

        perPrimary[{event,root}].push_back({face,u,v});
        ++fullFace[face];
        ++total;

        if(std::abs(u)<=half && std::abs(v)<=half)
            ++centralPatch[face];
    }

    int bestSide=0;
    for(int fID=1;fID<=3;++fID)
        if(centralPatch[fID]>centralPatch[bestSide]) bestSide=fID;

    std::vector<Scenario> scenarios = {
        {"bottom_only",{5}},
        {"bottom_plus_best_side",{5,bestSide}},
        {"bottom_plus_+X",{5,0}},
        {"bottom_plus_-X",{5,1}},
        {"bottom_plus_+Z",{5,2}},
        {"bottom_plus_-Z",{5,3}},
        {"bottom_plus_X_pair",{5,0,1}},
        {"bottom_plus_Z_pair",{5,2,3}},
        {"all_4_sides",{0,1,2,3}},
        {"bottom_plus_4_sides",{5,0,1,2,3}},
        {"all_6_faces",{0,1,2,3,4,5}}
    };

    const int thresholds[] = {1,10,50,100};

    std::cout<<std::fixed<<std::setprecision(2);
    std::cout<<"\nCENTRAL "<<sensorSize_mm<<" x "<<sensorSize_mm
             <<" mm VIRTUAL APERTURE STUDY\n";
    std::cout<<"Input: "<<inputName<<"\n";
    std::cout<<"Useful root primaries: "<<perPrimary.size()<<"\n";
    std::cout<<"Total firstExit photons: "<<total<<"\n";
    std::cout<<"Aperture cut: |LocalU| <= "<<half
             <<" mm and |LocalV| <= "<<half<<" mm\n";
    std::cout<<"Bad FaceID: "<<badFace<<"\n\n";

    std::cout<<"=== EACH FACE: FULL FACE -> CENTRAL PATCH ===\n";
    for(int fID=0;fID<6;++fID){
        std::cout<<std::setw(3)<<FaceName(fID)
                 <<" full="<<std::setw(8)<<fullFace[fID]
                 <<" patch="<<std::setw(8)<<centralPatch[fID]
                 <<" patch/full="<<std::setw(7)<<pct(centralPatch[fID],fullFace[fID])<<"%"
                 <<" patch/all="<<std::setw(7)<<pct(centralPatch[fID],total)<<"%\n";
    }

    std::cout<<"\nBest side for centered "<<sensorSize_mm<<" mm patch: "
             <<FaceName(bestSide)<<"\n";

    std::cout<<"\n=== CONFIGURATION COVERAGE ===\n";
    std::cout<<std::left
             <<std::setw(28)<<"scenario"
             <<std::setw(8)<<"Nface"
             <<std::setw(14)<<"patchPhotons"
             <<std::setw(14)<<"allPhoton%"
             <<std::setw(11)<<">=1 evt%"
             <<std::setw(11)<<">=10 evt%"
             <<std::setw(11)<<">=50 evt%"
             <<std::setw(11)<<">=100 evt%"
             <<"\n";

    for(const auto& s:scenarios){
        long long photons=0;
        long long evPass[4]={0,0,0,0};

        for(const auto& kv:perPrimary){
            long long n=0;
            for(const auto& c:kv.second){
                if(!s.faces.count(c.face)) continue;
                if(std::abs(c.u)<=half && std::abs(c.v)<=half)
                    ++n;
            }

            photons += n;
            for(int j=0;j<4;++j)
                if(n>=thresholds[j]) ++evPass[j];
        }

        std::cout<<std::left
                 <<std::setw(28)<<s.name
                 <<std::setw(8)<<s.faces.size()
                 <<std::setw(14)<<photons
                 <<std::setw(13)<<pct(photons,total)<<"%"
                 <<std::setw(10)<<pct(evPass[0],perPrimary.size())<<"%"
                 <<std::setw(10)<<pct(evPass[1],perPrimary.size())<<"%"
                 <<std::setw(10)<<pct(evPass[2],perPrimary.size())<<"%"
                 <<std::setw(10)<<pct(evPass[3],perPrimary.size())<<"%"
                 <<"\n";
    }

    std::cout<<"\nInterpretation:\n";
    std::cout<<"- This applies a centered square aperture on the Core boundary only.\n";
    std::cout<<"- It is NOT yet propagation through Acrylic/Coupling/concentrator/PDE.\n";
    std::cout<<"- Event thresholds are boundary photons inside the virtual aperture.\n";
    std::cout<<"- Use these numbers to decide whether full-face conclusions survive realistic area restriction.\n";

    const bool pass=(badFace==0 && total==cb->GetEntries());
    std::cout<<"\n25MM COVERAGE STRUCTURAL SANITY: "
             <<(pass?"PASS":"FAIL")<<"\n";
}
