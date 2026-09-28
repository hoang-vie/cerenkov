#include <TFile.h>
#include <TTree.h>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <utility>

namespace {
using Key = std::pair<int,int>;
const char* FN(int f){
    static const char* n[]={"+X","-X","+Z","-Z","+Y","-Y"};
    return (f>=0&&f<6)?n[f]:"?";
}
double pct(long long a,long long b){ return b ? 100.0*double(a)/double(b) : 0.0; }
struct ExitInfo{ int face=-1; float u=0,v=0; };
}

void analyze_core_exit_to_sensor_transport(
    const char* inputName="study_cry_core_boundary_fullall_1000_run0.root",
    double patchSize_mm=25.0)
{
    TFile f(inputName,"READ");
    auto* cb=(TTree*)f.Get("CoreBoundaryCrossing");
    auto* si=(TTree*)f.Get("SensorIncident");
    auto* ds=(TTree*)f.Get("DetectedSignal");
    if(!cb||!si||!ds){ std::cerr<<"Missing required tree(s)\n"; return; }

    const double h=patchSize_mm/2.0;

    int e=0,p=0,face=0; float u=0,v=0;
    cb->SetBranchAddress("EventID",&e);
    cb->SetBranchAddress("PhotonID",&p);
    cb->SetBranchAddress("FaceID",&face);
    cb->SetBranchAddress("LocalU_mm",&u);
    cb->SetBranchAddress("LocalV_mm",&v);

    std::map<Key,ExitInfo> exitMap;
    std::array<long long,6> exits{}, patchExits{};
    long long dup=0,badFace=0;

    for(Long64_t i=0;i<cb->GetEntries();++i){
        cb->GetEntry(i);
        if(face<0||face>=6){ ++badFace; continue; }
        Key k{e,p};
        if(exitMap.count(k)) ++dup;
        exitMap[k]={face,u,v};
        ++exits[face];
        if(std::abs(u)<=h && std::abs(v)<=h) ++patchExits[face];
    }

    int se=0,sp=0;
    si->SetBranchAddress("EventID",&se);
    si->SetBranchAddress("PhotonID",&sp);

    std::set<Key> sensorKeys;
    std::array<long long,6> sensorFromFace{}, sensorFromPatch{};
    long long sensorWithoutExit=0;

    for(Long64_t i=0;i<si->GetEntries();++i){
        si->GetEntry(i);
        Key k{se,sp};
        if(!sensorKeys.insert(k).second) continue;
        auto it=exitMap.find(k);
        if(it==exitMap.end()){ ++sensorWithoutExit; continue; }
        const auto& x=it->second;
        ++sensorFromFace[x.face];
        if(std::abs(x.u)<=h && std::abs(x.v)<=h) ++sensorFromPatch[x.face];
    }

    int de=0,dp=0;
    ds->SetBranchAddress("EventID",&de);
    ds->SetBranchAddress("PhotonID",&dp);

    std::set<Key> detKeys;
    std::array<long long,6> detFromFace{}, detFromPatch{};
    long long detWithoutExit=0;

    for(Long64_t i=0;i<ds->GetEntries();++i){
        ds->GetEntry(i);
        Key k{de,dp};
        if(!detKeys.insert(k).second) continue;
        auto it=exitMap.find(k);
        if(it==exitMap.end()){ ++detWithoutExit; continue; }
        const auto& x=it->second;
        ++detFromFace[x.face];
        if(std::abs(x.u)<=h && std::abs(x.v)<=h) ++detFromPatch[x.face];
    }

    std::cout<<std::fixed<<std::setprecision(2);
    std::cout<<"\nCORE EXIT -> SENSOR TRANSPORT BY FACE\n";
    std::cout<<"Input: "<<inputName<<"\n";
    std::cout<<"Virtual centered patch: "<<patchSize_mm<<" x "<<patchSize_mm<<" mm\n\n";

    std::cout<<"=== FULL-FACE FIRST EXIT ===\n";
    std::cout<<std::setw(4)<<"face"<<std::setw(12)<<"exit"<<std::setw(12)<<"sensor"
             <<std::setw(12)<<"det"<<std::setw(14)<<"exit->sens"
             <<std::setw(14)<<"sens->det"<<std::setw(14)<<"exit->det"<<"\n";
    for(int fID=0;fID<6;++fID){
        std::cout<<std::setw(4)<<FN(fID)<<std::setw(12)<<exits[fID]
                 <<std::setw(12)<<sensorFromFace[fID]<<std::setw(12)<<detFromFace[fID]
                 <<std::setw(13)<<pct(sensorFromFace[fID],exits[fID])<<"%"
                 <<std::setw(13)<<pct(detFromFace[fID],sensorFromFace[fID])<<"%"
                 <<std::setw(13)<<pct(detFromFace[fID],exits[fID])<<"%\n";
    }

    std::cout<<"\n=== CENTERED PATCH ONLY ===\n";
    std::cout<<std::setw(4)<<"face"<<std::setw(12)<<"patchExit"<<std::setw(12)<<"sensor"
             <<std::setw(12)<<"det"<<std::setw(14)<<"exit->sens"
             <<std::setw(14)<<"sens->det"<<std::setw(14)<<"exit->det"<<"\n";
    for(int fID=0;fID<6;++fID){
        std::cout<<std::setw(4)<<FN(fID)<<std::setw(12)<<patchExits[fID]
                 <<std::setw(12)<<sensorFromPatch[fID]<<std::setw(12)<<detFromPatch[fID]
                 <<std::setw(13)<<pct(sensorFromPatch[fID],patchExits[fID])<<"%"
                 <<std::setw(13)<<pct(detFromPatch[fID],sensorFromPatch[fID])<<"%"
                 <<std::setw(13)<<pct(detFromPatch[fID],patchExits[fID])<<"%\n";
    }

    std::cout<<"\n=== GLOBAL MATCHING ===\n";
    std::cout<<"Unique firstExit photons  : "<<exitMap.size()<<"\n";
    std::cout<<"Unique sensor photons     : "<<sensorKeys.size()<<"\n";
    std::cout<<"Unique detected photons   : "<<detKeys.size()<<"\n";
    std::cout<<"Sensor without firstExit  : "<<sensorWithoutExit<<"\n";
    std::cout<<"Detected without firstExit: "<<detWithoutExit<<"\n";
    std::cout<<"Duplicate firstExit keys  : "<<dup<<"\n";
    std::cout<<"Bad FaceID                : "<<badFace<<"\n";

    std::cout<<"\nEXIT-TO-SENSOR TRANSPORT SANITY: "
             <<((dup==0 && badFace==0)?"PASS":"FAIL")<<"\n";
}
