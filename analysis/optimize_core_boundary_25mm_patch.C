#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {
using Key = std::pair<int,int>;

struct Hit {
    Key key;
    float u;
    float v;
};

struct Opt {
    double u0=0, v0=0;
    long long photons=0;
    long long ev1=0, ev10=0, ev50=0, ev100=0;
};

const char* FN(int f){
    static const char* n[]={"+X","-X","+Z","-Z","+Y","-Y"};
    return (f>=0 && f<6) ? n[f] : "?";
}

double pct(long long a,long long b){
    return b ? 100.0*double(a)/double(b) : 0.0;
}
}

void optimize_core_boundary_25mm_patch(
    const char* inputName="study_cry_core_boundary_fullall_1000_run0.root",
    double patchSize_mm=25.0,
    double step_mm=1.0)
{
    TFile f(inputName,"READ");
    auto* t=(TTree*)f.Get("CoreBoundaryCrossing");
    if(!t){ std::cerr<<"Missing CoreBoundaryCrossing\n"; return; }

    int event=0,root=0,face=0;
    float u=0,v=0;
    t->SetBranchAddress("EventID",&event);
    t->SetBranchAddress("RootPrimaryTrackID",&root);
    t->SetBranchAddress("FaceID",&face);
    t->SetBranchAddress("LocalU_mm",&u);
    t->SetBranchAddress("LocalV_mm",&v);

    std::array<std::vector<Hit>,6> hits;
    std::array<double,6> minU,minV,maxU,maxV;
    for(int fID=0;fID<6;++fID){
        minU[fID]=minV[fID]= std::numeric_limits<double>::infinity();
        maxU[fID]=maxV[fID]=-std::numeric_limits<double>::infinity();
    }

    std::map<Key,int> useful;
    long long bad=0,total=0;

    for(Long64_t i=0;i<t->GetEntries();++i){
        t->GetEntry(i);
        if(face<0 || face>=6){ ++bad; continue; }

        hits[face].push_back({{event,root},u,v});
        ++useful[{event,root}];
        ++total;

        minU[face]=std::min(minU[face],double(u));
        maxU[face]=std::max(maxU[face],double(u));
        minV[face]=std::min(minV[face],double(v));
        maxV[face]=std::max(maxV[face],double(v));
    }

    const double h=patchSize_mm/2.0;

    auto eval=[&](int fID,double u0,double v0)->Opt{
        std::map<Key,long long> per;
        long long n=0;

        for(const auto& x:hits[fID]){
            if(std::abs(double(x.u)-u0)<=h &&
               std::abs(double(x.v)-v0)<=h){
                ++n;
                ++per[x.key];
            }
        }

        Opt o; o.u0=u0; o.v0=v0; o.photons=n;
        for(const auto& kv:per){
            if(kv.second>=1)   ++o.ev1;
            if(kv.second>=10)  ++o.ev10;
            if(kv.second>=50)  ++o.ev50;
            if(kv.second>=100) ++o.ev100;
        }
        return o;
    };

    std::cout<<std::fixed<<std::setprecision(2);
    std::cout<<"\nOPTIMIZE CENTRAL PATCH POSITION ON CORE BOUNDARY\n";
    std::cout<<"Input: "<<inputName<<"\n";
    std::cout<<"Patch: "<<patchSize_mm<<" x "<<patchSize_mm<<" mm\n";
    std::cout<<"Scan step: "<<step_mm<<" mm\n";
    std::cout<<"Useful root primaries: "<<useful.size()<<"\n";
    std::cout<<"Total firstExit photons: "<<total<<"\n";
    std::cout<<"NOTE: scan bounds are conservatively inferred from observed LocalU/V extrema.\n\n";

    std::array<Opt,6> best;
    std::array<Opt,6> center;

    for(int fID=0;fID<6;++fID){
        center[fID]=eval(fID,0.0,0.0);

        double uLo=minU[fID]+h, uHi=maxU[fID]-h;
        double vLo=minV[fID]+h, vHi=maxV[fID]-h;

        if(uLo>uHi || vLo>vHi){
            std::cout<<FN(fID)<<": face range too small for patch\n";
            continue;
        }

        Opt b;
        b.photons=-1;

        for(double u0=uLo; u0<=uHi+1e-9; u0+=step_mm){
            for(double v0=vLo; v0<=vHi+1e-9; v0+=step_mm){
                Opt o=eval(fID,u0,v0);
                if(o.photons>b.photons)
                    b=o;
            }
        }
        best[fID]=b;

        std::cout<<"=== "<<FN(fID)<<" ===\n";
        std::cout<<"Observed U range: ["<<minU[fID]<<", "<<maxU[fID]<<"] mm\n";
        std::cout<<"Observed V range: ["<<minV[fID]<<", "<<maxV[fID]<<"] mm\n";
        std::cout<<"Center patch : photons="<<center[fID].photons
                 <<"  all="<<pct(center[fID].photons,total)<<"%"
                 <<"  evt>=50="<<pct(center[fID].ev50,useful.size())<<"%"
                 <<"  evt>=100="<<pct(center[fID].ev100,useful.size())<<"%\n";
        std::cout<<"Best patch   : U0="<<b.u0<<" mm V0="<<b.v0<<" mm"
                 <<"  photons="<<b.photons
                 <<"  all="<<pct(b.photons,total)<<"%"
                 <<"  evt>=50="<<pct(b.ev50,useful.size())<<"%"
                 <<"  evt>=100="<<pct(b.ev100,useful.size())<<"%\n";
        std::cout<<"Gain vs center: "
                 <<pct(b.photons-center[fID].photons,center[fID].photons)
                 <<"%\n\n";
    }

    int bestSide=0;
    for(int fID=1;fID<=3;++fID)
        if(best[fID].photons>best[bestSide].photons) bestSide=fID;

    std::cout<<"=== SUMMARY ===\n";
    std::cout<<"Best optimized side: "<<FN(bestSide)<<"\n";
    std::cout<<"Bottom centered photons  : "<<center[5].photons<<"\n";
    std::cout<<"Bottom optimized photons : "<<best[5].photons
             <<" at ("<<best[5].u0<<","<<best[5].v0<<") mm\n";
    std::cout<<"Best side optimized      : "<<FN(bestSide)
             <<" photons="<<best[bestSide].photons
             <<" at ("<<best[bestSide].u0<<","<<best[bestSide].v0<<") mm\n";

    std::cout<<"\nPATCH POSITION OPTIMIZATION SANITY: "
             <<((bad==0 && total==t->GetEntries())?"PASS":"FAIL")<<"\n";
}
