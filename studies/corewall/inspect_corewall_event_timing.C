#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

namespace {
double Quantile(std::vector<double> v, double q) {
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(v.begin(),v.end());
  const double x=q*(v.size()-1);
  const size_t lo=size_t(std::floor(x)), hi=size_t(std::ceil(x));
  if (lo==hi) return v[lo];
  const double f=x-lo;
  return v[lo]*(1.0-f)+v[hi]*f;
}
}

void inspect_corewall_event_timing(
    const char* inputName="validation_corewall_muon_topology_run0.root",
    int eventID=71) {

  TFile* f=TFile::Open(inputName,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"Cannot open "<<inputName<<"\n"; return; }

  auto* tPho=dynamic_cast<TTree*>(f->Get("Photon"));
  auto* tTrack=dynamic_cast<TTree*>(f->Get("Track"));
  if(!tPho || !tTrack){ std::cerr<<"Missing Photon/Track\n"; f->Close(); return; }

  struct TI { int parent=0,pdg=0; double birth=0; };
  auto key=[](int e,int t)->unsigned long long {
    return (unsigned long long(unsigned(e))<<32) | unsigned(t);
  };
  std::unordered_map<unsigned long long,TI> tracks;
  int te=0,tid=0,parent=0,pdg=0; float vt=0;
  tTrack->SetBranchAddress("EventID",&te);
  tTrack->SetBranchAddress("TrackID",&tid);
  tTrack->SetBranchAddress("ParentID",&parent);
  tTrack->SetBranchAddress("PDG",&pdg);
  tTrack->SetBranchAddress("VertexTime_ns",&vt);
  for(Long64_t i=0;i<tTrack->GetEntries();++i){
    tTrack->GetEntry(i);
    if(te!=eventID) continue;
    tracks[key(te,tid)]={parent,pdg,vt};
  }

  int pe=0,parentPho=0,reached=0; float wt=0,bt=0;
  tPho->SetBranchAddress("EventID",&pe);
  tPho->SetBranchAddress("ParentTrackID",&parentPho);
  tPho->SetBranchAddress("ReachedCoreWall",&reached);
  tPho->SetBranchAddress("FirstWallTime_ns",&wt);
  tPho->SetBranchAddress("BirthTime_ns",&bt);

  std::vector<double> allT,muT,eT,flight,emitterDelay;
  for(Long64_t i=0;i<tPho->GetEntries();++i){
    tPho->GetEntry(i);
    if(pe!=eventID || !reached) continue;
    allT.push_back(wt);
    flight.push_back(double(wt)-double(bt));
    auto it=tracks.find(key(pe,parentPho));
    if(it!=tracks.end()){
      emitterDelay.push_back(double(wt)-it->second.birth);
      if(it->second.parent==0 && std::abs(it->second.pdg)==13) muT.push_back(wt);
      if(std::abs(it->second.pdg)==11) eT.push_back(wt);
    }
  }

  auto print=[&](const char* label,const std::vector<double>& v){
    if(v.empty()){ std::cout<<label<<": EMPTY\n"; return; }
    auto mm=std::minmax_element(v.begin(),v.end());
    std::cout<<std::left<<std::setw(17)<<label
             <<" N="<<std::setw(6)<<v.size()
             <<" min="<<std::fixed<<std::setprecision(5)<<*mm.first
             <<" q10="<<Quantile(v,.10)
             <<" q25="<<Quantile(v,.25)
             <<" q50="<<Quantile(v,.50)
             <<" q75="<<Quantile(v,.75)
             <<" q90="<<Quantile(v,.90)
             <<" max="<<*mm.second<<"\n";
  };

  std::cout<<"\n=== CORE-WALL STEP-3B TIMING INSPECTION ===\n";
  std::cout<<"EventID="<<eventID<<"\n";
  print("wall time all",allT);
  print("wall time mu",muT);
  print("wall time e",eT);
  print("photon flight",flight);
  print("emitter->wall",emitterDelay);

  if(!allT.empty()){
    const double q10=Quantile(allT,.10),q33=Quantile(allT,1.0/3.0),
                 q67=Quantile(allT,2.0/3.0),q90=Quantile(allT,.90);
    std::cout<<"\nSuggested global-time slices for morphology (data-derived):\n"
             <<"  early : ["<<q10<<", "<<q33<<") ns\n"
             <<"  middle: ["<<q33<<", "<<q67<<") ns\n"
             <<"  late  : ["<<q67<<", "<<q90<<"] ns\n"
             <<"These exclude the outer 10% tails so rare late transport does not dominate the display.\n";
  }
  std::cout<<"CORE-WALL STEP-3B TIMING INSPECTION: READY_FOR_MAPS\n";
  f->Close();
}
