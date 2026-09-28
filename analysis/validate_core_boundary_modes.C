#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <string>

namespace {
bool LoadMetadata(TTree* cfg, std::string& schema, std::string& mode, std::string& detail) {
  if (!cfg || cfg->GetEntries()!=1) return false;
  char s[128]={}, m[128]={}, d[128]={};
  if (!cfg->GetBranch("SchemaVersion") ||
      !cfg->GetBranch("CoreBoundaryMode") ||
      !cfg->GetBranch("CoreBoundaryDetail")) return false;
  cfg->SetBranchAddress("SchemaVersion",s);
  cfg->SetBranchAddress("CoreBoundaryMode",m);
  cfg->SetBranchAddress("CoreBoundaryDetail",d);
  cfg->GetEntry(0);
  schema=s; mode=m; detail=d;
  return true;
}
}

void validate_core_boundary_modes(
    const char* prefix="validation_core_boundary_modes") {
  bool allPass=true;

  for (int run=0; run<4; ++run) {
    const std::string fn=std::string(prefix)+"_run"+std::to_string(run)+".root";
    TFile f(fn.c_str(),"READ");
    auto* cb=dynamic_cast<TTree*>(f.Get("CoreBoundaryCrossing"));
    auto* ph=dynamic_cast<TTree*>(f.Get("Photon"));
    auto* cfg=dynamic_cast<TTree*>(f.Get("RunConfig"));
    if (!cb || !ph || !cfg) {
      std::cout << "run" << run << ": missing required tree : FAIL\n";
      allPass=false;
      continue;
    }

    std::string schema,mode,detail;
    const bool metaReadable=LoadMetadata(cfg,schema,mode,detail);
    bool pass=metaReadable && (schema=="general-v1.9" || schema=="general-v1.10");

    if (run==0) {
      pass &= mode=="off" && detail=="basic" && cb->GetEntries()==0;
      std::cout << "run0 OFF entries=" << cb->GetEntries()
                << " metadata=" << mode << "/" << detail
                << " : " << (pass?"PASS":"FAIL") << "\n";
    } else if (run==1 || run==2) {
      float px=0,py=0,pz=0,cosN=0;
      cb->SetBranchAddress("PolX",&px);
      cb->SetBranchAddress("PolY",&py);
      cb->SetBranchAddress("PolZ",&pz);
      cb->SetBranchAddress("CosThetaNormal",&cosN);
      long long bad=0;
      for (Long64_t i=0;i<cb->GetEntries();++i) {
        cb->GetEntry(i);
        if (cosN < 0.999999) ++bad;
        if (run==1) {
          if (!(std::isnan(px) && std::isnan(py) && std::isnan(pz))) ++bad;
        } else {
          const double p2=double(px)*px+double(py)*py+double(pz)*pz;
          if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
              std::abs(std::sqrt(p2)-1.0)>1e-5) ++bad;
        }
      }
      const std::string expectedDetail=(run==1?"basic":"phaseSpace");
      pass &= mode=="firstExit" && detail==expectedDetail &&
              cb->GetEntries()==100 && bad==0;
      std::cout << "run" << run << " " << expectedDetail
                << " entries=" << cb->GetEntries() << " bad=" << bad
                << " : " << (pass?"PASS":"FAIL") << "\n";
    } else {
      int reached=0, nTIR=0, killed=0;
      ph->SetBranchAddress("ReachedCoreWall",&reached);
      ph->SetBranchAddress("NTIR",&nTIR);
      ph->SetBranchAddress("KilledBySafetyCut",&killed);
      long long reachedCount=0,tirCount=0,killedCount=0;
      for (Long64_t i=0;i<ph->GetEntries();++i) {
        ph->GetEntry(i);
        if (reached) ++reachedCount;
        if (nTIR>0) ++tirCount;
        if (killed) ++killedCount;
      }
      pass &= mode=="firstExit" && detail=="basic" &&
              cb->GetEntries()==0 && ph->GetEntries()==100 &&
              reachedCount==100 && tirCount==100 && killedCount==100;
      std::cout << "run3 TIR-only crossings=" << cb->GetEntries()
                << " photonRows=" << ph->GetEntries()
                << " reachedWall=" << reachedCount
                << " hasTIR=" << tirCount
                << " safetyKilled=" << killedCount
                << " : " << (pass?"PASS":"FAIL") << "\n";
    }

    allPass &= pass;
  }

  std::cout << "CORE-BOUNDARY MODE/DETAIL/TIR VALIDATION: "
            << (allPass?"PASS":"FAIL") << "\n";
}
