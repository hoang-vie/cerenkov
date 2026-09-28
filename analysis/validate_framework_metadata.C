#include <TFile.h>
#include <TTree.h>

#include <cstring>
#include <iostream>
#include <string>

void validate_framework_metadata(
    const char* fileName="validation_output_mt_run0.root") {
  TFile f(fileName,"READ");
  auto* cfg=dynamic_cast<TTree*>(f.Get("RunConfig"));
  if (!cfg || cfg->GetEntries()!=1) {
    std::cout << "FRAMEWORK METADATA VALIDATION: FAIL (RunConfig missing/not unique)\n";
    return;
  }

  const char* required[] = {
      "SchemaVersion", "ProjectRevision", "GitCommit", "GitBranch",
      "GitDirty", "CoreProfile", "MaterialProfileRegistryVersion"};
  for (const auto* name : required) {
    if (!cfg->GetBranch(name)) {
      std::cout << "FRAMEWORK METADATA VALIDATION: FAIL (missing branch "
                << name << ")\n";
      return;
    }
  }

  char schema[128]={};
  char revision[256]={};
  char commit[256]={};
  char branch[256]={};
  char coreProfile[128]={};
  char registry[256]={};
  int dirty=-1;

  cfg->SetBranchAddress("SchemaVersion",schema);
  cfg->SetBranchAddress("ProjectRevision",revision);
  cfg->SetBranchAddress("GitCommit",commit);
  cfg->SetBranchAddress("GitBranch",branch);
  cfg->SetBranchAddress("GitDirty",&dirty);
  cfg->SetBranchAddress("CoreProfile",coreProfile);
  cfg->SetBranchAddress("MaterialProfileRegistryVersion",registry);
  cfg->GetEntry(0);

  bool pass=true;
  pass &= std::string(schema)=="general-v1.10";
  pass &= dirty==0 || dirty==1;
  pass &= std::strlen(revision)>0;
  pass &= std::strlen(commit)>0;
  pass &= std::strlen(branch)>0;
  pass &= std::string(coreProfile)=="Water_Default";
  pass &= std::strlen(registry)>0;

  const std::string c(commit), r(revision);
  if (c=="unknown") {
    pass &= r=="unknown" || r=="unknown-dirty";
  } else {
    const std::string shortCommit=c.substr(0, c.size()<12 ? c.size() : 12);
    const std::string expected=shortCommit + (dirty ? "-dirty" : "");
    pass &= r==expected;
  }

  std::cout << "SchemaVersion : " << schema << "\n";
  std::cout << "ProjectRevision: " << revision << "\n";
  std::cout << "GitCommit     : " << commit << "\n";
  std::cout << "GitBranch     : " << branch << "\n";
  std::cout << "GitDirty      : " << dirty << "\n";
  std::cout << "CoreProfile   : " << coreProfile << "\n";
  std::cout << "ProfileRegistry: " << registry << "\n";
  std::cout << "FRAMEWORK METADATA VALIDATION: " << (pass?"PASS":"FAIL") << "\n";
}
