#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <utility>

namespace {
struct Vec3 { double x=0, y=0, z=0; };

Vec3 FaceNormal(int face) {
  switch (face) {
    case 0: return { 1, 0, 0};
    case 1: return {-1, 0, 0};
    case 2: return { 0, 0, 1};
    case 3: return { 0, 0,-1};
    case 4: return { 0, 1, 0};
    case 5: return { 0,-1, 0};
    default: return {};
  }
}

double Dot(const Vec3& a, const Vec3& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

double Mag(const Vec3& a) { return std::sqrt(Dot(a,a)); }

double FaceResidualMm(int face, double x, double y, double z) {
  switch (face) {
    case 0: return std::abs(x - 50.0);
    case 1: return std::abs(x + 50.0);
    case 2: return std::abs(z - 50.0);
    case 3: return std::abs(z + 50.0);
    case 4: return std::abs(y - 50.0);
    case 5: return std::abs(y + 50.0);
    default: return 1e9;
  }
}
}

void validate_core_boundary_first_exit_6faces(
    const char* prefix="validation_core_boundary_first_exit_6faces") {
  const int expectedFace[6] = {0,1,2,3,4,5};
  const int oppositeFace[6] = {1,0,3,2,5,4};
  const char* label[6] = {"+X","-X","+Z","-Z","+Y","-Y"};
  bool allPass = true;

  for (int run=0; run<6; ++run) {
    const std::string fn = std::string(prefix) + "_run" + std::to_string(run) + ".root";
    TFile f(fn.c_str(), "READ");
    auto* t = dynamic_cast<TTree*>(f.Get("CoreBoundaryCrossing"));
    auto* cfg = dynamic_cast<TTree*>(f.Get("RunConfig"));
    if (!t || !cfg) {
      std::cerr << "run" << run << ": missing CoreBoundaryCrossing or RunConfig\n";
      allPass = false;
      continue;
    }

    char schema[128] = {}, mode[128] = {}, detail[128] = {};
    cfg->SetBranchAddress("SchemaVersion", schema);
    cfg->SetBranchAddress("CoreBoundaryMode", mode);
    cfg->SetBranchAddress("CoreBoundaryDetail", detail);
    cfg->GetEntry(0);

    int event=0, photon=0, face=-1;
    float x=0,y=0,z=0,u=0,v=0,dx=0,dy=0,dz=0,cosN=0,wl=0,time=0,px=0,py=0,pz=0;
    t->SetBranchAddress("EventID", &event);
    t->SetBranchAddress("PhotonID", &photon);
    t->SetBranchAddress("FaceID", &face);
    t->SetBranchAddress("X_mm", &x);
    t->SetBranchAddress("Y_mm", &y);
    t->SetBranchAddress("Z_mm", &z);
    t->SetBranchAddress("LocalU_mm", &u);
    t->SetBranchAddress("LocalV_mm", &v);
    t->SetBranchAddress("DirX", &dx);
    t->SetBranchAddress("DirY", &dy);
    t->SetBranchAddress("DirZ", &dz);
    t->SetBranchAddress("CosThetaNormal", &cosN);
    t->SetBranchAddress("Wavelength_nm", &wl);
    t->SetBranchAddress("Time_ns", &time);
    t->SetBranchAddress("PolX", &px);
    t->SetBranchAddress("PolY", &py);
    t->SetBranchAddress("PolZ", &pz);

    Long64_t nExpected=0, nOpposite=0, nOther=0, badPrimary=0, duplicate=0;
    std::set<std::pair<int,int>> keys;
    double maxFaceResidual=0, maxUResidual=0, maxVResidual=0;
    double minCos=1.0, maxPolNormResidual=0, maxDirPol=0;

    for (Long64_t i=0; i<t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!keys.insert({event,photon}).second) ++duplicate;

      if (face == expectedFace[run]) {
        ++nExpected;
        const Vec3 dir{dx,dy,dz};
        const Vec3 pol{px,py,pz};
        const Vec3 normal = FaceNormal(face);
        const double faceResidual = FaceResidualMm(face,x,y,z);
        const double uResidual = std::abs(double(u)-7.0);
        const double vResidual = std::abs(double(v)+5.0);
        const double polNormResidual = std::abs(Mag(pol)-1.0);
        const double dirPol = std::abs(Dot(dir,pol));
        maxFaceResidual = std::max(maxFaceResidual, faceResidual);
        maxUResidual = std::max(maxUResidual, uResidual);
        maxVResidual = std::max(maxVResidual, vResidual);
        minCos = std::min(minCos, double(cosN));
        maxPolNormResidual = std::max(maxPolNormResidual, polNormResidual);
        maxDirPol = std::max(maxDirPol, dirPol);

        const bool bad =
            faceResidual > 1e-3 ||
            uResidual > 1e-3 || vResidual > 1e-3 ||
            Dot(dir,normal) < 0.999999 || cosN < 0.999999 ||
            std::abs(double(wl)-450.0) > 0.02 ||
            !std::isfinite(time) || time <= 0.0 ||
            !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
            polNormResidual > 1e-5 || dirPol > 1e-5;
        if (bad) {
          if (badPrimary < 5) {
            std::cout << "  bad primary row " << i
                      << " face=" << face << " xyz=(" << x << "," << y << "," << z << ")"
                      << " uv=(" << u << "," << v << ")"
                      << " dir=(" << dx << "," << dy << "," << dz << ")"
                      << " cosN=" << cosN << " wl=" << wl
                      << " pol=(" << px << "," << py << "," << pz << ")\n";
          }
          ++badPrimary;
        }
      } else if (face == oppositeFace[run]) {
        ++nOpposite;
      } else {
        ++nOther;
      }
    }

    const bool metadataPass =
        (std::string(schema) == "general-v1.9" || std::string(schema) == "general-v1.10") &&
        std::string(mode) == "firstExit" &&
        std::string(detail) == "phaseSpace";
    const bool pass =
        t->GetEntries() == 1000 &&
        nExpected > 900 && nOther == 0 && badPrimary == 0 && duplicate == 0 &&
        metadataPass;
    allPass &= pass;

    std::cout << "run" << run << " " << label[run]
              << " entries=" << t->GetEntries()
              << " primary=" << nExpected
              << " opposite(Fresnel return)=" << nOpposite
              << " other=" << nOther
              << " duplicates=" << duplicate
              << " badPrimary=" << badPrimary
              << " maxFaceResidual=" << maxFaceResidual << " mm"
              << " max|U-7|=" << maxUResidual << " mm"
              << " max|V+5|=" << maxVResidual << " mm"
              << " minCos=" << minCos
              << " max|pol|-1=" << maxPolNormResidual
              << " max|dir.pol|=" << maxDirPol
              << " metadata=" << (metadataPass?"PASS":"FAIL")
              << " : " << (pass?"PASS":"FAIL") << "\n";
  }

  std::cout << "CORE-BOUNDARY FIRST-EXIT 6-FACE VALIDATION: "
            << (allPass?"PASS":"FAIL") << "\n";
}
