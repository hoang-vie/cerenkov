#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <cmath>
#include <iostream>

void validate_mcs_muon(
    const char* inputName="validation_dedx_muon_run0.root")
{
    TFile f(inputName, "READ");
    auto* t = dynamic_cast<TTree*>(f.Get("Track"));

    if (!t) {
        std::cerr << "Missing Track tree\n";
        return;
    }

    int trackID=0, parentID=0, pdg=0;
    int entrySet=0, exitSet=0;

    float vertexKE=0;
    float dirX=0, dirY=0, dirZ=0;
    float xIn=0, zIn=0;
    float xOut=0, zOut=0;

    t->SetBranchAddress("TrackID", &trackID);
    t->SetBranchAddress("ParentID", &parentID);
    t->SetBranchAddress("PDG", &pdg);
    t->SetBranchAddress("VertexKE_MeV", &vertexKE);

    t->SetBranchAddress("DirX", &dirX);
    t->SetBranchAddress("DirY", &dirY);
    t->SetBranchAddress("DirZ", &dirZ);

    t->SetBranchAddress("CoreEntrySet", &entrySet);
    t->SetBranchAddress("CoreEntryX_mm", &xIn);
    t->SetBranchAddress("CoreEntryZ_mm", &zIn);

    t->SetBranchAddress("CoreExitSet", &exitSet);
    t->SetBranchAddress("CoreExitX_mm", &xOut);
    t->SetBranchAddress("CoreExitZ_mm", &zOut);

    TH1D hx("hx",
           "Muon lateral displacement;#Delta x (mm);Tracks",
           200, -1.0, 1.0);

    TH1D hz("hz",
           "Muon lateral displacement;#Delta z (mm);Tracks",
           200, -1.0, 1.0);

    long long n = 0;
    double sumKE = 0.0;

    for (Long64_t i=0; i<t->GetEntries(); ++i) {
        t->GetEntry(i);

        if (parentID != 0) continue;
        if (std::abs(pdg) != 13) continue;
        if (!entrySet || !exitSet) continue;

        // Validation beam is nominally along -Y.
        if (dirY > -0.999) continue;

        hx.Fill(xOut - xIn);
        hz.Fill(zOut - zIn);

        sumKE += vertexKE;
        ++n;
    }

    if (n == 0) {
        std::cerr << "No usable primary through-going muons\n";
        return;
    }

    const double meanKE = sumKE / n;

    // muon rest mass, MeV
    const double m = 105.6583755;

    const double E = meanKE + m;
    const double p = std::sqrt(E*E - m*m);
    const double beta = p / E;

    // 10 cm liquid water.
    // X0 ~= 36.08 g/cm^2 and rho ~= 1 g/cm^3,
    // therefore x/X0 ~= 10/36.08.
    const double L_mm = 100.0;
    const double xOverX0 = 10.0 / 36.08;

    // Highland projected RMS scattering angle.
    const double theta0 =
        (13.6 / (beta*p))
        * std::sqrt(xOverX0)
        * (1.0 + 0.038*std::log(xOverX0));

    // Gaussian distributed-scattering approximation:
    // projected lateral RMS ~= L * theta0 / sqrt(3).
    const double sigmaY_mm =
        L_mm * theta0 / std::sqrt(3.0);

    const double rmsX = hx.GetStdDev();
    const double rmsZ = hz.GetStdDev();

    const double ratioX = rmsX / sigmaY_mm;
    const double ratioZ = rmsZ / sigmaY_mm;

    std::cout << "Selected through-going primary muons=" << n << "\n";
    std::cout << "Mean vertex KE=" << meanKE << " MeV\n";
    std::cout << "p=" << p << " MeV/c beta=" << beta << "\n";
    std::cout << "x/X0=" << xOverX0 << "\n";
    std::cout << "Highland theta0=" << theta0*1e3 << " mrad\n";
    std::cout << "Expected projected lateral RMS~"
              << sigmaY_mm << " mm\n";
    std::cout << "MC RMS dx=" << rmsX
              << " mm  ratio=" << ratioX << "\n";
    std::cout << "MC RMS dz=" << rmsZ
              << " mm  ratio=" << ratioZ << "\n";
    std::cout << "MC mean dx=" << hx.GetMean()
              << " mm\n";
    std::cout << "MC mean dz=" << hz.GetMean()
              << " mm\n";

    // Deliberately broad detector-level validation band.
    // Highland + L*theta/sqrt(3) is an approximation, not an exact
    // Geant4 reference prediction.
    const bool symmetry =
        std::abs(rmsX-rmsZ) / (0.5*(rmsX+rmsZ)) < 0.10;

    const bool scale =
        ratioX > 0.70 && ratioX < 1.30 &&
        ratioZ > 0.70 && ratioZ < 1.30;

    const bool centered =
        std::abs(hx.GetMean()) < 0.02 &&
        std::abs(hz.GetMean()) < 0.02;

    const bool pass = symmetry && scale && centered;

    std::cout << "XY symmetry: "
              << (symmetry ? "PASS" : "FAIL") << "\n";
    std::cout << "Highland-scale agreement: "
              << (scale ? "PASS" : "FAIL") << "\n";
    std::cout << "Zero-bias check: "
              << (centered ? "PASS" : "FAIL") << "\n";

    std::cout << "MUON MULTIPLE-SCATTERING VALIDATION: "
              << (pass ? "PASS" : "FAIL") << "\n";

    TCanvas c("cMCS","Muon multiple scattering",1200,500);
    c.Divide(2,1);
    c.cd(1); hx.Draw();
    c.cd(2); hz.Draw();

    c.Print("validation_plots/muon_mcs_validation.png");
    c.Print("validation_plots/muon_mcs_validation.pdf");
}
