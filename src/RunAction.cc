#include "RunAction.hh"

#include "BuildInfo.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SimulationMetadata.hh"

#include "G4AccumulableManager.hh"
#include "G4AnalysisManager.hh"
#include "G4Material.hh"
#include "G4OpticalParameters.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4VUserPhysicsList.hh"
#include "G4Version.hh"
#include "Randomize.hh"

#include <string>

namespace {
G4bool IsSingleWriterWorker() {
  return !G4Threading::IsMultithreadedApplication() ||
         (!G4Threading::IsMasterThread() && G4Threading::G4GetThreadId() == 0);
}
}

RunAction::RunAction() {
  auto* am = G4AnalysisManager::Instance();
  am->SetVerboseLevel(0);
  am->SetNtupleMerging(true);
  am->SetDefaultFileType("root");

  auto* accum = G4AccumulableManager::Instance();
  accum->Register(fGeneratedEvents);
  accum->Register(fTriggeredEvents);
  accum->Register(fDetailedEvents);
  accum->Register(fCRYSimulatedTime_s);

  am->CreateNtuple("TriggerEvent", "Thin event and trigger summary");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("TriggerPassed");
  am->CreateNtupleFColumn("EdepTop_MeV");
  am->CreateNtupleFColumn("EdepBottom_MeV");
  am->CreateNtupleFColumn("TimeTop_ns");
  am->CreateNtupleFColumn("TimeBottom_ns");
  am->CreateNtupleIColumn("NPrimaryNonOptical");
  am->CreateNtupleIColumn("NNonOpticalTracks");
  am->CreateNtupleIColumn("NOpticalTracks");
  am->CreateNtupleIColumn("NIncident");
  am->CreateNtupleIColumn("NDetected");
  am->CreateNtupleFColumn("EdepCoreTotal_MeV");
  am->FinishNtuple(kTriggerEvent);

  am->CreateNtuple("ScintDeposit", "Sparse raw scintillator step deposits for offline trigger reconstruction");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("PaddleID");
  am->CreateNtupleIColumn("TrackID");
  am->CreateNtupleFColumn("Time_ns");
  am->CreateNtupleFColumn("Edep_MeV");
  am->FinishNtuple(kScintDeposit);

  am->CreateNtuple("Track", "Generic non-optical track truth");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("TrackID");
  am->CreateNtupleIColumn("ParentID");
  am->CreateNtupleIColumn("RootPrimaryTrackID");
  am->CreateNtupleIColumn("PDG");
  am->CreateNtupleFColumn("Charge_e");
  am->CreateNtupleIColumn("CreatorProcessType");
  am->CreateNtupleIColumn("CreatorProcessSubType");
  am->CreateNtupleSColumn("CreatorProcessName");
  am->CreateNtupleIColumn("BirthVolumeID");
  am->CreateNtupleFColumn("VertexKE_MeV");
  am->CreateNtupleFColumn("VertexTime_ns");
  am->CreateNtupleFColumn("VertexX_mm");
  am->CreateNtupleFColumn("VertexY_mm");
  am->CreateNtupleFColumn("VertexZ_mm");
  am->CreateNtupleFColumn("DirX");
  am->CreateNtupleFColumn("DirY");
  am->CreateNtupleFColumn("DirZ");
  am->CreateNtupleIColumn("CoreEntrySet");
  am->CreateNtupleFColumn("CoreEntryX_mm");
  am->CreateNtupleFColumn("CoreEntryY_mm");
  am->CreateNtupleFColumn("CoreEntryZ_mm");
  am->CreateNtupleIColumn("CoreExitSet");
  am->CreateNtupleFColumn("CoreExitX_mm");
  am->CreateNtupleFColumn("CoreExitY_mm");
  am->CreateNtupleFColumn("CoreExitZ_mm");
  am->CreateNtupleFColumn("EndTime_ns");
  am->CreateNtupleFColumn("EndProperTime_ns");
  am->CreateNtupleFColumn("EndX_mm");
  am->CreateNtupleFColumn("EndY_mm");
  am->CreateNtupleFColumn("EndZ_mm");
  am->CreateNtupleIColumn("EndVolumeID");
  am->CreateNtupleIColumn("EndProcessType");
  am->CreateNtupleIColumn("EndProcessSubType");
  am->CreateNtupleSColumn("EndProcessName");
  am->CreateNtupleIColumn("EndTrackStatus");
  am->FinishNtuple(kTrack);

  am->CreateNtuple("TrackVolume", "Per-track per-volume deposited energy and path length");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("TrackID");
  am->CreateNtupleIColumn("VolumeID");
  am->CreateNtupleFColumn("Edep_MeV");
  am->CreateNtupleFColumn("DeltaKE_MeV");
  am->CreateNtupleFColumn("Length_mm");
  am->CreateNtupleIColumn("NSteps");
  am->FinishNtuple(kTrackVolume);

  am->CreateNtuple("Photon", "One compact summary row per optical photon");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("PhotonID");
  am->CreateNtupleIColumn("ParentTrackID");
  am->CreateNtupleIColumn("SourceTrackID");
  am->CreateNtupleIColumn("CreatorProcessType");
  am->CreateNtupleIColumn("CreatorProcessSubType");
  am->CreateNtupleIColumn("OriginVolumeID");
  am->CreateNtupleFColumn("Wavelength_nm");
  am->CreateNtupleFColumn("BirthTime_ns");
  am->CreateNtupleFColumn("EmissionAngle_deg");
  am->CreateNtupleFColumn("BirthX_mm");
  am->CreateNtupleFColumn("BirthY_mm");
  am->CreateNtupleFColumn("BirthZ_mm");
  am->CreateNtupleFColumn("BirthDirX");
  am->CreateNtupleFColumn("BirthDirY");
  am->CreateNtupleFColumn("BirthDirZ");
  am->CreateNtupleIColumn("ReachedCoreWall");
  am->CreateNtupleFColumn("FirstWallTime_ns");
  am->CreateNtupleFColumn("FirstWallX_mm");
  am->CreateNtupleFColumn("FirstWallY_mm");
  am->CreateNtupleFColumn("FirstWallZ_mm");
  am->CreateNtupleIColumn("ExtractedToReadout");
  am->CreateNtupleFColumn("ExtractionTime_ns");
  am->CreateNtupleIColumn("NBoundaryReflections");
  am->CreateNtupleIColumn("NTIR");
  am->CreateNtupleIColumn("NTyvekReflections");
  am->CreateNtupleIColumn("NRayleigh");
  am->CreateNtupleIColumn("EndReason");
  am->CreateNtupleIColumn("EndVolumeID");
  am->CreateNtupleFColumn("EndTime_ns");
  am->CreateNtupleIColumn("KilledBySafetyCut");
  am->CreateNtupleFColumn("FirstRayleighTrackLength_mm");
  am->CreateNtupleFColumn("EndTrackLength_mm");
  am->CreateNtupleIColumn("AcrylicAbsOutOfRange");
  // Appended for backward-compatible column numbering. For Cerenkov photons
  // this is the exact step-average beta used by G4Cerenkov angle sampling.
  am->CreateNtupleFColumn("EmitterBeta");
  am->FinishNtuple(kPhoton);

  am->CreateNtuple("SensorIncident", "Photon phase space at readout entrance");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("PhotonID");
  am->CreateNtupleIColumn("ChannelID");
  am->CreateNtupleFColumn("Wavelength_nm");
  am->CreateNtupleFColumn("Time_ns");
  am->CreateNtupleFColumn("LocalX_mm");
  am->CreateNtupleFColumn("LocalY_mm");
  am->CreateNtupleFColumn("LocalZ_mm");
  am->CreateNtupleFColumn("LocalDirX");
  am->CreateNtupleFColumn("LocalDirY");
  am->CreateNtupleFColumn("LocalDirZ");
  am->CreateNtupleFColumn("LocalPolX");
  am->CreateNtupleFColumn("LocalPolY");
  am->CreateNtupleFColumn("LocalPolZ");
  am->CreateNtupleFColumn("IncidentAngle_deg");
  am->CreateNtupleIColumn("NBoundaryReflections");
  am->CreateNtupleIColumn("NTIR");
  am->CreateNtupleIColumn("NTyvekReflections");
  am->CreateNtupleIColumn("NRayleigh");
  am->CreateNtupleFColumn("EntranceX_mm");
  am->CreateNtupleFColumn("EntranceY_mm");
  am->CreateNtupleFColumn("EntranceZ_mm");
  am->CreateNtupleIColumn("Forward");
  am->FinishNtuple(kSensorIncident);

  am->CreateNtuple("DetectedSignal", "Detected optical photons in DirectSiPM mode");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("PhotonID");
  am->CreateNtupleIColumn("ChannelID");
  am->CreateNtupleFColumn("Time_ns");
  am->FinishNtuple(kDetectedSignal);

  am->CreateNtuple("ChannelEvent", "Sparse per-event readout-channel aggregate");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("ChannelID");
  am->CreateNtupleIColumn("NIncident");
  am->CreateNtupleIColumn("NDetected");
  am->CreateNtupleFColumn("FirstIncidentTime_ns");
  am->CreateNtupleFColumn("FirstDetectedTime_ns");
  am->FinishNtuple(kChannelEvent);

  am->CreateNtuple("RunConfig", "Run configuration and provenance");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("GeometryVersion");
  am->CreateNtupleFColumn("CoreX_mm");
  am->CreateNtupleFColumn("CoreY_mm");
  am->CreateNtupleFColumn("CoreZ_mm");
  am->CreateNtupleFColumn("ScintX_mm");
  am->CreateNtupleFColumn("ScintY_mm");
  am->CreateNtupleFColumn("ScintZ_mm");
  am->CreateNtupleFColumn("ReadoutX_mm");
  am->CreateNtupleFColumn("ReadoutY_mm");
  am->CreateNtupleFColumn("ReadoutZ_mm");
  am->CreateNtupleFColumn("CouplingX_mm");
  am->CreateNtupleFColumn("CouplingY_mm");
  am->CreateNtupleFColumn("CouplingZ_mm");
  am->CreateNtupleFColumn("HoleX_mm");
  am->CreateNtupleFColumn("HoleY_mm");
  am->CreateNtupleFColumn("AcrylicThickness_mm");
  am->CreateNtupleFColumn("TyvekThickness_mm");
  am->CreateNtupleFColumn("ScintGap_mm");
  am->CreateNtupleFColumn("SourceDistance_mm");
  am->CreateNtupleSColumn("TyvekConfig");
  am->CreateNtupleSColumn("CoreMaterial");
  am->CreateNtupleIColumn("UseAcrylicShell");
  am->CreateNtupleSColumn("SensorMode");
  am->CreateNtupleFColumn("ScintThreshold_MeV");
  am->CreateNtupleFColumn("CoincidenceWindow_ns");
  am->CreateNtupleFColumn("IntegrationWindow_ns");
  am->CreateNtupleSColumn("OutputEventMode");
  am->CreateNtupleSColumn("IncidentDetail");
  am->CreateNtupleSColumn("PhotonBirthDetail");
  am->CreateNtupleIColumn("BoundaryDebug");
  am->CreateNtupleIColumn("MaxOpticalBounces");
  am->CreateNtupleIColumn("ChannelCount");
  am->CreateNtupleIColumn("CompressionLevel");
  am->CreateNtupleSColumn("AcrylicRindexDataset");
  am->CreateNtupleSColumn("AcrylicAbsDataset");
  am->CreateNtupleSColumn("TyvekDataset");
  am->CreateNtupleSColumn("GeneratorMode");
  am->CreateNtupleSColumn("CRYTimeMode");
  am->CreateNtupleSColumn("CRYConfig");
  am->CreateNtupleSColumn("Geant4Version");
  am->CreateNtupleSColumn("PhysicsList");
  am->CreateNtupleDColumn("RNGSeed");
  am->CreateNtupleSColumn("AirRindexDataset");
  am->CreateNtupleSColumn("WaterRindexDataset");
  am->CreateNtupleSColumn("WaterAbsDataset");
  am->CreateNtupleSColumn("WaterRayleighDataset");
  am->CreateNtupleSColumn("CouplingRindexDataset");
  am->CreateNtupleSColumn("SensorPDEDataset");
  am->CreateNtupleSColumn("GunParticle");
  am->CreateNtupleFColumn("GunEnergy_MeV");
  am->CreateNtupleFColumn("GunX_mm");
  am->CreateNtupleFColumn("GunY_mm");
  am->CreateNtupleFColumn("GunZ_mm");
  am->CreateNtupleFColumn("GunDirX");
  am->CreateNtupleFColumn("GunDirY");
  am->CreateNtupleFColumn("GunDirZ");
  am->CreateNtupleFColumn("GunTime_ns");
  am->CreateNtupleIColumn("GunUseAutoY");
  am->CreateNtupleFColumn("GunPolX");
  am->CreateNtupleFColumn("GunPolY");
  am->CreateNtupleFColumn("GunPolZ");
  am->CreateNtupleSColumn("SchemaVersion");
  am->CreateNtupleIColumn("NThreads");
  am->CreateNtupleSColumn("RNGEngine");
  am->CreateNtupleSColumn("CRYVersion");
  am->CreateNtupleFColumn("DefaultCut_mm");
  am->CreateNtupleIColumn("CerenkovEnabled");
  am->CreateNtupleIColumn("ScintillationEnabled");
  am->CreateNtupleIColumn("AbsorptionEnabled");
  am->CreateNtupleIColumn("RayleighEnabled");
  am->CreateNtupleIColumn("MieEnabled");
  am->CreateNtupleIColumn("BoundaryEnabled");
  am->CreateNtupleIColumn("WLSEnabled");
  am->CreateNtupleIColumn("WLS2Enabled");
  am->CreateNtupleFColumn("AcrylicAbsValidMin_nm");
  am->CreateNtupleFColumn("AcrylicAbsValidMax_nm");
  am->CreateNtupleSColumn("AcrylicAbsOutOfRangePolicy");
  am->CreateNtupleSColumn("CouplingAbsDataset");
  am->CreateNtupleSColumn("CouplingMaterialModel");
  am->CreateNtupleSColumn("SensorModel");
  am->CreateNtupleSColumn("ProjectRevision");
  am->CreateNtupleFColumn("AcrylicRindexValidMin_nm");
  am->CreateNtupleFColumn("AcrylicRindexValidMax_nm");
  am->CreateNtupleFColumn("CouplingRindexValidMin_nm");
  am->CreateNtupleFColumn("CouplingRindexValidMax_nm");
  am->CreateNtupleFColumn("CouplingAbsValidMin_nm");
  am->CreateNtupleFColumn("CouplingAbsValidMax_nm");
  am->CreateNtupleFColumn("WaterRindexValidMin_nm");
  am->CreateNtupleFColumn("WaterRindexValidMax_nm");
  am->CreateNtupleFColumn("WaterAbsValidMin_nm");
  am->CreateNtupleFColumn("WaterAbsValidMax_nm");
  am->CreateNtupleFColumn("WaterRayleighValidMin_nm");
  am->CreateNtupleFColumn("WaterRayleighValidMax_nm");
  am->CreateNtupleSColumn("WaterRayleighModelStatus");
  am->CreateNtupleSColumn("OpticalDataOutOfRangePolicy");
  am->CreateNtupleSColumn("CoreEntryExitSemantic");
  am->CreateNtupleSColumn("SensorPDEOperatingPoint");
  am->CreateNtupleSColumn("OpticalDataFingerprint_FNV1a64");

  // general-v1.8: additive profile-aware metadata.  Columns 0..98 above are
  // deliberately preserved byte-for-byte in name/order/type for v1.6 analysis
  // compatibility; generic material-role provenance starts at column 99.
  am->CreateNtupleSColumn("MaterialProfileRegistryVersion");             // 99
  am->CreateNtupleSColumn("LegacyRunConfigCompatibility");               // 100

  am->CreateNtupleSColumn("CoreProfile");                                // 101
  am->CreateNtupleSColumn("CoreProfileG4Material");                       // 102
  am->CreateNtupleSColumn("CoreProfileRindexDataset");                    // 103
  am->CreateNtupleSColumn("CoreProfileAbsDataset");                       // 104
  am->CreateNtupleSColumn("CoreProfileRayleighDataset");                  // 105
  am->CreateNtupleFColumn("CoreProfileRindexValidMin_nm");                // 106
  am->CreateNtupleFColumn("CoreProfileRindexValidMax_nm");                // 107
  am->CreateNtupleFColumn("CoreProfileAbsValidMin_nm");                   // 108
  am->CreateNtupleFColumn("CoreProfileAbsValidMax_nm");                   // 109
  am->CreateNtupleFColumn("CoreProfileRayleighValidMin_nm");              // 110
  am->CreateNtupleFColumn("CoreProfileRayleighValidMax_nm");              // 111
  am->CreateNtupleSColumn("CoreProfileRindexStatus");                     // 112
  am->CreateNtupleSColumn("CoreProfileAbsStatus");                        // 113
  am->CreateNtupleSColumn("CoreProfileRayleighStatus");                   // 114
  am->CreateNtupleSColumn("CoreProfileOutOfRangePolicy");                 // 115

  am->CreateNtupleIColumn("UseOuterShell");                               // 116
  am->CreateNtupleFColumn("OuterShellThickness_mm");                      // 117
  am->CreateNtupleSColumn("OuterShellProfile");                           // 118
  am->CreateNtupleSColumn("OuterShellProfileG4Material");                 // 119
  am->CreateNtupleSColumn("OuterShellProfileRindexDataset");              // 120
  am->CreateNtupleSColumn("OuterShellProfileAbsDataset");                 // 121
  am->CreateNtupleSColumn("OuterShellProfileRayleighDataset");            // 122
  am->CreateNtupleFColumn("OuterShellProfileRindexValidMin_nm");          // 123
  am->CreateNtupleFColumn("OuterShellProfileRindexValidMax_nm");          // 124
  am->CreateNtupleFColumn("OuterShellProfileAbsValidMin_nm");             // 125
  am->CreateNtupleFColumn("OuterShellProfileAbsValidMax_nm");             // 126
  am->CreateNtupleFColumn("OuterShellProfileRayleighValidMin_nm");        // 127
  am->CreateNtupleFColumn("OuterShellProfileRayleighValidMax_nm");        // 128
  am->CreateNtupleSColumn("OuterShellProfileRindexStatus");               // 129
  am->CreateNtupleSColumn("OuterShellProfileAbsStatus");                  // 130
  am->CreateNtupleSColumn("OuterShellProfileRayleighStatus");             // 131
  am->CreateNtupleSColumn("OuterShellProfileOutOfRangePolicy");           // 132

  am->CreateNtupleSColumn("CouplingProfile");                             // 133
  am->CreateNtupleSColumn("CouplingProfileG4Material");                   // 134
  am->CreateNtupleSColumn("CouplingProfileRindexDataset");                // 135
  am->CreateNtupleSColumn("CouplingProfileAbsDataset");                   // 136
  am->CreateNtupleSColumn("CouplingProfileRayleighDataset");              // 137
  am->CreateNtupleFColumn("CouplingProfileRindexValidMin_nm");            // 138
  am->CreateNtupleFColumn("CouplingProfileRindexValidMax_nm");            // 139
  am->CreateNtupleFColumn("CouplingProfileAbsValidMin_nm");               // 140
  am->CreateNtupleFColumn("CouplingProfileAbsValidMax_nm");               // 141
  am->CreateNtupleFColumn("CouplingProfileRayleighValidMin_nm");          // 142
  am->CreateNtupleFColumn("CouplingProfileRayleighValidMax_nm");          // 143
  am->CreateNtupleSColumn("CouplingProfileRindexStatus");                 // 144
  am->CreateNtupleSColumn("CouplingProfileAbsStatus");                    // 145
  am->CreateNtupleSColumn("CouplingProfileRayleighStatus");               // 146
  am->CreateNtupleSColumn("CouplingProfileOutOfRangePolicy");             // 147

  // general-v1.8: explicitly declared empirical ABS validity.
  // Existing *AbsValidMin/Max fields remain the numerical/tabulated
  // transport-support range of the active dataset.
  am->CreateNtupleFColumn("CoreProfileAbsMeasuredMin_nm");                  // 148
  am->CreateNtupleFColumn("CoreProfileAbsMeasuredMax_nm");                  // 149
  am->CreateNtupleFColumn("OuterShellProfileAbsMeasuredMin_nm");            // 150
  am->CreateNtupleFColumn("OuterShellProfileAbsMeasuredMax_nm");            // 151
  am->CreateNtupleFColumn("CouplingProfileAbsMeasuredMin_nm");              // 152
  am->CreateNtupleFColumn("CouplingProfileAbsMeasuredMax_nm");              // 153

  // general-v1.9: append-only Core-boundary recorder metadata.  The complete
  // v1.8 RunConfig prefix (columns 0..153) remains unchanged.
  am->CreateNtupleSColumn("CoreBoundaryMode");                              // 154
  am->CreateNtupleSColumn("CoreBoundaryDetail");                            // 155
  am->CreateNtupleSColumn("CoreBoundarySemantic");                          // 156

  // general-v1.10: append-only build/source provenance.  The complete v1.9
  // RunConfig prefix (columns 0..156) remains unchanged.
  am->CreateNtupleSColumn("GitCommit");                                     // 157
  am->CreateNtupleSColumn("GitBranch");                                     // 158
  am->CreateNtupleIColumn("GitDirty");                                      // 159

  am->FinishNtuple(kRunConfig);

  am->CreateNtuple("ChannelMap", "Readout channel geometry");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("ChannelID");
  am->CreateNtupleIColumn("FaceID");
  am->CreateNtupleSColumn("FaceName");
  am->CreateNtupleIColumn("GridU");
  am->CreateNtupleIColumn("GridV");
  am->CreateNtupleFColumn("LocalU_mm");
  am->CreateNtupleFColumn("LocalV_mm");
  am->CreateNtupleFColumn("CenterX_mm");
  am->CreateNtupleFColumn("CenterY_mm");
  am->CreateNtupleFColumn("CenterZ_mm");
  am->CreateNtupleFColumn("NormalX");
  am->CreateNtupleFColumn("NormalY");
  am->CreateNtupleFColumn("NormalZ");
  am->CreateNtupleFColumn("SizeU_mm");
  am->CreateNtupleFColumn("SizeV_mm");
  am->CreateNtupleFColumn("EntranceCenterX_mm");
  am->CreateNtupleFColumn("EntranceCenterY_mm");
  am->CreateNtupleFColumn("EntranceCenterZ_mm");
  am->FinishNtuple(kChannelMap);

  am->CreateNtuple("VolumeMap", "Stable volume category map for this schema");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("VolumeID");
  am->CreateNtupleSColumn("VolumeName");
  am->CreateNtupleSColumn("MaterialName");
  am->FinishNtuple(kVolumeMap);

  am->CreateNtuple("RunStatistics", "Per-worker run counters; sum rows for run totals");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("WorkerThreadID");
  am->CreateNtupleIColumn("GeneratedEvents");
  am->CreateNtupleIColumn("TriggeredEvents");
  am->CreateNtupleIColumn("DetailedEvents");
  am->CreateNtupleDColumn("CRYSimulatedTime_s");
  am->FinishNtuple(kRunStatistics);

  am->CreateNtuple("BoundaryDebug", "Optical boundary log plus tagged first-Rayleigh debug rows; debug only");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("PhotonID");
  am->CreateNtupleIColumn("StepNumber");
  am->CreateNtupleIColumn("PreVolumeID");
  am->CreateNtupleIColumn("PostVolumeID");
  am->CreateNtupleIColumn("BoundaryStatus");
  am->CreateNtupleFColumn("Time_ns");
  am->CreateNtupleFColumn("X_mm");
  am->CreateNtupleFColumn("Y_mm");
  am->CreateNtupleFColumn("Z_mm");
  am->CreateNtupleFColumn("PreDirX");
  am->CreateNtupleFColumn("PreDirY");
  am->CreateNtupleFColumn("PreDirZ");
  am->CreateNtupleFColumn("PostDirX");
  am->CreateNtupleFColumn("PostDirY");
  am->CreateNtupleFColumn("PostDirZ");
  am->CreateNtupleFColumn("PrePolX");
  am->CreateNtupleFColumn("PrePolY");
  am->CreateNtupleFColumn("PrePolZ");
  am->CreateNtupleFColumn("PostPolX");
  am->CreateNtupleFColumn("PostPolY");
  am->CreateNtupleFColumn("PostPolZ");
  am->CreateNtupleSColumn("BoundaryStatusName");
  am->FinishNtuple(kBoundaryDebug);

  am->CreateNtuple("ReadoutPlaneConfig", "Configured readout layout per detector face");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("FaceID");
  am->CreateNtupleSColumn("FaceName");
  am->CreateNtupleIColumn("Enabled");
  am->CreateNtupleIColumn("AutoFill");
  am->CreateNtupleIColumn("NU");
  am->CreateNtupleIColumn("NV");
  am->CreateNtupleFColumn("PitchU_mm");
  am->CreateNtupleFColumn("PitchV_mm");
  am->CreateNtupleFColumn("OffsetU_mm");
  am->CreateNtupleFColumn("OffsetV_mm");
  am->CreateNtupleFColumn("GapU_mm");
  am->CreateNtupleFColumn("GapV_mm");
  am->CreateNtupleFColumn("MarginU_mm");
  am->CreateNtupleFColumn("MarginV_mm");
  am->FinishNtuple(kReadoutPlaneConfig);

  am->CreateNtuple("PrimaryTruth", "Lightweight non-optical primary truth aligned with the thin event tier");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("TrackID");
  am->CreateNtupleIColumn("PDG");
  am->CreateNtupleFColumn("VertexKE_MeV");
  am->CreateNtupleFColumn("VertexTime_ns");
  am->CreateNtupleFColumn("VertexX_mm");
  am->CreateNtupleFColumn("VertexY_mm");
  am->CreateNtupleFColumn("VertexZ_mm");
  am->CreateNtupleFColumn("DirX");
  am->CreateNtupleFColumn("DirY");
  am->CreateNtupleFColumn("DirZ");
  am->FinishNtuple(kPrimaryTruth);

  am->CreateNtuple("CoreBoundaryCrossing",
                   "First successful optical-photon exit from the Core logical volume");
  am->CreateNtupleIColumn("RunID");
  am->CreateNtupleIColumn("EventID");
  am->CreateNtupleIColumn("PhotonID");
  am->CreateNtupleIColumn("ParentTrackID");
  am->CreateNtupleIColumn("SourceTrackID");
  am->CreateNtupleIColumn("RootPrimaryTrackID");
  am->CreateNtupleIColumn("OriginVolumeID");
  am->CreateNtupleIColumn("FaceID");
  am->CreateNtupleFColumn("X_mm");
  am->CreateNtupleFColumn("Y_mm");
  am->CreateNtupleFColumn("Z_mm");
  am->CreateNtupleFColumn("LocalU_mm");
  am->CreateNtupleFColumn("LocalV_mm");
  am->CreateNtupleFColumn("DirX");
  am->CreateNtupleFColumn("DirY");
  am->CreateNtupleFColumn("DirZ");
  am->CreateNtupleFColumn("CosThetaNormal");
  am->CreateNtupleFColumn("Wavelength_nm");
  am->CreateNtupleFColumn("Time_ns");
  am->CreateNtupleFColumn("PolX");
  am->CreateNtupleFColumn("PolY");
  am->CreateNtupleFColumn("PolZ");
  am->FinishNtuple(kCoreBoundaryCrossing);

}

RunAction::~RunAction() = default;

void RunAction::CountEvent(G4bool triggered, G4bool detailsWritten) {
  fGeneratedEvents += 1;
  if (triggered) fTriggeredEvents += 1;
  if (detailsWritten) fDetailedEvents += 1;
}

void RunAction::BeginOfRunAction(const G4Run* run) {
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  G4AccumulableManager::Instance()->Reset();

  auto* det = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  auto* am = G4AnalysisManager::Instance();

  if (det) {
    am->SetCompressionLevel(det->GetCompressionLevel());
    const G4String file = det->GetOutputFileName() + "_run" +
                          std::to_string(run->GetRunID());
    am->SetFileName(file);
  }
  am->OpenFile();

  const auto* gen = dynamic_cast<const PrimaryGeneratorAction*>(
      G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  fCRYTimeStart_s = (gen && gen->GetGeneratorMode() == "CRY")
      ? gen->GetCRYTimeSimulated() : 0.0;

  if (!det || !IsSingleWriterWorker()) return;

  const auto core = det->GetCoreFullSize();
  const auto scint = det->GetScintSize();
  const auto readout = det->GetReadoutSize();
  const auto coupling = det->GetCouplingSize();
  const auto hole = det->GetTyvekHoleSize();

  am->FillNtupleIColumn(kRunConfig,0,run->GetRunID());
  am->FillNtupleIColumn(kRunConfig,1,static_cast<G4int>(det->GetGeometryVersion()));
  am->FillNtupleFColumn(kRunConfig,2,core.x()/mm);
  am->FillNtupleFColumn(kRunConfig,3,core.y()/mm);
  am->FillNtupleFColumn(kRunConfig,4,core.z()/mm);
  am->FillNtupleFColumn(kRunConfig,5,scint.x()/mm);
  am->FillNtupleFColumn(kRunConfig,6,scint.y()/mm);
  am->FillNtupleFColumn(kRunConfig,7,scint.z()/mm);
  am->FillNtupleFColumn(kRunConfig,8,readout.x()/mm);
  am->FillNtupleFColumn(kRunConfig,9,readout.y()/mm);
  am->FillNtupleFColumn(kRunConfig,10,readout.z()/mm);
  am->FillNtupleFColumn(kRunConfig,11,coupling.x()/mm);
  am->FillNtupleFColumn(kRunConfig,12,coupling.y()/mm);
  am->FillNtupleFColumn(kRunConfig,13,coupling.z()/mm);
  am->FillNtupleFColumn(kRunConfig,14,hole.x()/mm);
  am->FillNtupleFColumn(kRunConfig,15,hole.y()/mm);
  // R1 keeps legacy RunConfig column names/indices unchanged for regression.
  // Values are sourced from the new generic outer-shell API.
  am->FillNtupleFColumn(kRunConfig,16,det->GetOuterShellThickness()/mm);
  am->FillNtupleFColumn(kRunConfig,17,det->GetTyvekThickness()/mm);
  am->FillNtupleFColumn(kRunConfig,18,det->GetScintGap()/mm);
  am->FillNtupleFColumn(kRunConfig,19,det->GetSourceDistance()/mm);
  am->FillNtupleSColumn(kRunConfig,20,det->GetTyvekConfig());
  am->FillNtupleSColumn(kRunConfig,21,det->GetCoreMaterialName());
  am->FillNtupleIColumn(kRunConfig,22,det->GetUseOuterShell()?1:0);
  am->FillNtupleSColumn(kRunConfig,23,det->GetSensorMode());
  am->FillNtupleFColumn(kRunConfig,24,det->GetScintThreshold()/MeV);
  am->FillNtupleFColumn(kRunConfig,25,det->GetCoincidenceWindow()/ns);
  am->FillNtupleFColumn(kRunConfig,26,det->GetScintIntegrationWindow()/ns);
  am->FillNtupleSColumn(kRunConfig,27,det->GetOutputEventMode());
  am->FillNtupleSColumn(kRunConfig,28,det->GetIncidentDetail());
  am->FillNtupleSColumn(kRunConfig,29,det->GetPhotonBirthDetail());
  am->FillNtupleIColumn(kRunConfig,30,det->GetBoundaryDebug()?1:0);
  am->FillNtupleIColumn(kRunConfig,31,det->GetMaxOpticalBounces());
  am->FillNtupleIColumn(kRunConfig,32,det->GetChannelCount());
  am->FillNtupleIColumn(kRunConfig,33,det->GetCompressionLevel());
  am->FillNtupleSColumn(kRunConfig,34,det->GetAcrylicRindexDataset());
  am->FillNtupleSColumn(kRunConfig,35,det->GetAcrylicAbsDataset());
  am->FillNtupleSColumn(kRunConfig,36,det->GetTyvekDataset());
  am->FillNtupleSColumn(kRunConfig,37,gen ? gen->GetGeneratorMode() : "Unknown");
  am->FillNtupleSColumn(kRunConfig,38,gen ? gen->GetTimeMode() : "Unknown");
  am->FillNtupleSColumn(kRunConfig,39,gen ? gen->GetCRYConfigString() : "");
  am->FillNtupleSColumn(kRunConfig,40,G4Version);
  am->FillNtupleSColumn(kRunConfig,41,"FTFP_BERT+G4EmStandardPhysics_option4+G4OpticalPhysics");
  am->FillNtupleDColumn(kRunConfig,42,static_cast<G4double>(CherenkovMetadata::InitialSeed.load()));
  am->FillNtupleSColumn(kRunConfig,43,det->GetAirRindexDataset());
  am->FillNtupleSColumn(kRunConfig,44,det->GetWaterRindexDataset());
  am->FillNtupleSColumn(kRunConfig,45,det->GetWaterAbsDataset());
  am->FillNtupleSColumn(kRunConfig,46,det->GetWaterRayleighDataset());
  am->FillNtupleSColumn(kRunConfig,47,det->GetCouplingRindexDataset());
  am->FillNtupleSColumn(kRunConfig,48,det->GetSensorPdeDataset());

  G4ThreeVector gunPos;
  G4ThreeVector gunDir(0.,-1.,0.);
  G4double gunEnergy = 0.0;
  G4double gunTime = 0.0;
  G4bool gunAutoY = false;
  G4ThreeVector gunPol(1.,0.,0.);
  G4String gunParticle;
  if (gen) {
    gunPos = gen->GetGunPosition();
    gunDir = gen->GetGunDirection();
    gunEnergy = gen->GetGunEnergy();
    gunTime = gen->GetGunTime();
    gunAutoY = gen->GetGunUseAutoY();
    gunPol = gen->GetGunPolarization();
    gunParticle = gen->GetGunParticleName();
    if (gunAutoY) gunPos.setY(det->GetGunPosY());
  }
  am->FillNtupleSColumn(kRunConfig,49,gunParticle);
  am->FillNtupleFColumn(kRunConfig,50,gunEnergy/MeV);
  am->FillNtupleFColumn(kRunConfig,51,gunPos.x()/mm);
  am->FillNtupleFColumn(kRunConfig,52,gunPos.y()/mm);
  am->FillNtupleFColumn(kRunConfig,53,gunPos.z()/mm);
  am->FillNtupleFColumn(kRunConfig,54,gunDir.x());
  am->FillNtupleFColumn(kRunConfig,55,gunDir.y());
  am->FillNtupleFColumn(kRunConfig,56,gunDir.z());
  am->FillNtupleFColumn(kRunConfig,57,gunTime/ns);
  am->FillNtupleIColumn(kRunConfig,58,gunAutoY?1:0);
  am->FillNtupleFColumn(kRunConfig,59,gunPol.x());
  am->FillNtupleFColumn(kRunConfig,60,gunPol.y());
  am->FillNtupleFColumn(kRunConfig,61,gunPol.z());

  auto* runManager = G4RunManager::GetRunManager();
  G4int nThreads = 1;
#ifdef G4MULTITHREADED
  if (G4Threading::IsMultithreadedApplication()) {
    if (auto* masterRunManager = G4MTRunManager::GetMasterRunManager()) {
      nThreads = masterRunManager->GetNumberOfThreads();
    }
  }
#endif
  const auto* physics = runManager->GetUserPhysicsList();
  const G4double defaultCut = physics ? physics->GetDefaultCutValue() : -1.0;
  const auto* optical = G4OpticalParameters::Instance();

  am->FillNtupleSColumn(kRunConfig,62,"general-v1.10");
  am->FillNtupleIColumn(kRunConfig,63,nThreads);
  am->FillNtupleSColumn(kRunConfig,64,"CLHEP::RanecuEngine");
  am->FillNtupleSColumn(kRunConfig,65,"1.7");
  am->FillNtupleFColumn(kRunConfig,66,defaultCut >= 0.0 ? defaultCut/mm : -1.0);
  am->FillNtupleIColumn(kRunConfig,67,optical->GetProcessActivation("Cerenkov") ? 1 : 0);
  am->FillNtupleIColumn(kRunConfig,68,optical->GetProcessActivation("Scintillation") ? 1 : 0);
  am->FillNtupleIColumn(kRunConfig,69,optical->GetProcessActivation("OpAbsorption") ? 1 : 0);
  am->FillNtupleIColumn(kRunConfig,70,optical->GetProcessActivation("OpRayleigh") ? 1 : 0);
  am->FillNtupleIColumn(kRunConfig,71,optical->GetProcessActivation("OpMieHG") ? 1 : 0);
  am->FillNtupleIColumn(kRunConfig,72,optical->GetProcessActivation("OpBoundary") ? 1 : 0);
  am->FillNtupleIColumn(kRunConfig,73,optical->GetProcessActivation("OpWLS") ? 1 : 0);
  am->FillNtupleIColumn(kRunConfig,74,optical->GetProcessActivation("OpWLS2") ? 1 : 0);
  const auto acrAbsRange = det->GetDatasetWavelengthRange(det->GetAcrylicAbsDataset());
  const auto acrRindexRange = det->GetDatasetWavelengthRange(det->GetAcrylicRindexDataset());
  const auto coupRindexRange = det->GetDatasetWavelengthRange(det->GetCouplingRindexDataset());
  const auto coupAbsRange = det->GetDatasetWavelengthRange(det->GetCouplingAbsDataset());
  const auto waterRindexRange = det->GetDatasetWavelengthRange(det->GetWaterRindexDataset());
  const auto waterAbsRange = det->GetDatasetWavelengthRange(det->GetWaterAbsDataset());
  const auto waterRayleighRange = det->GetDatasetWavelengthRange(det->GetWaterRayleighDataset());

  // Profile-aware role metadata.  Empty optional datasets are represented by
  // empty strings and validity range (-1,-1), rather than being passed to the
  // file-range loader.
  const auto* waterProfile = det->GetMaterialProfile("Water_Default");
  const auto* coreProfile = det->GetMaterialProfile(det->GetCoreProfileName());
  const auto* shellProfile = det->GetMaterialProfile(det->GetOuterShellProfileName());
  const auto* couplingProfile = det->GetMaterialProfile(det->GetCouplingProfileName());
  const auto rangeOrUnset = [&](const G4String& dataset) -> std::pair<G4double,G4double> {
    if (dataset.empty()) return {-1.0, -1.0};
    return det->GetDatasetWavelengthRange(dataset);
  };
  const auto coreRindexRange = rangeOrUnset(coreProfile ? coreProfile->rindexDataset : "");
  const auto coreAbsRange = rangeOrUnset(coreProfile ? coreProfile->absDataset : "");
  const auto coreRayleighRange = rangeOrUnset(coreProfile ? coreProfile->rayleighDataset : "");
  const auto shellRindexRange = rangeOrUnset(shellProfile ? shellProfile->rindexDataset : "");
  const auto shellAbsRange = rangeOrUnset(shellProfile ? shellProfile->absDataset : "");
  const auto shellRayleighRange = rangeOrUnset(shellProfile ? shellProfile->rayleighDataset : "");
  const auto couplingProfileRindexRange =
      rangeOrUnset(couplingProfile ? couplingProfile->rindexDataset : "");
  const auto couplingProfileAbsRange =
      rangeOrUnset(couplingProfile ? couplingProfile->absDataset : "");
  const auto couplingProfileRayleighRange =
      rangeOrUnset(couplingProfile ? couplingProfile->rayleighDataset : "");

  // Explicit empirical validity is distinct from numerical dataset support.
  // Do NOT fall back to the dataset range here: (-1,-1) means that the
  // profile has not declared a measured ABS domain.
  const auto measuredAbsRangeOrUnset =
      [](const MaterialProfile* profile) -> std::pair<G4double,G4double> {
    if (!profile ||
        profile->absMeasuredMinNm <= 0.0 ||
        profile->absMeasuredMaxNm <= profile->absMeasuredMinNm) {
      return {-1.0, -1.0};
    }
    return {profile->absMeasuredMinNm, profile->absMeasuredMaxNm};
  };

  const auto coreAbsMeasuredRange =
      measuredAbsRangeOrUnset(coreProfile);
  const auto shellAbsMeasuredRange =
      measuredAbsRangeOrUnset(shellProfile);
  const auto couplingAbsMeasuredRange =
      measuredAbsRangeOrUnset(couplingProfile);

  am->FillNtupleFColumn(kRunConfig,75,static_cast<G4float>(acrAbsRange.first));
  am->FillNtupleFColumn(kRunConfig,76,static_cast<G4float>(acrAbsRange.second));
  am->FillNtupleSColumn(kRunConfig,77,"Geant4EndpointClamp; validity stored explicitly");
  am->FillNtupleSColumn(kRunConfig,78,det->GetCouplingAbsDataset());
  am->FillNtupleSColumn(kRunConfig,79,"OpticalGrease_PDMS_Generic");
  am->FillNtupleSColumn(kRunConfig,80,"Hamamatsu_S13360-3050CS_effective_PDE");
  am->FillNtupleSColumn(kRunConfig,81,CherenkovBuildInfo::ProjectRevision);
  am->FillNtupleFColumn(kRunConfig,82,static_cast<G4float>(acrRindexRange.first));
  am->FillNtupleFColumn(kRunConfig,83,static_cast<G4float>(acrRindexRange.second));
  am->FillNtupleFColumn(kRunConfig,84,static_cast<G4float>(coupRindexRange.first));
  am->FillNtupleFColumn(kRunConfig,85,static_cast<G4float>(coupRindexRange.second));
  am->FillNtupleFColumn(kRunConfig,86,static_cast<G4float>(coupAbsRange.first));
  am->FillNtupleFColumn(kRunConfig,87,static_cast<G4float>(coupAbsRange.second));
  am->FillNtupleFColumn(kRunConfig,88,static_cast<G4float>(waterRindexRange.first));
  am->FillNtupleFColumn(kRunConfig,89,static_cast<G4float>(waterRindexRange.second));
  am->FillNtupleFColumn(kRunConfig,90,static_cast<G4float>(waterAbsRange.first));
  am->FillNtupleFColumn(kRunConfig,91,static_cast<G4float>(waterAbsRange.second));
  am->FillNtupleFColumn(kRunConfig,92,static_cast<G4float>(waterRayleighRange.first));
  am->FillNtupleFColumn(kRunConfig,93,static_cast<G4float>(waterRayleighRange.second));
  am->FillNtupleSColumn(kRunConfig,94,waterProfile ? waterProfile->rayleighStatus : "");
  am->FillNtupleSColumn(kRunConfig,95,"Geant4EndpointClamp; compare photon wavelengths to per-dataset valid ranges");
  am->FillNtupleSColumn(kRunConfig,96,"CoreLogicalVolumeBoundary");
  am->FillNtupleSColumn(kRunConfig,97,"Hamamatsu_typical_PDE_curve_Vover=3V");
  am->FillNtupleSColumn(kRunConfig,98,det->GetOpticalDataFingerprint());

  am->FillNtupleSColumn(kRunConfig,99,det->GetMaterialProfileRegistryVersion());
  am->FillNtupleSColumn(kRunConfig,100,"general-v1.6 RunConfig columns 0-98 preserved");

  am->FillNtupleSColumn(kRunConfig,101,coreProfile ? coreProfile->profileName : "");
  am->FillNtupleSColumn(kRunConfig,102,
      (coreProfile && coreProfile->material) ? coreProfile->material->GetName() : "");
  am->FillNtupleSColumn(kRunConfig,103,coreProfile ? coreProfile->rindexDataset : "");
  am->FillNtupleSColumn(kRunConfig,104,coreProfile ? coreProfile->absDataset : "");
  am->FillNtupleSColumn(kRunConfig,105,coreProfile ? coreProfile->rayleighDataset : "");
  am->FillNtupleFColumn(kRunConfig,106,static_cast<G4float>(coreRindexRange.first));
  am->FillNtupleFColumn(kRunConfig,107,static_cast<G4float>(coreRindexRange.second));
  am->FillNtupleFColumn(kRunConfig,108,static_cast<G4float>(coreAbsRange.first));
  am->FillNtupleFColumn(kRunConfig,109,static_cast<G4float>(coreAbsRange.second));
  am->FillNtupleFColumn(kRunConfig,110,static_cast<G4float>(coreRayleighRange.first));
  am->FillNtupleFColumn(kRunConfig,111,static_cast<G4float>(coreRayleighRange.second));
  am->FillNtupleSColumn(kRunConfig,112,coreProfile ? coreProfile->rindexStatus : "");
  am->FillNtupleSColumn(kRunConfig,113,coreProfile ? coreProfile->absStatus : "");
  am->FillNtupleSColumn(kRunConfig,114,coreProfile ? coreProfile->rayleighStatus : "");
  am->FillNtupleSColumn(kRunConfig,115,coreProfile ? coreProfile->outOfRangePolicy : "");

  am->FillNtupleIColumn(kRunConfig,116,det->GetUseOuterShell()?1:0);
  am->FillNtupleFColumn(kRunConfig,117,det->GetOuterShellThickness()/mm);
  am->FillNtupleSColumn(kRunConfig,118,shellProfile ? shellProfile->profileName : "");
  am->FillNtupleSColumn(kRunConfig,119,
      (shellProfile && shellProfile->material) ? shellProfile->material->GetName() : "");
  am->FillNtupleSColumn(kRunConfig,120,shellProfile ? shellProfile->rindexDataset : "");
  am->FillNtupleSColumn(kRunConfig,121,shellProfile ? shellProfile->absDataset : "");
  am->FillNtupleSColumn(kRunConfig,122,shellProfile ? shellProfile->rayleighDataset : "");
  am->FillNtupleFColumn(kRunConfig,123,static_cast<G4float>(shellRindexRange.first));
  am->FillNtupleFColumn(kRunConfig,124,static_cast<G4float>(shellRindexRange.second));
  am->FillNtupleFColumn(kRunConfig,125,static_cast<G4float>(shellAbsRange.first));
  am->FillNtupleFColumn(kRunConfig,126,static_cast<G4float>(shellAbsRange.second));
  am->FillNtupleFColumn(kRunConfig,127,static_cast<G4float>(shellRayleighRange.first));
  am->FillNtupleFColumn(kRunConfig,128,static_cast<G4float>(shellRayleighRange.second));
  am->FillNtupleSColumn(kRunConfig,129,shellProfile ? shellProfile->rindexStatus : "");
  am->FillNtupleSColumn(kRunConfig,130,shellProfile ? shellProfile->absStatus : "");
  am->FillNtupleSColumn(kRunConfig,131,shellProfile ? shellProfile->rayleighStatus : "");
  am->FillNtupleSColumn(kRunConfig,132,shellProfile ? shellProfile->outOfRangePolicy : "");

  am->FillNtupleSColumn(kRunConfig,133,couplingProfile ? couplingProfile->profileName : "");
  am->FillNtupleSColumn(kRunConfig,134,
      (couplingProfile && couplingProfile->material) ? couplingProfile->material->GetName() : "");
  am->FillNtupleSColumn(kRunConfig,135,couplingProfile ? couplingProfile->rindexDataset : "");
  am->FillNtupleSColumn(kRunConfig,136,couplingProfile ? couplingProfile->absDataset : "");
  am->FillNtupleSColumn(kRunConfig,137,couplingProfile ? couplingProfile->rayleighDataset : "");
  am->FillNtupleFColumn(kRunConfig,138,static_cast<G4float>(couplingProfileRindexRange.first));
  am->FillNtupleFColumn(kRunConfig,139,static_cast<G4float>(couplingProfileRindexRange.second));
  am->FillNtupleFColumn(kRunConfig,140,static_cast<G4float>(couplingProfileAbsRange.first));
  am->FillNtupleFColumn(kRunConfig,141,static_cast<G4float>(couplingProfileAbsRange.second));
  am->FillNtupleFColumn(kRunConfig,142,static_cast<G4float>(couplingProfileRayleighRange.first));
  am->FillNtupleFColumn(kRunConfig,143,static_cast<G4float>(couplingProfileRayleighRange.second));
  am->FillNtupleSColumn(kRunConfig,144,couplingProfile ? couplingProfile->rindexStatus : "");
  am->FillNtupleSColumn(kRunConfig,145,couplingProfile ? couplingProfile->absStatus : "");
  am->FillNtupleSColumn(kRunConfig,146,couplingProfile ? couplingProfile->rayleighStatus : "");
  am->FillNtupleSColumn(kRunConfig,147,couplingProfile ? couplingProfile->outOfRangePolicy : "");

  am->FillNtupleFColumn(
      kRunConfig,148,static_cast<G4float>(coreAbsMeasuredRange.first));
  am->FillNtupleFColumn(
      kRunConfig,149,static_cast<G4float>(coreAbsMeasuredRange.second));

  am->FillNtupleFColumn(
      kRunConfig,150,static_cast<G4float>(shellAbsMeasuredRange.first));
  am->FillNtupleFColumn(
      kRunConfig,151,static_cast<G4float>(shellAbsMeasuredRange.second));

  am->FillNtupleFColumn(
      kRunConfig,152,static_cast<G4float>(couplingAbsMeasuredRange.first));
  am->FillNtupleFColumn(
      kRunConfig,153,static_cast<G4float>(couplingAbsMeasuredRange.second));

  am->FillNtupleSColumn(kRunConfig,154,det->GetCoreBoundaryMode());
  am->FillNtupleSColumn(kRunConfig,155,det->GetCoreBoundaryDetail());
  am->FillNtupleSColumn(
      kRunConfig,156,
      "FirstSuccessfulCoreExit; post-boundary position/direction/polarization; "
      "FaceID convention +X,-X,+Z,-Z,+Y,-Y");
  am->FillNtupleSColumn(kRunConfig,157,CherenkovBuildInfo::GitCommit);
  am->FillNtupleSColumn(kRunConfig,158,CherenkovBuildInfo::GitBranch);
  am->FillNtupleIColumn(kRunConfig,159,CherenkovBuildInfo::GitDirty);

  am->AddNtupleRow(kRunConfig);

  const auto& planes = det->GetPlaneConfigs();
  for (G4int face=0; face<6; ++face) {
    const auto& p = planes[face];
    am->FillNtupleIColumn(kReadoutPlaneConfig,0,run->GetRunID());
    am->FillNtupleIColumn(kReadoutPlaneConfig,1,face);
    am->FillNtupleSColumn(kReadoutPlaneConfig,2,det->FaceName(face));
    am->FillNtupleIColumn(kReadoutPlaneConfig,3,p.enabled?1:0);
    am->FillNtupleIColumn(kReadoutPlaneConfig,4,p.autoFill?1:0);
    am->FillNtupleIColumn(kReadoutPlaneConfig,5,p.nU);
    am->FillNtupleIColumn(kReadoutPlaneConfig,6,p.nV);
    am->FillNtupleFColumn(kReadoutPlaneConfig,7,p.pitchU/mm);
    am->FillNtupleFColumn(kReadoutPlaneConfig,8,p.pitchV/mm);
    am->FillNtupleFColumn(kReadoutPlaneConfig,9,p.offsetU/mm);
    am->FillNtupleFColumn(kReadoutPlaneConfig,10,p.offsetV/mm);
    am->FillNtupleFColumn(kReadoutPlaneConfig,11,p.gapU/mm);
    am->FillNtupleFColumn(kReadoutPlaneConfig,12,p.gapV/mm);
    am->FillNtupleFColumn(kReadoutPlaneConfig,13,p.marginU/mm);
    am->FillNtupleFColumn(kReadoutPlaneConfig,14,p.marginV/mm);
    am->AddNtupleRow(kReadoutPlaneConfig);
  }

  for (const auto& ch : det->GetChannelMap()) {
    am->FillNtupleIColumn(kChannelMap,0,run->GetRunID());
    am->FillNtupleIColumn(kChannelMap,1,ch.channelID);
    am->FillNtupleIColumn(kChannelMap,2,ch.faceID);
    am->FillNtupleSColumn(kChannelMap,3,det->FaceName(ch.faceID));
    am->FillNtupleIColumn(kChannelMap,4,ch.gridU);
    am->FillNtupleIColumn(kChannelMap,5,ch.gridV);
    am->FillNtupleFColumn(kChannelMap,6,ch.localU/mm);
    am->FillNtupleFColumn(kChannelMap,7,ch.localV/mm);
    am->FillNtupleFColumn(kChannelMap,8,ch.center.x()/mm);
    am->FillNtupleFColumn(kChannelMap,9,ch.center.y()/mm);
    am->FillNtupleFColumn(kChannelMap,10,ch.center.z()/mm);
    am->FillNtupleFColumn(kChannelMap,11,ch.normal.x());
    am->FillNtupleFColumn(kChannelMap,12,ch.normal.y());
    am->FillNtupleFColumn(kChannelMap,13,ch.normal.z());
    am->FillNtupleFColumn(kChannelMap,14,ch.sizeU/mm);
    am->FillNtupleFColumn(kChannelMap,15,ch.sizeV/mm);
    am->FillNtupleFColumn(kChannelMap,16,ch.entranceCenter.x()/mm);
    am->FillNtupleFColumn(kChannelMap,17,ch.entranceCenter.y()/mm);
    am->FillNtupleFColumn(kChannelMap,18,ch.entranceCenter.z()/mm);
    am->AddNtupleRow(kChannelMap);
  }

  for (const auto& v : det->GetVolumeMap()) {
    am->FillNtupleIColumn(kVolumeMap,0,run->GetRunID());
    am->FillNtupleIColumn(kVolumeMap,1,v.volumeID);
    am->FillNtupleSColumn(kVolumeMap,2,v.name);
    am->FillNtupleSColumn(kVolumeMap,3,v.materialName);
    am->AddNtupleRow(kVolumeMap);
  }
}

void RunAction::EndOfRunAction(const G4Run* run) {
  auto* accum = G4AccumulableManager::Instance();
  auto* am = G4AnalysisManager::Instance();

  const G4bool eventProcessingInstance =
      !G4Threading::IsMultithreadedApplication() || !IsMaster();

  // ROOT ntuple rows are produced on event-processing workers and merged by
  // G4AnalysisManager.  Do not FillNtuple... on the master: in ROOT ntuple
  // merging mode the master owns the merge target rather than a fillable
  // worker ntuple, which otherwise produces Analysis_W001 at end of run.
  if (eventProcessingInstance) {
    const auto* gen = dynamic_cast<const PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    if (gen && gen->GetGeneratorMode() == "CRY") {
      const G4double elapsed = gen->GetCRYTimeSimulated() - fCRYTimeStart_s;
      if (elapsed >= 0.0) fCRYSimulatedTime_s += elapsed;
    }

    const G4int workerID = G4Threading::IsMultithreadedApplication()
        ? G4Threading::G4GetThreadId() : 0;
    am->FillNtupleIColumn(kRunStatistics,0,run->GetRunID());
    am->FillNtupleIColumn(kRunStatistics,1,workerID);
    am->FillNtupleIColumn(kRunStatistics,2,fGeneratedEvents.GetValue());
    am->FillNtupleIColumn(kRunStatistics,3,fTriggeredEvents.GetValue());
    am->FillNtupleIColumn(kRunStatistics,4,fDetailedEvents.GetValue());
    am->FillNtupleDColumn(kRunStatistics,5,fCRYSimulatedTime_s.GetValue());
    am->AddNtupleRow(kRunStatistics);
  }

  // Keep the accumulables merged on the master for a correct global state and
  // future run-level summaries, but the persisted RunStatistics totals are
  // reconstructed by summing the per-worker rows above.
  if (IsMaster()) accum->Merge();

  am->Write();
  am->CloseFile();
}
