#ifndef EventAction_h
#define EventAction_h 1

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class G4Event;
class G4Step;
class G4Track;
class G4VProcess;
class RunAction;

struct ScintDeposit {
  G4int trackID = 0;
  G4double edepMeV = 0.0;
  G4double timeNs = 0.0;
};

struct TrackRecord {
  G4int trackID = 0;
  G4int parentID = 0;
  G4int rootPrimaryTrackID = 0;
  G4int pdg = 0;
  G4float chargeE = 0.0f;
  G4int creatorProcessType = -1;
  G4int creatorProcessSubType = -1;
  G4String creatorProcessName = "Primary";
  G4int birthVolumeID = -1;

  G4float vertexKE_MeV = 0.0f;
  G4float vertexTime_ns = 0.0f;
  G4float vertexX_mm = 0.0f;
  G4float vertexY_mm = 0.0f;
  G4float vertexZ_mm = 0.0f;
  G4float dirX = 0.0f;
  G4float dirY = 0.0f;
  G4float dirZ = 0.0f;

  G4int coreEntrySet = 0;
  G4float coreEntryX_mm = 0.0f;
  G4float coreEntryY_mm = 0.0f;
  G4float coreEntryZ_mm = 0.0f;
  G4int coreExitSet = 0;
  G4float coreExitX_mm = 0.0f;
  G4float coreExitY_mm = 0.0f;
  G4float coreExitZ_mm = 0.0f;

  G4float endTime_ns = std::numeric_limits<G4float>::quiet_NaN();
  G4float endProperTime_ns = std::numeric_limits<G4float>::quiet_NaN();
  G4float endX_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float endY_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float endZ_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4int endVolumeID = -1;
  G4int endProcessType = -1;
  G4int endProcessSubType = -1;
  G4String endProcessName = "Unknown";
  G4int endTrackStatus = 0;
};

struct TrackVolumeRecord {
  G4int trackID = 0;
  G4int volumeID = -1;
  G4float edep_MeV = 0.0f;
  G4float deltaKE_MeV = 0.0f;
  G4float length_mm = 0.0f;
  G4int nSteps = 0;
};

struct PhotonRecord {
  G4int photonID = 0;
  G4int parentTrackID = 0;
  G4int sourceTrackID = 0;
  G4int creatorProcessType = -1;
  G4int creatorProcessSubType = -1;
  G4int originVolumeID = -1;

  G4float wavelength_nm = 0.0f;
  G4float birthTime_ns = 0.0f;
  G4float emissionAngle_deg = std::numeric_limits<G4float>::quiet_NaN();
  G4float emitterBeta = std::numeric_limits<G4float>::quiet_NaN();

  G4float birthX_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float birthY_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float birthZ_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float birthDirX = std::numeric_limits<G4float>::quiet_NaN();
  G4float birthDirY = std::numeric_limits<G4float>::quiet_NaN();
  G4float birthDirZ = std::numeric_limits<G4float>::quiet_NaN();

  G4int reachedCoreWall = 0;
  G4float firstWallTime_ns = std::numeric_limits<G4float>::quiet_NaN();
  G4float firstWallX_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float firstWallY_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float firstWallZ_mm = std::numeric_limits<G4float>::quiet_NaN();

  G4int extractedToReadout = 0;
  G4float extractionTime_ns = std::numeric_limits<G4float>::quiet_NaN();

  G4int nBoundaryReflections = 0;
  G4int nTIR = 0;
  G4int nTyvekReflections = 0;
  G4int nRayleigh = 0;

  G4int endReason = 0;
  G4int endVolumeID = -1;
  G4float endTime_ns = std::numeric_limits<G4float>::quiet_NaN();
  G4int killedBySafetyCut = 0;
  G4float firstRayleighTrackLength_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float endTrackLength_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4int acrylicAbsOutOfRange = 0;
};

struct SensorIncidentRecord {
  G4int photonID = 0;
  G4int channelID = -1;
  G4float wavelength_nm = 0.0f;
  G4float time_ns = 0.0f;

  G4float localX_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float localY_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float localZ_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float localDirX = std::numeric_limits<G4float>::quiet_NaN();
  G4float localDirY = std::numeric_limits<G4float>::quiet_NaN();
  G4float localDirZ = std::numeric_limits<G4float>::quiet_NaN();
  G4float localPolX = std::numeric_limits<G4float>::quiet_NaN();
  G4float localPolY = std::numeric_limits<G4float>::quiet_NaN();
  G4float localPolZ = std::numeric_limits<G4float>::quiet_NaN();
  G4float incidentAngle_deg = std::numeric_limits<G4float>::quiet_NaN();

  // Entrance-plane coordinates: same local x/y axes as SensorBox, but z=0
  // is the inward-facing physical/virtual readout plane rather than box center.
  G4float entranceX_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float entranceY_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float entranceZ_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4int forward = 0;  // 1 iff signed localDirZ > 0 (toward detector exterior)

  G4int nBoundaryReflections = 0;
  G4int nTIR = 0;
  G4int nTyvekReflections = 0;
  G4int nRayleigh = 0;
};

struct DetectionRecord {
  G4int photonID = 0;
  G4int channelID = -1;
  G4float time_ns = 0.0f;
};

struct ChannelAccumulator {
  G4int nIncident = 0;
  G4int nDetected = 0;
  G4float firstIncidentTime_ns = std::numeric_limits<G4float>::quiet_NaN();
  G4float firstDetectedTime_ns = std::numeric_limits<G4float>::quiet_NaN();
};

struct BoundaryDebugRecord {
  G4int photonID = 0;
  G4int stepNumber = 0;
  G4int preVolumeID = -1;
  G4int postVolumeID = -1;
  G4int boundaryStatus = -1;
  G4String boundaryStatusName = "Unknown";
  G4float time_ns = 0.0f;
  G4float x_mm = 0.0f;
  G4float y_mm = 0.0f;
  G4float z_mm = 0.0f;
  G4float preDirX = std::numeric_limits<G4float>::quiet_NaN();
  G4float preDirY = std::numeric_limits<G4float>::quiet_NaN();
  G4float preDirZ = std::numeric_limits<G4float>::quiet_NaN();
  G4float postDirX = std::numeric_limits<G4float>::quiet_NaN();
  G4float postDirY = std::numeric_limits<G4float>::quiet_NaN();
  G4float postDirZ = std::numeric_limits<G4float>::quiet_NaN();
  G4float prePolX = std::numeric_limits<G4float>::quiet_NaN();
  G4float prePolY = std::numeric_limits<G4float>::quiet_NaN();
  G4float prePolZ = std::numeric_limits<G4float>::quiet_NaN();
  G4float postPolX = std::numeric_limits<G4float>::quiet_NaN();
  G4float postPolY = std::numeric_limits<G4float>::quiet_NaN();
  G4float postPolZ = std::numeric_limits<G4float>::quiet_NaN();
};

struct CoreBoundaryCrossingRecord {
  G4int photonID = 0;
  G4int parentTrackID = 0;
  G4int sourceTrackID = 0;
  G4int rootPrimaryTrackID = 0;
  G4int originVolumeID = -1;
  G4int faceID = -1;

  G4float x_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float y_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float z_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float localU_mm = std::numeric_limits<G4float>::quiet_NaN();
  G4float localV_mm = std::numeric_limits<G4float>::quiet_NaN();

  G4float dirX = std::numeric_limits<G4float>::quiet_NaN();
  G4float dirY = std::numeric_limits<G4float>::quiet_NaN();
  G4float dirZ = std::numeric_limits<G4float>::quiet_NaN();
  G4float cosThetaNormal = std::numeric_limits<G4float>::quiet_NaN();

  G4float wavelength_nm = std::numeric_limits<G4float>::quiet_NaN();
  G4float time_ns = std::numeric_limits<G4float>::quiet_NaN();

  G4float polX = std::numeric_limits<G4float>::quiet_NaN();
  G4float polY = std::numeric_limits<G4float>::quiet_NaN();
  G4float polZ = std::numeric_limits<G4float>::quiet_NaN();
};

class EventAction : public G4UserEventAction {
public:
  explicit EventAction(RunAction* runAction);
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  TrackRecord& RegisterTrack(const G4Track* track, G4int birthVolumeID);
  void UpdateTrackEnd(const G4Track* track, G4int endVolumeID, const G4VProcess* endProcess);
  void AddTrackVolume(G4int trackID, G4int volumeID, G4double edep, G4double deltaKE, G4double length);
  void SetCoreEntry(G4int trackID, const G4ThreeVector& pos);
  void SetCoreExit(G4int trackID, const G4ThreeVector& pos);

  PhotonRecord& RegisterPhoton(const G4Track* track, G4int originVolumeID,
                               G4double wavelengthNm, G4double emissionAngleDeg,
                               G4double emitterBeta,
                               G4double birthGlobalTime,
                               const G4ThreeVector& birthPosition,
                               const G4ThreeVector& birthDirection);
  PhotonRecord* GetPhoton(G4int photonID);
  G4int ResolveRootPrimary(G4int trackID) const;
  G4int ResolveSourceTrack(G4int trackID) const;

  void AddSensorIncident(const SensorIncidentRecord& hit);
  void AddDetection(const DetectionRecord& hit);
  void AddBoundaryDebug(const BoundaryDebugRecord& hit);
  G4bool AddCoreBoundaryCrossing(const CoreBoundaryCrossingRecord& hit);

  void AddScintTop(G4int trackID, G4double edepMeV, G4double timeNs);
  void AddScintBottom(G4int trackID, G4double edepMeV, G4double timeNs);

private:
  RunAction* fRunAction = nullptr;

  std::unordered_map<G4int, TrackRecord> fTracks;
  std::map<std::pair<G4int,G4int>, TrackVolumeRecord> fTrackVolumes;
  std::unordered_map<G4int, PhotonRecord> fPhotons;
  std::unordered_map<G4int, G4int> fRootPrimaryByAnyTrack;
  std::unordered_map<G4int, G4int> fSourceTrackByAnyTrack;

  std::vector<SensorIncidentRecord> fSensorIncidents;
  std::vector<DetectionRecord> fDetections;
  std::vector<BoundaryDebugRecord> fBoundaryDebug;
  std::vector<CoreBoundaryCrossingRecord> fCoreBoundaryCrossings;
  std::unordered_set<G4int> fCoreBoundaryRecorded;
  std::unordered_map<G4int, ChannelAccumulator> fChannels;

  std::vector<ScintDeposit> fScintTopDeposits;
  std::vector<ScintDeposit> fScintBottomDeposits;

  G4double fEdepTopMeV = 0.0;
  G4double fEdepBottomMeV = 0.0;
  G4double fTriggerTopTimeNs = -1.0;
  G4double fTriggerBottomTimeNs = -1.0;
  G4bool fTriggerPassed = false;
};

#endif
