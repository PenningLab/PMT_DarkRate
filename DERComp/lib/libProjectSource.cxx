namespace std {}
using namespace std;
#include "libProjectHeaders.h"

#include "libLinkDef.h"

#include "libProjectDict.cxx"

struct DeleteObjectFunctor {
   template <typename T>
   void operator()(const T *ptr) const {
      delete ptr;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T,Q> &) const {
      // Do nothing
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T,Q*> &ptr) const {
      delete ptr.second;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T*,Q> &ptr) const {
      delete ptr.first;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T*,Q*> &ptr) const {
      delete ptr.first;
      delete ptr.second;
   }
};

#ifndef DetectorMCTruthEvent_cxx
#define DetectorMCTruthEvent_cxx
DetectorMCTruthEvent::DetectorMCTruthEvent() {
}
DetectorMCTruthEvent::DetectorMCTruthEvent(const DetectorMCTruthEvent & rhs)
   : TObject(const_cast<DetectorMCTruthEvent &>( rhs ))
   , pulses(const_cast<DetectorMCTruthEvent &>( rhs ).pulses)
   , vertices(const_cast<DetectorMCTruthEvent &>( rhs ).vertices)
   , darkCounts(const_cast<DetectorMCTruthEvent &>( rhs ).darkCounts)
   , fParentPosition_mm(const_cast<DetectorMCTruthEvent &>( rhs ).fParentPosition_mm)
   , fParentDirection(const_cast<DetectorMCTruthEvent &>( rhs ).fParentDirection)
   , fParentTime_ns(const_cast<DetectorMCTruthEvent &>( rhs ).fParentTime_ns)
   , fParentEnergy_keV(const_cast<DetectorMCTruthEvent &>( rhs ).fParentEnergy_keV)
   , iBaccEventNumber(const_cast<DetectorMCTruthEvent &>( rhs ).iBaccEventNumber)
   , iRunNumber(const_cast<DetectorMCTruthEvent &>( rhs ).iRunNumber)
   , iDEREventNumber(const_cast<DetectorMCTruthEvent &>( rhs ).iDEREventNumber)
   , sParentParticleName(const_cast<DetectorMCTruthEvent &>( rhs ).sParentParticleName)
   , sParentVolumeName(const_cast<DetectorMCTruthEvent &>( rhs ).sParentVolumeName)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   DetectorMCTruthEvent &modrhs = const_cast<DetectorMCTruthEvent &>( rhs );
   modrhs.pulses.clear();
   modrhs.vertices.clear();
   modrhs.darkCounts.clear();
   modrhs.sParentParticleName.clear();
   modrhs.sParentVolumeName.clear();
}
DetectorMCTruthEvent::~DetectorMCTruthEvent() {
}
#endif // DetectorMCTruthEvent_cxx

#ifndef DetectorPulseMCTruth_cxx
#define DetectorPulseMCTruth_cxx
DetectorPulseMCTruth::DetectorPulseMCTruth() {
}
DetectorPulseMCTruth::DetectorPulseMCTruth(const DetectorPulseMCTruth & rhs)
   : iPulseIdentifier(const_cast<DetectorPulseMCTruth &>( rhs ).iPulseIdentifier)
   , iVertexNumber(const_cast<DetectorPulseMCTruth &>( rhs ).iVertexNumber)
   , iPheCount(const_cast<DetectorPulseMCTruth &>( rhs ).iPheCount)
   , fFirstPheTime_ns(const_cast<DetectorPulseMCTruth &>( rhs ).fFirstPheTime_ns)
   , fLastPheTime_ns(const_cast<DetectorPulseMCTruth &>( rhs ).fLastPheTime_ns)
   , iPMTIndices(const_cast<DetectorPulseMCTruth &>( rhs ).iPMTIndices)
   , iPMTHits(const_cast<DetectorPulseMCTruth &>( rhs ).iPMTHits)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   DetectorPulseMCTruth &modrhs = const_cast<DetectorPulseMCTruth &>( rhs );
   modrhs.iPMTIndices.clear();
   modrhs.iPMTHits.clear();
}
DetectorPulseMCTruth::~DetectorPulseMCTruth() {
}
#endif // DetectorPulseMCTruth_cxx

#ifndef DetectorVertexMCTruth_cxx
#define DetectorVertexMCTruth_cxx
DetectorVertexMCTruth::DetectorVertexMCTruth() {
}
DetectorVertexMCTruth::DetectorVertexMCTruth(const DetectorVertexMCTruth & rhs)
   : fPosition_mm(const_cast<DetectorVertexMCTruth &>( rhs ).fPosition_mm)
   , fElectricFieldDirection(const_cast<DetectorVertexMCTruth &>( rhs ).fElectricFieldDirection)
   , fTime_ns(const_cast<DetectorVertexMCTruth &>( rhs ).fTime_ns)
   , fEnergyDep_keV(const_cast<DetectorVertexMCTruth &>( rhs ).fEnergyDep_keV)
   , fElectricFieldStrength_Vcm(const_cast<DetectorVertexMCTruth &>( rhs ).fElectricFieldStrength_Vcm)
   , iRawS1Photons(const_cast<DetectorVertexMCTruth &>( rhs ).iRawS1Photons)
   , iRawS2Photons(const_cast<DetectorVertexMCTruth &>( rhs ).iRawS2Photons)
   , iRawScintPhotons(const_cast<DetectorVertexMCTruth &>( rhs ).iRawScintPhotons)
   , iS1PhotonHits(const_cast<DetectorVertexMCTruth &>( rhs ).iS1PhotonHits)
   , iS2PhotonHits(const_cast<DetectorVertexMCTruth &>( rhs ).iS2PhotonHits)
   , iScintPhotonHits(const_cast<DetectorVertexMCTruth &>( rhs ).iScintPhotonHits)
   , iDetectedS1Photons(const_cast<DetectorVertexMCTruth &>( rhs ).iDetectedS1Photons)
   , iDetectedS2Photons(const_cast<DetectorVertexMCTruth &>( rhs ).iDetectedS2Photons)
   , iDetectedScintPhotons(const_cast<DetectorVertexMCTruth &>( rhs ).iDetectedScintPhotons)
   , iS1PulseIndex(const_cast<DetectorVertexMCTruth &>( rhs ).iS1PulseIndex)
   , iS2PulseIndex(const_cast<DetectorVertexMCTruth &>( rhs ).iS2PulseIndex)
   , iScintillationPulseIndex(const_cast<DetectorVertexMCTruth &>( rhs ).iScintillationPulseIndex)
   , iArtifactPulseIndices(const_cast<DetectorVertexMCTruth &>( rhs ).iArtifactPulseIndices)
   , sParticleName(const_cast<DetectorVertexMCTruth &>( rhs ).sParticleName)
   , sVolumeName(const_cast<DetectorVertexMCTruth &>( rhs ).sVolumeName)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   DetectorVertexMCTruth &modrhs = const_cast<DetectorVertexMCTruth &>( rhs );
   modrhs.iArtifactPulseIndices.clear();
   modrhs.sParticleName.clear();
   modrhs.sVolumeName.clear();
}
DetectorVertexMCTruth::~DetectorVertexMCTruth() {
}
#endif // DetectorVertexMCTruth_cxx

#ifndef DarkCountMCTruth_cxx
#define DarkCountMCTruth_cxx
DarkCountMCTruth::DarkCountMCTruth() {
}
DarkCountMCTruth::DarkCountMCTruth(const DarkCountMCTruth & rhs)
   : iPMTIndex(const_cast<DarkCountMCTruth &>( rhs ).iPMTIndex)
   , fTime_ns(const_cast<DarkCountMCTruth &>( rhs ).fTime_ns)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
}
DarkCountMCTruth::~DarkCountMCTruth() {
}
#endif // DarkCountMCTruth_cxx

