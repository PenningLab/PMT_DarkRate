//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Jun 21 14:10:48 2019 by ROOT version 6.16/00)
//      from the StreamerInfo in file lz_201906192219_000010_000000_raw.root
//////////////////////////////////////////////////////////


#ifndef DetectorMCTruthEvent_h
#define DetectorMCTruthEvent_h
class DetectorMCTruthEvent;

#include "Rtypes.h"
#include "TObject.h"
#include "Riostream.h"
#include <vector>
#include "DetectorPulseMCTruth.h"
#include "DetectorVertexMCTruth.h"
#include "DarkCountMCTruth.h"
#include "TVector3.h"
#include <string>

class DetectorMCTruthEvent : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<DetectorPulseMCTruth> pulses;      // (DetectorPulseMCTruth)
   vector<DetectorVertexMCTruth> vertices;    // (DetectorVertexMCTruth)
   vector<DarkCountMCTruth>      darkCounts;    // (DarkCountMCTruth)
   TVector3                      fParentPosition_mm;    //
   TVector3                      fParentDirection;      //
   double                        fParentTime_ns;        //
   float                         fParentEnergy_keV;     //
   unsigned int                  iBaccEventNumber;      //
   unsigned int                  iRunNumber;            //
   unsigned short                iDEREventNumber;       //
   string                        sParentParticleName;    //
   string                        sParentVolumeName;      //

   DetectorMCTruthEvent();
   DetectorMCTruthEvent(const DetectorMCTruthEvent & );
   virtual ~DetectorMCTruthEvent();

   ClassDef(DetectorMCTruthEvent,7); // Generated by MakeProject.
};
#endif