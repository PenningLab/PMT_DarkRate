//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Jun 21 14:10:48 2019 by ROOT version 6.16/00)
//      from the StreamerInfo in file lz_201906192219_000010_000000_raw.root
//////////////////////////////////////////////////////////


#ifndef DetectorPulseMCTruth_h
#define DetectorPulseMCTruth_h
class DetectorPulseMCTruth;

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>

class DetectorPulseMCTruth {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned short iPulseIdentifier;    //1 - S1, 2 - S2, 3 - other
   int            iVertexNumber;       //corresponding vertex index, -1 if no vertex is associated with it
   unsigned int   iPheCount;           //number of phes in the pulse
   double         fFirstPheTime_ns;    //time of the first phe in the pulse
   double         fLastPheTime_ns;     //time of the last phe in the pulse
   vector<unsigned short> iPMTIndices;         //Indices of PMTs with photon hits
   vector<unsigned int>   iPMTHits;            //Number of photons on a pmt

   DetectorPulseMCTruth();
   DetectorPulseMCTruth(const DetectorPulseMCTruth & );
   virtual ~DetectorPulseMCTruth();

};
#endif
