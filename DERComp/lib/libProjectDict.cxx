// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME libProjectDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "libProjectHeaders.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *DetectorPulseMCTruth_Dictionary();
   static void DetectorPulseMCTruth_TClassManip(TClass*);
   static void *new_DetectorPulseMCTruth(void *p = 0);
   static void *newArray_DetectorPulseMCTruth(Long_t size, void *p);
   static void delete_DetectorPulseMCTruth(void *p);
   static void deleteArray_DetectorPulseMCTruth(void *p);
   static void destruct_DetectorPulseMCTruth(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DetectorPulseMCTruth*)
   {
      ::DetectorPulseMCTruth *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DetectorPulseMCTruth));
      static ::ROOT::TGenericClassInfo 
         instance("DetectorPulseMCTruth", "DetectorPulseMCTruth.h", 16,
                  typeid(::DetectorPulseMCTruth), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DetectorPulseMCTruth_Dictionary, isa_proxy, 4,
                  sizeof(::DetectorPulseMCTruth) );
      instance.SetNew(&new_DetectorPulseMCTruth);
      instance.SetNewArray(&newArray_DetectorPulseMCTruth);
      instance.SetDelete(&delete_DetectorPulseMCTruth);
      instance.SetDeleteArray(&deleteArray_DetectorPulseMCTruth);
      instance.SetDestructor(&destruct_DetectorPulseMCTruth);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DetectorPulseMCTruth*)
   {
      return GenerateInitInstanceLocal((::DetectorPulseMCTruth*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DetectorPulseMCTruth*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DetectorPulseMCTruth_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::DetectorPulseMCTruth*)0x0)->GetClass();
      DetectorPulseMCTruth_TClassManip(theClass);
   return theClass;
   }

   static void DetectorPulseMCTruth_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *DetectorVertexMCTruth_Dictionary();
   static void DetectorVertexMCTruth_TClassManip(TClass*);
   static void *new_DetectorVertexMCTruth(void *p = 0);
   static void *newArray_DetectorVertexMCTruth(Long_t size, void *p);
   static void delete_DetectorVertexMCTruth(void *p);
   static void deleteArray_DetectorVertexMCTruth(void *p);
   static void destruct_DetectorVertexMCTruth(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DetectorVertexMCTruth*)
   {
      ::DetectorVertexMCTruth *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DetectorVertexMCTruth));
      static ::ROOT::TGenericClassInfo 
         instance("DetectorVertexMCTruth", "DetectorVertexMCTruth.h", 18,
                  typeid(::DetectorVertexMCTruth), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DetectorVertexMCTruth_Dictionary, isa_proxy, 4,
                  sizeof(::DetectorVertexMCTruth) );
      instance.SetNew(&new_DetectorVertexMCTruth);
      instance.SetNewArray(&newArray_DetectorVertexMCTruth);
      instance.SetDelete(&delete_DetectorVertexMCTruth);
      instance.SetDeleteArray(&deleteArray_DetectorVertexMCTruth);
      instance.SetDestructor(&destruct_DetectorVertexMCTruth);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DetectorVertexMCTruth*)
   {
      return GenerateInitInstanceLocal((::DetectorVertexMCTruth*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DetectorVertexMCTruth*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DetectorVertexMCTruth_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::DetectorVertexMCTruth*)0x0)->GetClass();
      DetectorVertexMCTruth_TClassManip(theClass);
   return theClass;
   }

   static void DetectorVertexMCTruth_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *DarkCountMCTruth_Dictionary();
   static void DarkCountMCTruth_TClassManip(TClass*);
   static void *new_DarkCountMCTruth(void *p = 0);
   static void *newArray_DarkCountMCTruth(Long_t size, void *p);
   static void delete_DarkCountMCTruth(void *p);
   static void deleteArray_DarkCountMCTruth(void *p);
   static void destruct_DarkCountMCTruth(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DarkCountMCTruth*)
   {
      ::DarkCountMCTruth *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DarkCountMCTruth));
      static ::ROOT::TGenericClassInfo 
         instance("DarkCountMCTruth", "DarkCountMCTruth.h", 14,
                  typeid(::DarkCountMCTruth), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DarkCountMCTruth_Dictionary, isa_proxy, 4,
                  sizeof(::DarkCountMCTruth) );
      instance.SetNew(&new_DarkCountMCTruth);
      instance.SetNewArray(&newArray_DarkCountMCTruth);
      instance.SetDelete(&delete_DarkCountMCTruth);
      instance.SetDeleteArray(&deleteArray_DarkCountMCTruth);
      instance.SetDestructor(&destruct_DarkCountMCTruth);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DarkCountMCTruth*)
   {
      return GenerateInitInstanceLocal((::DarkCountMCTruth*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DarkCountMCTruth*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DarkCountMCTruth_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::DarkCountMCTruth*)0x0)->GetClass();
      DarkCountMCTruth_TClassManip(theClass);
   return theClass;
   }

   static void DarkCountMCTruth_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_DetectorMCTruthEvent(void *p = 0);
   static void *newArray_DetectorMCTruthEvent(Long_t size, void *p);
   static void delete_DetectorMCTruthEvent(void *p);
   static void deleteArray_DetectorMCTruthEvent(void *p);
   static void destruct_DetectorMCTruthEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DetectorMCTruthEvent*)
   {
      ::DetectorMCTruthEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DetectorMCTruthEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DetectorMCTruthEvent", ::DetectorMCTruthEvent::Class_Version(), "DetectorMCTruthEvent.h", 22,
                  typeid(::DetectorMCTruthEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::DetectorMCTruthEvent::Dictionary, isa_proxy, 4,
                  sizeof(::DetectorMCTruthEvent) );
      instance.SetNew(&new_DetectorMCTruthEvent);
      instance.SetNewArray(&newArray_DetectorMCTruthEvent);
      instance.SetDelete(&delete_DetectorMCTruthEvent);
      instance.SetDeleteArray(&deleteArray_DetectorMCTruthEvent);
      instance.SetDestructor(&destruct_DetectorMCTruthEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DetectorMCTruthEvent*)
   {
      return GenerateInitInstanceLocal((::DetectorMCTruthEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DetectorMCTruthEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr DetectorMCTruthEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DetectorMCTruthEvent::Class_Name()
{
   return "DetectorMCTruthEvent";
}

//______________________________________________________________________________
const char *DetectorMCTruthEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DetectorMCTruthEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DetectorMCTruthEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DetectorMCTruthEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DetectorMCTruthEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DetectorMCTruthEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DetectorMCTruthEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DetectorMCTruthEvent*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DetectorPulseMCTruth(void *p) {
      return  p ? new(p) ::DetectorPulseMCTruth : new ::DetectorPulseMCTruth;
   }
   static void *newArray_DetectorPulseMCTruth(Long_t nElements, void *p) {
      return p ? new(p) ::DetectorPulseMCTruth[nElements] : new ::DetectorPulseMCTruth[nElements];
   }
   // Wrapper around operator delete
   static void delete_DetectorPulseMCTruth(void *p) {
      delete ((::DetectorPulseMCTruth*)p);
   }
   static void deleteArray_DetectorPulseMCTruth(void *p) {
      delete [] ((::DetectorPulseMCTruth*)p);
   }
   static void destruct_DetectorPulseMCTruth(void *p) {
      typedef ::DetectorPulseMCTruth current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DetectorPulseMCTruth

namespace ROOT {
   // Wrappers around operator new
   static void *new_DetectorVertexMCTruth(void *p) {
      return  p ? new(p) ::DetectorVertexMCTruth : new ::DetectorVertexMCTruth;
   }
   static void *newArray_DetectorVertexMCTruth(Long_t nElements, void *p) {
      return p ? new(p) ::DetectorVertexMCTruth[nElements] : new ::DetectorVertexMCTruth[nElements];
   }
   // Wrapper around operator delete
   static void delete_DetectorVertexMCTruth(void *p) {
      delete ((::DetectorVertexMCTruth*)p);
   }
   static void deleteArray_DetectorVertexMCTruth(void *p) {
      delete [] ((::DetectorVertexMCTruth*)p);
   }
   static void destruct_DetectorVertexMCTruth(void *p) {
      typedef ::DetectorVertexMCTruth current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DetectorVertexMCTruth

namespace ROOT {
   // Wrappers around operator new
   static void *new_DarkCountMCTruth(void *p) {
      return  p ? new(p) ::DarkCountMCTruth : new ::DarkCountMCTruth;
   }
   static void *newArray_DarkCountMCTruth(Long_t nElements, void *p) {
      return p ? new(p) ::DarkCountMCTruth[nElements] : new ::DarkCountMCTruth[nElements];
   }
   // Wrapper around operator delete
   static void delete_DarkCountMCTruth(void *p) {
      delete ((::DarkCountMCTruth*)p);
   }
   static void deleteArray_DarkCountMCTruth(void *p) {
      delete [] ((::DarkCountMCTruth*)p);
   }
   static void destruct_DarkCountMCTruth(void *p) {
      typedef ::DarkCountMCTruth current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DarkCountMCTruth

//______________________________________________________________________________
void DetectorMCTruthEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class DetectorMCTruthEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(DetectorMCTruthEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(DetectorMCTruthEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DetectorMCTruthEvent(void *p) {
      return  p ? new(p) ::DetectorMCTruthEvent : new ::DetectorMCTruthEvent;
   }
   static void *newArray_DetectorMCTruthEvent(Long_t nElements, void *p) {
      return p ? new(p) ::DetectorMCTruthEvent[nElements] : new ::DetectorMCTruthEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_DetectorMCTruthEvent(void *p) {
      delete ((::DetectorMCTruthEvent*)p);
   }
   static void deleteArray_DetectorMCTruthEvent(void *p) {
      delete [] ((::DetectorMCTruthEvent*)p);
   }
   static void destruct_DetectorMCTruthEvent(void *p) {
      typedef ::DetectorMCTruthEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DetectorMCTruthEvent

namespace ROOT {
   static TClass *vectorlEunsignedsPshortgR_Dictionary();
   static void vectorlEunsignedsPshortgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPshortgR(void *p = 0);
   static void *newArray_vectorlEunsignedsPshortgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPshortgR(void *p);
   static void deleteArray_vectorlEunsignedsPshortgR(void *p);
   static void destruct_vectorlEunsignedsPshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned short>*)
   {
      vector<unsigned short> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned short>", -2, "vector", 339,
                  typeid(vector<unsigned short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPshortgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned short>) );
      instance.SetNew(&new_vectorlEunsignedsPshortgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPshortgR);
      instance.SetDelete(&delete_vectorlEunsignedsPshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPshortgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned short> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned short>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned short>*)0x0)->GetClass();
      vectorlEunsignedsPshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPshortgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short> : new vector<unsigned short>;
   }
   static void *newArray_vectorlEunsignedsPshortgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short>[nElements] : new vector<unsigned short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPshortgR(void *p) {
      delete ((vector<unsigned short>*)p);
   }
   static void deleteArray_vectorlEunsignedsPshortgR(void *p) {
      delete [] ((vector<unsigned short>*)p);
   }
   static void destruct_vectorlEunsignedsPshortgR(void *p) {
      typedef vector<unsigned short> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned short>

namespace ROOT {
   static TClass *vectorlEunsignedsPintgR_Dictionary();
   static void vectorlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPintgR(void *p = 0);
   static void *newArray_vectorlEunsignedsPintgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPintgR(void *p);
   static void deleteArray_vectorlEunsignedsPintgR(void *p);
   static void destruct_vectorlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned int>*)
   {
      vector<unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned int>", -2, "vector", 339,
                  typeid(vector<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned int>) );
      instance.SetNew(&new_vectorlEunsignedsPintgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPintgR);
      instance.SetDelete(&delete_vectorlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPintgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned int>*)0x0)->GetClass();
      vectorlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int> : new vector<unsigned int>;
   }
   static void *newArray_vectorlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int>[nElements] : new vector<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPintgR(void *p) {
      delete ((vector<unsigned int>*)p);
   }
   static void deleteArray_vectorlEunsignedsPintgR(void *p) {
      delete [] ((vector<unsigned int>*)p);
   }
   static void destruct_vectorlEunsignedsPintgR(void *p) {
      typedef vector<unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned int>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 339,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEDetectorVertexMCTruthgR_Dictionary();
   static void vectorlEDetectorVertexMCTruthgR_TClassManip(TClass*);
   static void *new_vectorlEDetectorVertexMCTruthgR(void *p = 0);
   static void *newArray_vectorlEDetectorVertexMCTruthgR(Long_t size, void *p);
   static void delete_vectorlEDetectorVertexMCTruthgR(void *p);
   static void deleteArray_vectorlEDetectorVertexMCTruthgR(void *p);
   static void destruct_vectorlEDetectorVertexMCTruthgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<DetectorVertexMCTruth>*)
   {
      vector<DetectorVertexMCTruth> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<DetectorVertexMCTruth>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<DetectorVertexMCTruth>", -2, "vector", 339,
                  typeid(vector<DetectorVertexMCTruth>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEDetectorVertexMCTruthgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<DetectorVertexMCTruth>) );
      instance.SetNew(&new_vectorlEDetectorVertexMCTruthgR);
      instance.SetNewArray(&newArray_vectorlEDetectorVertexMCTruthgR);
      instance.SetDelete(&delete_vectorlEDetectorVertexMCTruthgR);
      instance.SetDeleteArray(&deleteArray_vectorlEDetectorVertexMCTruthgR);
      instance.SetDestructor(&destruct_vectorlEDetectorVertexMCTruthgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<DetectorVertexMCTruth> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<DetectorVertexMCTruth>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEDetectorVertexMCTruthgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<DetectorVertexMCTruth>*)0x0)->GetClass();
      vectorlEDetectorVertexMCTruthgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEDetectorVertexMCTruthgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEDetectorVertexMCTruthgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DetectorVertexMCTruth> : new vector<DetectorVertexMCTruth>;
   }
   static void *newArray_vectorlEDetectorVertexMCTruthgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DetectorVertexMCTruth>[nElements] : new vector<DetectorVertexMCTruth>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEDetectorVertexMCTruthgR(void *p) {
      delete ((vector<DetectorVertexMCTruth>*)p);
   }
   static void deleteArray_vectorlEDetectorVertexMCTruthgR(void *p) {
      delete [] ((vector<DetectorVertexMCTruth>*)p);
   }
   static void destruct_vectorlEDetectorVertexMCTruthgR(void *p) {
      typedef vector<DetectorVertexMCTruth> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<DetectorVertexMCTruth>

namespace ROOT {
   static TClass *vectorlEDetectorPulseMCTruthgR_Dictionary();
   static void vectorlEDetectorPulseMCTruthgR_TClassManip(TClass*);
   static void *new_vectorlEDetectorPulseMCTruthgR(void *p = 0);
   static void *newArray_vectorlEDetectorPulseMCTruthgR(Long_t size, void *p);
   static void delete_vectorlEDetectorPulseMCTruthgR(void *p);
   static void deleteArray_vectorlEDetectorPulseMCTruthgR(void *p);
   static void destruct_vectorlEDetectorPulseMCTruthgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<DetectorPulseMCTruth>*)
   {
      vector<DetectorPulseMCTruth> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<DetectorPulseMCTruth>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<DetectorPulseMCTruth>", -2, "vector", 339,
                  typeid(vector<DetectorPulseMCTruth>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEDetectorPulseMCTruthgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<DetectorPulseMCTruth>) );
      instance.SetNew(&new_vectorlEDetectorPulseMCTruthgR);
      instance.SetNewArray(&newArray_vectorlEDetectorPulseMCTruthgR);
      instance.SetDelete(&delete_vectorlEDetectorPulseMCTruthgR);
      instance.SetDeleteArray(&deleteArray_vectorlEDetectorPulseMCTruthgR);
      instance.SetDestructor(&destruct_vectorlEDetectorPulseMCTruthgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<DetectorPulseMCTruth> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<DetectorPulseMCTruth>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEDetectorPulseMCTruthgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<DetectorPulseMCTruth>*)0x0)->GetClass();
      vectorlEDetectorPulseMCTruthgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEDetectorPulseMCTruthgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEDetectorPulseMCTruthgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DetectorPulseMCTruth> : new vector<DetectorPulseMCTruth>;
   }
   static void *newArray_vectorlEDetectorPulseMCTruthgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DetectorPulseMCTruth>[nElements] : new vector<DetectorPulseMCTruth>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEDetectorPulseMCTruthgR(void *p) {
      delete ((vector<DetectorPulseMCTruth>*)p);
   }
   static void deleteArray_vectorlEDetectorPulseMCTruthgR(void *p) {
      delete [] ((vector<DetectorPulseMCTruth>*)p);
   }
   static void destruct_vectorlEDetectorPulseMCTruthgR(void *p) {
      typedef vector<DetectorPulseMCTruth> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<DetectorPulseMCTruth>

namespace ROOT {
   static TClass *vectorlEDarkCountMCTruthgR_Dictionary();
   static void vectorlEDarkCountMCTruthgR_TClassManip(TClass*);
   static void *new_vectorlEDarkCountMCTruthgR(void *p = 0);
   static void *newArray_vectorlEDarkCountMCTruthgR(Long_t size, void *p);
   static void delete_vectorlEDarkCountMCTruthgR(void *p);
   static void deleteArray_vectorlEDarkCountMCTruthgR(void *p);
   static void destruct_vectorlEDarkCountMCTruthgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<DarkCountMCTruth>*)
   {
      vector<DarkCountMCTruth> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<DarkCountMCTruth>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<DarkCountMCTruth>", -2, "vector", 339,
                  typeid(vector<DarkCountMCTruth>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEDarkCountMCTruthgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<DarkCountMCTruth>) );
      instance.SetNew(&new_vectorlEDarkCountMCTruthgR);
      instance.SetNewArray(&newArray_vectorlEDarkCountMCTruthgR);
      instance.SetDelete(&delete_vectorlEDarkCountMCTruthgR);
      instance.SetDeleteArray(&deleteArray_vectorlEDarkCountMCTruthgR);
      instance.SetDestructor(&destruct_vectorlEDarkCountMCTruthgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<DarkCountMCTruth> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<DarkCountMCTruth>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEDarkCountMCTruthgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<DarkCountMCTruth>*)0x0)->GetClass();
      vectorlEDarkCountMCTruthgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEDarkCountMCTruthgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEDarkCountMCTruthgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DarkCountMCTruth> : new vector<DarkCountMCTruth>;
   }
   static void *newArray_vectorlEDarkCountMCTruthgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DarkCountMCTruth>[nElements] : new vector<DarkCountMCTruth>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEDarkCountMCTruthgR(void *p) {
      delete ((vector<DarkCountMCTruth>*)p);
   }
   static void deleteArray_vectorlEDarkCountMCTruthgR(void *p) {
      delete [] ((vector<DarkCountMCTruth>*)p);
   }
   static void destruct_vectorlEDarkCountMCTruthgR(void *p) {
      typedef vector<DarkCountMCTruth> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<DarkCountMCTruth>

namespace {
  void TriggerDictionaryInitialization_libProjectDict_Impl() {
    static const char* headers[] = {
"libProjectHeaders.h",
0
    };
    static const char* includePaths[] = {
"/yorktown/root6/include",
"/yorktown/root6/etc",
"/yorktown/root6/etc/cling",
"/yorktown/root6/include",
"/yorktown/root6/include",
"/yorktown/DERPulses/lib/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libProjectDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$DetectorPulseMCTruth.h")))  __attribute__((annotate("$clingAutoload$libProjectHeaders.h")))  DetectorPulseMCTruth;
class __attribute__((annotate("$clingAutoload$DetectorVertexMCTruth.h")))  __attribute__((annotate("$clingAutoload$libProjectHeaders.h")))  DetectorVertexMCTruth;
class __attribute__((annotate("$clingAutoload$DarkCountMCTruth.h")))  __attribute__((annotate("$clingAutoload$libProjectHeaders.h")))  DarkCountMCTruth;
class __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$DetectorMCTruthEvent.h")))  __attribute__((annotate("$clingAutoload$libProjectHeaders.h")))  DetectorMCTruthEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libProjectDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "libProjectHeaders.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"DarkCountMCTruth", payloadCode, "@",
"DetectorMCTruthEvent", payloadCode, "@",
"DetectorPulseMCTruth", payloadCode, "@",
"DetectorVertexMCTruth", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libProjectDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libProjectDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libProjectDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libProjectDict() {
  TriggerDictionaryInitialization_libProjectDict_Impl();
}
