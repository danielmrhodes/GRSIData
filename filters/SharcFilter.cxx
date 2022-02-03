// To make a new filter, copy this file under a new name in the "filter" directory.
// The "FilterCondition" function should return a boolean value.
// The boolean indicates whether the event should be kept or not.

#include "TRuntimeObjects.h"
#include "TSharc.h"

extern "C"
bool FilterCondition(TRuntimeObjects& obj) {
 
  std::shared_ptr<TSharc> sharc = obj.GetDetector<TSharc>();
  if(!sharc)
    return false;
  
  if(sharc->GetMultiplicity() == 0)
    return false;

  return true;
}
