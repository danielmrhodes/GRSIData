#include "TRuntimeObjects.h"
#include "TGenericDetector.h"

void MakeGenericDetector(TRuntimeObjects& obj, TGenericDetector& gendet) {

  std::string dirname = "GenericDetector";
  for(int i=0;i<gendet.GetMultiplicity();i++) {

    TGRSIDetectorHit* hit = gendet.GetGenericDetectorHit(i);
    
    int det = hit->GetDetector();
    double en = hit->GetEnergy();
    double tm = hit->GetTime()/1e9;

    obj.FillHistogram(dirname,"Energy",8192,0,2048,en);
    obj.FillHistogram(dirname,Form("Energy_Det%02d",det),8192,0,2048,en);
    obj.FillHistogram(dirname,"Energy_v_Det",20,0,20,det,4000,0,4000,en);

    obj.FillHistogram(dirname,"Energy_v_Time",1800,0,1800,tm,4000,0,4000,en);
    obj.FillHistogram(dirname,Form("Energy_v_Time_Det%02d",det),1800,0,1800,tm,4000,0,4000,en);

  }

  return;
  
}

extern "C" 
void MakeAnalysisHistograms(TRuntimeObjects& obj) {

  std::shared_ptr<TGenericDetector> gendet = obj.GetDetector<TGenericDetector>();
  if(gendet) {
    MakeGenericDetector(obj,*gendet);
  }

  return;

}
