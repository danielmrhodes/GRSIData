#include "TRuntimeObjects.h"
#include "TTigress.h"
#include "TSharc.h"
#include "TTrific.h"

#include "GCutG.h"

void MakeTigress(TRuntimeObjects& obj, TTigress& tigress) {

  std::string dirname = "Tigress";
  for(int i=0;i<tigress.GetMultiplicity();i++) {
    TTigressHit* hit = tigress.GetTigressHit(i);
    
    int det = hit->GetDetector();
    double en = hit->GetEnergy();

    obj.FillHistogram(dirname,"Energy",4000,0,4000,en);
    obj.FillHistogram(dirname,"Energy_v_Det",20,0,20,det,4000,0,4000,en);

  }

  for(int i=0;i<tigress.GetAddbackMultiplicity();i++) {
    TTigressHit* hit = tigress.GetAddbackHit(i);
    
    int det = hit->GetDetector();
    double en = hit->GetEnergy();

    obj.FillHistogram(dirname,"AddbackEnergy",4000,0,4000,en);
    obj.FillHistogram(dirname,"AddbackEnergy_v_Det",20,0,20,det,4000,0,4000,en);

  }

  return;
  
}

void MakeSharc(TRuntimeObjects& obj, TSharc& sharc) {

  std::string dirname = "Sharc";
  for(int i=0;i<sharc.GetMultiplicity();i++) {
    TSharcHit* hit = sharc.GetSharcHit(i);
    
    int det = hit->GetDetector();
    int numF = 100*det + hit->GetFrontStrip();
    int numB = 100*det + hit->GetBackStrip();

    double enF = hit->GetDeltaFrontE();
    double enB = hit->GetDeltaBackE();

    double theta = hit->GetPosition().Theta()*TMath::RadToDeg();

    obj.FillHistogram(dirname,"dE_Front",10000,0,10000,enF);
    obj.FillHistogram(dirname,"dE_Front_v_FrontStrip",1700,0,1700,numF,5000,0,10000,enF);
    obj.FillHistogram(dirname,"dE_Front_v_Theta",180,0,180,theta,5000,0,10000,enF);

    obj.FillHistogram(dirname,"dE_Back",10000,0,10000,enB);
    obj.FillHistogram(dirname,"dE_Back_v_BackStrip",1700,0,1700,numB,5000,0,10000,enB);
    obj.FillHistogram(dirname,"dE_Back_v_Theta",180,0,180,theta,5000,0,10000,enB);

    if(hit->GetPadAddress() == -1)
      continue;

    double padE = hit->GetPad().GetEnergy();
    obj.FillHistogram(dirname,"dE_Front_v_padE",2000,0,16000,padE,2000,0,8000,enF);

  }

  return;
  
}

void MakeTrific(TRuntimeObjects& obj, TTrific& trific) {

  double sumE = 0;
  double sumdE = 0;
  double vetoE = 0;

  std::string dirname = "Trific";
  for(int i=0;i<trific.GetMultiplicity();i++) {
    TTrificHit* hit = trific.GetTrificHit(i);

    int det = hit->GetDetector();
    if(det == 99)
      continue;
    
    double en = hit->GetEnergy();
    if(det > 17)
      vetoE += en;
    if(det > 2 && det < 20)
      sumE += en;
    if(det > 2 && det < 10)
      sumdE += en;

    obj.FillHistogram(dirname,"Energy",4000,0,4000,en);
    obj.FillHistogram(dirname,"Energy_v_Det",25,0,25,det,2000,0,2500,en);

  }

  if(sumE >0 && sumdE > 0 && sumE != sumdE && vetoE < 2000)
    obj.FillHistogram(dirname,"dE_v_E",2000,2000,15000,sumdE,2000,7000,24000,sumE);

  return;
  
}

void MakeSharcTigress(TRuntimeObjects& obj, TSharc& sharc, TTigress& tigress) {

  std::string dirname = "SharcTigress";
  for(int i=0;i<sharc.GetMultiplicity();i++) {
    TSharcHit* shc_hit = sharc.GetSharcHit(i);
    
    long s_time = shc_hit->GetTimeStamp();

    for(int j=0;j<tigress.GetAddbackMultiplicity();j++) {
      TTigressHit* tig_hit = tigress.GetAddbackHit(j);

      long t_time = tig_hit->GetTimeStamp();
      long tdiff = t_time - s_time;

      double en = tig_hit->GetEnergy();

      obj.FillHistogram(dirname,"TDiff",4000,-4000,4000,tdiff);
      obj.FillHistogram(dirname,"TigE_v_TDiff",1000,-2000,2000,tdiff,2000,0,4000,en);

    }

  }

  return;
  
}

void MakeTrificTigress(TRuntimeObjects& obj, TTrific& trific, TTigress& tigress) {

  std::string dirname = "TrificTigress";
  for(int i=0;i<trific.GetMultiplicity();i++) {
    TTrificHit* trif_hit = trific.GetTrificHit(i);
    
    long s_time = trif_hit->GetTimeStamp();

    for(int j=0;j<tigress.GetAddbackMultiplicity();j++) {
      TTigressHit* tig_hit = tigress.GetAddbackHit(j);

      long t_time = tig_hit->GetTimeStamp();
      long tdiff = t_time - s_time;

      double en = tig_hit->GetEnergy();
      
      obj.FillHistogram(dirname,"TDiff",4000,-4000,4000,tdiff);
      obj.FillHistogram(dirname,"TigE_v_TDiff",1000,-2000,2000,tdiff,2000,0,4000,en);

    }

  }

  return;
  
}

void MakeSharcGatedTigress(TRuntimeObjects& obj, TSharc& sharc, TTigress& tigress, GCutG* gate) {

  if(!gate)
    return;

  std::string dirname = Form("SharcGatedTigress_%s",gate->GetName());

  std::vector<TSharcHit*> good_sharc_hits;
  for(int i=0;i<sharc.GetMultiplicity();i++) {
    TSharcHit* hit = sharc.GetSharcHit(i);

    double enF = hit->GetDeltaFrontE();
    double theta = hit->GetPosition().Theta()*TMath::RadToDeg();
    
    if(gate->IsInside(theta,enF)) {
      good_sharc_hits.push_back(hit);
    }

  }
  
  for(int i=0;i<tigress.GetAddbackMultiplicity();i++) {
    TTigressHit* tig_hit = tigress.GetAddbackHit(i);

    long tig_time = tig_hit->GetTimeStamp();
    
    int partner = -1;
    for(unsigned int j=0;j<good_sharc_hits.size();j++) {
      TSharcHit* shc_hit = good_sharc_hits.at(j);

      long shc_time = shc_hit->GetTimeStamp();
      long tdiff = tig_time - shc_time;

      if(tdiff > 40 || tdiff < 90) {
	partner = j;
	break;
      }
    }

    if(partner < 0)
      continue;
    
    good_sharc_hits.erase(good_sharc_hits.begin()+partner);
      
    int det = tig_hit->GetDetector();
    double en = tig_hit->GetDoppler(0.115);

    obj.FillHistogram(dirname,"AddbackDopEn",4000,0,4000,en);
    obj.FillHistogram(dirname,"AddbackDopEn_v_Det",20,0,20,det,4000,0,4000,en);
    
  }

  return;
  
}

void MakeTrificGatedTigress(TRuntimeObjects& obj, TTrific& trific, TTigress& tigress, GCutG* gate) {

  if(!gate)
    return;

  std::string dirname = Form("TrificGatedTigress_%s",gate->GetName());

  std::vector<TTrificHit*> good_trific_hits;
  for(int i=0;i<trific.GetMultiplicity();i++) {
    TTrificHit* hit = trific.GetTrificHit(i);

    int det = hit->GetDetector();
    if(det == 99)
      continue;

    double en = hit->GetEnergy();    
    if(gate->IsInside(det,en)) {
      good_trific_hits.push_back(hit);
    }

  }

  for(int i=0;i<tigress.GetAddbackMultiplicity();i++) {
    TTigressHit* tig_hit = tigress.GetAddbackHit(i);

    long tig_time = tig_hit->GetTimeStamp();
    
    int partner = -1;
    for(unsigned int j=0;j<good_trific_hits.size();j++) {
      TTrificHit* trf_hit = good_trific_hits.at(j);

      long trf_time = trf_hit->GetTimeStamp();
      long tdiff = tig_time - trf_time;

      if(tdiff > 33 || tdiff < 52) {
	partner = j;
	break;
      }
    }

    if(partner < 0)
      continue;
    
    good_trific_hits.erase(good_trific_hits.begin()+partner);
      
    int det = tig_hit->GetDetector();
    double en = tig_hit->GetDoppler(0.115);

    obj.FillHistogram(dirname,"AddbackDopEn",4000,0,4000,en);
    obj.FillHistogram(dirname,"AddbackDopEn_v_Det",20,0,20,det,4000,0,4000,en);
    
  }

  return;
  
}

int gates_loaded = 0;
std::vector<GCutG*> sharc_gates;
std::vector<GCutG*> trific_gates;

extern "C" 
void MakeAnalysisHistograms(TRuntimeObjects& obj) {

  TList *gates = &(obj.GetGates());
  if(gates_loaded!=gates->GetSize()) {
    
    TIter iter(gates);
    while(TObject *obj1 = iter.Next()) {
      
      GCutG *gate = (GCutG*)obj1;
      std::string tag = gate->GetTag();
      std::string name = gate->GetName();
      
      if(!tag.compare("TE")) {
	gates_loaded++;
        std::cout << "Time Gates not implemented: " << name << std::endl;
      }
      else if(!tag.compare("Sharc")) {
        sharc_gates.push_back(gate);
	gates_loaded++;
        std::cout << "Sharc Gate: " << tag << " " << name << std::endl;
      }
      else if(!tag.compare("Trific")) {
        trific_gates.push_back(gate);
	gates_loaded++;
        std::cout << "Trific Gate: " << tag << " " << name << std::endl;
      }
      else {
	std::cout << "Unknown Gate: Name = " << name << ", Tag = " << tag << std::endl;
	gates_loaded++;
      }
      
    }
  }

  std::shared_ptr<TTigress> tigress = obj.GetDetector<TTigress>();
  std::shared_ptr<TSharc> sharc = obj.GetDetector<TSharc>();
  std::shared_ptr<TTrific> trific = obj.GetDetector<TTrific>();
  
  if(tigress) {
    MakeTigress(obj,*tigress);
  }

  if(sharc) {
    MakeSharc(obj,*sharc);
  }

  if(trific) {
    MakeTrific(obj,*trific);
  }

  if(sharc && tigress) {
    MakeSharcTigress(obj,*sharc,*tigress);
    
    for(GCutG* gate : sharc_gates) {
      MakeSharcGatedTigress(obj,*sharc,*tigress,gate);
    }
  }
  
  if(trific && tigress) {
    MakeTrificTigress(obj,*trific,*tigress);

    for(GCutG* gate : trific_gates) {
      MakeTrificGatedTigress(obj,*trific,*tigress,gate);
    }
  }



  return;

}
