#include <iostream>
#include <numeric>	
#include "TEmma.h"
#include "TRandom.h"
#include "TMath.h"

////////////////////////////////////////////////////////////
//
// TEmma
//
// The TEmma class defines the observables and algorithms used
// when analyzing EMMA data.
//
////////////////////////////////////////////////////////////

/// \cond CLASSIMP
ClassImp(TEmma)
/// \endcond

double TEmma::fAnodeTrigger = 0;
double TEmma::fICEnergy = 0;
double TEmma::fXdiff = 0;
double TEmma::fXsum = 0;
double TEmma::fYdiff = 0;
double TEmma::fYsum = 0;
double TEmma::fXlength = 80.; //Size of X focal plane in mm
double TEmma::fYlength = 30.; //Size of Y focal plane in mm
short TEmma::fFail = 0;

TEmma::TEmma() : TDetector()
{
#if MAJOR_ROOT_VERSION < 6
	Class()->IgnoreTObjectStreamer(kTRUE);
#endif
	Clear();
}

TEmma::TEmma(const TEmma& rhs) : TDetector()
{
#if MAJOR_ROOT_VERSION < 6
	Class()->IgnoreTObjectStreamer(kTRUE);
#endif
	rhs.Copy(*this);
}

void TEmma::Copy(TObject& rhs) const
{
	TDetector::Copy(rhs);
	static_cast<TEmma&>(rhs).fEmmaICHits   = fEmmaICHits;
	static_cast<TEmma&>(rhs).fEmmaSiHits   = fEmmaSiHits;
	static_cast<TEmma&>(rhs).fEmmaSSBHits   = fEmmaSSBHits;
	static_cast<TEmma&>(rhs).fEmmaTdcHits   = fEmmaTdcHits;
	static_cast<TEmma&>(rhs).fEmmaAnodeHits   = fEmmaAnodeHits;
	static_cast<TEmma&>(rhs).fEmmaTriggerHits   = fEmmaTriggerHits;
}

TEmma::~TEmma()
{
	// Default Destructor
}

void TEmma::Print(Option_t*) const
{
	// Prints out TEmma members, currently does nothing.
	std::cout<<fHits.size()<<" fHits"<<std::endl;
}

TEmma& TEmma::operator=(const TEmma& rhs)
{
	rhs.Copy(*this);
	return *this;
}

void TEmma::AddFragment(const std::shared_ptr<const TFragment>& frag, TChannel* chan)
{
	if(frag == nullptr || chan == nullptr) {
		return;
	}
	TEmmaHit dethit(*frag);
	if(chan->GetMnemonic()->SubSystem() == TMnemonic::EMnemonic::kI) { // IC ADC Data 
		fEmmaICHits.push_back(std::move(dethit));
	} else if(chan->GetMnemonic()->SubSystem() == TMnemonic::EMnemonic::kS) { // Si at focal plane
		fEmmaSiHits.push_back(std::move(dethit));
	} else if(chan->GetMnemonic()->SubSystem() == TMnemonic::EMnemonic::kT) { // EMMA focal plane Trigger
		fEmmaTriggerHits.push_back(std::move(dethit));
	} else if(chan->GetMnemonic()->SubSystem() == TMnemonic::EMnemonic::kP) { // PGAC
		if(chan->GetMnemonic()->CollectedCharge() == TMnemonic::EMnemonic::kN) { // Anode data 
			switch (chan->GetMnemonic()->OutputSensor()) {
				case TMnemonic::EMnemonic::kX:
					fEmmaAnodeHits.push_back(std::move(dethit)); // Anode ADC data  
					break;

				case TMnemonic::EMnemonic::kT:
					dethit.SetTdcNumber(frag->GetSegment());
					fEmmaTdcHits.push_back(std::move(dethit)); // PGAC Anode TDC data
					break;

				default:
					break;
			}
		} else if(chan->GetMnemonic()->CollectedCharge() == TMnemonic::EMnemonic::kP) { // PGAC Cathode TDC data
			switch (chan->GetMnemonic()->OutputSensor()) {
				case TMnemonic::EMnemonic::kL:
					dethit.SetTdcNumber(10);
					break;

				case TMnemonic::EMnemonic::kR:
					dethit.SetTdcNumber(11);
					break;

				case TMnemonic::EMnemonic::kU:
					dethit.SetTdcNumber(12);
					break;

				case TMnemonic::EMnemonic::kD:
					dethit.SetTdcNumber(13);
					break;

				default:
					break;
			};
			fEmmaTdcHits.push_back(std::move(dethit));
		}
	} else if(chan->GetMnemonic()->SubSystem() == TMnemonic::EMnemonic::kO) { // ORTEC SSBs at target position
		fEmmaSSBHits.push_back(std::move(dethit));
	} else return;
}

TVector3 TEmma::GetPosition(double left, double right, double top, double bottom, double delayL, double delayR, double delayT, double delayB )
{
	// Calculates recoil position from PGAC TDC data including delays
	double Xdiff = (left+delayL) - (right+delayR);
	double Xsum = (left) + (right);
	double Ydiff = (bottom+delayB) - (top+delayT);
	double Ysum = (bottom) + (top);

	double Xpos = ( Xdiff / Xsum )*fXlength;
	double Ypos = ( Ydiff / Ysum )*fYlength;

	return TVector3(Xpos, Ypos, 1);
}

TEmmaHit* TEmma::GetICHit(const int& i)
{
	if(i < GetICMultiplicity()) {
		return &fEmmaICHits.at(i);
	}
	std::cerr<<"EMMA IC hits are out of range"<<std::endl;
	throw grsi::exit_exception(1);
	return nullptr;
}

TEmmaHit* TEmma::GetAnodeHit(const int& i)
{
	if(i < GetAnodeMultiplicity()) {
		return &fEmmaAnodeHits.at(i);
	}
	std::cerr<<"EMMA PGAC Anode hits are out of range"<<std::endl;
	throw grsi::exit_exception(1);
	return nullptr;
}

TEmmaHit* TEmma::GetTdcHit(const int& i)
{
	if(i < GetTdcMultiplicity()) {
		return &fEmmaTdcHits.at(i);
	}
	std::cerr<<"EMMA TDC hits are out of range"<<std::endl;
	throw grsi::exit_exception(1);
	return nullptr;
}

TEmmaHit* TEmma::GetSiHit(const int& i)
{
	if(i < GetSiMultiplicity()) {
		return &fEmmaSiHits.at(i);
	}
	std::cerr<<"EMMA Si hits are out of range"<<std::endl;
	throw grsi::exit_exception(1);
	return nullptr;
}

TEmmaHit* TEmma::GetTriggerHit(const int& i)
{
	if(i < GetTriggerMultiplicity()) {
		return &fEmmaTriggerHits.at(i);
	}
	std::cerr<<"EMMA Trigger hits are out of range"<<std::endl;
	throw grsi::exit_exception(1);
	return nullptr;
}

TEmmaHit* TEmma::GetSSBHit(const int& i)
{
	if(i < GetSSBMultiplicity()) {
		return &fEmmaSSBHits.at(i);
	}
	std::cerr<<"EMMA SSB hits are out of range"<<std::endl;
	throw grsi::exit_exception(1);
	return nullptr;
}

void TEmma::BuildHits() 
{ 
	// Everything below is reproduced from the EMMA sort code given to me by Nick Esker
	// The EMMA PGAC has all cathode signals chained together in the X and Y direction with readouts at each end fed into a TDC.  
	// Build hits subtracts the trigger time (from an anode wire) and returns a left/right/up/down value which is used in GetPosition()

	std::vector <double> tdcArray;
	std::vector <double> icArray;
	if(fEmmaTdcHits.size() > 4 ) {   // Require a Good hit to contain only the PGAC TDC signals
		TEmmaHit * hit = new TEmmaHit();
		for(size_t i = 0; i < fEmmaTdcHits.size(); ++i) { 
			hit->SetTimeStamp(fEmmaTdcHits[i].GetTimeStamp());
			hit->SetAddress(fEmmaTdcHits[i].GetAddress());
			if(fEmmaTdcHits[i].GetTdcNumber() < 10) tdcArray.push_back(fEmmaTdcHits[i].GetEnergy());
			if(fEmmaTdcHits[i].GetTdcNumber() == 10) hit->SetLeft(fEmmaTdcHits[i].GetEnergy());
			if(fEmmaTdcHits[i].GetTdcNumber() == 11) hit->SetRight(fEmmaTdcHits[i].GetEnergy());
			if(fEmmaTdcHits[i].GetTdcNumber() == 12) hit->SetTop(fEmmaTdcHits[i].GetEnergy());
			if(fEmmaTdcHits[i].GetTdcNumber() == 13) hit->SetBottom(fEmmaTdcHits[i].GetEnergy());
		}

		if(tdcArray.size() != 0) {
			fAnodeTrigger = * std::min_element(tdcArray.begin(), tdcArray.end());
			if(hit->GetLeft() != 0 && hit->GetRight() != 0 && hit->GetTop() != 0 && hit->GetBottom() != 0 && fAnodeTrigger != 0) {
				hit->SetLeft((hit->GetLeft() - fAnodeTrigger));
				hit->SetRight((hit->GetRight() - fAnodeTrigger));
				hit->SetTop((hit->GetTop() - fAnodeTrigger));
				hit->SetBottom((hit->GetBottom() - fAnodeTrigger));
				hit->SetAnodeTrigger(fAnodeTrigger);
				fHits.push_back(hit);
			} else {
				//std::cout<<"TDC Array Failed"<<std::endl;
				fFail = 0;
				if(hit->GetLeft() == 0) fFail++;
				if(hit->GetRight() == 0) fFail++;
				if(hit->GetTop() == 0) fFail++;
				if(hit->GetBottom() == 0) fFail++;
				hit->SetFailedFill(fFail);
				fHits.push_back(hit);
			}
		} else {
			return;
		}
	}
	return;
}
void TEmma::Clear(Option_t* opt)
{
	TDetector::Clear(opt);
	fHits.clear();
	fEmmaICHits.clear();
	fEmmaAnodeHits.clear();
	fEmmaTdcHits.clear();
	fEmmaSiHits.clear();
	fEmmaTriggerHits.clear();
	fEmmaSSBHits.clear();
}
