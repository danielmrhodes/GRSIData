#ifndef TGRSIMNEMONIC_H
#define TGRSIMNEMONIC_H

#include <string>
#include "TMnemonic.h"
#include "Globals.h"
#include "TClass.h"

enum class EDigitizer : char { kDefault, kGRF16, kGRF4G, kTIG10, kTIG64, kCaen, kMadc, kV1190 };

class TGRSIMnemonic : public TMnemonic {
public:
   TGRSIMnemonic() : TMnemonic() { Clear(); }
   TGRSIMnemonic(const char* name) : TGRSIMnemonic() { TMnemonic::Parse(name); }
   ~TGRSIMnemonic() override = default;

   // standard C++ makes these enumerations global to the class. ie, the name of the enumeration
   // EMnemonic or ESystem has no effect on the clashing of enumerated variable names.
   // These separations exist only to easily see the difference when looking at the code here.
   enum class ESystem {
      kTigress,         //0
      kSharc,
      kTriFoil,
      kRF,
      kCSM,
      kSiLi,
      kSiLiS3,
      kGeneric,
      kS3,
      kBambino,
      kTip,             //10
      kGriffin,
      kSceptar,
      kPaces,
      kLaBr,
      kTAC,
      kZeroDegree,
      kDescant,
		kGriffinBgo,
		kLaBrBgo,
      kFipps,           //20
		kBgo,
		kTdrClover,
		kTdrCloverBgo,
		kTdrTigress,
		kTdrTigressBgo,
		kTdrSiLi,
		kTdrPlastic,
		kEmma,
		kEmmaS3,
		kTrific,
		kClear            //30
   };

   ESystem   System() const { return fSystem; }

   void Parse(std::string* name) override;

   void EnumerateDigitizer(TPriorityValue<std::string>& digitizerName, TPriorityValue<EDigitizer>& digitizerType, TPriorityValue<int>& timeStampUnit) override;

	TClass* GetClassType() const override;
   void Print(Option_t* opt = "") const override;
   void Clear(Option_t* opt = "") override;

private:
   ESystem fSystem;

   void EnumerateSystem();

   /// \cond CLASSIMP
   ClassDefOverride(TGRSIMnemonic, 1)
   /// \endcond
};

#endif
