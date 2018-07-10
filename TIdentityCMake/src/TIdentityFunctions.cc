//Author A. Rustamov
//An interface for PID functions

#include "TIdentityFunctions.h"

using namespace std;

ClassImp(TIdentityFunctions)

TIdentityFunctions::TIdentityFunctions()
{
  funs = NULL;
}

TIdentityFunctions::~TIdentityFunctions()
{
  ;
}

Double_t TIdentityFunctions::GetValue(Int_t i, Double_t x)
{
  return funs(i,x);
}


