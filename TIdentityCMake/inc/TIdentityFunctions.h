#ifndef TIDENTITYFUNCTIONS_H
#define TIDENTITYFUNCTIONS_H
#include "TNamed.h"
#include "fstream"
#include "iostream"

typedef Double_t (*fptr)(Int_t, Double_t);
class TIdentityFunctions:public TNamed
{
 public:
  TIdentityFunctions();
  ~TIdentityFunctions();
  Double_t GetValue(Int_t i, Double_t x);
  //private:
  fptr funs;
  ClassDef(TIdentityFunctions,0)
    };
#endif
