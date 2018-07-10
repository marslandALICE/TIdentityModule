#ifndef TIDENTITYBASE_H
#define TIDENTITYBASE_H
#include "TNamed.h"
#include "fstream"
#include "iostream"
#include "TString.h"

class TIdentityBase:public TNamed
{
 protected:
  virtual void Run()       = 0;
  virtual void SetFunctionPointers(Double_t (*fun)(Int_t, Double_t)) = 0;
 public:
  ClassDef(TIdentityBase,0)
    };   
#endif
    
    
