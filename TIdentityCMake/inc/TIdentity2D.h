#ifndef TIDENTITY2D_H
#define TIDENTITY2D_H
#include "TNamed.h"
#include "TBranch.h"
#include "fstream"
#include "iostream"
#include "TString.h"
#include "TIdentityBase.h"
#include "stdlib.h"
#include "TIdentityFunctions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include "TString.h"
#include "TH1.h"
#include "TMatrixD.h"
class TFile;
class TTree;
class TH1D;

class TIdentity2D : public TIdentityBase
{
public:
  TIdentity2D();
  TIdentity2D(Int_t size);
  virtual ~TIdentity2D();
  //
  //
  static Double_t GetValue(Double_t *, Double_t *);
  static Double_t GetFunctions(Double_t *, Double_t *);
  static Double_t GetFunctionsMix(Double_t *, Double_t *);
  static TIdentityFunctions *fFunctions;
  //
  //
  virtual void SetFunctionPointers(fptr);
  virtual void Run();
  //
  //
  void InitIden2D(Int_t);
  void CalcMoments();
  void GetMoments();
  void SetFileName(TString _fileName) { fTIdenFileName = _fileName; }
  void GetTree(Long_t &n, TString idenTreeName);
  void Finalize();
  void GetBins(const Int_t nExtraBins, Double_t *);
  void InitFunctions();
  void ResetValues();
  void AddParticles();
  void AddIntegrals(Int_t);
  void SetLimits(Float_t, Float_t, Double_t, Double_t, Int_t);
  void SetUseSign(Int_t _useSign) { fUseSign = _useSign; }
  void SetSeparateSign(Bool_t _useSeparateSign) { fSeparateSign = _useSeparateSign; }
  void Reset();
  void SetBranchNames(const Int_t tmpNBranches, TString tmpBranchNameArr[])
  {
    fNBranches = tmpNBranches;
    fBranchNames = new TString[fNBranches];
    fBranchVariables = new Float_t[fNBranches];
    for (Int_t i = 0; i < fNBranches; i++)
    {
      fBranchNames[i] = tmpBranchNameArr[i];
      fBranchVariables[i] = 0.;
    }
  }
  //
  //
  Double_t GetSecondMoment(Int_t);
  Double_t GetMixedMoment(Int_t, Int_t);
  Double_t GetNuDyn(Int_t, Int_t);
  Double_t GetMean(Int_t);
  Double_t GetMeanI(Int_t);
  Double_t MyIntegral(TF1 *);
  Double_t GetIntegral(Int_t, Int_t, Int_t);
  Double_t GetIntegralMix(Int_t, Int_t);
  Double_t GetWI(Int_t, Int_t, Int_t);
  Double_t GetDeDx() { return fMyDeDx; }
  Double_t GetAverCount(Int_t i) { return fAverCount[i]; }
  Int_t GetIndex(Int_t, Int_t &, Int_t &);
  Int_t AddEntry();
  Int_t GetNEvents() const { return fCountVeto; }
  Bool_t GetEntry(Int_t);
  TTree *GetTreeFromChain(TString treeList, TString treeName);
  //
  //

private:
  static Int_t fNParticles;
  static Int_t fNMixParticles;
  //
  TMatrixD *fTidenMatrixA;
  Double_t *fTidenVectorB;
  Float_t fMindEdx;
  Float_t fMaxdEdx;
  Double_t fDedxBinWidth;
  //
  Double_t **fWI;
  Double_t **fWI2;
  Double_t **fWIMix;
  Double_t *fRecMoments;
  Double_t *fW_sum;
  Double_t *fAver;
  Double_t *fAverMixed;
  Double_t *fAver2;
  Double_t *fAverI;

  Double_t *fMValue;

  Float_t *fBranchVariables;
  TBranch *fMyBinBrach; // just to check which kind of tree is used
  TString *fBranchNames;
  TString fTIdenFileName;
  Char_t fTFunctionsName[255];
  //
  Double_t fAverCount[3];
  Double_t fMyDeDx;
  Float_t fDEdx;
  Int_t fNBranches;
  Int_t fEventNumOldVersion;
  Int_t fCountVeto;
  Int_t fMyBin[3];
  Int_t fUseSign;
  Bool_t fSeparateSign;
  Int_t fSign;
  Int_t fSize_size;
  Int_t fTSize;
  Int_t fSizeMatrix;
  Int_t fTSizeMixed;
  Int_t fCount;
  Int_t fCountPart;
  Int_t fCountPartNeg;
  Int_t fCountPartPos;
  //
  ULong64_t fPrevEvt;
  ULong64_t fPrevEvtVeto;
  ULong64_t fEventNum;
  Long_t fTreeEntries;
  UInt_t fCutBit;
  //
  TFile *fTIdentityFile;
  TFile *fDebugFile;
  TTree *fTIdentityTree;
  TH1D **fHistWs;
  TH1D **fHistOmegas;
  TF1 *fTFunctions[10];
  TF1 *fIFunctions[50][50];
  TF1 *fIFunctions2[50][50];
  TF1 *fIFunctionsMix[50][50];
  //
  std::vector<double> *fW;
  std::vector<double> *fW2;
  std::vector<double> *fWmixed;
  std::vector<double> fCountVec;
  std::vector<double> fCountVec2;
  std::vector<double> fCountVecMix;

  ClassDef(TIdentity2D, 0)
};
#endif
