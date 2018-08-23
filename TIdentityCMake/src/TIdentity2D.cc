//Author A. Rustamov
//TIdentity class for calculation of all second moments
#include "TIdentity2D.h"
#include "TIdentityFunctions.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"

using namespace std;
ClassImp(TIdentity2D)

Int_t TIdentity2D::fNParticles;
Int_t TIdentity2D::fNMixParticles;

TIdentityFunctions *TIdentity2D::fFunctions = NULL;

TIdentity2D::TIdentity2D()
{
  InitIden2D(4);
}

TIdentity2D::TIdentity2D(Int_t size)
{

  InitIden2D(size);
}

void TIdentity2D::InitIden2D(Int_t size)
{
  cout << " _________________________________________________________________________________" << endl;
  cout << "|                                                                                 |" << endl;
  cout << "|  *   ******     *******    *      *    *********   *   *********   *         *  |" << endl;
  cout << "|  *   *      *   *          * *    *        *       *       *        *       *   |" << endl;
  cout << "|  *   *      *   *          *  *   *        *       *       *         *     *    |" << endl;
  cout << "|  *   *      *   * *****    *   *  *        *       *       *          * * *     |" << endl;
  cout << "|  *   *      *   *          *    * *        *       *       *             *      |" << endl;
  cout << "|  *   *      *   *          *     **        *       *       *            *       |" << endl;
  cout << "|  *   ******     ********   *      *        *       *       *         * *        |" << endl;
  cout << "|                                                                                 |" << endl;
  cout << "| based on:                                                                       |" << endl;
  cout << "| 1. PRC 84, 024902 (2011)                                                        |" << endl;
  cout << "| 2. PRC 86, 044906 (2012)                                                        |" << endl;
  cout << "| see also: PRC 83, 054907 (2011)                                                 |" << endl;
  cout << "|                                                                                 |" << endl;
  cout << "|contact: a.rustamov@cern.ch                                                      |" << endl;
  cout << "|_________________________________________________________________________________|" << endl;
  cout << " " << endl;

  fDebugFile = new TFile("TIdenDebug.root", "recreate");
  fCountPart = 0;
  fCountPartNeg = 0;
  fCountPartPos = 0;

  fUseSign = -1000;
  fSeparateSign=kFALSE;
  fSign = -1000;
  fMyBinBrach = 0x0;
  fCountVeto = 0;
  fMyDeDx = 0.;
  fDEdx = 0.;

  fPrevEvtVeto = -1;
  fTSize = size;
  fTSizeMixed = size * (size - 1) / 2;
  fSizeMatrix = fTSize + fTSizeMixed;
  fW = new vector<double>[fTSize];
  fW2 = new vector<double>[fTSize];
  fWmixed = new vector<double>[fTSizeMixed];
  if (!fFunctions) fFunctions = new TIdentityFunctions();
  fW_sum = new Double_t[fTSize];
  fPrevEvt = -1000;
  fAver = new Double_t[fTSize];
  fAver2 = new Double_t[fTSize];
  fAverI = new Double_t[fTSize];
  fAverMixed = new Double_t[fTSizeMixed];

  fMValue = new Double_t[fTSize];

  fNParticles = fTSize;
  fNMixParticles = fTSizeMixed;

  fSize_size = 3000; // step size for integral calculation

  fWI = new Double_t *[fTSize];
  fWI2 = new Double_t *[fTSize];
  for (Int_t i = 0; i < fTSize; i++)
  {
    fWI[i] = new Double_t[fTSize];
    fWI2[i] = new Double_t[fTSize];
    for (Int_t j = 0; j < fTSize; j++)
    {
      fWI[i][j] = 0.;
      fWI2[i][j] = 0.;
    }
  }
  fWIMix = new Double_t *[fTSizeMixed];
  for (Int_t i = 0; i < fTSizeMixed; i++)
  {
    fWIMix[i] = new Double_t[fTSize];
    for (Int_t j = 0; j < fTSize; j++)
    {
      fWIMix[i][j] = 0.;
    }
  }
  for (Int_t i = 0; i < fTSize; i++)
  {
    fAverI[i] = 0;
  }
  fTidenMatrixA = new TMatrixD(fSizeMatrix, fSizeMatrix);
  fTidenVectorB = new Double_t[fSizeMatrix];
  fRecMoments = new Double_t[fSizeMatrix];
  for (Int_t i = 0; i < fSizeMatrix; i++)
  {
    fRecMoments[i] = 0.;
  }
}

void TIdentity2D::Reset()
{
  fCountVec.clear();
  fCountVec2.clear();
  fCountVecMix.clear();
  fCountPart = fCountPartNeg = fCountPartPos = 0;

  for (Int_t i = 0; i < fTSize; i++)
  {
    // cout << " TIdentity2D::Reset.Info: resetting vectors" << endl;
    fW[i].clear();
    fW2[i].clear();
    fAver[i] = 0;
    fAver2[i] = 0;
    fAverI[i] = 0;
    fW_sum[i] = 0;
  }

  for (Int_t i = 0; i < fTSizeMixed; i++)
  {
    fWmixed[i].clear();
    fAverMixed[i] = 0;
  }

  fCountVeto = 0;
  fPrevEvtVeto = -1;

  fPrevEvt = -1000;
  for (Int_t i = 0; i < fTSize; i++)
  {
    for (Int_t j = 0; j < fTSize; j++)
    {
      fWI[i][j] = 0.;
      fWI2[i][j] = 0.;
    }
  }
  for (Int_t i = 0; i < fTSizeMixed; i++)
  {
    for (Int_t j = 0; j < fTSize; j++)
    {
      fWIMix[i][j] = 0.;
    }
  }
  for (Int_t i = 0; i < fTSize; i++)
  {
    fAverI[i] = 0;
  }

  for (Int_t i = 0; i < fSizeMatrix; i++)
  {
    fRecMoments[i] = 0.;
    fTidenVectorB[i] = 0;
  }
}

TIdentity2D::~TIdentity2D()
{
  if (fW)
  {
    delete[] fW;
    fW = 0;
  }
  if (fW2)
  {
    delete[] fW2;
    fW2 = 0;
  }
  if (fWmixed)
  {
    delete[] fWmixed;
    fWmixed = 0;
  }
  if (fW_sum)
  delete[] fW_sum;
  if (fAver)
  delete[] fAver;
  if (fAver2)
  delete[] fAver2;
  if (fAverI)
  delete[] fAverI;
  if (fAverMixed)
  delete[] fAverMixed;
  if (fMValue)
  delete[] fMValue;
  if (fRecMoments)
  delete[] fRecMoments;
  if (fTidenMatrixA)
  delete fTidenMatrixA;
  if (fTidenVectorB)
  delete[] fTidenVectorB;
}

void TIdentity2D::GetTree(Long_t &nent, TString idenTreeName)
{

  if (fTIdenFileName.Contains(".root"))
  {
    fTIdentityFile = new TFile(fTIdenFileName);
    cout << " TIdentity2D::GetTree.Info: We are reading the file " << fTIdenFileName << endl;
    fTIdentityTree = (TTree *)fTIdentityFile->Get(idenTreeName);
  }
  else
  {
    fTIdentityTree = GetTreeFromChain(fTIdenFileName, idenTreeName);
  }
  if (!fTIdentityTree)
  cout << " TIdentity2D::GetTree.Error: tree could not be read" << endl;
  cout << " TIdentity2D::GetTree.Info: ======================== " << endl;
  fMyBinBrach = (TBranch *)fTIdentityTree->FindBranch("myBin");
  cout << " TIdentity2D::GetTree.Info: ======================== " << endl;
  if (!fMyBinBrach)
  { // new version of tree format
    fTIdentityTree->SetBranchAddress("gid", &fEventNum);
    fTIdentityTree->SetBranchAddress("dEdx", &fDEdx);
    fTIdentityTree->SetBranchAddress("sign", &fSign);
    fTIdentityTree->SetBranchAddress("cutBit", &fCutBit);
    for (Int_t i = 0; i < fNBranches; i++)
    {
      fTIdentityTree->SetBranchAddress(fBranchNames[i], &fBranchVariables[i]);
    }
  }
  else
  { // old version of tree format
    fTIdentityTree->SetBranchAddress("sign", &fSign);
    fTIdentityTree->SetBranchAddress("myBin", fMyBin);
    fTIdentityTree->SetBranchAddress("myDeDx", &fMyDeDx);
    fTIdentityTree->SetBranchAddress("evtNum", &fEventNumOldVersion);
  }
  fTreeEntries = (Long_t)fTIdentityTree->GetEntries();
  nent = fTreeEntries;
  InitFunctions();
}

TTree *TIdentity2D::GetTreeFromChain(TString treeList, TString treeName)
{
  cout << " TIdentity2D::GetTreeFromChain.Info: Files added to the chain" << endl;
  ifstream fileTmp(treeList);
  char file[255];
  TChain *chain = NULL;
  chain = new TChain(treeName);
  while (fileTmp)
  {
    fileTmp.getline(file, sizeof(file)); // delim defaults to '\n'
    if (fileTmp)
    cout << " TIdentity2D::GetTreeFromChain.Info: " << file << endl;
    chain->Add(file);
  }
  fileTmp.close();
  return (TTree *)chain;
}

Double_t TIdentity2D::GetValue(Double_t *xval, Double_t *par)
{
  Double_t xx = xval[0];
  Int_t index = (Int_t)par[0];
  return fFunctions->GetValue(index, xx);
}

Double_t TIdentity2D::GetFunctions(Double_t *xx, Double_t *par)
{
  Int_t j = (Int_t)par[0];
  Int_t i = (Int_t)par[1];
  Int_t k = (Int_t)par[2];
  Double_t val[fNParticles];
  Double_t sumVal = 0;
  for (Int_t ii = 0; ii < fNParticles; ii++)
  {
    val[ii] = fFunctions->GetValue(ii, xx[0]);
    sumVal += val[ii];
  }
  Double_t relVal[fNParticles];
  if (sumVal < 1e-15)
  return 0.;
  for (Int_t m = 0; m < fNParticles; m++)
  {
    relVal[m] = val[m] / sumVal;
    if (k == 2)
    relVal[m] *= relVal[m];
  }
  return relVal[j] * val[i];
}

Double_t TIdentity2D::GetFunctionsMix(Double_t *xx, Double_t *par)
{
  Int_t j = (Int_t)par[0];
  Int_t k = (Int_t)par[1];
  Double_t val[fNParticles];
  Double_t sumVal = 0;
  for (Int_t ii = 0; ii < fNParticles; ii++)
  {
    val[ii] = fFunctions->GetValue(ii, xx[0]);
    sumVal += val[ii];
  }
  if (sumVal < 1e-15)
  return 0.;
  Double_t myVal[fNMixParticles];
  Int_t t = 0;
  for (Int_t m = 0; m < fNParticles - 1; m++)
  {
    for (Int_t n = m + 1; n < fNParticles; n++)
    {
      myVal[t] = val[m] * val[n] / sumVal / sumVal;
      t++;
    }
  }
  return myVal[j] * val[k];
}

void TIdentity2D::InitFunctions()
{
  for (Int_t i = 0; i < fTSize + 1; i++)
  {
    sprintf(fTFunctionsName, "func[%d]", i); // ??? function names ???
    fTFunctions[i] = new TF1(fTFunctionsName, GetValue, fMindEdx, fMaxdEdx, 1);
    fTFunctions[i]->SetParameter(0, i);
  }

  for (Int_t i = 0; i < fTSize; i++)
  {
    for (Int_t j = 0; j < fTSize; j++)
    {
      sprintf(fTFunctionsName, "Ifunc[%d][%d]", i, j);
      fIFunctions[i][j] = new TF1(fTFunctionsName, GetFunctions, fMindEdx, fMaxdEdx, 3);
      fIFunctions[i][j]->SetParameter(0, i);
      fIFunctions[i][j]->SetParameter(1, j);
      fIFunctions[i][j]->SetParameter(2, 1);

      sprintf(fTFunctionsName, "Ifunc2[%d][%d]", i, j);
      fIFunctions2[i][j] = new TF1(fTFunctionsName, GetFunctions, fMindEdx, fMaxdEdx, 3);
      fIFunctions2[i][j]->SetParameter(0, i);
      fIFunctions2[i][j]->SetParameter(1, j);
      fIFunctions2[i][j]->SetParameter(2, 2);
    }
  }
  for (Int_t i = 0; i < fTSizeMixed; i++)
  {
    for (Int_t j = 0; j < fTSize; j++)
    {
      sprintf(fTFunctionsName, "IfuncMix[%d][%d]", i, j);
      fIFunctionsMix[i][j] = new TF1(fTFunctionsName, GetFunctionsMix, fMindEdx, fMaxdEdx, 2);
      fIFunctionsMix[i][j]->SetParameter(0, i);
      fIFunctionsMix[i][j]->SetParameter(1, j);
    }
  }
}

void TIdentity2D::SetFunctionPointers(fptr fun)
{
  fFunctions->funs = fun;
}

void TIdentity2D::Run()
{
  ;
}

Bool_t TIdentity2D::GetEntry(Int_t i)
{
  fTIdentityTree->GetEntry(i);
  if (fMyBinBrach) fEventNum = fEventNumOldVersion;
  if (fMyBinBrach && fSeparateSign) fMyDeDx = fMyDeDx*fSign; // treat particle and antiparticle separately

  if (fEventNum != fPrevEvtVeto && fPrevEvtVeto > 0)
  {
    fCountVeto++;
  }
  fPrevEvtVeto = fEventNum;
  if (fDEdx != 0) fMyDeDx = fDEdx;
  if ((fMyDeDx < fMindEdx || fMyDeDx > fMaxdEdx) || fMyDeDx == 0) return kFALSE;
  //
  // secure the usage of sign=0 which is sum of + and - particles
  //
  if (!(fSign == fUseSign || fUseSign == 0)) return kFALSE;
  return kTRUE;
}

Int_t TIdentity2D::AddEntry()
{

  if (fEventNum == fPrevEvt)
  {
    AddParticles();
  }
  else
  {
    if (fCount != 0)
    {
      fCountVec.push_back(fCountPart);
      fCountVec2.push_back(fCountPart * fCountPart);
      fCountVecMix.push_back(fCountPartNeg * fCountPartPos);

      for (Int_t m = 0; m < fTSize; m++)
      {
        fW[m].push_back(fW_sum[m]);
        fW2[m].push_back(fW_sum[m] * fW_sum[m]);
      }
      //
      // Debug hists
      for (Int_t i = 0; i < fTSize; i++)
      {
        fHistWs[i]->Fill(fW[i][fW[i].size() - 1]);
      }
      //
      Int_t t = 0;
      for (Int_t m = 0; m < fTSize - 1; m++)
      {
        for (Int_t n = m + 1; n < fTSize; n++)
        {
          fWmixed[t].push_back(fW_sum[m] * fW_sum[n]);
          t++;
        }
      }
    }
    ResetValues();
    fCount = 0;
    AddParticles();
  }
  fPrevEvt = fEventNum;
  return 1;
}

void TIdentity2D::ResetValues()
{
  for (Int_t s = 0; s < fTSize; s++)
  {
    fW_sum[s] = 0.;
  }
  fCountPart = 0;
  fCountPartNeg = 0;
  fCountPartPos = 0;
}

void TIdentity2D::AddParticles()
{
  //Double_t fMValue[fTSize] = {0.};
  //Double_t fMValue[1000] = {0.};

  Double_t sumfMValue = 0;

  for (Int_t i = 0; i < fTSize; i++)
  {

    if (fMyDeDx == 0)
    {
      fMValue[i] = 0;
      fCount++;
      continue;
    }
    fMValue[i] = fTFunctions[i]->Eval(fMyDeDx);
    sumfMValue += fMValue[i];
  }
  //
  // Debug hists for omega values
  for (Int_t i = 0; i < fTSize; i++)
  {
    fHistOmegas[i]->Fill(fMValue[i] / sumfMValue);
  }
  //
  fCountPart += 1;
  if (fSign == -1)
  {
    fCountPartNeg += 1;
  }
  if (fSign == 1)
  {
    fCountPartPos += 1;
  }

  if (sumfMValue > 1e-10)
  {
    fCount++;
    for (Int_t i = 0; i < fTSize; i++)
    {
      fW_sum[i] += fMValue[i] / sumfMValue;
      ;
    }
  }
}

void TIdentity2D::Finalize()
{
  cout << " " << endl;
  cout << " TIdentity2D::Finalize.Info: ***************************************************" << endl;
  cout << " TIdentity2D::Finalize.Info: ************ number of analyzed events: " << fCountVeto + 1 << " ******" << endl;
  cout << " TIdentity2D::Finalize.Info: ***************************************************" << endl;
  cout << " " << endl;

  for (Int_t m = 0; m < fTSize; m++)
  {
    fAver[m] = accumulate(fW[m].begin(), fW[m].end(), 0.0) / fCountVeto;
    fAver2[m] = accumulate(fW2[m].begin(), fW2[m].end(), 0.0) / fCountVeto;
  }
  Int_t t = 0;
  for (Int_t m = 0; m < fTSize - 1; m++)
  {
    for (Int_t n = m + 1; n < fTSize; n++)
    {
      fAverMixed[t] = accumulate(fWmixed[t].begin(), fWmixed[t].end(), 0.0) / fCountVeto;
      t++;
    }
  }

  fAverCount[0] = accumulate(fCountVec.begin(), fCountVec.end(), 0.0) / fCountVeto;
  fAverCount[1] = accumulate(fCountVec2.begin(), fCountVec2.end(), 0.0) / fCountVeto;
  fAverCount[2] = accumulate(fCountVecMix.begin(), fCountVecMix.end(), 0.0) / fCountVeto;

  //
  // Write some output to data
  fDebugFile->cd();
  for (Int_t i = 0; i < fTSize; i++)
  fHistWs[i]->Write();
  for (Int_t i = 0; i < fTSize; i++)
  fHistOmegas[i]->Write();
  fDebugFile->Close();
  delete fDebugFile;
  //
}

void TIdentity2D::GetBins(const Int_t nExtraBins, Double_t *bins)
{
  // TString branchNames[nBranches]={"eta","cent","ptot","cRows","tpcchi2"};
  if (!fMyBinBrach)
  {
    bins[0] = fEventNum;
    bins[1] = fDEdx;
    bins[2] = fSign;
    bins[3] = fCutBit;
    for (Int_t i = 0; i < nExtraBins; i++)
    bins[i + 4] = fBranchVariables[i];
  }
  else
  {
    // 0 --> eta, 1 --> cent, 2 --> ptot
    bins[0] = fMyBin[0];
    bins[1] = fMyBin[1];
    bins[2] = fMyBin[2];
    bins[3] = fSign;
    bins[4] = fMyDeDx;
  }
}

Double_t TIdentity2D::GetIntegral(Int_t i, Int_t j, Int_t k)
{
  if (k == 1)
  {
    return MyIntegral(fIFunctions[i][j]);
  }
  else
  return MyIntegral(fIFunctions2[i][j]);
}

Double_t TIdentity2D::GetIntegralMix(Int_t i, Int_t j)
{
  return MyIntegral(fIFunctionsMix[i][j]);
}

void TIdentity2D::SetLimits(Float_t min, Float_t max, Double_t BW, Double_t nParticlesMax, Int_t wDistResol)
{
  fMindEdx = min;
  fMaxdEdx = max;
  fDedxBinWidth = (max - min) / BW;
  //
  //
  Int_t nBinsNparticles = nParticlesMax * wDistResol;
  cout << " TIdentity2D::SetLimits.Info: dEdx histogram bin range is being set for " << fTSize << " particle " << endl;
  fHistWs = new TH1D *[fTSize];
  fHistOmegas = new TH1D *[fTSize];
  for (Int_t i = 0; i < fTSize; i++)
  {
    fHistWs[i] = new TH1D(Form("hW_%d", i), Form("hW_%d", i), nBinsNparticles, 0., nParticlesMax);
    fHistOmegas[i] = new TH1D(Form("hOmega_%d", i), Form("hOmega_%d", i), 100, 0., 1.);
  }
}

Double_t TIdentity2D::MyIntegral(TF1 *Fun)
{

  Double_t xx, sum = 0;
  Double_t step = (fMaxdEdx - fMindEdx) / (2 * fSize_size);
  for (Int_t i = 1; i < 2 * fSize_size + 1; i++)
  {
    xx = fMindEdx + step * i;
    sum += Fun->Eval(xx);
  }
  return sum * step;
}

void TIdentity2D::AddIntegrals(Int_t lsign)
{
  if (fUseSign != lsign)
  return;

  for (Int_t i = 0; i < fTSize; i++)
  {
    for (Int_t j = 0; j < fTSize; j++)
    {
      fWI[i][j] += GetIntegral(i, j, 1);
      fWI2[i][j] += GetIntegral(i, j, 2);
    }
  }

  for (Int_t t = 0; t < fTSizeMixed; t++)
  {
    for (Int_t j = 0; j < fTSize; j++)
    {
      fWIMix[t][j] += GetIntegralMix(t, j);
    }
  }

  for (Int_t i = 0; i < fTSize; i++)
  {
    fAverI[i] += MyIntegral(fTFunctions[i]) / fDedxBinWidth;
  }
}

Double_t TIdentity2D::GetWI(Int_t i, Int_t j, Int_t k)
{
  if (k == 1)
  return fWI[i][j] / fAverI[j] / fDedxBinWidth;
  else if (k == 2)
  return fWI2[i][j] / fAverI[j] / fDedxBinWidth;
  else
  return fWIMix[i][j] / fAverI[j] / fDedxBinWidth;
}

void TIdentity2D::CalcMoments()
{
  Int_t t = fTSize;
  Int_t tt = fTSize;
  for (Int_t i = 0; i < fTSize; i++)
  {
    t = fTSize;
    for (Int_t j = 0; j < fTSize; j++)
    {
      (*fTidenMatrixA)(i, j) = GetWI(i, j, 1) * GetWI(i, j, 1);
      if (j == fTSize - 1)
      continue;
      for (Int_t k = j + 1; k < fTSize; k++)
      {
        (*fTidenMatrixA)(i, t) = 2. * GetWI(i, j, 1) * GetWI(i, k, 1);
        t++;
      }
    }
  }

  t = fTSize;
  for (Int_t i = 0; i < fTSize - 1; i++)
  {
    for (Int_t j = i + 1; j < fTSize; j++)
    {
      tt = fTSize;
      for (Int_t k = 0; k < fTSize; k++)
      {
        (*fTidenMatrixA)(t, k) = GetWI(i, k, 1) * GetWI(j, k, 1);
        if (k == fTSize - 1)
        continue;
        for (Int_t kk = k + 1; kk < fTSize; kk++)
        {
          (*fTidenMatrixA)(t, tt) = GetWI(i, k, 1) * GetWI(j, kk, 1) + GetWI(i, kk, 1) * GetWI(j, k, 1);
          tt++;
        }
      }
      t++;
    }
  }

  TMatrixD invA = (*fTidenMatrixA).Invert();

  for (Int_t i = 0; i < fTSize; i++)
  {
    fTidenVectorB[i] = fAver2[i];
    for (Int_t j = 0; j < fTSize; j++)
    {
      fTidenVectorB[i] -= fAver[j] * (GetWI(i, j, 2) - GetWI(i, j, 1) * GetWI(i, j, 1));
    }
  }
  Int_t indA, indB;
  for (Int_t kk = 0; kk < fTSizeMixed; kk++)
  {
    fTidenVectorB[kk + fTSize] = fAverMixed[kk];
    GetIndex(kk, indA, indB);
    for (Int_t m = 0; m < fTSize; m++)
    {
      fTidenVectorB[kk + fTSize] -= fAver[m] * (GetWI(kk, m, 0) - GetWI(indA, m, 1) * GetWI(indB, m, 1));
    }
  }

  for (Int_t k = 0; k < fSizeMatrix; k++)
  {
    for (Int_t tt = 0; tt < fSizeMatrix; tt++)
    {
      fRecMoments[k] += invA(k, tt) * fTidenVectorB[tt];
    }
  }
}

void TIdentity2D::GetMoments()
{
  ;
}

Double_t TIdentity2D::GetSecondMoment(Int_t i)
{
  if (i >= fSizeMatrix)
  {
    cout << " TIdentity2D::GetSecondMoment.Info: out of bound" << endl;
    return -1000.;
  }
  return fRecMoments[i];
}

Double_t TIdentity2D::GetMixedMoment(Int_t i, Int_t j)
{
  Int_t tmp = i;
  Int_t ind = fTSize - 1;
  if (i > j)
  {
    i = j;
    j = tmp;
  }
  for (Int_t k = 0; k < fTSize - 1; k++)
  {
    for (Int_t kk = k + 1; kk < fTSize; kk++)
    {
      ind++;
      if (i == k && j == kk)
      break;
    }
    if (i == k)
    break;
  }
  return fRecMoments[ind];
}

Double_t TIdentity2D::GetMean(Int_t i)
{
  if (i >= fTSize)
  {
    cout << " TIdentity2D::GetMean.Info: out of bound" << endl;
    return -1000.;
  }
  return fAver[i];
}

Double_t TIdentity2D::GetMeanI(Int_t i)
{
  if (i >= fTSize)
  {
    cout << " TIdentity2D::GetMeanI.Info: out of bound" << endl;
    return -1000.;
  }
  return fAverI[i];
}

Double_t TIdentity2D::GetNuDyn(Int_t i, Int_t j)
{
  Double_t nydyn = GetSecondMoment(i) / GetMean(i) / GetMean(i);
  nydyn += GetSecondMoment(j) / GetMean(j) / GetMean(j);
  nydyn -= 2. * GetMixedMoment(i, j) / GetMean(i) / GetMean(j);
  nydyn -= (1. / GetMean(i) + 1. / GetMean(j));
  return nydyn;
}

Int_t TIdentity2D::GetIndex(Int_t k, Int_t &a, Int_t &b)
{
  Int_t t = 0;
  for (Int_t m = 0; m < fTSize - 1; m++)
  {
    for (Int_t n = m + 1; n < fTSize; n++)
    {
      if (t == k)
      {
        a = m;
        b = n;
        return 1;
      }
      t++;
    }
  }
  return 1;
}
