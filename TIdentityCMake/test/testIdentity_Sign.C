#include "TIdentity2D.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TObject.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include <iomanip>
#include "iostream"
#include "string"
using namespace std;
using std::cout;
using std::setw;

// =======================================================================================================
// =======================================================================================================
// Helper Functions
void      InitializeObjects();
void      ReadFitParamsFromLineShapes(TString paramTreeName);
void      RetrieveMoments(TIdentity2D *tidenObj);
Double_t  EvalFitValue(Int_t particle, Double_t x);
// =======================================================================================================
// =======================================================================================================
//
// ======= Modification Part =============================================================================
const Int_t fNParticles          = 5;                       // number of particles species
const Int_t fnSignBins           = 3;                       // process cases --> 0: sum of particles and anti-particles, 1: particles, -1: antiparticles
TString     treeIdentity         = "tracks";                // name of the input TTree
TString     lookUpCloneArrayName = "funcLineShapesCArr";    // name of the lookup table which contains line shapes
const Int_t nBinsLineShape       = 1000;                    // number of bins in case of the line shapes are used as histograms
Int_t       fnTestEntries        = 0;                       // 0 for full sample processing and n>0 for testing
Int_t       lookUpTableLineMode  = 0;                       // 0 for hist (faster but limited precision due to number of bins) and 1 for func
//
// fixed tree branches --> [0]=event; [1]=dEdx; [2]=sign; [3]=cutBit; ||||  [4]=cent;
Double_t fTreeVariablesArray[5];                            // 4 default and + 1 extra branch
const Int_t nBranches = 1;                                  // number of extra branches
TString branchNames[nBranches]={"cent"};                    // names of the extra branches
// =======================================================================================================
//
// Inputs
Char_t  inputfileNameDataTree[255];     //file name of tree
Char_t  inputfileNameLineShapes[255];   // file name for fit function
TString fileNameDataTree = "";
TString fileNameLineShapes = "";
// =======================================================================================================
//
Double_t nEvents = 0;
Double_t nnorm   = 1.;
Int_t fUsedSign;
Int_t fSignBin;
TFile *fLineShapesLookUpTable = NULL;
TClonesArray *cloneArrFunc=NULL;
TH1D ***hLineShape;
TF1 ***fLineShape;
TFile *outFile;
TH1D *hDedxDebug[2];
static TH1D *hDedxDebugLineShapes[2][fNParticles];
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kDe=4,
  kElEl=5,kPiPi=6,kKaKa=7,kPrPr=8,kDeDe=9,
  kElPi=10,kElKa=11,kElPr=12,kElDe=13,
  kPiKa=14,kPiPr=15,kPiDe=16,
  kKaPr=17,kKaDe=18,
};

//
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
//
int main(int argc, char *argv[])
{

  //
  //  Arguments: $dataTree $lineShapes $sign
  //

  cout << " main.Info: NUMBER OF ARGUMENTS "<<argc<<endl;
  if(argc == 4)
  {
    sprintf(inputfileNameDataTree,"%s",argv[1]);
    sprintf(inputfileNameLineShapes,"%s",argv[2]);
    fUsedSign = atoi(argv[3]);
    cout<<" main.Info: read file names from input "<<endl;
  }
  else
  {
    cout<<" main.Error: wrong input list"<<endl;
  }
  //
  InitializeObjects();
  fileNameDataTree   = inputfileNameDataTree;
  fileNameLineShapes = inputfileNameLineShapes;
  //
  // Initialize objects and get the bin information
  TROOT IdentityMethod("IdentityMethod","compiled identity method");
  ReadFitParamsFromLineShapes(fileNameLineShapes);
  //
  // Create the TIdentity2D object and start analysis
  TIdentity2D *iden4 = new TIdentity2D(fNParticles);      // Set the number of particles to 4
  iden4 -> SetFileName(fileNameDataTree);
  iden4 -> SetBranchNames(nBranches,branchNames);
  iden4 -> SetFunctionPointers(EvalFitValue);
  iden4 -> SetLimits(0.,50.,1500.,30.,50); // --> (dEdxMin,dEdxMax,nBinsUsed in dEdx), if slice histograms are scaled wrt binwidth, then binwidth=1
  iden4 -> SetUseSign(fUsedSign);  // pass input sign value to TIdentity module
  Long_t nEntries;
  iden4 -> GetTree(nEntries,treeIdentity);
  iden4 -> Reset();
  //
  // track by track loop --> read all track info  and add tracks to the iden4 object
  if (fnTestEntries>0) nEntries = fnTestEntries;
  Int_t fUsedBins[fnSignBins]={0};
  Int_t countEntry = 0;
  for( Int_t i = 0; i < nEntries; i++ )
  {
    //
    // Read the entries and corresponding bins
    if( !iden4 ->  GetEntry(i) ) continue;
    iden4 -> GetBins(nBranches, fTreeVariablesArray);    // reads identity tree and retrives binned info
    //
    if (fTreeVariablesArray[2]==-1) { fSignBin=0; hDedxDebug[0]->Fill(fTreeVariablesArray[1]); }
    if (fTreeVariablesArray[2]== 1) { fSignBin=1; hDedxDebug[1]->Fill(fTreeVariablesArray[1]); }
    //
    // tag the bins an read only required bins
    if (fUsedSign==0) { fUsedBins[fSignBin] = 1; iden4 -> AddEntry(); countEntry++; }
    else { fUsedBins[fSignBin] = 1; iden4 -> AddEntry(); countEntry++; }
  }
  //
  cout << " main.Info: Total number of tracks processed = " << countEntry << endl;
  iden4 -> Finalize();
  //
  // Combine different phase space bins
  for(Int_t isign = 0; isign < fnSignBins-1; isign++) {
    if(fUsedBins[isign] != 1) continue;
    fSignBin = isign;  // to be used in retrival of obj from the lookup table
    iden4  -> AddIntegrals(fUsedSign); // real sign information passed
  }
  //
  iden4 -> CalcMoments();
  RetrieveMoments(iden4);
  delete iden4;
  //
  // Dump some debug output
  outFile->cd();
  for (Int_t isign=0;isign<2;isign++){
    hDedxDebug[isign]->Write();
    for (Int_t i=0;i<fNParticles;i++){
      hDedxDebugLineShapes[isign][i] = (TH1D*)hLineShape[i][isign]->Clone();
      hDedxDebugLineShapes[isign][i]->SetName(Form("LineShape_sign_%d_part_%d",isign,i));
      hDedxDebugLineShapes[isign][i]->Write();
    }
  }

  outFile -> Close();
  delete outFile;
  return 1;
}
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
void ReadFitParamsFromLineShapes(TString paramTreeName)
{

  //
  // Read the lineShapes from the lookup table and fill them into pointer arrays
  //

  fLineShapesLookUpTable = new TFile(paramTreeName);
  cloneArrFunc   = (TClonesArray*)fLineShapesLookUpTable->Get(lookUpCloneArrayName);
  if (!cloneArrFunc) cout << " ReadFitParamsFromLineShapes.Error: cloneArrFunc is empty " << endl;
  for (Int_t ipart = 0; ipart<fNParticles; ipart++){
    for (Int_t isign = 0; isign<fnSignBins; isign++){
      TString objName = Form("particle_%d_bin_%d",ipart,isign);
      fLineShape[ipart][isign] = (TF1*)cloneArrFunc->FindObject(objName);
      fLineShape[ipart][isign]->SetName(objName);
      fLineShape[ipart][isign]->SetNpx(nBinsLineShape);
      hLineShape[ipart][isign] = (TH1D*)fLineShape[ipart][isign]->GetHistogram();
      hLineShape[ipart][isign]->SetName(objName);
    }
  }

}
// =======================================================================================================
Double_t EvalFitValue(Int_t particle, Double_t x)
{

  //
  // Calculates rho_j_x quantity for a given track
  //

  Double_t rho_j_x = 0.;
  Int_t bin = hLineShape[particle][fSignBin]->FindBin(x);
  if (lookUpTableLineMode==0) rho_j_x = hLineShape[particle][fSignBin]->GetBinContent(bin);
  if (lookUpTableLineMode==1) rho_j_x = fLineShape[particle][fSignBin]->Eval(x);
  return rho_j_x;

}
// =======================================================================================================
void InitializeObjects()
{

  //
  // Initialize output debug file and the pointers to hold the lineShapes
  //

  cout << " ================================================================================= " << endl;
  cout << " InitializeObjects.Info: Input sign            = " << fUsedSign                   << endl;
  cout << " InitializeObjects.Info: treeIdentity          = " << treeIdentity                << endl;
  cout << " InitializeObjects.Info: data Tree             = " << inputfileNameDataTree       << endl;
  cout << " InitializeObjects.Info: Line Shapes           = " << inputfileNameLineShapes     << endl;
  cout << " ================================================================================= " << endl;
  //
  outFile = new TFile(Form("TIdentity_Moments_%d.root",fUsedSign),"recreate");
  for (Int_t i=0; i<2; i++) hDedxDebug[i] = new TH1D(Form("hDedxDebug_%d",i),Form("hDedxDebug_%d",i),1500,0.,50.);
  //
  // Initialize pointers to lookup table
  fLineShape = new TF1 **[fNParticles];
  for (Int_t ipart = 0; ipart<fNParticles; ipart++){
    fLineShape[ipart] = new TF1*[fnSignBins];
    for (Int_t isign = 0; isign<fnSignBins; isign++){
      fLineShape[ipart][isign] = NULL;
    }
  }
  hLineShape = new TH1D **[fNParticles];
  for (Int_t ipart = 0; ipart<fNParticles; ipart++){
    hLineShape[ipart] = new TH1D*[fnSignBins];
    for (Int_t isign = 0; isign<fnSignBins; isign++){
      hLineShape[ipart][isign] = NULL;
    }
  }

}
// =======================================================================================================
void RetrieveMoments(TIdentity2D *tidenObj)
{

  //
  // Retrive moments and print them
  //
  const Int_t nMoments = 19;
  Double_t momArray[nMoments]={0.}, intArray[nMoments]={0.};
  //
  momArray[kEl] = tidenObj -> GetMean(kEl);
  momArray[kPi] = tidenObj -> GetMean(kPi);
  momArray[kKa] = tidenObj -> GetMean(kKa);
  momArray[kPr] = tidenObj -> GetMean(kPr);
  momArray[kDe] = tidenObj -> GetMean(kDe);
  //
  // Second Moments
  momArray[kElEl] = tidenObj -> GetSecondMoment(kEl);
  momArray[kPiPi] = tidenObj -> GetSecondMoment(kPi);
  momArray[kKaKa] = tidenObj -> GetSecondMoment(kKa);
  momArray[kPrPr] = tidenObj -> GetSecondMoment(kPr);
  momArray[kDeDe] = tidenObj -> GetSecondMoment(kDe);
  //
  // Mixed Moments
  momArray[kElPi] = tidenObj -> GetMixedMoment(kEl,kPi);
  momArray[kElKa] = tidenObj -> GetMixedMoment(kEl,kKa);
  momArray[kElPr] = tidenObj -> GetMixedMoment(kEl,kPr);
  momArray[kElDe] = tidenObj -> GetMixedMoment(kEl,kDe);
  momArray[kPiKa] = tidenObj -> GetMixedMoment(kPi,kKa);
  momArray[kPiPr] = tidenObj -> GetMixedMoment(kPi,kPr);
  momArray[kPiDe] = tidenObj -> GetMixedMoment(kPi,kDe);
  momArray[kKaPr] = tidenObj -> GetMixedMoment(kKa,kPr);
  momArray[kKaDe] = tidenObj -> GetMixedMoment(kKa,kDe);
  //
  //Integrals:
  intArray[kEl] = tidenObj -> GetMeanI(kEl);
  intArray[kPi] = tidenObj -> GetMeanI(kPi);
  intArray[kKa] = tidenObj -> GetMeanI(kKa);
  intArray[kPr] = tidenObj -> GetMeanI(kPr);
  intArray[kDe] = tidenObj -> GetMeanI(kDe);
  //
  // Printing
  nnorm   = momArray[kPi]/intArray[kPi];
  nEvents = tidenObj -> GetNEvents();
  cout << " =============================== Summary of Moments =============================== "<<endl;
  cout << " #Particles  : " << fNParticles << endl;
  cout << " events      : " << nEvents+1 << endl;
  cout << " ================================================================================== "<<endl;
  cout << " electron         : "<< momArray[kEl]   <<" int: "<< intArray[kEl]*nnorm << "  ratio: " << momArray[kEl]/(intArray[kEl]*nnorm) << endl;
  cout << " pion             : "<< momArray[kPi]   <<" int: "<< intArray[kPi]*nnorm << "  ratio: " << momArray[kPi]/(intArray[kPi]*nnorm) << endl;
  cout << " kaon             : "<< momArray[kKa]   <<" int: "<< intArray[kKa]*nnorm << "  ratio: " << momArray[kKa]/(intArray[kKa]*nnorm) << endl;
  cout << " proton           : "<< momArray[kPr]   <<" int: "<< intArray[kPr]*nnorm << "  ratio: " << momArray[kPr]/(intArray[kPr]*nnorm) << endl;
  cout << " deuteron         : "<< momArray[kDe]   <<" int: "<< intArray[kDe]*nnorm << "  ratio: " << momArray[kDe]/(intArray[kDe]*nnorm) << endl;
  cout << " electron2        : "<< momArray[kElEl]  <<endl;
  cout << " pion2            : "<< momArray[kPiPi]  <<endl;
  cout << " kaon2            : "<< momArray[kKaKa]  <<endl;
  cout << " proton2          : "<< momArray[kPrPr]  <<endl;
  cout << " deuteron2        : "<< momArray[kDeDe]  <<endl;
  cout << " electronPion     : "<< momArray[kElPi]  <<endl;
  cout << " electronKaon     : "<< momArray[kElKa]  <<endl;
  cout << " electronProton   : "<< momArray[kElPr]  <<endl;
  cout << " electronDeuteron : "<< momArray[kElDe]  <<endl;
  cout << " kaonProton       : "<< momArray[kKaPr]  <<endl;
  cout << " kaonDeuteron     : "<< momArray[kKaDe]  <<endl;
  cout << " ================================================================================== "<<endl;

}
