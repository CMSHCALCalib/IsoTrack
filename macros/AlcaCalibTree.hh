//////////////////////////////////////////////////////////
// Header with AlcaCalibTree class
//////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TString.h>
#include <TF1.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>


//**********************************************************
// Class with TTree containing parameters of selected events
//**********************************************************
class AlcaCalibTree {
public :
  TChain          *fChain;   //!pointer to the analyzed TTree
  //TChain          *inChain;   //!pointer to the analyzed TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
   Int_t           t_Run;
   Int_t           t_Event;
   Int_t           t_DataType;
   Int_t           t_ieta;
   Int_t           t_iphi;
   Double_t        t_EventWeight;
   Int_t           t_nVtx;
   Int_t           t_nTrk;
   Int_t           t_goodPV;
   Double_t        t_l1pt;
   Double_t        t_l1eta;
   Double_t        t_l1phi;
   Double_t        t_l3pt;
   Double_t        t_l3eta;
   Double_t        t_l3phi;
   Double_t        t_p;
   Double_t        t_pt;
   Double_t        t_phi;
   Double_t        t_mindR1;
   Double_t        t_mindR2;
   Double_t        t_eMipDR;
   Double_t        t_eHcal;
   Double_t        t_eHcal10;
   Double_t        t_eHcal30;
   Double_t        t_hmaxNearP;
   Double_t        t_rhoh;
   Bool_t          t_selectTk;
   Bool_t          t_qltyFlag;
   Bool_t          t_qltyMissFlag;
   Bool_t          t_qltyPVFlag;
   Double_t        t_gentrackP;
   vector<unsigned int> *t_DetIds;
   vector<double>  *t_HitEnergies;
   vector<bool>    *t_trgbits;
   vector<unsigned int> *t_DetIds1;
   vector<unsigned int> *t_DetIds3;
   vector<double>  *t_HitEnergies1;
   vector<double>  *t_HitEnergies3;
  
  // List of branches
   TBranch        *b_t_Run;   //!
   TBranch        *b_t_Event;   //!
   TBranch        *b_t_DataType;   //!
   TBranch        *b_t_ieta;   //!
   TBranch        *b_t_iphi;   //!
   TBranch        *b_t_EventWeight;   //!
   TBranch        *b_t_nVtx;   //!
   TBranch        *b_t_nTrk;   //!
   TBranch        *b_t_goodPV;   //!
   TBranch        *b_t_l1pt;   //!
   TBranch        *b_t_l1eta;   //!
   TBranch        *b_t_l1phi;   //!
   TBranch        *b_t_l3pt;   //!
   TBranch        *b_t_l3eta;   //!
   TBranch        *b_t_l3phi;   //!
   TBranch        *b_t_p;   //!
   TBranch        *b_t_pt;   //!
   TBranch        *b_t_phi;   //!
   TBranch        *b_t_mindR1;   //!
   TBranch        *b_t_mindR2;   //!
   TBranch        *b_t_eMipDR;   //!
   TBranch        *b_t_eHcal;   //!
   TBranch        *b_t_eHcal10;   //!
   TBranch        *b_t_eHcal30;   //!
   TBranch        *b_t_hmaxNearP;   //!
   TBranch        *b_t_rhoh;   //!
   TBranch        *b_t_selectTk;   //!
   TBranch        *b_t_qltyFlag;   //!
   TBranch        *b_t_qltyMissFlag;   //!
   TBranch        *b_t_qltyPVFlag;   //!
   TBranch        *b_t_gentrackP;   //!
   TBranch        *b_t_DetIds;   //!
   TBranch        *b_t_HitEnergies;   //!
   TBranch        *b_t_trgbits;   //!
   TBranch        *b_t_DetIds1;   //!
   TBranch        *b_t_DetIds3;   //!
   TBranch        *b_t_HitEnergies1;   //!
   TBranch        *b_t_HitEnergies3;   //!

  //--- constructor & destructor
  //AlcaCalibTree(TTree *tree=0);
  AlcaCalibTree(TChain *tree);
  virtual ~AlcaCalibTree();
  
  //--- functions
  virtual Int_t      GetEntry(Long64_t entry);
  virtual Long64_t   LoadTree(Long64_t entry);
  //virtual void     Init(TTree *tree);
  virtual void       Init(TChain *tree);
  virtual Bool_t     Notify();  
  };

//**********************************************************
// AlcaCalibTree constructor
//**********************************************************
//AlcaCalibTree::AlcaCalibTree(TTree *tree) : fChain(0) {
AlcaCalibTree::AlcaCalibTree(TChain *tree)
{ //: fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("output.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("HcalIsoTrkAnalyzer");
    dir->GetObject("CalibTree",tree);
  }
  Init(tree);
}

//**********************************************************
// AlcaCalibTree destructor
//**********************************************************
AlcaCalibTree::~AlcaCalibTree() {

  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

//**********************************************************
// Get entry function
//**********************************************************
Int_t AlcaCalibTree::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

//**********************************************************
// Load tree function
//**********************************************************
Long64_t AlcaCalibTree::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

//**********************************************************
// Initialisation of TTree
//**********************************************************
void AlcaCalibTree::Init(TChain *tree) {
  // Set object pointer
   t_DetIds = 0;
   t_HitEnergies = 0;
   t_trgbits = 0;
   t_DetIds1 = 0;
   t_DetIds3 = 0;
   t_HitEnergies1 = 0;
   t_HitEnergies3 = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
   fChain->SetBranchAddress("t_Run", &t_Run, &b_t_Run);
   fChain->SetBranchAddress("t_Event", &t_Event, &b_t_Event);
   fChain->SetBranchAddress("t_DataType", &t_DataType, &b_t_DataType);
   fChain->SetBranchAddress("t_ieta", &t_ieta, &b_t_ieta);
   fChain->SetBranchAddress("t_iphi", &t_iphi, &b_t_iphi);
   fChain->SetBranchAddress("t_EventWeight", &t_EventWeight, &b_t_EventWeight);
   fChain->SetBranchAddress("t_nVtx", &t_nVtx, &b_t_nVtx);
   fChain->SetBranchAddress("t_nTrk", &t_nTrk, &b_t_nTrk);
   fChain->SetBranchAddress("t_goodPV", &t_goodPV, &b_t_goodPV);
   fChain->SetBranchAddress("t_l1pt", &t_l1pt, &b_t_l1pt);
   fChain->SetBranchAddress("t_l1eta", &t_l1eta, &b_t_l1eta);
   fChain->SetBranchAddress("t_l1phi", &t_l1phi, &b_t_l1phi);
   fChain->SetBranchAddress("t_l3pt", &t_l3pt, &b_t_l3pt);
   fChain->SetBranchAddress("t_l3eta", &t_l3eta, &b_t_l3eta);
   fChain->SetBranchAddress("t_l3phi", &t_l3phi, &b_t_l3phi);
   fChain->SetBranchAddress("t_p", &t_p, &b_t_p);
   fChain->SetBranchAddress("t_pt", &t_pt, &b_t_pt);
   fChain->SetBranchAddress("t_phi", &t_phi, &b_t_phi);
   fChain->SetBranchAddress("t_mindR1", &t_mindR1, &b_t_mindR1);
   fChain->SetBranchAddress("t_mindR2", &t_mindR2, &b_t_mindR2);
   fChain->SetBranchAddress("t_eMipDR", &t_eMipDR, &b_t_eMipDR);
   fChain->SetBranchAddress("t_eHcal", &t_eHcal, &b_t_eHcal);
   fChain->SetBranchAddress("t_eHcal10", &t_eHcal10, &b_t_eHcal10);
   fChain->SetBranchAddress("t_eHcal30", &t_eHcal30, &b_t_eHcal30);
   fChain->SetBranchAddress("t_hmaxNearP", &t_hmaxNearP, &b_t_hmaxNearP);
   fChain->SetBranchAddress("t_rhoh", &t_rhoh, &b_t_rhoh);
   fChain->SetBranchAddress("t_selectTk", &t_selectTk, &b_t_selectTk);
   fChain->SetBranchAddress("t_qltyFlag", &t_qltyFlag, &b_t_qltyFlag);
   fChain->SetBranchAddress("t_qltyMissFlag", &t_qltyMissFlag, &b_t_qltyMissFlag);
   fChain->SetBranchAddress("t_qltyPVFlag", &t_qltyPVFlag, &b_t_qltyPVFlag);
   fChain->SetBranchAddress("t_gentrackP", &t_gentrackP, &b_t_gentrackP);
   fChain->SetBranchAddress("t_DetIds", &t_DetIds, &b_t_DetIds);
   fChain->SetBranchAddress("t_HitEnergies", &t_HitEnergies, &b_t_HitEnergies);
   fChain->SetBranchAddress("t_trgbits", &t_trgbits, &b_t_trgbits);
   fChain->SetBranchAddress("t_DetIds1", &t_DetIds1, &b_t_DetIds1);
   fChain->SetBranchAddress("t_DetIds3", &t_DetIds3, &b_t_DetIds3);
   fChain->SetBranchAddress("t_HitEnergies1", &t_HitEnergies1, &b_t_HitEnergies1);
   fChain->SetBranchAddress("t_HitEnergies3", &t_HitEnergies3, &b_t_HitEnergies3);
  Notify();
}

//**********************************************************
// Notification when opening new file
//**********************************************************
Bool_t AlcaCalibTree::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}


