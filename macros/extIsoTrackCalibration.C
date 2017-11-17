//////////////////////////////////////////////////////////
// Extended L3 iterative procedure 
// for IsoTrack calibration
// requires header file CalibTree.hh
// (based on CalibTree.C from CMSSW_7_4)
// CalibTree class contains ROOT-tree
// generated with IsoTrackCalibration plugin
//
// Version 1.0, October 2017
// based on version 6.2 of L3_IsoTrackCalibration.C
// with added option to inter(extra)polate factors
// to uncalibrated cells from neighbor calibrated cells  
//
// Version 2.0, October 2017
// add HEP, HEM hists for final 
//
// Author: M. Chadeeva
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
//#include <TProfile.h>
#include <TLegend.h>
#include <TString.h>
#include <TF1.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>

//**********************************************************
// Constants
//**********************************************************
const int NVERTEX_MIN = 1;
const int NVERTEX_MAX = 100;
const int MAX_CALIBRATED_IETA = 23;
const int MAX_TRACK_IETA = 23;
const int MIN_RUNNUM = 297494; // min run number 
const int MAX_RUNNUM = 900000; // max run number 
//const char *l3prefix5 = "dat2017";

const unsigned int MAXNUM_SUBDET = 100;
const bool SINGLE_REFERENCE_RESPONSE = false;
const int FIRST_IETA_TR = 15;
const int FIRST_IETA_HE = 21;
const unsigned int N_DEPTHS = 5;

// for SiPM studies
const int N_MODULE_GROUPS = 18;
const int ABSIETA_MIN = 17;
const int ABSIETA_MAX = 23;
const double ZERO_ANGLE_IN_DEGREE = 10;

/*
// old detID format for CMSSW_7_X 
const unsigned int PHI_MASK       = 0x7F;
const unsigned int ETA_OFFSET     = 7;
const unsigned int ETA_MASK       = 0x3F;
const unsigned int ZSIDE_MASK     = 0x2000;
const unsigned int DEPTH_OFFSET   = 14;
const unsigned int DEPTH_MASK     = 0x1F;
const unsigned int DEPTH_SET      = 0x1C000;
unsigned int MASK(0xFF80); //merge phi 
*/
// new detID format for CMSSW_8_X
const unsigned int PHI_MASK       = 0x3FF;
const unsigned int ETA_OFFSET     = 10;
const unsigned int ETA_MASK       = 0x1FF;
const unsigned int ZSIDE_MASK     = 0x80000;
const unsigned int DEPTH_OFFSET   = 20;
const unsigned int DEPTH_MASK     = 0xF;
const unsigned int DEPTH_SET      = 0xF00000;

//--------------------------------------------
// merging depths

const unsigned int MERGE_PHI_AND_DEPTHS = 1;
const unsigned int MASK(0xFFC00); //merge phi and depth
/*
const unsigned int MERGE_PHI_AND_DEPTHS = 0;
const unsigned int MASK(0xFFFC00); //merge phi
*/
//-------------------------------------------
//-------------------------------------------
// individual ieta rings
const unsigned int MASK2(0); // no second mask
const int N_ETA_RINGS_PER_BIN = 1;
/*
// twin (even+odd) ieta rings
const unsigned int MASK2(0x80);
const int N_ETA_RINGS_PER_BIN = 2;

// 4-fold ieta rings
const unsigned int MASK2(0x180);
const int N_ETA_RINGS_PER_BIN = 4;
*/
//-------------------------------------------
const int MAX_ONESIDE_ETA_RINGS = 30;
const int HALF_NUM_ETA_BINS =
  (MAX_ONESIDE_ETA_RINGS + 1*(N_ETA_RINGS_PER_BIN>1))/N_ETA_RINGS_PER_BIN;
const int NUM_ETA_BINS = 2*HALF_NUM_ETA_BINS + 1;

const int MIN_N_TRACKS_PER_CELL = 50;
const int MIN_N_ENTRIES_FOR_FIT = 150;
const double MAX_REL_UNC_FACTOR = 0.2;
const unsigned UNCALIB_FACT_TO_ONE = 0; // 0: extrapolate behind MAX_CALIBRATED_IETA
                                        // 1: set factor=1 for uncalibrated cells

const bool APPLY_CORRECTION_FOR_PU = true;
/*
//--- base correction for PU as of 2015 MC studies
const char *l3prefix1 = "base";
const double DELTA_CUT = 0.0; 
const double LINEAR_COR_COEF[5] = { -0.375, -0.375, -0.375, -0.375, -0.375 };
const double SQUARE_COR_COEF[5] = { -0.450, -0.450, -0.450, -0.450, -0.450 };
const int FIRST_IETA_FWD_1 = 25;
const int FIRST_IETA_FWD_2 = 26;
//const double UPPER_LIMIT_DELTA_PU_COR = 2.0;
*/
//--- optimized correction for PU as of 2016 MC studies
// tuned to get MPV after correction
// at the level of that of single pion response w/o PU and w/o correction
const char *l3prefix1 = "opt";
const double DELTA_CUT = 0.02; 
const double LINEAR_COR_COEF[5] = { -0.35, -0.35, -0.35, -0.35, -0.45 };
const double SQUARE_COR_COEF[5] = { -0.65, -0.65, -0.65, -0.30, -0.10 };
const int FIRST_IETA_FWD_1 = 25;
const int FIRST_IETA_FWD_2 = 26;

//const double UPPER_LIMIT_DELTA_PU_COR = 1.5;

const double UPPER_LIMIT_RESPONSE_BEFORE_COR = 3.0;
const double FLEX_SEL_FIRST_CONST = 20.0;  // 20.(for exp); or 16*2
const double FLEX_SEL_SECOND_CONST = 18.0; // 18.(for exp);

const double MIN_RESPONSE_HIST = 0.0;
const double MAX_RESPONSE_HIST = UPPER_LIMIT_RESPONSE_BEFORE_COR;
const int NBIN_RESPONSE_HIST = 480;
const int NBIN_RESPONSE_HIST_IND = 120;
const double FIT_RMS_INTERVAL = 1.5;
const double RESOLUTION_HCAL = 0.3;

//std::cout.precision(3);

//**********************************************************
// Header with CalibTree class definition
//**********************************************************

#include "CalibTree.hh"

//**********************************************************
// Description of function to run iteration
//**********************************************************
unsigned int runIterations(const char *inFileDir = ".", 
			   const char *inFileNamePrefix = "outputFromAnalyzer",
			   const int firstInputFileEnum = 0,
			   const int lastInputFileEnum = 1,
			   const bool preselectedSample = false,
			   const unsigned maxNumberOfIterations = 1,
			   const char *l3prefix5 = "Dat17",
			   const double minHcalEnergy = 10.0,
			   const double minPt = 7.0,
			   const double limitForChargeIsolation = 10.0,
			   const double minTrackMomentum = 40.0,
			   const double maxTrackMomentum = 60.0,
			   const double limitForMipInEcal = 1.0,
			   const bool shiftResponse = 1,
			   const unsigned int subSample = 2,
			   const bool isCrosscheck = false,
			   const char *inTxtFilePrefix = "test",
			   const char *treeDirName = "IsoTrackCalibration", 
			   const char *treeName = "CalibTree",
			   unsigned int Debug = 0)
{
  // Debug:  0-no debugging; 1-short debug; >1 - number of events to be shown in detail
  // subSample: extract factors from odd (0), even(1) or all(2) events
  // limitForChargeIsolation:  <0 - flex. sel.
  // and corr. for PU
  
  if ( isCrosscheck )
    std::cout << "Test with previously extracted factors..." << std::endl;
  else
    std::cout << "Extracting factors using L3 algorithm and isolated tracks..." << std::endl;

  char l3prefix0[10] = "_ref1";
  if ( !SINGLE_REFERENCE_RESPONSE ) sprintf(l3prefix0,"_ref3");
  char l3prefix2[20] = "_noCor";
  if ( APPLY_CORRECTION_FOR_PU ) sprintf(l3prefix2,"_%s%02d",
					 l3prefix1, int(100*DELTA_CUT)
					 );				
  char l3prefix3[10] = "_mean";
  if ( shiftResponse ) sprintf(l3prefix3,"_mpv");
  char l3prefix4[8] = "_3dep";
  if ( MERGE_PHI_AND_DEPTHS ) sprintf(l3prefix4,"_merged");
  char l3prefix[100];
  sprintf(l3prefix,"%s%1d_vtx%03d-%03d_cal%02d_tr%02d%s%s%s%s",
	  l3prefix5,
	  UNCALIB_FACT_TO_ONE,
	  NVERTEX_MIN,
	  NVERTEX_MAX,
	  MAX_CALIBRATED_IETA,
	  MAX_TRACK_IETA,
	  l3prefix0, l3prefix2, l3prefix3, l3prefix4);
    
  char fnameInput[120];
  char fnameOutRoot[120];
  char fnameOutTxt[120] = "dummy";
  char fnameInTxt[120]  = "dummy";
  char tname[100];
  
  TGraph *g_converge1 = new TGraph(maxNumberOfIterations);
  TGraph *g_converge2 = new TGraph(maxNumberOfIterations);
  TGraph *g_converge3 = new TGraph(maxNumberOfIterations);

  if ( preselectedSample ) sprintf(tname, "%s", treeName );
  else sprintf(tname, "%s/%s", treeDirName, treeName );
  
  TChain tree(tname);

  //--- combine tree from several enumerated files with the same prefix
  //    or one file w/o number (firstInputFileEnum = lastInputFileEnum < 0 )

  for ( int ik = firstInputFileEnum; ik <= lastInputFileEnum; ik++ ) {
    if ( ik < 0 ) 
      sprintf(fnameInput, "%s/%s.root", inFileDir, inFileNamePrefix);
    else if (ik < 10 )
      sprintf(fnameInput, "%s/%s_%1d.root", inFileDir, inFileNamePrefix, ik);
    else if (ik < 100 )
      sprintf(fnameInput, "%s/%s_%2d.root", inFileDir, inFileNamePrefix, ik);
    else if (ik < 1000 )
      sprintf(fnameInput, "%s/%s_%3d.root", inFileDir, inFileNamePrefix, ik);
    else
      sprintf(fnameInput, "%s/%s_%4d.root", inFileDir, inFileNamePrefix, ik);

    if ( !gSystem->Which("./", fnameInput ) ) { // check file availability
      std::cout << "File " << fnameInput << " doesn't exist." << std::endl;
    }
    else {
      tree.Add(fnameInput);
      std::cout << "Add tree from " << fnameInput 
	        << "   total number of entries (tracks): "
		<< tree.GetEntries() << std::endl;
    }
  }
  if ( tree.GetEntries() == 0 ) {
    std:: cout << "Tree is empty." << std::endl;
    return -2;
  }

  //--- Initialize tree
  CalibTree t(&tree,
	      minHcalEnergy, minPt,
	      limitForMipInEcal, limitForChargeIsolation,
	      minTrackMomentum, maxTrackMomentum);
  
  char isoPrefix[14];
  if ( limitForChargeIsolation < 0 )
    sprintf(isoPrefix, "flex%02d-%02d-%02d",
	    int(t.limCharIso),
	    int(abs(FLEX_SEL_FIRST_CONST)),
	    int(abs(FLEX_SEL_SECOND_CONST))
	    );
  else
    sprintf(isoPrefix, "const%03d",
	    int(t.limCharIso*10)
	    );

    
  //--- Define files
  if ( isCrosscheck ) {
    sprintf(fnameInTxt, "%s_%s_p%02d-%02d_pt%02d_eh%02d_ee%1d_step%1d.txt",
	    inTxtFilePrefix,
	    isoPrefix, 
	    int(minTrackMomentum), int(maxTrackMomentum), int(minPt),
	    int(minHcalEnergy), int(limitForMipInEcal),
	    N_ETA_RINGS_PER_BIN
	    );
    sprintf(fnameOutRoot,
	    "test_%1d_%s_by_%s_%s_p%02d-%02d_pt%02d_eh%02d_ee%1d_step%1d.root",
	    subSample,
	    inFileNamePrefix,
	    inTxtFilePrefix,
	    isoPrefix,
	    int(minTrackMomentum), int(maxTrackMomentum), int(minPt),
	    int(minHcalEnergy), int(limitForMipInEcal),
	    N_ETA_RINGS_PER_BIN
	    );
  }
  else {    
    sprintf(fnameOutTxt,
	    "%s_%1d_%s_i%02d_%s_p%02d-%02d_pt%02d_eh%02d_ee%1d_step%1d.txt",
	    l3prefix,
	    subSample,
	    inFileNamePrefix,
	    maxNumberOfIterations,
	    isoPrefix,
	    int(minTrackMomentum), int(maxTrackMomentum), int(minPt),
	    int(minHcalEnergy), int(limitForMipInEcal),
	    N_ETA_RINGS_PER_BIN
	    );
    sprintf(fnameOutRoot,
	    "%s_%1d_%s_i%02d_%s_p%02d-%02d_pt%02d_eh%02d_ee%1d_step%1d.root",
	    l3prefix,
	    subSample,
	    inFileNamePrefix,
	    maxNumberOfIterations,
	    isoPrefix,
	    int(minTrackMomentum), int(maxTrackMomentum), int(minPt),
	    int(minHcalEnergy), int(limitForMipInEcal),
	    N_ETA_RINGS_PER_BIN
	    );
  }  
  if ( !t.openOutputRootFile(fnameOutRoot) ) {
    std::cout << "Problems with booking output file " << fnameOutRoot << std::endl;
    return -1;
  }
  std::cout << "Correction for PU: ";
  if ( APPLY_CORRECTION_FOR_PU ) {
    std::cout << " applied for delta > " << DELTA_CUT << std::endl;
    std::cout << " for HB(ieta<" << FIRST_IETA_TR << "): "<< LINEAR_COR_COEF[0]
	      << " ; " << SQUARE_COR_COEF[0]
	      << std::endl;
    std::cout << " for TR(ieta<" << FIRST_IETA_HE << "): " << LINEAR_COR_COEF[1]
	      << " ; " << SQUARE_COR_COEF[1]
	      << std::endl;
    std::cout << " for HE(ieta<" << FIRST_IETA_FWD_1 << "): " << LINEAR_COR_COEF[2]
	      << " ; " << SQUARE_COR_COEF[2]
	      << std::endl;
    std::cout << " for HE(ieta<" << FIRST_IETA_FWD_2 << "): " << LINEAR_COR_COEF[3]
	      << " ; " << SQUARE_COR_COEF[3]
	      << std::endl;
    std::cout << " for HE(ieta>=" << FIRST_IETA_FWD_2 << "): " << LINEAR_COR_COEF[4]
	      << " ; " << SQUARE_COR_COEF[4]
	      << std::endl;
  }
  else
    std::cout << " no " << std::endl;
    
  /*
  std::cout << "Constant coefficient from charge isolation: "
	    << t.constForFlexSel << std::endl; 
  */

  unsigned int numOfSavedFactors(0);
  int nEventsWithGoodTrack(0);
  double MPVfromLastFit(0);
  
  if ( isCrosscheck ) {
    // open txt file and fill map with factors
    if ( t.getFactorsFromFile(fnameInTxt, Debug) ) {
      nEventsWithGoodTrack = t.firstLoop(subSample, false, true, Debug);
      std::cout << "Number of events with good track = "
		<< nEventsWithGoodTrack << std::endl;
      MPVfromLastFit = t.lastLoop(subSample, maxNumberOfIterations, true, Debug);
      std::cout << "Finish testing " << t.factors.size() << " factors from file "
		<< fnameInTxt << std::endl;
      std::cout << "MPV from fit after last iteration = "
		<< MPVfromLastFit << std::endl;
      std::cout << "Test plots saved in " << fnameOutRoot << std::endl;
    }
    else {
      std::cout << "File " << fnameInTxt << " doesn't exist." << std::endl;
    }
  }
  else {
    //--- Prepare initial histograms and count good track
    nEventsWithGoodTrack = t.firstLoop(subSample, shiftResponse, false, Debug);
    std::cout << "Number of events with good track = "
	      << nEventsWithGoodTrack << std::endl;
    //--- Iterate
    for ( unsigned int k = 0; k < maxNumberOfIterations; ++k ) {
      g_converge1->SetPoint( k, k+1, t.loopForIteration(subSample, k+1, Debug) );
      g_converge2->SetPoint( k, k+1, t.maxZtestFromWeights );
      g_converge3->SetPoint( k, k+1, t.maxSys2StatRatio );
    }
    //--- Finish
    MPVfromLastFit = t.lastLoop(subSample, maxNumberOfIterations, false, Debug);
    numOfSavedFactors = t.saveFactorsInFile(fnameOutTxt);

    sprintf(tname,"Mean deviation for subdetectors with Ntrack>%d",
	    MIN_N_TRACKS_PER_CELL);
    g_converge1->SetTitle(tname);
    g_converge1->GetXaxis()->SetTitle("iteration");
    t.foutRootFile->WriteTObject(g_converge1, "g_cvgD");
    sprintf(tname,"Max abs(Z-test) for factors");
    g_converge2->SetTitle(tname);
    g_converge2->GetXaxis()->SetTitle("iteration");
    t.foutRootFile->WriteTObject(g_converge2, "g_cvgW");
    sprintf(tname,"Max ratio of syst. to stat. uncertainty");
    g_converge3->SetTitle(tname);
    g_converge3->GetXaxis()->SetTitle("iteration");
    t.foutRootFile->WriteTObject(g_converge3, "g_cvgR");
    
    std::cout << "Finish adjusting factors after "
	      << maxNumberOfIterations << " iterations" << std::endl;
    std::cout << "MPV from fit after last iteration = "
	      << MPVfromLastFit << std::endl;
    std::cout << "Table with " << numOfSavedFactors << " factors"
	      << " with more than " << MIN_N_TRACKS_PER_CELL << " tracks/subdetector"
	      << " (from " << t.factors.size() << " available)"
	      << " is written in file " << fnameOutTxt << std::endl;
    std::cout << "Plots saved in " << fnameOutRoot << std::endl;
  }

  return numOfSavedFactors;
}

//**********************************************************
// Initial loop over events in the tree
//**********************************************************
Int_t CalibTree::firstLoop(unsigned int subsample,
			   bool shiftResp,
			   bool istest,
			   unsigned int debug)
{
  char name[100];
  unsigned int ndebug(0);
  double maxRespForGoodTrack(0);
  double minRespForGoodTrack(1000);
  int nRespOverHistLimit(0);
  int minRunNumber(1000000);
  int maxRunNumber(0);
  
  int ntrk_ieta[NUM_ETA_BINS];
  for ( int j = 0; j < NUM_ETA_BINS; j++ ) {
    ntrk_ieta[j] = 0;
  }
  
  char scorr[80] = "correction for PU";
  char sxlabel[80] ="E^{cor}_{hcal}/(p_{track} - E_{ecal})"; 
  if ( !APPLY_CORRECTION_FOR_PU ) {
    sprintf(scorr,"no correction for PU");
    sprintf(sxlabel,"E_{hcal}/(p_{track} - E_{ecal})");
  }
  
  TF1* f1 = new TF1("f1","gaus", MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  TH1F* e2p[NUM_ETA_BINS];

  h_nvtx = new TH1D("h_nvtx","Number of vertices in selected events",100,0,100);
  h_nvtx->GetXaxis()->SetTitle("N_{vtx}");
  h_cluster = new TProfile("h_cluster","Number of subdetectors in cluster",
			   2*MAX_ONESIDE_ETA_RINGS,
			   -MAX_ONESIDE_ETA_RINGS, MAX_ONESIDE_ETA_RINGS); 
  h_cluster->GetXaxis()->SetTitle("i#eta of track");
  h_cluster->GetYaxis()->SetTitle("<N_{subdet}>");

  //--------- initialize histograms for response -----------------------------------------
  sprintf(name,"Initial HB+HE: %s", scorr);
  e2p_init = new TH1F("e2p_init", name,
		      NBIN_RESPONSE_HIST, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2p_init->Sumw2();
  e2p_init->GetXaxis()->SetTitle(sxlabel);
  
  sprintf(name,"Initial HB: %s", scorr);
  e2pHB_init = new TH1F("e2pHB_init", name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHB_init->Sumw2();
  e2pHB_init->GetXaxis()->SetTitle(sxlabel);

  sprintf(name,"Initial TR: %s", scorr);
  e2pTR_init = new TH1F("e2pTR_init", name,
			NBIN_RESPONSE_HIST/10, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pTR_init->Sumw2();
  e2pTR_init->GetXaxis()->SetTitle(sxlabel);

  sprintf(name,"Initial HE: %s", scorr);
  e2pHE_init = new TH1F("e2pHE_init", name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHE_init->Sumw2();
  e2pHE_init->GetXaxis()->SetTitle(sxlabel);

  //-- for SiPM module studies
  char snm[20];
  sprintf(snm,"e2pHEP[00]");
  sprintf(name,"Initial HEP, iphi 1,2,71,72");
  e2pHEP[0] = new TH1F(snm, name,
		       NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHEP[0]->Sumw2();
  e2pHEP[0]->GetXaxis()->SetTitle(sxlabel);

  sprintf(snm,"e2pHEM[00]");
  sprintf(name,"Initial HEM, iphi 1,2,71,72");
  e2pHEM[0] = new TH1F(snm, name,
		       NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHEM[0]->Sumw2();
  e2pHEM[0]->GetXaxis()->SetTitle(sxlabel);
  
  for ( int imd = 1; imd < N_MODULE_GROUPS; imd++ ) {
    int imin = 3 + 4*(imd - 1); //(10 + 20*imd)*TMath::Pi()/180; 
    int imax = imin + 3; //20*TMath::Pi()/180; 
    sprintf(snm,"e2pHEP[%02d]",imd);
    sprintf(name,"Initial HEP, iphi %d-%d", imin, imax);
    e2pHEP[imd] = new TH1F(snm, name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
    e2pHEP[imd]->Sumw2();
    e2pHEP[imd]->GetXaxis()->SetTitle(sxlabel);

    sprintf(snm,"e2pHEM[%02d]",imd);
    sprintf(name,"Initial HEM, iphi %d-%d", imin, imax);
    e2pHEM[imd] = new TH1F(snm, name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
    e2pHEM[imd]->Sumw2();
    e2pHEM[imd]->GetXaxis()->SetTitle(sxlabel);
  }
  
//--- initialize chain ----------------------------------------
  if (fChain == 0) return 0;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nb = 0;
  
  int nSelectedEvents(0);

  if ( debug > 0 ) { 
    std::cout << "---------- First loop -------------------------- " << std::endl;
  }
// ----------------------- loop over events -------------------------------------  
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if ( ientry < 0 || ndebug > debug ) break;   
    nb = fChain->GetEntry(jentry);   //nbytes += nb;
    
    if ( (jentry%2 == subsample) ) continue;   // only odd or even events
 
// --------------- selection of good track --------------------    

    if ( t_Run < minRunNumber ) minRunNumber = t_Run;
    if ( t_Run > maxRunNumber ) maxRunNumber = t_Run;
    
    if ( t_nVtx < NVERTEX_MIN  ||  t_nVtx > NVERTEX_MAX 
	 || t_Run < MIN_RUNNUM  || t_Run > MAX_RUNNUM ) continue;
    if ( t_ieta < -MAX_TRACK_IETA || t_ieta > MAX_TRACK_IETA ) continue;
    
    if ( !goodTrack() ) continue;
    
    h_nvtx->Fill(t_nVtx,1.0);

    if ( debug > 1 ) {
      ndebug++;
      std::cout << "***Entry (Track) Number : " << ientry << "(" << jentry << ")"
		<< " p/eHCal/eMipDR/nDets : " << t_p << "/" << t_eHcal 
		<< "/" << t_eMipDR << "/" << (*t_DetIds).size() 
		<< std::endl;
    }
    
    double eTotal(0.0);
    //double eTotalWithEcal(0.0);
      
    // ---- loop over active subdetectors in the event for total energy ---
    unsigned int nDets = (*t_DetIds).size();
    h_cluster->Fill(t_ieta, nDets);
      
    for (unsigned int idet = 0; idet < nDets; idet++) { 
      eTotal += (*t_HitEnergies)[idet];
    }
    //eTotalWithEcal = eTotal + t_eMipDR;    

// --- Correction for PU  --------
    double eTotalCor(eTotal);
    //double eTotalWithEcalCor(eTotalWithEcal);
    double correctionForPU(1.0);
    double de2p(0.0);
    int abs_t_ieta = abs(t_ieta);
    
    if ( APPLY_CORRECTION_FOR_PU ) { 
      de2p = (t_eHcal30 - t_eHcal10)/t_p;

      if ( de2p > DELTA_CUT ) {
	int icor = int(abs_t_ieta >= FIRST_IETA_TR) + int(abs_t_ieta >= FIRST_IETA_HE)
	  + int(abs_t_ieta >= FIRST_IETA_FWD_1) + int(abs_t_ieta >= FIRST_IETA_FWD_2);
	correctionForPU = (1 + LINEAR_COR_COEF[icor]*(t_eHcal/t_p)*de2p
			   *(1 + SQUARE_COR_COEF[icor]*de2p));
      }
    }    

    // check for possibility to correct for PU
    if ( correctionForPU <= 0 || correctionForPU > 1 ) continue;
    nSelectedEvents++;

    eTotalCor = eTotal*correctionForPU;
    //eTotalWithEcalCor = eTotalCor + t_eMipDR;

    //double response = eTotalWithEcalCor/t_p;
    double response = eTotalCor/(t_p - t_eMipDR);

    std::map<unsigned int, bool> sameSubdet;
    sameSubdet.clear();
    double resp2 = response*response;

    for (unsigned int idet = 0; idet < nDets; idet++) { 
      unsigned int detId = ( (*t_DetIds)[idet] & MASK ) | MASK2 ;

      if ( debug > 1 ) {
	unsigned int detId0 = ( (*t_DetIds)[idet] & MASK ) ;
	std::cout << "jentry/idet/detId :: ieta/z/depth ::: "
		  << std::dec
		  << jentry << " / "
		  << ((*t_DetIds)[idet]) << " / "
		  << detId0 << "(" << detId << ")" << " :: "
		  << ((detId0>>ETA_OFFSET) & ETA_MASK)
		  << "(" << ((detId>>ETA_OFFSET) & ETA_MASK) << ")" << " / "
		  << ((detId0&ZSIDE_MASK) ? 1 : -1)
		  << "(" << ((detId&ZSIDE_MASK) ? 1 : -1) << ")" << " / "
		  << ((detId0>>DEPTH_OFFSET)&DEPTH_MASK)
		  << "(" << ((detId>>DEPTH_OFFSET)&DEPTH_MASK) << ")"
		  << std::endl;
      }
      if (nPhiMergedInEvent.find(detId) != nPhiMergedInEvent.end()) 
	nPhiMergedInEvent[detId]++;
      else 
	nPhiMergedInEvent.insert(std::pair<unsigned int,int>(detId, 1));
		
      if (nTrks.find(detId) != nTrks.end()) {
	if ( sameSubdet.find(detId) == sameSubdet.end() ) {
	  nTrks[detId]++;
	  nSubdetInEvent[detId] += nDets;
	  sumOfResponse[detId] += response;
	  sumOfResponseSquared[detId] += resp2;
	  sameSubdet.insert(std::pair<unsigned int,bool>(detId, true));
	}
      }
      else {
	nTrks.insert(std::pair<unsigned int,int>(detId, 1));
	nSubdetInEvent.insert(std::pair<unsigned int,int>(detId, nDets));
	sumOfResponse.insert(std::pair<unsigned int,double>(detId,response));
	sumOfResponseSquared.insert(std::pair<unsigned int,double>(detId,resp2));
	sameSubdet.insert(std::pair<unsigned int,bool>(detId, true));
	subDetector_trk.insert(std::pair<unsigned int,
			       int>( detId,((*t_DetIds)[idet] &0xe000000) / 0x2000000 ));
      }
      
    }

// --- Fill initial histograms ---------------------------      
    e2p_init->Fill(response ,1.0);

    if ( abs_t_ieta < FIRST_IETA_TR )
      e2pHB_init->Fill(response ,1.0);
    else if ( abs_t_ieta < FIRST_IETA_HE )
      e2pTR_init->Fill(response ,1.0);
    else
      e2pHE_init->Fill(response ,1.0);

    // --- For SiPM module crosscheck ----
    double abs_phi = abs(t_phi)*360/TMath::TwoPi();
    if ( abs_t_ieta <= ABSIETA_MAX && abs_t_ieta >= ABSIETA_MIN
	 && abs_phi <= ZERO_ANGLE_IN_DEGREE ) {
	if ( t_ieta < 0 ) e2pHEM[0]->Fill(response ,1.0);
	else e2pHEP[0]->Fill(response ,1.0);
      }
    double phi = 72*(t_phi + (t_phi<0)*TMath::TwoPi())/TMath::TwoPi();    
    for ( int imd = 1; imd < N_MODULE_GROUPS; imd++ ) {
      int imin = 3 + 4*(imd - 1); 
      int imax = imin + 3;
      if ( abs_t_ieta <= ABSIETA_MAX && abs_t_ieta >= ABSIETA_MIN
	   && phi <= imax && phi >= imin ) {
	if ( t_ieta < 0 ) e2pHEM[imd]->Fill(response ,1.0);
	else e2pHEP[imd]->Fill(response ,1.0);
      }
    }

    if ( debug > 1 ) {
      std::cout << "***Entry : " << ientry
		<< " ***ieta/p/Ecal/nDet : "
		<< t_ieta << "/" << t_p
		<< "/" << t_eMipDR << "/" << (*t_DetIds).size() 
		<< " ***Etot/E10/E30/Ecor/cPU : " << t_eHcal
		<< "/" << t_eHcal10 << "/" << t_eHcal30
		<< "/" << eTotalCor << "/" << correctionForPU
		<< "(" << de2p << ")"
		<< std::endl;
    }
    if ( response > maxRespForGoodTrack  )
      maxRespForGoodTrack = response;
    if ( response < minRespForGoodTrack )
      minRespForGoodTrack = response;
    if ( response > MAX_RESPONSE_HIST )
      nRespOverHistLimit++;

    int jj = HALF_NUM_ETA_BINS + int(t_ieta/N_ETA_RINGS_PER_BIN);
    ntrk_ieta[jj]++;
    
  } // ------------------- end of loop over events -------------------------------------

  for ( int j = 0; j < NUM_ETA_BINS; j++ ) {
    if ( maxNumOfTracksForIeta < ntrk_ieta[j] ) maxNumOfTracksForIeta = ntrk_ieta[j];
  }

//---------------------- Fill individual e2p[ieta] ----------------------
  int n_ieta_bins = 0;
  
  if ( istest ) {
    n_ieta_bins = 2*2.5*pow(maxNumOfTracksForIeta,1/3.0);
    for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
      sprintf(name,"e2p[%02d]", i);
      e2p[i] = new TH1F(name, "",
			n_ieta_bins, //NBIN_RESPONSE_HIST_IND,
			MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
      e2p[i]->Sumw2();
    }
// ----------------------- second loop over events -------------------------------------
// --------- to fill e2p[i] ---------------------------------
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if ( ientry < 0 || ndebug > debug ) break;   
      nb = fChain->GetEntry(jentry);   //nbytes += nb;
    
      if ( (jentry%2 == subsample) ) continue;   // only odd or even events
 
// --------------- selection of good track --------------------    

      if ( t_nVtx < NVERTEX_MIN  ||  t_nVtx > NVERTEX_MAX 
	   || t_Run < MIN_RUNNUM  || t_Run > MAX_RUNNUM ) continue;
      if ( t_ieta < -MAX_TRACK_IETA || t_ieta > MAX_TRACK_IETA ) continue;
    
      if ( !goodTrack() ) continue;

      double eTotal(0.0);
      
    // ---- loop over active subdetectors in the event for total energy ---
      unsigned int nDets = (*t_DetIds).size();
      
      for (unsigned int idet = 0; idet < nDets; idet++) { 
	eTotal += (*t_HitEnergies)[idet];
      }

// --- Correction for PU  --------
      double eTotalCor(eTotal);
      double correctionForPU(1.0);
      double de2p(0.0);
      int abs_t_ieta = abs(t_ieta);
    
      if ( APPLY_CORRECTION_FOR_PU ) { 
	de2p = (t_eHcal30 - t_eHcal10)/t_p;

	if ( de2p > DELTA_CUT ) {
	  int icor = int(abs_t_ieta >= FIRST_IETA_TR) + int(abs_t_ieta >= FIRST_IETA_HE)
	    + int(abs_t_ieta >= FIRST_IETA_FWD_1) + int(abs_t_ieta >= FIRST_IETA_FWD_2);
	  correctionForPU = (1 + LINEAR_COR_COEF[icor]*(t_eHcal/t_p)*de2p
			     *(1 + SQUARE_COR_COEF[icor]*de2p));
	}
      }

    // check for possibility to correct for PU
      if ( correctionForPU <= 0 || correctionForPU > 1 ) continue;

      eTotalCor = eTotal*correctionForPU;
      double response = eTotalCor/(t_p - t_eMipDR);

      int jeta = HALF_NUM_ETA_BINS + int(t_ieta/N_ETA_RINGS_PER_BIN);
      e2p[jeta]->Fill(response,1.0);
    
    } // ------------------- end of second loop over events
    for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
      int inum = (i - HALF_NUM_ETA_BINS)*N_ETA_RINGS_PER_BIN;
      sprintf(name,"e2pI_%d", inum);
      if ( e2p[i]->GetEntries() > 0 ) e2p[i]->Clone(name);
    }
  } //------ end istest

  double jeta[N_DEPTHS][MAXNUM_SUBDET];
  double nTrk[N_DEPTHS][MAXNUM_SUBDET];
  double nSub[N_DEPTHS][MAXNUM_SUBDET];
  double nPhi[N_DEPTHS][MAXNUM_SUBDET];
  double rms[N_DEPTHS][MAXNUM_SUBDET];
  unsigned int kdep[N_DEPTHS];
  for ( unsigned ik = 0; ik < N_DEPTHS; ik++ ) { kdep[ik] = 0; }

  // fill number of tracks
  std::map <unsigned int,int>::iterator nTrksItr = nTrks.begin();
  for (nTrksItr = nTrks.begin(); nTrksItr != nTrks.end(); nTrksItr++ ) {
    unsigned int detId = nTrksItr->first;
    int depth= ((detId>>DEPTH_OFFSET) & DEPTH_MASK) + int(MERGE_PHI_AND_DEPTHS);
    int zside= (detId&ZSIDE_MASK) ? 1 : -1;
    unsigned int kcur = kdep[depth-1];
    
    jeta[depth-1][kcur] = int((detId>>ETA_OFFSET) & ETA_MASK)*zside;
    nTrk[depth-1][kcur] = nTrksItr->second;
    nSub[depth-1][kcur] = double(nSubdetInEvent[detId])/double(nTrksItr->second);
    nPhi[depth-1][kcur] = double(nPhiMergedInEvent[detId])/double(nTrksItr->second);
    if ( nTrk[depth-1][kcur] > 1 ) 
      rms[depth-1][kcur] = sqrt((sumOfResponseSquared[detId] -
				 pow(sumOfResponse[detId],2)/nTrk[depth-1][kcur])
				/(nTrk[depth-1][kcur] - 1));
    else rms[depth-1][kcur] = RESOLUTION_HCAL;
    kdep[depth-1]++;
  }
  for ( unsigned ik = 0; ik < N_DEPTHS; ik++ ) {
    double x[MAXNUM_SUBDET];
    double ytrk[MAXNUM_SUBDET], ysub[MAXNUM_SUBDET];
    double yphi[MAXNUM_SUBDET], yrms[MAXNUM_SUBDET];
    for ( unsigned im = 0; im < MAXNUM_SUBDET; im++ ) {
      x[im] = jeta[ik][im];
      ytrk[im] = nTrk[ik][im];
      ysub[im] = nSub[ik][im];
      yphi[im] = nPhi[ik][im];
      yrms[im] = rms[ik][im];
    }
    TGraph*  g_ntrk = new TGraph(kdep[ik], x, ytrk);
    sprintf(name, "Number of tracks for depth %1d", ik+1);
    g_ntrk->SetTitle(name);
    sprintf(name, "nTrk_depth%1d", ik+1);
    foutRootFile->WriteTObject(g_ntrk, name);

    TGraph*  g_nsub = new TGraph(kdep[ik], x, ysub);
    sprintf(name, "Mean number of active subdetectors, depth %1d", ik+1);
    g_nsub->SetTitle(name);
    sprintf(name, "nSub_depth%1d", ik+1);
    foutRootFile->WriteTObject(g_nsub, name);

    TGraph*  g_nphi = new TGraph(kdep[ik], x, yphi);
    sprintf(name, "Mean number of phi-merged subdetectors, depth %1d", ik+1);
    g_nphi->SetTitle(name);
    sprintf(name, "nPhi_depth%1d", ik+1);
    foutRootFile->WriteTObject(g_nphi, name);

    TGraph*  g_rms = new TGraph(kdep[ik], x, yrms);
    sprintf(name, "RMS of samples, depth %1d", ik+1);
    g_rms->SetTitle(name);
    sprintf(name, "rms_depth%1d", ik+1);
    foutRootFile->WriteTObject(g_rms, name);
  }

  //--- estimate ratio mean/MPV
  double xl = e2p_init->GetMean() - FIT_RMS_INTERVAL*e2p_init->GetRMS();
  double xr = e2p_init->GetMean() + FIT_RMS_INTERVAL*e2p_init->GetRMS();
  e2p_init->Fit("f1","QN", "R", xl, xr);
  xl = f1->GetParameter(1) - FIT_RMS_INTERVAL*f1->GetParameter(2);
  xr = f1->GetParameter(1) + FIT_RMS_INTERVAL*f1->GetParameter(2);
  e2p_init->Fit("f1","QN", "R", xl, xr);

  if ( shiftResp && (f1->GetParameter(1) != 0) ) {
    referenceResponse = e2p_init->GetMean()/f1->GetParameter(1);
    std::cout << "Use reference response=<mean from sample>/<mpv from fit>:"
	      << e2p_init->GetMean() << "/" << f1->GetParameter(1)
	      << " = " << referenceResponse //<< std::endl
	      << " (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	      << std::endl;
  }
  else {
    referenceResponse = 1;
    std::cout << "Use reference response = 1" << std::endl
	      << "<mean from sample>/<mpv from fit> = "
	      << e2p_init->GetMean()/f1->GetParameter(1)
	      << "  (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	      << std::endl;
  }
  //---- for HB
  xl = e2pHB_init->GetMean() - FIT_RMS_INTERVAL*e2pHB_init->GetRMS();
  xr = e2pHB_init->GetMean() + FIT_RMS_INTERVAL*e2pHB_init->GetRMS();
  e2pHB_init->Fit("f1","QN", "R", xl, xr);
  xl = f1->GetParameter(1) - FIT_RMS_INTERVAL*f1->GetParameter(2);
  xr = f1->GetParameter(1) + FIT_RMS_INTERVAL*f1->GetParameter(2);
  e2pHB_init->Fit("f1","QN", "R", xl, xr);

  if ( shiftResp && (f1->GetParameter(1) != 0) ) {
    referenceResponseHB = e2pHB_init->GetMean()/f1->GetParameter(1);
    std::cout << "In HB <mean from sample>/<mpv from fit> = "
	      << e2pHB_init->GetMean() << "/" << f1->GetParameter(1)
	      << " = " << referenceResponseHB //<< std::endl
	      << " (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	      << std::endl;
  }
  else {
    referenceResponseHB = 1;
    std::cout << "Use reference response in HB = 1" << std::endl
	      << "<mean from sample>/<mpv from fit> = "
	      << e2pHB_init->GetMean()/f1->GetParameter(1)
	      << "  (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	      << std::endl;
  }
  //---- for TR
  xl = e2pTR_init->GetMean() - FIT_RMS_INTERVAL*e2pTR_init->GetRMS();
  xr = e2pTR_init->GetMean() + FIT_RMS_INTERVAL*e2pTR_init->GetRMS();
  e2pTR_init->Fit("f1","QN", "R", xl, xr);
  xl = f1->GetParameter(1) - FIT_RMS_INTERVAL*f1->GetParameter(2);
  xr = f1->GetParameter(1) + FIT_RMS_INTERVAL*f1->GetParameter(2);
  e2pTR_init->Fit("f1","QN", "R", xl, xr);

  if ( shiftResp && (f1->GetParameter(1) != 0) ) {
    referenceResponseTR = e2pTR_init->GetMean()/f1->GetParameter(1);
    std::cout << "In TR <mean from sample>/<mpv from fit> = "
	      << e2pTR_init->GetMean() << "/" << f1->GetParameter(1)
	      << " = " << referenceResponseTR //<< std::endl
	      << " (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	      << std::endl;
  }
  else {
    referenceResponseTR = 1;
    std::cout << "Use reference response in TR = 1" << std::endl
	      << "<mean from sample>/<mpv from fit> = "
	      << e2pTR_init->GetMean()/f1->GetParameter(1)
	      << "  (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	      << std::endl;
  }
  //---- for HE
  xl = e2pHE_init->GetMean() - FIT_RMS_INTERVAL*e2pHE_init->GetRMS();
  xr = e2pHE_init->GetMean() + FIT_RMS_INTERVAL*e2pHE_init->GetRMS();
  e2pHE_init->Fit("f1","QN", "R", xl, xr);
  xl = f1->GetParameter(1) - FIT_RMS_INTERVAL*f1->GetParameter(2);
  xr = f1->GetParameter(1) + FIT_RMS_INTERVAL*f1->GetParameter(2);
  e2pHE_init->Fit("f1","QN", "R", xl, xr);

  if ( shiftResp && (f1->GetParameter(1) != 0) ) {
    referenceResponseHE = e2pHE_init->GetMean()/f1->GetParameter(1);
    std::cout << "In HE <mean from sample>/<mpv from fit> = "
	      << e2pHE_init->GetMean() << "/" << f1->GetParameter(1)
	      << " = " << referenceResponseHE //<< std::endl
	      << " (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	      << std::endl;
  }
  else {
    referenceResponseHE = 1;
    std::cout << "Use reference response in HE = 1" << std::endl
	      << "<mean from sample>/<mpv from fit> = "
	      << e2pHE_init->GetMean()/f1->GetParameter(1)
	      << "  (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	      << std::endl;
  }

  //----- print additional info
  std::cout << "Maximal response for good tracks = " 
	    << maxRespForGoodTrack << std::endl
	    << nRespOverHistLimit
	    << " events with response > " << MAX_RESPONSE_HIST
	    << "(hist limit for mean estimate)"
	    << std::endl;
  std::cout << "Minimal response for good tracks = " 
	    << minRespForGoodTrack
	    << std::endl;
  std::cout << "Maximum number of selected tracks per ieta bin = " 
	    << maxNumOfTracksForIeta
	    << std::endl;
  std::cout << "Number of selected tracks in HB = "
	    << e2pHB_init->GetEntries()
	    << std::endl;
  std::cout << "Number of selected tracks in TR = "
	    << e2pTR_init->GetEntries()
	    << std::endl;
  std::cout << "Number of selected tracks in HE = "
	    << e2pHE_init->GetEntries()
	    << std::endl;

  if ( istest ) {
    TGraph *g_chi = new TGraph(NUM_ETA_BINS);  
    TGraphErrors* g_e2pFit = new TGraphErrors(NUM_ETA_BINS);
    TGraphErrors* g_e2pMean = new TGraphErrors(NUM_ETA_BINS);
  
    int ipointF(0);
    int ipointM(0);
    for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
      int ieta = (i - HALF_NUM_ETA_BINS)*N_ETA_RINGS_PER_BIN;
      if ( N_ETA_RINGS_PER_BIN > 1 ) {
	ieta = (i > HALF_NUM_ETA_BINS) ? ieta+1 : ieta-1;
      }
      if ( abs(ieta) > MAX_TRACK_IETA ) continue;

      int nhistentries = e2p[i]->GetEntries();
      if ( nhistentries < 1 ) continue;
      else {
	g_e2pMean->SetPoint(ipointM, ieta, e2p[i]->GetMean());
	g_e2pMean->SetPointError(ipointM, 0, e2p[i]->GetMeanError());
	ipointM++;

	if ( nhistentries > MIN_N_ENTRIES_FOR_FIT ) {
	  //int nrebin = n_ieta_bins/(2*2.5*pow(nhistentries,1/3.0));
	  //if ( nrebin > 2 ) e2p[i]->Rebin(nrebin);
	
	  double xl = e2p[i]->GetMean() - FIT_RMS_INTERVAL*e2p[i]->GetRMS();
	  double xr = e2p[i]->GetMean() + FIT_RMS_INTERVAL*e2p[i]->GetRMS();
	  e2p[i]->Fit("f1","QN", "R", xl, xr);
	  xl = f1->GetParameter(1) - FIT_RMS_INTERVAL*f1->GetParameter(2);
	  xr = f1->GetParameter(1) + FIT_RMS_INTERVAL*f1->GetParameter(2);
	  e2p[i]->Fit("f1","QN", "R", xl, xr);
	  g_e2pFit->SetPoint(ipointF, ieta, f1->GetParameter(1));
	  g_e2pFit->SetPointError(ipointF, 0, f1->GetParError(1));
	  g_chi->SetPoint(ipointF, ieta, f1->GetChisquare()/f1->GetNDF());
	  ipointF++;
	}
      }
    }
    // fill number of tracks per ieta

    for ( int k = ipointF; k < NUM_ETA_BINS; k++ ) {
      g_e2pFit->RemovePoint(ipointF);
    }
    for ( int k = ipointM; k < NUM_ETA_BINS; k++ ) {
      g_e2pMean->RemovePoint(ipointM);
    }
    sprintf(name, "Test: initial response from fit");
    g_e2pFit->SetTitle(name);
    g_e2pFit->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "respFit_0");
    foutRootFile->WriteTObject(g_e2pFit, name);

    sprintf(name, "Test: initial mean response");
    g_e2pMean->SetTitle(name);
    g_e2pMean->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "respMean_0");
    foutRootFile->WriteTObject(g_e2pMean, name);

    sprintf(name, "Test: initial chi2/NDF");
    g_chi->SetTitle(name);
    g_chi->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "chi2ndf_0");
    foutRootFile->WriteTObject(g_chi, name); 

//--- delete hists ---------------------------
    for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
      delete e2p[i];
    }
  }

  std::cout << "Run range (all sample): " << minRunNumber << " -- " << maxRunNumber << std::endl;

  return nSelectedEvents;
}

//**********************************************************
// Loop over events in the tree for current iteration
//**********************************************************
Double_t CalibTree::loopForIteration(unsigned int subsample,
				     unsigned int nIter,
				     unsigned int debug )
{
  char name[500];
  double meanDeviation = 0;
  unsigned int ndebug(0);
    
  TF1* f1 = new TF1("f1","gaus", MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  TH1F* e2p[NUM_ETA_BINS];

  int n_ieta_bins = 2*2.5*pow(maxNumOfTracksForIeta,1/3.0);
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    sprintf(name,"e2p[%02d]", i);
    e2p[i] = new TH1F(name, "",
		      n_ieta_bins, //NBIN_RESPONSE_HIST_IND,
		      MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
    e2p[i]->Sumw2();
  }

  std::map<unsigned int, std::pair<double,double> > sumsForFactorCorrection;
  std::map<unsigned int, double> sumOfWeightsSquared;

  if ( debug > 0 ) {
    std::cout.precision(3);
    std::cout << "-------------------------------------------- nIter = "
	      << nIter << std::endl;
  }
//--- initialize chain ----------------------------------------
  if (fChain == 0) return 0;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nb = 0;
  
// ----------------------- loop over events -------------------------------------  
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if ( ientry < 0 || ndebug > debug ) break;   
    nb = fChain->GetEntry(jentry);   //nbytes += nb;
    
    if ( (jentry%2 == subsample) ) continue;   // only odd or even events

    // --------------- selection of good track --------------------
    if ( t_nVtx < NVERTEX_MIN  ||  t_nVtx > NVERTEX_MAX
	 || t_Run < MIN_RUNNUM  || t_Run > MAX_RUNNUM ) continue;
    if ( t_ieta	< -MAX_TRACK_IETA || t_ieta > MAX_TRACK_IETA ) continue;
    
    if ( !goodTrack() ) continue;

    if ( debug > 1 ) {
      ndebug++;
      std::cout << "***Entry (Track) Number : " << ientry 
		<< " p/eHCal/eMipDR/nDets : " << t_p << "/" << t_eHcal 
		<< "/" << t_eMipDR << "/" << (*t_DetIds).size() 
		<< std::endl;
    }
    
    double eTotal(0.0);
    //double eTotalWithEcal(0.0);
      
    // ---- first loop over active subdetectors in the event for total energy ---

    for (unsigned int idet = 0; idet < (*t_DetIds).size(); idet++) { 
      double hitEnergy(0);	
      unsigned int detId = ( (*t_DetIds)[idet] & MASK ) | MASK2 ;
	
      if (factors.find(detId) != factors.end()) 
	hitEnergy = factors[detId] * (*t_HitEnergies)[idet];
      else 
	hitEnergy = (*t_HitEnergies)[idet];

      eTotal += hitEnergy;
    }

    //eTotalWithEcal = eTotal + t_eMipDR;    

// --- Correction for PU   --------      
    double eTotalCor(eTotal);
    //double eTotalWithEcalCor(eTotalWithEcal);
    double correctionForPU(1.0);

    if ( APPLY_CORRECTION_FOR_PU ) {
      double de2p = (t_eHcal30 - t_eHcal10)/t_p;
      if ( de2p > DELTA_CUT ) {
	int abs_t_ieta = abs(t_ieta);
	int icor = int(abs_t_ieta >= FIRST_IETA_TR) + int(abs_t_ieta >= FIRST_IETA_HE)
	  + int(abs_t_ieta >= FIRST_IETA_FWD_1) + int(abs_t_ieta >= FIRST_IETA_FWD_2);
	correctionForPU = (1 + LINEAR_COR_COEF[icor]*(t_eHcal/t_p)*de2p
			   *(1 + SQUARE_COR_COEF[icor]*de2p));
      }
    }    

    // check for possibility to correct for PU
    if ( correctionForPU <= 0 || correctionForPU > 1 ) continue;

    eTotalCor = eTotal*correctionForPU;
    //eTotalWithEcalCor = eTotalCor + t_eMipDR;
     
    //double response = eTotalWithEcalCor/t_p; // - referenceResponse;
    double response = eTotalCor/(t_p - t_eMipDR); // - referenceResponse;

    int jeta = HALF_NUM_ETA_BINS + int(t_ieta/N_ETA_RINGS_PER_BIN);
    e2p[jeta]->Fill(response,1.0);
      
// ---- second loop over active subdetectors in the event  -----------------
    for (unsigned int idet = 0; idet < (*t_DetIds).size(); idet++) {
      double hitEnergy(0);
      unsigned int detId = ( (*t_DetIds)[idet] & MASK ) | MASK2 ;
		 
      if (factors.find(detId) != factors.end())
	hitEnergy = factors[detId] * (*t_HitEnergies)[idet];
      else
	hitEnergy = (*t_HitEnergies)[idet];

      double cellWeight = hitEnergy/eTotal;   
      //double trackWeight = (cellWeight * t_p) / eTotalWithEcalCor; // old method
      double trackWeight = cellWeight*response;   // new method
      double cellweight2 = cellWeight*cellWeight;
      
      if( sumsForFactorCorrection.find(detId) != sumsForFactorCorrection.end() ) {
	cellWeight  += sumsForFactorCorrection[detId].first;
	trackWeight += sumsForFactorCorrection[detId].second;
	sumsForFactorCorrection[detId] = std::pair<double,double>(cellWeight,trackWeight);
	sumOfWeightsSquared[detId] += cellweight2;
      }
      else {
	sumsForFactorCorrection.insert(std::pair<unsigned int,
				       std::pair<double,double> >(detId,
								  std::pair<double,double>(cellWeight,
											   trackWeight)));
	sumOfWeightsSquared.insert(std::pair<unsigned int,double>(detId, cellweight2));
     }
	
      if ( debug > 1 ) { //|| hitEnergy < -0.5) {
	double f = 1;
	int zside= (detId&ZSIDE_MASK) ? 1 : -1;
	if (factors.find(detId) != factors.end()) f = factors[detId];
	std::cout << jentry << "::: "
		  << " Ncells: " << (*t_DetIds).size()
		  << " !! detId(ieta)/e/f : " 
	    //    << std::hex << (*t_DetIds)[idet] << ":"
		  << detId << "(" << int((detId>>ETA_OFFSET) & ETA_MASK)*zside << ")"
		  << "/" << hitEnergy
		  << "/" << f
		  << " ||| cellW/trW : " << cellWeight << " / " << trackWeight
		  << " ||| E/Ecor/p : " << eTotal
		  << " / " << eTotalCor
		  << " / " << t_p
		  << " || e10/e30/cF : " << t_eHcal10
		  << " / " << t_eHcal30
		  << " / " << correctionForPU
		  << std::endl;
      }
    }  // --------------- end of second loop over cells ----------
  } // ------------------- end of loop over events -------------------------------------

//----- Save initial histograms for each ieta -------------
  if ( nIter == 1 ) {
    for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
      int inum = (i - HALF_NUM_ETA_BINS)*N_ETA_RINGS_PER_BIN;
      sprintf(name,"e2pI_%d", inum);
      if ( e2p[i]->GetEntries() > 0 ) e2p[i]->Clone(name);
    }
  }

//----- Graphs to be saved in root file ----------------
  if ( debug > 0 ) {
    std::cout << "Fit and calculate means..." << std::endl;
    std::cout << "Number of plots (ieta bins) = " << NUM_ETA_BINS << std::endl;
  }
  TGraph *g_chi = new TGraph(NUM_ETA_BINS);  
  TGraphErrors* g_e2pFit = new TGraphErrors(NUM_ETA_BINS);
  TGraphErrors* g_e2pMean = new TGraphErrors(NUM_ETA_BINS);
  TGraph *g_nhistentries = new TGraph(NUM_ETA_BINS);
  
  int ipointF(0);
  int ipointM(0);
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    int ieta = (i - HALF_NUM_ETA_BINS)*N_ETA_RINGS_PER_BIN;
    if ( N_ETA_RINGS_PER_BIN > 1 ) {
      ieta = (i > HALF_NUM_ETA_BINS) ? ieta+1 : ieta-1;
    }
    int nhistentries = e2p[i]->GetEntries();
    /*
      if ( debug > 0 ) {
	std::cout << "i / entries / ieta :::"
		  << i
		  << " / " << nhistentries
		  << " / " << ieta
		  << std::endl;
      }
    */
     if ( nIter == 1 ) {
       g_nhistentries->SetPoint(i, ieta, nhistentries);
     }

    if ( nhistentries < 1 ) continue;
    else {
      g_e2pMean->SetPoint(ipointM, ieta, e2p[i]->GetMean());
      g_e2pMean->SetPointError(ipointM, 0, e2p[i]->GetMeanError());
      ipointM++;

      if ( nhistentries > MIN_N_ENTRIES_FOR_FIT ) {
	//int nrebin = n_ieta_bins/(2*2.5*pow(nhistentries,1/3.0));
	//if ( nrebin > 2 ) e2p[i]->Rebin(nrebin);
	
	double xl = e2p[i]->GetMean() - FIT_RMS_INTERVAL*e2p[i]->GetRMS();
	double xr = e2p[i]->GetMean() + FIT_RMS_INTERVAL*e2p[i]->GetRMS();
	e2p[i]->Fit("f1","QN", "R", xl, xr);
	xl = f1->GetParameter(1) - FIT_RMS_INTERVAL*f1->GetParameter(2);
	xr = f1->GetParameter(1) + FIT_RMS_INTERVAL*f1->GetParameter(2);
	e2p[i]->Fit("f1","QN", "R", xl, xr);
	g_e2pFit->SetPoint(ipointF, ieta, f1->GetParameter(1));
	g_e2pFit->SetPointError(ipointF, 0, f1->GetParError(1));
	g_chi->SetPoint(ipointF, ieta, f1->GetChisquare()/f1->GetNDF());
	ipointF++;
      }
    }
  }
  // fill number of tracks per ieta
  if ( nIter == 1 ) {
    sprintf(name, "Number of selected tracks");
    g_nhistentries->SetTitle(name);
    g_nhistentries->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "selTrks");
    foutRootFile->WriteTObject(g_nhistentries, name);
  }

  for ( int k = ipointF; k < NUM_ETA_BINS; k++ ) {
    g_e2pFit->RemovePoint(ipointF);
  }
  for ( int k = ipointM; k < NUM_ETA_BINS; k++ ) {
    g_e2pMean->RemovePoint(ipointM);
  }
  sprintf(name, "Response from fit, iteration %d", nIter-1);
  g_e2pFit->SetTitle(name);
  g_e2pFit->GetXaxis()->SetTitle("i#eta");
  sprintf(name, "respFit_%d", nIter-1);
  foutRootFile->WriteTObject(g_e2pFit, name);

  sprintf(name, "Mean response, iteration %d", nIter-1);
  g_e2pMean->SetTitle(name);
  g_e2pMean->GetXaxis()->SetTitle("i#eta");
  sprintf(name, "respMean_%d", nIter-1);
  foutRootFile->WriteTObject(g_e2pMean, name);

  sprintf(name, "Chi2/NDF, iteration %d", nIter-1);
  g_chi->SetTitle(name);
  g_chi->GetXaxis()->SetTitle("i#eta");
  sprintf(name, "chi2ndf_%d", nIter-1);
  foutRootFile->WriteTObject(g_chi, name);

// --- convergence criteria and correction factors -----------------------------------

  double MeanConvergenceDelta(0),  MaxRelDeviationWeights(0), MaxRatioUncertainties(0);
  double dets[MAXNUM_SUBDET];
  double ztest[MAXNUM_SUBDET], sys2statRatio[MAXNUM_SUBDET];

  if ( debug > 0 ) std::cout << "Calculate correction factors..." << std::endl;

  unsigned int kount(0), mkount(0);
  unsigned int maxKountW(0), maxKountR(0);

  double fac[N_DEPTHS][MAXNUM_SUBDET];
  double dfac[N_DEPTHS][MAXNUM_SUBDET];
  double ieta[N_DEPTHS][MAXNUM_SUBDET];
  double dieta[N_DEPTHS][MAXNUM_SUBDET];
  unsigned iscalib[N_DEPTHS][MAXNUM_SUBDET];

  unsigned int kdep[N_DEPTHS];
  unsigned int kdepU[N_DEPTHS];
  for ( unsigned ik = 0; ik < N_DEPTHS; ik++ ) { kdep[ik] = 0; kdepU[ik] = 0; }

//-------------- loop over all cells, fill calibrated ---------
  for (std::map <unsigned int,
	 std::pair<double,double> >::iterator sumsForFactorCorrectionItr
	 = sumsForFactorCorrection.begin();
       sumsForFactorCorrectionItr != sumsForFactorCorrection.end();
       sumsForFactorCorrectionItr++) {

    unsigned int detId = sumsForFactorCorrectionItr->first;
    int zside = (detId&ZSIDE_MASK) ? 1 : -1;
    int depth = ((detId>>DEPTH_OFFSET)&DEPTH_MASK) + int(MERGE_PHI_AND_DEPTHS);
    unsigned int kcur = kdep[depth-1];    
    ieta[depth-1][kcur] = int((detId>>ETA_OFFSET) & ETA_MASK)*zside;
    dieta[depth-1][kcur] = 0;
    
    double sumOfWeights = (sumsForFactorCorrectionItr->second).first;
    int nSubDetTracks(0);
    double subdetRMS(RESOLUTION_HCAL);
    if ( nTrks.find(detId) != nTrks.end() ) {
      nSubDetTracks = nTrks[detId];
      if ( nSubDetTracks > 1 ) 
	subdetRMS = sqrt((sumOfResponseSquared[detId] 
			  - pow(sumOfResponse[detId],2)/double(nSubDetTracks))
			 /double(nSubDetTracks - 1));
    }
    else {
      std::cout << "!!!!!!! No tracks for subdetector " << detId << std::endl;
      continue;
    }
    double NcellMean = double(nSubdetInEvent[detId])/double(nSubDetTracks);

    double ratioWeights(1);
    if ( abs(sumOfWeights) > 0 )
      ratioWeights = sqrt(sumOfWeightsSquared[detId])/sumOfWeights;
    double correctionRMS = subdetRMS*ratioWeights*sqrt(NcellMean);
    
    double absErrorW(0);
    double absErrorWprevious(0);
    double factorPrevious(1);
    double factorCorrection(1);

    double refR = referenceResponse;
    if ( !SINGLE_REFERENCE_RESPONSE ) {
      if ( abs(ieta[depth-1][kcur]) < FIRST_IETA_TR )
	refR = referenceResponseHB;
      else if ( abs(ieta[depth-1][kcur]) < FIRST_IETA_HE )
	refR = referenceResponseTR;
      else
	refR = referenceResponseHE;
    }

    if ( abs(sumOfWeights) > 0 )  
      factorCorrection = 1 + refR
	- (sumsForFactorCorrectionItr->second).second / sumOfWeights;
    
    //------- set factor=1 to uncalibrated -----------------------
    if ( correctionRMS/factorCorrection > MAX_REL_UNC_FACTOR
	 || nSubDetTracks < MIN_N_TRACKS_PER_CELL
	 || abs(ieta[depth-1][kcur]) >  MAX_CALIBRATED_IETA
	 ) {
      correctionRMS = sqrt(pow(correctionRMS,2) + pow((factorCorrection - 1),2));
      factorCorrection = 1;
      iscalib[depth-1][kcur] = 0;
    }
    else {
      iscalib[depth-1][kcur] = 1;
      if (factorCorrection > 1) MeanConvergenceDelta += (1 - 1/factorCorrection);
      else                      MeanConvergenceDelta += (1 - factorCorrection);
      mkount++;
    }
    //---- fill map with calibrated and uncalibrated factors
    if (factors.find(detId) != factors.end()) {
      factorPrevious = factors[detId];
      factors[detId] *= factorCorrection;
      absErrorWprevious = uncFromWeights[detId];
      absErrorW = factorPrevious*correctionRMS;
      uncFromWeights[detId] = absErrorW;
      uncFromDeviation[detId] = factorPrevious*abs(factorCorrection - 1);
    }
    else {
      factorPrevious = 1;
      factors.insert(std::pair<unsigned int, double>(detId, factorCorrection));
      subDetector_final.insert(std::pair<unsigned int, double>(detId,
							       subDetector_trk[detId]));
      absErrorW = correctionRMS;
      absErrorWprevious = 0;
      uncFromWeights.insert(std::pair<unsigned int, double>(detId, absErrorW));
      uncFromDeviation.insert(std::pair<unsigned int, double>(detId,
							      abs(factorCorrection - 1)));
    }

    if ( debug > 0 ) {
      //if ( ieta[depth-1][kcur] == 27 && depth == 2 ) {
      std::cout.precision(3);
      std::cout << detId // << " (" << mkount << ")"
	        << " *** ieta/depth | rw | cw | tw | fCor | nTrk | Ncell | C |::: "
	        << ieta[depth-1][kcur] << "/" << depth << " | "
		<< ratioWeights << " | "
	        << sumOfWeights << " | "
	        << (sumsForFactorCorrectionItr->second).second << " | "
	        << factorCorrection << " | "
		<< nSubDetTracks << " | "
		<< NcellMean << " | "
		<< correctionRMS << " |"
	        << std::endl;
    }

    dets[kount] = detId;

    fac[depth-1][kcur] = factors[detId];
    dfac[depth-1][kcur] =
      sqrt(pow(uncFromWeights[detId],2) + pow(uncFromDeviation[detId],2));
    
    sys2statRatio[kount] = abs(factorPrevious*(factorCorrection - 1))/absErrorW;
    if ( sys2statRatio[kount] > MaxRatioUncertainties ) {
      MaxRatioUncertainties = sys2statRatio[kount];
      maxKountR = kount;
    }
    ztest[kount] = factorPrevious*(factorCorrection - 1)
      /sqrt(pow(absErrorWprevious,2) + pow(absErrorW,2));
    if ( abs(ztest[kount]) > MaxRelDeviationWeights ) {
      MaxRelDeviationWeights = abs(ztest[kount]);
      maxKountW = kount;
    } 
    kount++;
    kdep[depth-1]++;
  } // ---- end of loop over all cells
  
//-------------- loop over all cells, fill uncalibrated ---------------------------------
  if ( UNCALIB_FACT_TO_ONE == 0 ) {
    for (std::map <unsigned int,
	   std::pair<double,double> >::iterator sumsForFactorCorrectionItr
	   = sumsForFactorCorrection.begin();
	 sumsForFactorCorrectionItr != sumsForFactorCorrection.end();
	 sumsForFactorCorrectionItr++) {

      unsigned int detId = sumsForFactorCorrectionItr->first;
      int zside = (detId&ZSIDE_MASK) ? 1 : -1;
      int depth = ((detId>>DEPTH_OFFSET)&DEPTH_MASK) + int(MERGE_PHI_AND_DEPTHS);
      unsigned int kcur = kdepU[depth-1];
      unsigned jk = depth - 1;
      double ieta_cur = int((detId>>ETA_OFFSET) & ETA_MASK)*zside;
    
      if ( iscalib[jk][kcur] == 0) {
	double fL(1), sL(0);
	double fR(1), sR(0);
	bool iscalibL(0), iscalibR(0);
	if ( abs(ieta_cur) <= MAX_CALIBRATED_IETA ) { // try to find calibrated neighbors
	  for ( unsigned im = 0; im < MAXNUM_SUBDET; im++ ) {
	    if ( (ieta[jk][im] - ieta_cur) == 1 ) {
	      fR = fac[jk][im]; sR = dfac[jk][im];
	      iscalibR = iscalib[jk][im];
	    }
	    if ( (ieta[jk][im] - ieta_cur) == -1 ) {
	      fL = fac[jk][im]; sL = dfac[jk][im];
	      iscalibL = iscalib[jk][im];
	    }
	  }
	  if ( iscalibR == 1 || iscalibL == 1 ) {
	    fac[jk][kcur] = 0.5*(fL + fR);
	    dfac[jk][kcur] = sqrt((pow(sR,2) + pow(sL,2))/2);
	    factors[detId] = fac[jk][kcur];
	    uncFromWeights[detId] = dfac[jk][kcur];
	    uncFromDeviation[detId] = 0;
	  }
	}
	else { // take factor from MAX_CALIBRATED_IETA
	  for ( unsigned im = 0; im < MAXNUM_SUBDET; im++ ) {
	    if ( abs(ieta[jk][im]) == MAX_CALIBRATED_IETA
		 && iscalib[jk][im] == 1
		 && (ieta[jk][im]*ieta_cur) > 0 )
	      fR = fac[jk][im];
	  }
	  fac[jk][kcur] = fR;
	  factors[detId] = fR;
	  dfac[jk][kcur] = uncFromWeights[detId];
	  //uncFromWeights[detId] = dfac[jk][kcur];
	  //uncFromDeviation[detId] = 0;
	}
      }
      kdepU[jk]++;
    }
  }
//---- write current plots -----------------------------
  if ( debug > 0 ) std::cout << "Write graphs..." << std::endl;

  for ( unsigned ik = 0; ik < N_DEPTHS; ik++ ) {
    double x[MAXNUM_SUBDET], dx[MAXNUM_SUBDET], y[MAXNUM_SUBDET], dy[MAXNUM_SUBDET];
    for ( unsigned im = 0; im < MAXNUM_SUBDET; im++ ) {
      x[im] = ieta[ik][im];
      dx[im] = dieta[ik][im];
      y[im] = fac[ik][im];
      dy[im] = dfac[ik][im];
    }
    TGraphErrors*  g_fac = new TGraphErrors(kdep[ik], x, y, dx, dy);
    sprintf(name, "Extracted correction factors, depth %d", ik+1);
    g_fac->SetTitle(name);
    g_fac->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "Cfacs_depth%1d_%d", ik+1, nIter-1);
    foutRootFile->WriteTObject(g_fac, name);
  }

  TGraph  *g_ztest, *g_sys2stat;

  g_ztest = new TGraph(kount, dets, ztest); 
  sprintf(name, "Z-test (unc. from weights) vs detId for iter %d", nIter-1);
  g_ztest->SetTitle(name);
  sprintf(name, "Ztest_detId_%d", nIter-1);
  foutRootFile->WriteTObject(g_ztest, name);

  g_sys2stat = new TGraph(kount, dets, sys2statRatio); 
  sprintf(name, "Ratio of syst. to stat. unc. vs detId for iter %d", nIter-1);
  g_sys2stat->SetTitle(name);
  sprintf(name, "Sys2stat_detId_%d", nIter-1);
  foutRootFile->WriteTObject(g_sys2stat, name);

  std::cout << "----------Iteration " << nIter << "--------------------" << std::endl;
  maxZtestFromWeights = MaxRelDeviationWeights; 
  std::cout << "Max abs(Z-test) with stat errors from weights = "
	    << maxZtestFromWeights << " for subdetector " << maxKountW << std::endl;
  maxSys2StatRatio = MaxRatioUncertainties; 
  std::cout << "Max ratio of syst.(f_cur - f_prev) to stat. uncertainty = "
	    << maxSys2StatRatio << " for subdetector " << maxKountR << std::endl;

  meanDeviation = (mkount > 0) ? (MeanConvergenceDelta/mkount) : 0;
  std::cout << "Mean absolute deviation from previous iteration = " << meanDeviation
	    << " for " << mkount
	    << " from " << kount << " DetIds" << std::endl;

//--- delete hists ---------------------------
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    delete e2p[i];
  }

  return meanDeviation;
}
//**********************************************************
// Last loop over events in the tree
//**********************************************************
Double_t CalibTree::lastLoop(unsigned int subsample,
			     unsigned int maxIter,
			     bool isTest,
			     unsigned int debug)
{
  char name[100];
  unsigned int ndebug(0);
  
  char stest[80] = "test";
  if ( !isTest )
    sprintf(stest,"after %2d iterations", maxIter);
  char scorr[80] = "correction for PU";
  char sxlabel[80] ="(E^{cor}_{hcal} + E_{ecal})/p_{track}"; 
  if ( !APPLY_CORRECTION_FOR_PU ) {
    sprintf(scorr,"no correction for PU");
    sprintf(sxlabel,"(E_{hcal} + E_{ecal})/p_{track}");
  } 
    
  TF1* f1 = new TF1("f1","gaus", MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  TH1F* e2p[NUM_ETA_BINS];

  int n_ieta_bins = 2*2.5*pow(maxNumOfTracksForIeta,1/3.0);
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    sprintf(name,"e2p[%02d]", i);
    e2p[i] = new TH1F(name, "",
		      n_ieta_bins, //NBIN_RESPONSE_HIST_IND,
		      MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
    e2p[i]->Sumw2();
  }

  sprintf(name,"HB+HE: %s, %s", stest, scorr);
  e2p_last = new TH1F("e2p_last", name,
		      NBIN_RESPONSE_HIST, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2p_last->Sumw2();
  e2p_last->GetXaxis()->SetTitle(sxlabel);
  
  sprintf(name,"HB: %s, %s", stest, scorr);
  e2pHB_last = new TH1F("e2pHB_last", name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHB_last->Sumw2();
  e2pHB_last->GetXaxis()->SetTitle(sxlabel);

  sprintf(name,"Initial TR: %s", scorr);
  e2pTR_last = new TH1F("e2pTR_last", name,
			NBIN_RESPONSE_HIST/10, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pTR_last->Sumw2();
  e2pTR_last->GetXaxis()->SetTitle(sxlabel);
  
  sprintf(name,"HE: %s, %s", stest, scorr);
  e2pHE_last = new TH1F("e2pHE_last", name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHE_last->Sumw2();
  e2pHE_last->GetXaxis()->SetTitle(sxlabel);
 
  //-- for SiPM module studies
  char snm[20];
  sprintf(snm,"e2pHEPf[00]");
  sprintf(name,"Final HEP, iphi 1,2,71,72");
  e2pHEPf[0] = new TH1F(snm, name,
		       NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHEPf[0]->Sumw2();
  e2pHEPf[0]->GetXaxis()->SetTitle(sxlabel);

  sprintf(snm,"e2pHEMf[00]");
  sprintf(name,"Final HEM, iphi 1,2,71,72");
  e2pHEMf[0] = new TH1F(snm, name,
		       NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHEMf[0]->Sumw2();
  e2pHEMf[0]->GetXaxis()->SetTitle(sxlabel);
  
  for ( int imd = 1; imd < N_MODULE_GROUPS; imd++ ) {
    int imin = 3 + 4*(imd - 1); //(10 + 20*imd)*TMath::Pi()/180; 
    int imax = imin + 3; //20*TMath::Pi()/180; 
    sprintf(snm,"e2pHEPf[%02d]",imd);
    sprintf(name,"Final HEP, iphi %d-%d", imin, imax);
    e2pHEPf[imd] = new TH1F(snm, name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
    e2pHEPf[imd]->Sumw2();
    e2pHEPf[imd]->GetXaxis()->SetTitle(sxlabel);

    sprintf(snm,"e2pHEMf[%02d]",imd);
    sprintf(name,"Final HEM, iphi %d-%d", imin, imax);
    e2pHEMf[imd] = new TH1F(snm, name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
    e2pHEMf[imd]->Sumw2();
    e2pHEMf[imd]->GetXaxis()->SetTitle(sxlabel);
  }
//--- initialize chain ----------------------------------------
  if (fChain == 0) return 0;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nb = 0;
  
  if ( debug > 0 ) { 
    std::cout << "------------- Last loop after " << maxIter << " iterations"
	      << std::endl;
  }
// ----------------------- loop over events -------------------------------------  
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if ( ientry < 0 || ndebug > debug ) break;   
    nb = fChain->GetEntry(jentry);   //nbytes += nb;
    
    if ( (jentry%2 == subsample) ) continue;   // only odd or even events

    // --------------- selection of good track --------------------
    
    if ( t_nVtx < NVERTEX_MIN  ||  t_nVtx > NVERTEX_MAX
	 || t_Run < MIN_RUNNUM  || t_Run > MAX_RUNNUM ) continue;
    if ( t_ieta	< -MAX_TRACK_IETA || t_ieta > MAX_TRACK_IETA ) continue;

    if ( !goodTrack() ) continue;
    
    if ( debug > 1 ) {
      ndebug++;
      std::cout << "***Entry (Track) Number : " << ientry 
		<< " p/eHCal/eMipDR/nDets : " << t_p << "/" << t_eHcal 
		<< "/" << t_eMipDR << "/" << (*t_DetIds).size() 
		<< std::endl;
    }
    
    double eTotal(0.0);
    //double eTotalWithEcal(0.0);
      
    // ---- loop over active cells in the event for total energy ---

    for (unsigned int idet = 0; idet < (*t_DetIds).size(); idet++) { 
      double hitEnergy(0);	
      unsigned int detId = ( (*t_DetIds)[idet] & MASK ) | MASK2 ;
	
      if (factors.find(detId) != factors.end()) 
	hitEnergy = factors[detId] * (*t_HitEnergies)[idet];
      else 
	hitEnergy = (*t_HitEnergies)[idet];

      eTotal += hitEnergy;
    }

    //eTotalWithEcal = eTotal + t_eMipDR;    

// --- Correction for PU   --------      
    double eTotalCor(eTotal);
    //double eTotalWithEcalCor(eTotalWithEcal);
    double correctionForPU(1.0);
    int abs_t_ieta = abs(t_ieta);

    if ( APPLY_CORRECTION_FOR_PU ) {       
      double de2p = (t_eHcal30 - t_eHcal10)/t_p;
      if ( de2p > DELTA_CUT ) {
	int icor = int(abs_t_ieta >= FIRST_IETA_TR) + int(abs_t_ieta >= FIRST_IETA_HE)
	  + int(abs_t_ieta >= FIRST_IETA_FWD_1) + int(abs_t_ieta >= FIRST_IETA_FWD_2);
	correctionForPU = (1 + LINEAR_COR_COEF[icor]*(t_eHcal/t_p)*de2p
			   *(1 + SQUARE_COR_COEF[icor]*de2p));
      }
    }      
    // check for possibility to correct for PU
    if ( correctionForPU <= 0 || correctionForPU > 1 ) continue;

    eTotalCor = eTotal*correctionForPU;
    //eTotalWithEcalCor = eTotalCor + t_eMipDR;
    double response = eTotalCor/(t_p - t_eMipDR);

    int jeta = HALF_NUM_ETA_BINS + int(t_ieta/N_ETA_RINGS_PER_BIN);
    e2p[jeta]->Fill(response,1.0);
          
// --- Fill final histograms ---------------------------      
    e2p_last->Fill(response, 1.0);

    if ( abs_t_ieta < FIRST_IETA_TR )
      e2pHB_last->Fill(response, 1.0);
    else if ( abs_t_ieta < FIRST_IETA_HE )
      e2pTR_last->Fill(response, 1.0);
    else
      e2pHE_last->Fill(response, 1.0);

    // --- For SiPM module crosscheck ----
    double abs_phi = abs(t_phi)*360/TMath::TwoPi();
    if ( abs_t_ieta <= ABSIETA_MAX && abs_t_ieta >= ABSIETA_MIN
	 && abs_phi <= ZERO_ANGLE_IN_DEGREE ) {
	if ( t_ieta < 0 ) e2pHEMf[0]->Fill(response ,1.0);
	else e2pHEPf[0]->Fill(response ,1.0);
      }
    double phi = 72*(t_phi + (t_phi<0)*TMath::TwoPi())/TMath::TwoPi();    
    for ( int imd = 1; imd < N_MODULE_GROUPS; imd++ ) {
      int imin = 3 + 4*(imd - 1); 
      int imax = imin + 3;
      if ( abs_t_ieta <= ABSIETA_MAX && abs_t_ieta >= ABSIETA_MIN
	   && phi <= imax && phi >= imin ) {
	if ( t_ieta < 0 ) e2pHEMf[imd]->Fill(response ,1.0);
	else e2pHEPf[imd]->Fill(response ,1.0);
      }
    }
      
  } // ------------------- end of loop over events -------------------------------------
//----- Save initial histograms for each ieta -------------
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    int inum = (i - HALF_NUM_ETA_BINS)*N_ETA_RINGS_PER_BIN;
    sprintf(name,"e2pL_%d", inum);
    if ( e2p[i]->GetEntries() > 0 ) e2p[i]->Clone(name);
  }

  if ( isTest ) {
    double fac[N_DEPTHS][MAXNUM_SUBDET]; 
    double dfac[N_DEPTHS][MAXNUM_SUBDET];
    double ieta[N_DEPTHS][MAXNUM_SUBDET];
    double dieta[N_DEPTHS][MAXNUM_SUBDET];
    unsigned int kdep[N_DEPTHS];
    for ( unsigned ik = 0; ik < N_DEPTHS; ik++ ) { kdep[ik] = 0; }
  
    std::map<unsigned int, double>::iterator factorsItr = factors.begin();
    for (factorsItr=factors.begin(); factorsItr != factors.end(); factorsItr++){

      unsigned int detId = factorsItr->first;
      int zside = (detId&ZSIDE_MASK) ? 1 : -1;
      int depth = ((detId>>DEPTH_OFFSET)&DEPTH_MASK) + int(MERGE_PHI_AND_DEPTHS);

      unsigned int kcur = kdep[depth-1];
      ieta[depth-1][kcur] = int((detId>>ETA_OFFSET) & ETA_MASK)*zside;
      dieta[depth-1][kcur] = 0;
      fac[depth-1][kcur] = factorsItr->second;
      dfac[depth-1][kcur] = 0;
      kdep[depth-1]++;
    }
    for ( unsigned ik = 0; ik < N_DEPTHS; ik++ ) {
      double x[MAXNUM_SUBDET], dx[MAXNUM_SUBDET], y[MAXNUM_SUBDET], dy[MAXNUM_SUBDET];
      for ( unsigned im = 0; im < MAXNUM_SUBDET; im++ ) {
	x[im] = ieta[ik][im];
	dx[im] = dieta[ik][im];
	y[im] = fac[ik][im];
	dy[im] = dfac[ik][im];
      }
      TGraphErrors*  g_fac = new TGraphErrors(kdep[ik], x, y, dx, dy);
      sprintf(name, "Applied correction factors, depth %1d", ik+1);
      g_fac->SetTitle(name);
      g_fac->GetXaxis()->SetTitle("i#eta");
      sprintf(name, "Cfacs_depth%1d", ik+1);
      foutRootFile->WriteTObject(g_fac, name);
    }
  }
  
  TGraph *g_chi = new TGraph(NUM_ETA_BINS);  
  TGraphErrors* g_e2pFit = new TGraphErrors(NUM_ETA_BINS);
  TGraphErrors* g_e2pMean = new TGraphErrors(NUM_ETA_BINS);
  
  int ipointF(0);
  int ipointM(0);
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    int ieta = (i - HALF_NUM_ETA_BINS)*N_ETA_RINGS_PER_BIN;
    if ( N_ETA_RINGS_PER_BIN > 1 ) {
      ieta = (i > HALF_NUM_ETA_BINS) ? ieta+1 : ieta-1;
    }
    if ( abs(ieta) > MAX_TRACK_IETA ) continue;

    int nhistentries = e2p[i]->GetEntries();
    if ( nhistentries < 1 ) continue;
    else {
      g_e2pMean->SetPoint(ipointM, ieta, e2p[i]->GetMean());
      g_e2pMean->SetPointError(ipointM, 0, e2p[i]->GetMeanError());
      ipointM++;

      if ( nhistentries > MIN_N_ENTRIES_FOR_FIT ) {
	//int nrebin = n_ieta_bins/(2*2.5*pow(nhistentries,1/3.0));
	//if ( nrebin > 2 ) e2p[i]->Rebin(nrebin);
	
	double xl = e2p[i]->GetMean() - FIT_RMS_INTERVAL*e2p[i]->GetRMS();
	double xr = e2p[i]->GetMean() + FIT_RMS_INTERVAL*e2p[i]->GetRMS();
	e2p[i]->Fit("f1","QN", "R", xl, xr);
	xl = f1->GetParameter(1) - FIT_RMS_INTERVAL*f1->GetParameter(2);
	xr = f1->GetParameter(1) + FIT_RMS_INTERVAL*f1->GetParameter(2);
	e2p[i]->Fit("f1","QN", "R", xl, xr);
	g_e2pFit->SetPoint(ipointF, ieta, f1->GetParameter(1));
	g_e2pFit->SetPointError(ipointF, 0, f1->GetParError(1));
	g_chi->SetPoint(ipointF, ieta, f1->GetChisquare()/f1->GetNDF());
	ipointF++;
      }
    }
  }
  // fill number of tracks per ieta

  for ( int k = ipointF; k < NUM_ETA_BINS; k++ ) {
    g_e2pFit->RemovePoint(ipointF);
  }
  for ( int k = ipointM; k < NUM_ETA_BINS; k++ ) {
    g_e2pMean->RemovePoint(ipointM);
  }
  if ( isTest ) {
    sprintf(name, "Test: response from fit");
    g_e2pFit->SetTitle(name);
    g_e2pFit->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "respFit_1");
    foutRootFile->WriteTObject(g_e2pFit, name);

    sprintf(name, "Test: mean response");
    g_e2pMean->SetTitle(name);
    g_e2pMean->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "respMean_1");
    foutRootFile->WriteTObject(g_e2pMean, name);

    sprintf(name, "Test: Chi2/NDF");
    g_chi->SetTitle(name);
    g_chi->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "chi2ndf_1");
    foutRootFile->WriteTObject(g_chi, name);
  }
  else {
    sprintf(name,"Response from fit: %d iterations, %s", maxIter, scorr);
    g_e2pFit->SetTitle(name);
    g_e2pFit->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "respFit_%d", maxIter);
    foutRootFile->WriteTObject(g_e2pFit, name);

    sprintf(name,"Mean response: %d iterations, %s", maxIter, scorr);
    g_e2pMean->SetTitle(name);
    g_e2pMean->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "respMean_%d", maxIter);
    foutRootFile->WriteTObject(g_e2pMean, name);

    sprintf(name,"Chi2/NDF: %d iterations, %s", maxIter, scorr);
    g_chi->SetTitle(name);
    g_chi->GetXaxis()->SetTitle("i#eta");
    sprintf(name, "chi2ndf_%d", maxIter);
    foutRootFile->WriteTObject(g_chi, name);
  }

  //--- fit response distributions ---------------------------------

  double xl = e2p_last->GetMean() - FIT_RMS_INTERVAL*e2p_last->GetRMS();
  double xr = e2p_last->GetMean() + FIT_RMS_INTERVAL*e2p_last->GetRMS();
  e2p_last->Fit("f1","QN", "R", xl, xr);
  xl = f1->GetParameter(1) - FIT_RMS_INTERVAL*f1->GetParameter(2);
  xr = f1->GetParameter(1) + FIT_RMS_INTERVAL*f1->GetParameter(2);
  e2p_last->Fit("f1","QN", "R", xl, xr);

  double fitMPV = f1->GetParameter(1);
  /*
  xl = e2pHB_last->GetMean() - FIT_RMS_INTERVAL*e2pHB_last->GetRMS();
  xr = e2pHB_last->GetMean() + FIT_RMS_INTERVAL*e2pHB_last->GetRMS();
  e2pHB_last->Fit("f1","QN", "R", xl, xr);

  xl = e2pHE_last->GetMean() - FIT_RMS_INTERVAL*e2pHE_last->GetRMS();
  xr = e2pHE_last->GetMean() + FIT_RMS_INTERVAL*e2pHE_last->GetRMS();
  e2pHE_last->Fit("f1","QN", "R", xl, xr);
  */
  
//--- delete hists ---------------------------
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    delete e2p[i];
  }

  return fitMPV;
}

//**********************************************************
// Isolated track selection
//**********************************************************
Bool_t CalibTree::goodTrack()
{

  //double maxCharIso = limCharIso*exp(abs(ieta)*constForFlexSel);

  bool ok = (    (t_selectTk)
	      && (t_qltyMissFlag)
		 && (t_hmaxNearP < limCharIso) //maxCharIso)
	      && (t_eMipDR < limMipEcal) 
	      && (t_p > minTrackMom) && (t_p < maxTrackMom)
	      && (t_pt >= minTrackPt)               // constraint on track pt
	      && (t_eHcal >= minEnrHcal)            // constraint on Hcal energy
	      && (t_eHcal/t_p < UPPER_LIMIT_RESPONSE_BEFORE_COR)
		 // reject events with too big cluster energy
		 //&& ((t_eHcal30 - t_eHcal10)/t_p < UPPER_LIMIT_DELTA_PU_COR)
		 // reject events with too high PU in the ring around cluster
	     );
  return ok;
}
//**********************************************************
// Save txt file with calculated factors
//**********************************************************
unsigned int CalibTree::saveFactorsInFile(std::string txtFileName)
{
  char sprnt[100];
  
  FILE* foutTxtFile = fopen(txtFileName.c_str(),"w+");
  fprintf(foutTxtFile,
	  "%1s%16s%16s%16s%9s%11s\n","#", "eta", "depth", "det", "value", "DetId");

  std::cout << "New factors:" << std::endl;
  std::map<unsigned int, double>::iterator factorsItr = factors.begin();
  unsigned int indx(0);
  unsigned int isave(0);
  
  for (factorsItr=factors.begin(); factorsItr != factors.end(); factorsItr++, indx++){
    unsigned int detId = factorsItr->first;
    int ieta = (detId>>ETA_OFFSET) & ETA_MASK;
    int zside= (detId&ZSIDE_MASK) ? 1 : -1;
    int depth= ((detId>>DEPTH_OFFSET)&DEPTH_MASK) + int(MERGE_PHI_AND_DEPTHS);

    double erWeight = 100*uncFromWeights[detId]/factorsItr->second;
    double erDev = 100*uncFromDeviation[detId]/factorsItr->second;
    double erTotal = 100*sqrt(pow(uncFromWeights[detId],2)
			  + pow(uncFromDeviation[detId],2))/factorsItr->second;
    
    if ( N_ETA_RINGS_PER_BIN < 2 ) { 
      sprintf(sprnt,
	      "DetId[%3d] %x (%3d,%1d)  %6.4f  : %6d  [%8.3f%% + %8.3f%% = %8.3f%%]",
	      indx, detId, ieta*zside, depth,
	      factorsItr->second, nTrks[detId],
	      erWeight, erDev, erTotal);
      std::cout << sprnt << std::endl;
    }
    else {
      int ieta_min = ieta - (N_ETA_RINGS_PER_BIN - 1);
      sprintf(sprnt,
	      "DetId[%3d] %x (%3d:%3d,%1d)  %6.4f  : %6d  [%8.3f%% + %8.3f%% = %8.3f%%]",
	      indx, detId, ieta_min*zside, ieta*zside, depth,
	      factorsItr->second, nTrks[detId],
	      erWeight, erDev, erTotal);
      std::cout << sprnt << std::endl;
    }
	    /*
    std::cout << "DetId[" << indx << "] " << std::hex  << (detId) << std::dec 
	      << "(" << ieta*zside << "," << depth << ") ( nTrks:" 
	      << nTrks[detId] << ") : " << factorsItr->second
	      << ""
	      << std::endl;
	    */
    
    const char* subDetector[2] = {"HB","HE"};
    if ( nTrks[detId] < MIN_N_TRACKS_PER_CELL ) continue;
    isave++;
    fprintf(foutTxtFile, "%17i%16i%16s%9.5f%11X\n", 
	    ieta*zside, depth, subDetector[subDetector_final[detId]-1],
	    factorsItr->second, detId);
  }
  fclose(foutTxtFile);
  foutTxtFile = NULL;
  return isave;
}
//**********************************************************
// Get factors from txt file
//**********************************************************
Bool_t CalibTree::getFactorsFromFile(std::string txtFileName,
				     unsigned int dbg)
{

  if ( !gSystem->Which("./", txtFileName.c_str() ) ) return false;

  FILE* finTxtFile = fopen(txtFileName.c_str(),"r");
  int flag;

  char header[80]; 
  for ( unsigned int i = 0; i < 6; i++ ) { 
    flag = fscanf(finTxtFile, "%7s", header);
  }

  int eta;
  int depth;
  char det[2]; 
  double cellFactor;
  unsigned int detId;
  unsigned int nReadFactors(0);
  
  while ( fscanf(finTxtFile, "%3d", &eta) != EOF )
    {
      flag = fscanf(finTxtFile, "%2d", &depth);
      flag = fscanf(finTxtFile, "%10s", det);
      flag = fscanf(finTxtFile, "%lf", &cellFactor);
      flag = fscanf(finTxtFile, "%x", &detId);
      factors.insert( std::pair<unsigned int, double>(detId, cellFactor) );
      nReadFactors++;
      if ( dbg > 0 ) 
	std::cout << "  " << std::dec << cellFactor
		  << "  " << std::hex << detId << std::endl; 
    }

  std::cout << std::dec << nReadFactors << " factors read from file "
	    << txtFileName
	    << std::endl;
  
  return true;
}
//**********************************************************
