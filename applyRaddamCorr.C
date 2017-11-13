//********************************************************
// Get preselected tree, apply raddam corrections
// and loose selections, write new tree into new file
//********************************************************

// new detID format since CMSSW_8_X
const unsigned int PHI_MASK       = 0x3FF;
const unsigned int ETA_OFFSET     = 10;
const unsigned int ETA_MASK       = 0x1FF;
const unsigned int ZSIDE_MASK     = 0x80000;
const unsigned int DEPTH_OFFSET   = 20;
const unsigned int DEPTH_MASK     = 0xF;
const unsigned int DEPTH_SET      = 0xF00000;

// Loose selections
const double CHARISO = 10;
const double MINMOM = 40;
const double MAXMOM = 60;
const double MINPT = 7;
const double MINEHCAL = 10;
const double MAXEECAL = 1;
const double MAXRESP = 3;

//**********************************************************
// Header with CalibTree class definition
//**********************************************************

#include "ModCalibTree.hh"
#include "AlcaCalibTree.hh"

//**********************************************************
//**********************************************************
int modifyTree( bool applyCorrections = true,
	        const char* txtFileName = "CorrFactor_rechit_HE_SeptChByCh.txt",
		const unsigned nblock = 13,
		bool applyCorrectionForHEP17 = false, //additional correction factor for HEP17
	        const char* newFileNamePrefix = "dat2017BCrd",
	        unsigned dataSetNumber = 1 // 1 - Bv1; 2 - Bv2 (from Alcareco)
		                           // 3 - Cv1; 4 - Cv2; 5 - Cv3 (from Alcareco)
	       )
{
  double factorHEP17 = 1;
  if ( applyCorrectionForHEP17 ) factorHEP17 = 1/1.07;
  
  int minRunNumber = 0;
  int maxRunNumber = 900000;
  bool isAlca = false;
  if ( dataSetNumber == 2 || dataSetNumber == 5 ) isAlca = true;

  const unsigned nphi = 72;
  const unsigned neta = 14;
  const unsigned ndepth = 3;
  double factP[neta][nphi][ndepth][nblock];
  double factM[neta][nphi][ndepth][nblock];
  int runBlocks[nblock];
  for ( unsigned jblock = 0; jblock < nblock; jblock++ ) {
    runBlocks[jblock] = 0;
    for ( unsigned jeta = 0; jeta < neta; jeta++ ) {
      for ( unsigned jdepth = 0; jdepth < ndepth; jdepth++ ) {
	for ( unsigned jphi = 0; jphi < nphi; jphi++ ) {
	  factP[jeta][jphi][jdepth][jblock] = 1;
	  factM[jeta][jphi][jdepth][jblock] = 1;
	}
      }
    }
  }
  
//**********************************************************
// Get correction factors from txt file
//**********************************************************
  char header[80]; 
  int flag;
  int runnum;
  int ieta;
  int iphi;
  int depth;
  double cellFactor;
  int k(0);
  char newFileName[80];

  if ( applyCorrections ) {
  if ( !gSystem->Which("./", txtFileName ) )
       return -1;

  FILE* finTxtFile = fopen(txtFileName,"r");

  // 1 header
  flag = fscanf(finTxtFile, "%s", header);
  // run block numbers
  for ( unsigned int i = 0; i < nblock; i++ ) { 
    flag = fscanf(finTxtFile, "%d", &runnum);
    runBlocks[i] = runnum;
    //std::cout << runBlocks[i] << std::endl;
  }
  // 3 headers
  flag = fscanf(finTxtFile, "%s", header);
  flag = fscanf(finTxtFile, "%s", header);
  flag = fscanf(finTxtFile, "%s", header);
 
  while ( fscanf(finTxtFile, "%3d", &ieta) != EOF )
    {
      flag = fscanf(finTxtFile, "%2d", &iphi);
      flag = fscanf(finTxtFile, "%1d", &depth);
      //std::cout << "  " << ieta << "  " << iphi << "  " << depth;
      for ( unsigned int j = 0; j < nblock; j++ ) {
	flag = fscanf(finTxtFile, "%lf", &cellFactor);
	if ( ieta < 0 ) {
	  k = -ieta - 16;
	  factM[k][iphi-1][depth-1][j] = cellFactor;
	  //std::cout << "  " << factM[k][iphi][depth][j];
	}
	else {
	  k = ieta - 16;
	  if ( iphi >= 63 && iphi <= 66 )
	    factP[k][iphi-1][depth-1][j] = cellFactor*factorHEP17;
	  else
	    factP[k][iphi-1][depth-1][j] = cellFactor;
	  //std::cout << "  " << factP[k][iphi][depth][j];
	}
	//std::cout << "  " << cellFactor;
      }
      //std::cout << std::endl; 
    }
  }
  else {
    runBlocks[0] = 0;
    runBlocks[1] = 900000;
  }
//---------------------------------------------------------
// Read old tree from file
//---------------------------------------------------------

  char fnameInput[120];
  char treeName[120];
  if ( isAlca ) sprintf(treeName,"HcalIsoTrkAnalyzer/CalibTree");
  else sprintf(treeName,"IsoTrackCalibration/CalibTree");
  TChain chain(treeName);

  if ( dataSetNumber == 1 ) {
    for ( int ik = 1; ik <= 9; ik++ ) { // 2017Bv1
      sprintf(fnameInput, "/afs/cern.ch/user/o/omarkin/public/data2017B/data2017B_%1d.root", ik);
      //sprintf(fnameInput, "./data2017B_%1d.root", ik);
      if ( !gSystem->Which("./", fnameInput ) ) { // check file availability
      std::cout << "File " << fnameInput << " doesn't exist." << std::endl;
      }
      else {
	chain.Add(fnameInput);
	std::cout << "Add tree from " << fnameInput 
		  << "   total number of entries (tracks): "
		  << chain.GetEntries() << std::endl;
      }
    }
    minRunNumber = 0;      // 297047;
    maxRunNumber = 900000; // 297723;
  }
  else if ( dataSetNumber == 2 ) {  // 2017B (v1+v2) from Shamik    
    for ( int ik = 1; ik <= 1; ik++ ) {
      sprintf(fnameInput,
	      "/afs/cern.ch/work/d/dbhowmik/public/Isotrack/Results_2017/2017B_C/JETHT2017BUSERAWFALSE.root"
	      );
      if ( !gSystem->Which("./", fnameInput ) ) { // check file availability
      std::cout << "File " << fnameInput << " doesn't exist." << std::endl;
      }
      else {
	chain.Add(fnameInput);
	std::cout << "Add tree from " << fnameInput 
		  << "   total number of entries (tracks): "
		  << chain.GetEntries() << std::endl;
      }
    }
    minRunNumber = 298641;
    maxRunNumber = 299329;
  } 
  else if ( dataSetNumber == 3 ) {  // 2017 Cv1
    for ( int ik = 1; ik <= 9; ik++ ) {
      sprintf(fnameInput, "/afs/cern.ch/user/o/omarkin/public/data2017C/v1/data2017C_%1d.root", ik);
      //sprintf(fnameInput, "./data2017C_%1d.root", ik);
      if ( !gSystem->Which("./", fnameInput ) ) { // check file availability
      std::cout << "File " << fnameInput << " doesn't exist." << std::endl;
      }
      else {
	chain.Add(fnameInput);
	std::cout << "Add tree from " << fnameInput 
		  << "   total number of entries (tracks): "
		  << chain.GetEntries() << std::endl;
      }
    }
    minRunNumber = 0;      // 299368;
    maxRunNumber = 900000; // 299649;
  }
  else if ( dataSetNumber == 4 ) {  // 2017 Cv2
    for ( int ik = 11; ik <= 19; ik++ ) {
      sprintf(fnameInput, "/afs/cern.ch/user/o/omarkin/public/data2017C/v2/data2017C_%2d.root", ik);
      //sprintf(fnameInput, "./data2017C_%1d.root", ik);
      if ( !gSystem->Which("./", fnameInput ) ) { // check file availability
      std::cout << "File " << fnameInput << " doesn't exist." << std::endl;
      }
      else {
	chain.Add(fnameInput);
	std::cout << "Add tree from " << fnameInput 
		  << "   total number of entries (tracks): "
		  << chain.GetEntries() << std::endl;
      }
    }
    minRunNumber = 0;      // 300079;
    maxRunNumber = 900000; // 300676; 
  }
  else if ( dataSetNumber == 5 ) {  // 2017C (v1+v2+v3) from Shamik    
    for ( int ik = 1; ik <= 1; ik++ ) {
      sprintf(fnameInput,
	      "/afs/cern.ch/work/d/dbhowmik/public/Isotrack/Results_2017/2017B_C/JETHT2017CUSERAWFALSE.root"
	      );
      if ( !gSystem->Which("./", fnameInput ) ) { // check file availability
      std::cout << "File " << fnameInput << " doesn't exist." << std::endl;
      }
      else {
	chain.Add(fnameInput);
	std::cout << "Add tree from " << fnameInput 
		  << "   total number of entries (tracks): "
		  << chain.GetEntries() << std::endl;
      }
    }
    minRunNumber = 300742;
    maxRunNumber = 302029; 
  }
  else {
    std::cout << "No data set available." << std::endl;
  }

  if ( applyCorrections )
    sprintf(newFileName,"%s_%1d.root", newFileNamePrefix, dataSetNumber);
  else
    sprintf(newFileName,"%sno_%1d.root", newFileNamePrefix, dataSetNumber);

  
  Long64_t nevt = chain.GetEntries();
  
  if ( nevt == 0 ) {
    std:: cout << "Tree is empty." << std::endl;
    return -2;
  }
  
  //--- Initialize tree
  if ( isAlca )
    AlcaCalibTree oldtree(&chain);
  else
    ModCalibTree oldtree(&chain);
  
  // new variables and tree
  // Declaration of leaf types
  Int_t           t_Run;
  Int_t           t_Event;
  Int_t           t_nVtx;
  Int_t           t_nTrk;
  Double_t        t_EventWeight;
  Double_t        t_p;
  Double_t        t_pt;
  Int_t           t_ieta;
  Double_t        t_phi;
  Double_t        t_eMipDR;
  Double_t        t_eHcal;
  Double_t        t_eHcal10;
  Double_t        t_eHcal30;
  Double_t        t_hmaxNearP;
  Bool_t          t_selectTk;
  Bool_t          t_qltyMissFlag;
  Bool_t          t_qltyPVFlag;
  std::vector<unsigned int> *t_DetIds;
  std::vector<double>  *t_HitEnergies;
 
  TFile *newfile = new TFile(newFileName,"recreate");
  TTree *newtree = new TTree("CalibTree", "CalibTree");
  newtree->Branch("t_Run",         &t_Run,         "t_Run/I");
  newtree->Branch("t_Event",       &t_Event,       "t_Event/I");
  newtree->Branch("t_nVtx",        &t_nVtx,        "t_nVtx/I");
  newtree->Branch("t_nTrk",        &t_nTrk,        "t_nTrk/I");
  newtree->Branch("t_EventWeight", &t_EventWeight, "t_EventWeight/D");
  newtree->Branch("t_p",           &t_p,           "t_p/D");
  newtree->Branch("t_pt",          &t_pt,          "t_pt/D");
  newtree->Branch("t_ieta",        &t_ieta,        "t_ieta/I");
  newtree->Branch("t_phi",         &t_phi,         "t_phi/D");
  newtree->Branch("t_eMipDR",      &t_eMipDR,      "t_eMipDR/D");
  newtree->Branch("t_eHcal",       &t_eHcal,       "t_eHcal/D");
  newtree->Branch("t_eHcal10",     &t_eHcal10,     "t_eHcal10/D");
  newtree->Branch("t_eHcal30",     &t_eHcal30,     "t_eHcal30/D");
  newtree->Branch("t_hmaxNearP",   &t_hmaxNearP,   "t_hmaxNearP/D");
  newtree->Branch("t_selectTk",    &t_selectTk,    "t_selectTk/O");
  newtree->Branch("t_qltyMissFlag",&t_qltyMissFlag,"t_qltyMissFlag/O");
  newtree->Branch("t_qltyPVFlag",  &t_qltyPVFlag,  "t_qltyPVFlag/O)");
  newtree->Branch("t_DetIds",       "std::vector<unsigned int>", &t_DetIds);
  newtree->Branch("t_HitEnergies",  "std::vector<double>",       &t_HitEnergies);

  
  Long64_t nentries = oldtree.fChain->GetEntriesFast();
  Long64_t nb = 0;
  for (Long64_t ij = 0; ij < nentries; ij++) {
    Long64_t ientry = oldtree.LoadTree(ij);
    if ( ientry < 0 ) break;   
    nb = oldtree.fChain->GetEntry(ij);   //nbytes += nb;
    
     if ( !((oldtree.t_selectTk)
	    && (oldtree.t_qltyMissFlag)
	    && (oldtree.t_hmaxNearP < CHARISO) //maxCharIso)
	    && (oldtree.t_eMipDR < MAXEECAL) 
	    && (oldtree.t_p > MINMOM) && (oldtree.t_p < MAXMOM)
	    && (oldtree.t_pt >= MINPT)               // constraint on track pt
	    && (oldtree.t_eHcal >= MINEHCAL)            // constraint on Hcal energy
	    && (oldtree.t_eHcal/oldtree.t_p < MAXRESP)
	    )
	  ) continue;

     t_DetIds->clear();
     t_HitEnergies->clear();

     if ( oldtree.t_Run < minRunNumber || oldtree.t_Run > maxRunNumber ) continue;
     
     t_Run = oldtree.t_Run;
     t_Event = oldtree.t_Event;
     t_nVtx = oldtree.t_nVtx;
     t_nTrk = oldtree.t_nTrk;
     t_EventWeight = oldtree.t_EventWeight;
     t_p = oldtree.t_p;
     t_pt = oldtree.t_pt;
     t_ieta = oldtree.t_ieta;
     t_phi = oldtree.t_phi;
     t_eMipDR = oldtree.t_eMipDR;
     //t_eHcal10 = oldtree.t_eHcal10;
     //t_eHcal30 = oldtree.t_eHcal30;
     t_hmaxNearP = oldtree.t_hmaxNearP;
     t_selectTk = oldtree.t_selectTk;
     t_qltyMissFlag = oldtree.t_qltyMissFlag;
     t_qltyPVFlag = oldtree.t_qltyPVFlag;
          
     unsigned int irun;
     if ( oldtree.t_Run < runBlocks[0] ) continue;
     else if ( oldtree.t_Run < runBlocks[1] ) irun = 0;
     else if ( oldtree.t_Run < runBlocks[2] ) irun = 1;
     else if ( oldtree.t_Run < runBlocks[3] ) irun = 2;
     else if ( oldtree.t_Run < runBlocks[4] ) irun = 3;
     else if ( oldtree.t_Run < runBlocks[5] ) irun = 4;
     else if ( oldtree.t_Run < runBlocks[6] ) irun = 5;
     else if ( oldtree.t_Run < runBlocks[7] ) irun = 6;
     else if ( oldtree.t_Run < runBlocks[8] ) irun = 7;
     else irun = 8;

    int id(0);
    double eTotal(0);
    double eTotalNew(0);
      
      unsigned int nDets = (*(oldtree.t_DetIds)).size();      
      for (unsigned int idet = 0; idet < nDets; idet++) {

	unsigned int detId = (*(oldtree.t_DetIds))[idet];
	double e = (*(oldtree.t_HitEnergies))[idet];
	eTotal += e;
	double enew(e);
	double corfact(1);

	iphi = detId & PHI_MASK;
	ieta = (detId>>ETA_OFFSET) & ETA_MASK;
	int zside= (detId&ZSIDE_MASK) ? 1 : -1;
	depth= (detId>>DEPTH_OFFSET) & DEPTH_MASK;
	
	if ( ieta > 15 ) {
	  if ( zside < 0 ) {
	    id = ieta - 16;
	    corfact = factM[id][iphi-1][depth-1][irun];
	  }
	  else {
	    id = ieta - 16;
	    corfact = factP[id][iphi-1][depth-1][irun];
	  }	  
	}
	enew = e*corfact;
	/*
	if ( corfact < 0.9 || corfact > 2 ) {
	  std::cout << zside << "/" << ieta << "/" << iphi << "/" << depth << "/" << irun
		    << "  " << e << "/" << enew << "/" << corfact
		    << std::endl;
	}
	*/
	eTotalNew += enew;
	t_DetIds->push_back(detId);
	t_HitEnergies->push_back(enew);
      }
      t_eHcal = eTotalNew;
      double ctotal(1);
      if ( eTotal > 0 ) ctotal = eTotalNew/eTotal;
      t_eHcal10 = oldtree.t_eHcal10*ctotal;
      t_eHcal30 = oldtree.t_eHcal30*ctotal;
      
      newtree->Fill();   
   }
   //newtree->Print();
  if ( isAlca ) 
    oldtree.~AlcaCalibTree();
  else
    oldtree.~ModCalibTree();
    
   newfile->cd();
   newtree->Write();
   //newtree->AutoSave();
   //delete newfile;
      
   return 0;
}

//-----------------------------------------------------------

