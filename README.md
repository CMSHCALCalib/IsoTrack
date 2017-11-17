
Isotrack calibration workflow

### Prepare the environment
a) setup CMSSW

    cmsrel CMSSW_9_4_0
    cd CMSSW_9_4_0/src
    cmsenv
    git cms-init
    git clone https://github.com/CMSHCALCalib/IsoTrack CMSHCALCalib/IsoTrack
    scram b -j4

### Preselection from RECO samples:
the isotrack code makes use of some of the utilities available in the CMSSW package `Calibration/IsolatedParticles` and it is organized in:

a) ntuplizer

b) macros

### Test the ntuplizer:

    cd CMSHCALCalib/IsoTrack/test
    #modify the GT, input and output parameters
    cmsRun proto_runIsoTrackCalibration_cfg.py

### Apply (raddam) corrections and loose selection conditions

a) use function `modifyTree()` from ROOT macro `applyRaddamCorrections.C`

b) loose selection conditions are set as constants:

    CHARISO = 10;    //charge isolation
    MINMOM = 40;     //min track momentum
    MAXMOM = 60;     //max track momentum
    MINPT = 7;       //min track pt
    MINEHCAL = 10;   //min Hcal energy   
    MAXEECAL = 1;    //max Ecal energy
    MAXRESP = 3;     //max response

c) two header files are necessary for two different inputs
    ModCalibTree.hh for RECO
    AlcaCalibTree.hh for ALCARECO

d) function arguments (defaults are indicated):

    applyCorrections = true                             /// one can apply correction or not
    txtFileName = "CorrFactor_rechit_HE_SeptChByCh.txt" /// name of input file with corrections
    nblock = 13                                         /// number of run blocks with corrections
    applyCorrectionForHEP17 = false                     /// apply additional correction factor for HEP17
    newFileNamePrefix = "dat2017BCrd"                   /// name prefix of new file
    dataSetNumber = 1                                   /// hardcoded paths to datasets:
   		   				       /// 1 - Bv1 (Reco); 2 - Bv2 (Alcareco); 3 - Cv1 (Reco);
						       /// 4 - Cv2 (Reco); 5 - Cv3 (Alcareco)


e) output files contain ROOT tree named `CalibTree`
 
### Iterative procedure to extract correction factors

a) use function `runIterations()` from ROOT macro `extIsoTrackCalibration.C`

b) the following values can be set as constants:

    NVERTEX_MIN = 1;           // min number of vertices in event
    NVERTEX_MAX = 100;         // max number of vertices in event
    MAX_CALIBRATED_IETA = 23;  // max calibrated cell ieta
    MAX_TRACK_IETA = 23;       // max number of track ieta
    MIN_RUNNUM = 297494;       // min run number 
    MAX_RUNNUM = 900000;       // max run number 
    SINGLE_REFERENCE_RESPONSE = false;  // 3-range or single reference response 
                                       // first range starts from 1
    FIRST_IETA_TR = 15;                 // start of the second range
    FIRST_IETA_HE = 21;                 // start of the third range
	      
one can merge depths by commenting and uncommenting the corresponding lines
one can merge ieta rings by commenting and uncommenting the corresponding lines

    MIN_N_TRACKS_PER_CELL = 50;     // min number of tracks per cell to get calibration factor
    MIN_N_ENTRIES_FOR_FIT = 150;    // min number of entries in histogram for fitting
    MAX_REL_UNC_FACTOR = 0.2;       // max relative uncertainty to claim the cell calibrated
    UNCALIB_FACT_TO_ONE = 0;        // treatment of uncalibrated cell factors
    			  	    // 0: extrapolate behind MAX_CALIBRATED_IETA
                                    // 1: set factor=1 for uncalibrated cells
    APPLY_CORRECTION_FOR_PU = true;  // apply correction for pileup

c) function arguments:

    inFileDir = “.”                         // input file directory
    inFileNamePrefix = "outputFromAnalyzer" // input file name prefix
                                           // input file numbering schema is prefix_<n>; nmin<=n<=nmax
    firstInputFileEnum = 0                  // nmin
    lastInputFileEnum = 1                   // nmax
    preselectedSample = false               // (true) use preselected (step 1 above)
                                           (false) use after loose selections (step 2 above)  
    maxNumberOfIterations = 1               // max number of iterations
    l3prefix5 = "Dat17"                     // output file prefix     
    minHcalEnergy = 10.0                    // min energy in Hcal (to reject mips)
    minPt = 7.0                             // min track pt
    limitForChargeIsolation = 10.0          // charge isolation
                                            <0 - flexible selections
    minTrackMomentum = 40.0                 // min track momentum 
    maxTrackMomentum = 60.0                 // max track momentum
    limitForMipInEcal = 1.0                 // max energy in Ecal
    shiftResponse = 1                       // take reference response form mpv (1) or sample mean (0) 
    subSample = 2                           // can take subsample: odd (0), even(1) or all(2) events
    isCrosscheck = false                    // no iterations but crosscheck
    inTxtFilePrefix = "test"                // text file name with factors for crosscheck 
    treeDirName = "IsoTrackCalibration"     // tree directory (used for preselected)
    treeName = "CalibTree"                  // tree name
    Debug = 0                               // 0-no debugging; 1-short debug; 
                                           >1 - number of events to be shown in detail


d) output root file contains histograms for initial and final distributions, including individual histograms at each ieta, and graphs of mpv and factors vs ieta for each iteration;
output txt file contains correction factors in the format required for DB.   

