//#############################################################################
//
// UserInputs.cpp - 
//
//
// D. Schmitz  June 18, 2005
//
//#############################################################################


#include "UserInputs.h"


//=============================================================================
// Constructors
//=============================================================================


UserInputs::UserInputs( char* jobOpts ){

  cout << endl;
  cout << "***************************************" << endl;
  cout << "*                                     *" << endl;
  cout << "*    Reading JobOptions File Inputs   *" << endl;
  cout << "*                                     *" << endl;
  cout << "***************************************" << endl;
  cout << endl;

  initialize( );

  ifstream *jobOptionsFileList = new ifstream( jobOpts );
  if( !jobOptionsFileList->is_open() ){
    cout << "Could not locate jobOptions file!!!!" << endl << endl;
    exit(0);
  }
    
  ifstream *jobOptionsFile = new ifstream( );

  char nextfile[256];
  jobOptionsFileList->getline( nextfile, 256 );

  do
  {   
    
    if( strcmp( nextfile, "" ) != 0 )
    {

      jobOptionsFile->open( nextfile );
      
      if( jobOptionsFile->is_open() ){
	cout << "Attempting to read jobOptions input from " << nextfile << endl << endl;
	readJobOptions( jobOptionsFile );
	jobOptionsFile->close();
      }
      else
	jobOptionsFile->clear();
    }
    
    jobOptionsFileList->getline( nextfile, 256 );

  }while( !jobOptionsFileList->eof() );



  //== read the given file last, so that the inputs overwrite inputs from files
  jobOptionsFileList->clear();
  readJobOptions( jobOptionsFileList );

  jobOptionsFileList->close();

  printUserInputs( );
    
  cout << endl;
  cout << "***************************************" << endl;
  cout << "*     Finished reading JobOptions     *" << endl;
  cout << "***************************************" << endl;
  cout << endl;


}


UserInputs::UserInputs( ){

}


//=============================================================================
// initialize( ) -- set all booleans to false
//=============================================================================
void UserInputs::initialize( ){

  inputFilesSet = false;
  correctionFileSet = false;
  rootOutputFileSet = false;
  psOutputFileSet = false;
  beamCharFileSet = false;
  haloCharFileSet = false;
  outputFileNameSet = false;
  SelEventListSet = false;
  MultiMatchEventListSet = false;
  modelSet = false;
  plotIndividualSet = false;
  RawWireVar = 0;

  printMod = 50000;

  isData = false;
  isMC = false;
  numEventsToProcessSet = false;
  numEventsToProcess = 10;
  numCenteringEventsSet = false;
  numCenteringEvents = 10;

  zSlabSize = 0.5;
  zSlabSizeSet = false;
  branchMaxDist = 99.;
  branchMaxDistSet = false;
  dedxNoBraggMax = 13.;
  dedxNoBraggMaxSet = false;
  clusterMaxDist =  99.;
  clusterMaxDistSet = false;

  zBeamCutoffSet = false;
  zBeamCutoff = 2;
  zTPCCutoffSet = false;
  zTPCCutoff = 2;
  zProjCutSet = false;
  zProjCut = 99;
  rCircleCutSet = false;
  rCircleCut = 5;
  MassCutMinSet = false;
  MassCutMin = 0;
  MassCutMaxSet = false;
  MassCutMax = 10000;

  pionMassCutMinSet = false;
  pionMassCutMin = 0;
  pionMassCutMaxSet = false;
  pionMassCutMax = 10000;

  pionCuts = 0;
  alphaCut = 99;
  alphaCutSet = false;
  pcPileupDist = 14;
  pcPileupDistSet = false;
  numPileupCut = 4;
  numPileupCutSet = false;
  numTracksShower = 2;
  numTracksShowerSet = false;
  lenTracksShower = 5;
  lenTracksShowerSet = false;
  tofMin = 0;
  tofMinSet = false;
  tofMax = 200;
  tofMaxSet = false;

  PhiCut = 10;
  PhiCutSet = false;
  ThetaCut = 10;
  ThetaCutSet = false;

  applyMassCut = 1;
  skipEventSelection = 0;
  verbose = 0;
  uniqueMatch = 0;
  BendDirFilter = 0;
  pickyTracksWC = 1;
  qualityTracksWC = 0;

  beamLength = 665.2; //cm
  beamLengthSet = false;
  tofOffset = 0;
  tofOffsetSet = false;

  beamPlanesSet = true;

  isProton = false;
  isPion = false;
  isKaon = false;

}



//=============================================================================
// readJobOptions()
//=============================================================================
void UserInputs::readJobOptions( ifstream *jobOptionsFile ){
  
  readIoFiles( jobOptionsFile );
  readDataSetParams( jobOptionsFile );
  readAnalysisCuts( jobOptionsFile );
  readCorrectionFiles( jobOptionsFile );
  readToyMcFlags( jobOptionsFile );
  readErrorAnalysisParams( jobOptionsFile );

}


//=============================================================================
// readIoFiles( )
//=============================================================================
void UserInputs::readIoFiles( ifstream *jobOptionsFile ){

  if( paramLookUp( jobOptionsFile, (char*)"inputFiles" ) ){
    inputFiles = new char[1000];
    *jobOptionsFile >> inputFiles;
    inputFilesSet = true;
  } 

  if( paramLookUp( jobOptionsFile, (char*)"correctionFile" ) ){
    correctionFile = new char[1000];
    *jobOptionsFile >> correctionFile;
    correctionFileSet = true;
  } 
 
  if( paramLookUp( jobOptionsFile, (char*)"rootOutputFile" ) ){
    rootOutputFile = new char[1000];
    *jobOptionsFile >> rootOutputFile;
    rootOutputFileSet = true;
  } 

  if( paramLookUp( jobOptionsFile, (char*)"SelEventList" ) ){
    SelEventList = new char[1000];
    *jobOptionsFile >> SelEventList;
    SelEventListSet = true;
  } 
  if( paramLookUp( jobOptionsFile, (char*)"MultiMatchEventList" ) ){
    MultiMatchEventList = new char[1000];
    *jobOptionsFile >> MultiMatchEventList;
    MultiMatchEventListSet = true;
  } 

  if( paramLookUp( jobOptionsFile, (char*)"psOutputFile" ) ){
    psOutputFile = new char[1000];
    *jobOptionsFile >> psOutputFile;
    psOutputFileSet = true;
  }
  if( paramLookUp( jobOptionsFile, (char*)"beamCharFile" ) ){
    beamCharFile = new char[1500];
    *jobOptionsFile >> beamCharFile;
    beamCharFileSet = true;
  }

    if( paramLookUp( jobOptionsFile, (char*)"haloCharFile" ) ){
    haloCharFile = new char[1500];
    *jobOptionsFile >> haloCharFile;
    haloCharFileSet = true;
  }


  if( paramLookUp( jobOptionsFile, (char*)"outputFileName" ) ){
    outputFileName = new char[1000];
    *jobOptionsFile >> outputFileName;
    outputFileNameSet = true;
  }

    if( paramLookUp( jobOptionsFile, (char*)"plotIndividual" ) ){
    plotIndividual = new char[1000];
    *jobOptionsFile >> plotIndividual;
    plotIndividualSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"printMod" ) ){
    *jobOptionsFile >> printMod;
  }

}


//=============================================================================
// readDataSetParams( )
//=============================================================================
void UserInputs::readDataSetParams( ifstream *jobOptionsFile ){

  char buffer[1000];

  if( paramLookUp( jobOptionsFile, (char*)"DataMC" ) ){
    *jobOptionsFile >> buffer;
    if( strcmp( buffer, "DATA" ) == 0 || strcmp( buffer, "Data" ) == 0 || strcmp( buffer, "data" ) == 0 )
      isData = true;
    if( strcmp( buffer, "MC" ) == 0 || strcmp( buffer, "Mc" ) == 0 || strcmp( buffer, "mc" ) == 0 )
      isMC = true;
  }
  char PIDbuffer[1000];
  if( paramLookUp( jobOptionsFile, (char*)"ParticleType" ) ){
    *jobOptionsFile >> PIDbuffer;
    if( strcmp( PIDbuffer, "PROTON" ) == 0 || strcmp( PIDbuffer, "Proton" ) == 0 || strcmp( PIDbuffer, "proton" ) == 0 )
      isProton = true;
    if( strcmp( PIDbuffer, "KAON" ) == 0 || strcmp( PIDbuffer, "Kaon" ) == 0 || strcmp( PIDbuffer, "kaon" ) == 0 )
      isKaon = true;
    if( strcmp( PIDbuffer, "PION" ) == 0 || strcmp( PIDbuffer, "Pion" ) == 0 || strcmp( PIDbuffer, "pion" ) == 0 )
      isPion = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"events" ) ){
    *jobOptionsFile >> numEventsToProcess;
    numEventsToProcessSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"centeringEvents" ) ){
    *jobOptionsFile >> numCenteringEvents;
    numCenteringEventsSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"Model" ) ){
    Model = new char[1000];
    *jobOptionsFile >> Model;
    modelSet = true;
  } 

  if( paramLookUp( jobOptionsFile, (char*)"RawWireVar" ) ){
    *jobOptionsFile >> RawWireVar;
  }



  //if( paramLookUp( jobOptionsFile, (char*)"target" ) ){
  //  target = new char[20];
  //  *jobOptionsFile >> target;
  //  targetSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"targetThickness" ) ){
  //  *jobOptionsFile >> targetThickness;
  //  targetThicknessSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"targetDensity" ) ){
  //  *jobOptionsFile >> targetDensity;
  //  targetDensitySet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"targetStartZ" ) ){
  //  *jobOptionsFile >> targetStartZ;
  //  targetStartZSet = true;
  //}
  //
  //if( paramLookUp( jobOptionsFile, (char*)"targetStopZ" ) ){
  //  *jobOptionsFile >> targetStopZ;
  //  targetStopZSet = true;
  //}  

  //if( paramLookUp( jobOptionsFile, (char*)"targetRadius" ) ){
  //  *jobOptionsFile >> targetRadius;
  //  targetRadiusSet = true;
  //}
 
  //if( paramLookUp( jobOptionsFile, (char*)"potTarget" ) ){
  //  *jobOptionsFile >> potTarget;
  //  potTargetSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"potEmpty" ) ){
  //  *jobOptionsFile >> potEmpty;
  //  potEmptySet = true;
  //}

  if( paramLookUp( jobOptionsFile, (char*)"verbose" ) ){
    *jobOptionsFile >> verbose;
  }
  if( paramLookUp( jobOptionsFile, (char*)"pionCuts" ) ){
    *jobOptionsFile >> pionCuts;
  }
  if( paramLookUp( jobOptionsFile, (char*)"applyMassCut" ) ){
    *jobOptionsFile >> applyMassCut;
  }
  if( paramLookUp( jobOptionsFile, (char*)"skipEventSelection" ) ){
    *jobOptionsFile >> skipEventSelection;
  }
  if( paramLookUp( jobOptionsFile, (char*)"pickyTracksWC" ) ){
    *jobOptionsFile >> pickyTracksWC;
  }
  if( paramLookUp( jobOptionsFile, (char*)"qualityTracksWC" ) ){
    *jobOptionsFile >> qualityTracksWC;
  }
  if( paramLookUp( jobOptionsFile, (char*)"BendDirFilter" ) ){
    *jobOptionsFile >> BendDirFilter;
  }

  if (!numCenteringEventsSet){ numCenteringEvents = numEventsToProcess;}


}


//=============================================================================
// readAnalysisCuts( )
//=============================================================================
void UserInputs::readAnalysisCuts( ifstream *jobOptionsFile ){

  //if( paramLookUp( jobOptionsFile, (char*)"numberIterations" ) )
  //  *jobOptionsFile >> numberIterations;

  //if( paramLookUp( jobOptionsFile, (char*)"multiplePass" ) )
  //  *jobOptionsFile >> multiplePass;

  //if( paramLookUp( jobOptionsFile, (char*)"correctionMethod" ) ){
  //  correctionMethod = new char[20];
  //  *jobOptionsFile >> correctionMethod;
  //  correctionMethodSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"unsmearingMethod" ) ){
  //  unsmearingMethod = new char[20];
  //  *jobOptionsFile >> unsmearingMethod;
  //  unsmearingMethodSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"beamidCut" ) ){
  //  *jobOptionsFile >> beamidCut;
  //  beamidCutSet = true;
  //}
  //  
  //if( paramLookUp( jobOptionsFile, (char*)"mwpcTargRadiusCut" ) ){
  //  *jobOptionsFile >> mwpcTargRadiusCut;
  //  mwpcTargRadiusCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"mwpcTargAngleCut" ) ){
  //  *jobOptionsFile >> mwpcTargAngleCut;
  //  mwpcTargAngleCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"t0Cut" ) ){
  //  *jobOptionsFile >> t0Cut;
  //  t0CutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"thxMinCut" ) ){
  //  *jobOptionsFile >> thxMinCut;
  //  thxMinCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"thxMaxCut" ) ){
  //  *jobOptionsFile >> thxMaxCut;
  //  thxMaxCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"thyMinCut" ) ){
  //  *jobOptionsFile >> thyMinCut;
  //  thyMinCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"thyMaxCut" ) ){
  //  *jobOptionsFile >> thyMaxCut;
  //  thyMaxCutSet = true;
  //}


  //if( paramLookUp( jobOptionsFile, (char*)"trackType" ) ){
  //  trackType = new char[20];
  //  *jobOptionsFile >> trackType;
  //  trackTypeSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"momScaleShift" ) ){
  //  *jobOptionsFile >> momScaleShift;
  //}


  //if( paramLookUp( jobOptionsFile, (char*)"hitsRoadNdc2Cut" ) ){
  //  *jobOptionsFile >> hitsRoadNdc2Cut;
  //  hitsRoadNdc2CutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"hitsRoadNdc1Cut" ) ){
  //  *jobOptionsFile >> hitsRoadNdc1Cut;
  //  hitsRoadNdc1CutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"ndc1_c1dCut" ) ){
  //  *jobOptionsFile >> ndc1_c1dCut;
  //  ndc1_c1dCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tofSelect" ) ){
  //  *jobOptionsFile >> tofSelect;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"betaCorrect" ) ){
  //  *jobOptionsFile >> betaCorrect;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"betaCut" ) ){
  //  *jobOptionsFile >> betaCut;
  //  betaCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"ltofCut" ) ){
  //  *jobOptionsFile >> ltofCut;
  //  ltofCutSet = true;
  //}
  //
  //if( paramLookUp( jobOptionsFile, (char*)"tofChi2Cut" ) ){
  //  *jobOptionsFile >> tofChi2Cut;
  //  tofChi2CutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tofPulseHeightCut" ) ){
  //  *jobOptionsFile >> tofPulseHeightCut;
  //  tofPulseHeightCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"chargeCut" ) ){
  //  *jobOptionsFile >> chargeCut;
  //  chargeCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"pidMethod" ) ){
  //  pidMethod = new char[20];  
  //  *jobOptionsFile >> pidMethod;
  //  pidMethodSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"pionProbCut" ) ){
  //  *jobOptionsFile >> pionProbCut;
  //  pionProbCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"protonProbCut" ) ){
  //  *jobOptionsFile >> protonProbCut;
  //  protonProbCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"includeKaons" ) ){
  //  *jobOptionsFile >> includeKaons;
  //}

  //if( includeKaons )
  //{
  //  if( paramLookUp( jobOptionsFile, (char*)"kaonMinP" ) ){
  //    *jobOptionsFile >> kaonMinP;
  //    kaonMinPSet = true;
  //  }
  //  if( paramLookUp( jobOptionsFile, (char*)"kaonMaxP" ) ){
  //    *jobOptionsFile >> kaonMaxP;
  //    kaonMaxPSet = true;
  //  }
  //  if( paramLookUp( jobOptionsFile, (char*)"kaonProbCut" ) ){
  //    *jobOptionsFile >> kaonProbCut;
  //    kaonProbCutSet = true;
  //  }
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"pionPriorFile" ) ){
  //  pionPriorFile = new char[1000];
  //  *jobOptionsFile >> pionPriorFile;
  //  pionPriorFileSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"protonPriorFile" ) ){
  //  protonPriorFile = new char[1000];
  //  *jobOptionsFile >> protonPriorFile;
  //  protonPriorFileSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"kaonPriorFile" ) ){
  //  kaonPriorFile = new char[1000];
  //  *jobOptionsFile >> kaonPriorFile;
  //  kaonPriorFileSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"useElectronVeto" ) ){
  //  *jobOptionsFile >> useElectronVeto;
  //}

  //if( useElectronVeto ){
  //  if( paramLookUp( jobOptionsFile, (char*)"electronVetoThreshold" ) ){
  //    *jobOptionsFile >> electronVetoThreshold;
  //    electronVetoThresholdSet = true;
  //  }
  // 
  //  if( paramLookUp( jobOptionsFile, (char*)"electronCkovCut" ) ){
  //    *jobOptionsFile >> electronCkovCut;
  //    electronCkovCutSet = true;
  //  }
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"useCkovThreshold" ) ){
  //  *jobOptionsFile >> useCkovThreshold;
  //  useCkovThresholdSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"NpeCut" ) ){
  //  *jobOptionsFile >> NpeCut;
  //  NpeCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"pdfReadFormat" ) ){
  //  pdfReadFormat = new char[20];
  //  *jobOptionsFile >> pdfReadFormat;
  //  pdfReadFormatSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"pionNumGaus" ) ){
  //  *jobOptionsFile >> pionNumGaus;
  //  pionNumGausSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"protonNumGaus" ) ){
  //  *jobOptionsFile >> protonNumGaus;
  //  protonNumGausSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"kaonNumGaus" ) ){
  //  *jobOptionsFile >> kaonNumGaus;
  //  kaonNumGausSet = true;
  //}

  if( paramLookUp( jobOptionsFile, (char*)"zBeamCutoff" ) ){
    *jobOptionsFile >> zBeamCutoff;
    zBeamCutoffSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"zTPCCutoff" ) ){
    *jobOptionsFile >> zTPCCutoff;
    zTPCCutoffSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"zProjCut" ) ){
    *jobOptionsFile >> zProjCut;
    zProjCutSet = true;
  }


  if( paramLookUp( jobOptionsFile, (char*)"rCircleCut" ) ){
    *jobOptionsFile >> rCircleCut;
    rCircleCutSet = true;
  }
  if( paramLookUp( jobOptionsFile, (char*)"MassCutMin" ) ){
    *jobOptionsFile >> MassCutMin;
    MassCutMinSet = true;
  }
  if( paramLookUp( jobOptionsFile, (char*)"MassCutMax" ) ){
    *jobOptionsFile >> MassCutMax;
    MassCutMaxSet = true;
  }
  if( paramLookUp( jobOptionsFile, (char*)"tofMin" ) ){
    *jobOptionsFile >> tofMin;
    tofMinSet = true;
  }
  if( paramLookUp( jobOptionsFile, (char*)"tofMax" ) ){
    *jobOptionsFile >> tofMax;
    tofMaxSet = true;
  }


  if( paramLookUp( jobOptionsFile, (char*)"PhiCut" ) ){
    *jobOptionsFile >> PhiCut;
    PhiCutSet = true;
  }
  if( paramLookUp( jobOptionsFile, (char*)"ThetaCut" ) ){
    *jobOptionsFile >> ThetaCut;
    ThetaCutSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"alphaCut" ) ){
    *jobOptionsFile >> alphaCut;
    alphaCutSet = true;
  }
  
  if( paramLookUp( jobOptionsFile, (char*)"pcPileupDist" ) ){
    *jobOptionsFile >> pcPileupDist;
    pcPileupDistSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"numPileupCut" ) ){
    *jobOptionsFile >> numPileupCut;
    numPileupCutSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"numTracksShower" ) ){
    *jobOptionsFile >> numTracksShower;
    numTracksShowerSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"lenTracksShower" ) ){
    *jobOptionsFile >> lenTracksShower;
    lenTracksShowerSet = true;
  }
    if( paramLookUp( jobOptionsFile, (char*)"zSlabSize" ) ){
    *jobOptionsFile >> zSlabSize;
    zSlabSizeSet = true;
  }
  
  if( paramLookUp( jobOptionsFile, (char*)"dedxNoBraggMax" ) ){
    *jobOptionsFile >> dedxNoBraggMax;
    dedxNoBraggMaxSet = true;

  }

  if( paramLookUp( jobOptionsFile, (char*)"clusterMaxDist" ) ){
    *jobOptionsFile >> clusterMaxDist;
    clusterMaxDistSet = true;

  }

  if( paramLookUp( jobOptionsFile, (char*)"branchMaxDist" ) ){
    *jobOptionsFile >> branchMaxDist;
    branchMaxDistSet = true;

  }

  if( paramLookUp( jobOptionsFile, (char*)"beamLength" ) ){
    *jobOptionsFile >> beamLength;
    beamLengthSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"tofOffset" ) ){
    *jobOptionsFile >> tofOffset;
    tofOffsetSet = true;
  }

  if( paramLookUp( jobOptionsFile, (char*)"mg1Angle" ) ){
    *jobOptionsFile >> mg1Angle;
  }
  if( paramLookUp( jobOptionsFile, (char*)"mg2Angle" ) ){
    *jobOptionsFile >> mg2Angle;
  }
  else{beamPlanesSet = false;}
  if( paramLookUp( jobOptionsFile, (char*)"mg1Const" ) ){
    *jobOptionsFile >> mg1Const;
  }
  else{beamPlanesSet = false;}

  if( paramLookUp( jobOptionsFile, (char*)"mg2Const" ) ){
    *jobOptionsFile >> mg2Const;
  }
  else{beamPlanesSet = false;}

  if( paramLookUp( jobOptionsFile, (char*)"uniqueMatch" ) ){
    *jobOptionsFile >> uniqueMatch;
    
  }
  
  //if( paramLookUp( jobOptionsFile, (char*)"chi2Vertex4MatchCut" ) ){
  //  *jobOptionsFile >> chi2Vertex4MatchCut;
  //  chi2Vertex4MatchCutSet = true;
  //} 

  //if( paramLookUp( jobOptionsFile, (char*)"vertex4RadiusCut" ) ){
  //  *jobOptionsFile >> vertex4RadiusCut;
  //  vertex4RadiusCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"vertex1RadiusCut" ) ){
  //  *jobOptionsFile >> vertex1RadiusCut;
  //  vertex1RadiusCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"hitsRoadNdc5Cut" ) ){
  //  *jobOptionsFile >> hitsRoadNdc5Cut;
  //  hitsRoadNdc5CutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"chi2Vertex2MatchCut" ) ){
  //  *jobOptionsFile >> chi2Vertex2MatchCut;
  //  chi2Vertex2MatchCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"ndc_c1dCut" ) ){
  //  *jobOptionsFile >> ndc_c1dCut;
  //  ndc_c1dCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"addNonGaussianTailsToMC" ) ){
  //  *jobOptionsFile >> addNonGaussianTailsToMC;
  //  addNonGaussianTailsToMCSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"fractionToTails" ) ){
  //  *jobOptionsFile >> fractionToTails;
  //  fractionToTailsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"increasePrecPtrueDiffScaleFactor" ) ){
  //  *jobOptionsFile >> increasePrecPtrueDiffScaleFactor;
  //  increasePrecPtrueDiffScaleFactorSet = true;
  //}
  //
  //if( paramLookUp( jobOptionsFile, (char*)"fractionCoincidencesCut" ) ){
  //  *jobOptionsFile >> fractionCoincidencesCut;
  //  fractionCoincidencesCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"mndc1nCut" ) ){
  //  *jobOptionsFile >> mndc1nCut;
  //  mndc1nCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"mndc2nCut" ) ){
  //  *jobOptionsFile >> mndc2nCut;
  //  mndc2nCutSet = true;
  //}
  //
  //if( paramLookUp( jobOptionsFile, (char*)"mpCut" ) ){
  //  *jobOptionsFile >> mpCut;
  //  mpCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"mthxCut" ) ){
  //  *jobOptionsFile >> mthxCut;
  //  mthxCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"mthyCut" ) ){
  //  *jobOptionsFile >> mthyCut;
  //  mthyCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"pDiffCut" ) ){
  //  *jobOptionsFile >> pDiffCut;
  //  pDiffCutSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"thDiffCut" ) ){
  //  *jobOptionsFile >> thDiffCut;
  //  thDiffCutSet = true;
  //}

}


//=============================================================================
// readCorrectionFiles
//=============================================================================
void UserInputs::readCorrectionFiles( ifstream *jobOptionsFile ){

  //if( paramLookUp( jobOptionsFile, (char*)"reconEffFile" ) ){
  //  reconEffFile = new char[1000];
  //  *jobOptionsFile >> reconEffFile;
  //  reconEffFileSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"absorptionFilePions" ) ){
  //  absorptionFilePions = new char[1000];
  //  *jobOptionsFile >> absorptionFilePions;
  //  absorptionFilePionsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"absorptionFileProtons" ) ){
  //  absorptionFileProtons = new char[1000];
  //  *jobOptionsFile >> absorptionFileProtons;
  //  absorptionFileProtonsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"absorptionFileKaons" ) ){
  //  absorptionFileKaons = new char[1000];
  //  *jobOptionsFile >> absorptionFileKaons;
  //  absorptionFileKaonsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tertiariesFilePions" ) ){
  //  tertiariesFilePions = new char[1000];
  //  *jobOptionsFile >> tertiariesFilePions;
  //  tertiariesFilePionsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tertiariesFileMuons" ) ){
  //  tertiariesFileMuons = new char[1000];
  //  *jobOptionsFile >> tertiariesFileMuons;
  //  tertiariesFileMuonsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tertiariesFileProtons" ) ){
  //  tertiariesFileProtons = new char[1000];
  //  *jobOptionsFile >> tertiariesFileProtons;
  //  tertiariesFileProtonsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"unfoldingFile1" ) ){
  //  unfoldingFile1 = new char[1000];
  //  *jobOptionsFile >> unfoldingFile1;
  //  unfoldingFile1Set = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"unfoldingFile2" ) ){
  //  unfoldingFile2 = new char[1000];
  //  *jobOptionsFile >> unfoldingFile2;
  //  unfoldingFile2Set = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"unfoldingAngleFile" ) ){
  //  unfoldingAngleFile = new char[1000];
  //  *jobOptionsFile >> unfoldingAngleFile;
  //  unfoldingAngleFileSet = true;
  //}
  //
  //if( paramLookUp( jobOptionsFile, (char*)"pidMatrixFile" ) ){
  //  pidMatrixFile = new char[1000];
  //  *jobOptionsFile >> pidMatrixFile;
  //  pidMatrixFileSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"electronVetoFilePions" ) ){
  //  electronVetoFilePions = new char[1000];
  //  *jobOptionsFile >> electronVetoFilePions;
  //  electronVetoFilePionsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"electronVetoFileProtons" ) ){
  //  electronVetoFileProtons = new char[1000];
  //  *jobOptionsFile >> electronVetoFileProtons;
  //  electronVetoFileProtonsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"kaonMatrixFile" ) ){
  //  kaonMatrixFile = new char[1000];
  //  *jobOptionsFile >> kaonMatrixFile;
  //  kaonMatrixFileSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tofPdfFile" ) ){
  //  tofPdfFile = new char[1000];
  //  *jobOptionsFile >> tofPdfFile;
  //  tofPdfFileSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tofPdfFileBins" ) ){
  //  tofPdfFileBins = new char[1000];
  //  *jobOptionsFile >> tofPdfFileBins;
  //  tofPdfFileBinsSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"ckovPdfFile" ) ){
  //  ckovPdfFile = new char[1000];
  //  *jobOptionsFile >> ckovPdfFile;
  //  ckovPdfFileSet = true;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"piOutFile" ) ){
  //  piOutFile = new char[200];
  //  *jobOptionsFile >> piOutFile;
  //  piOutFileSet = true;
  //} 

  //  if( paramLookUp( jobOptionsFile, (char*)"prOutFile" ) ){
  //  prOutFile = new char[200];
  //  *jobOptionsFile >> prOutFile;
  //  prOutFileSet = true;
  //} 

  //if( paramLookUp( jobOptionsFile, (char*)"piOutInPrFile" ) ){
  //  piOutInPrFile = new char[200];
  //  *jobOptionsFile >> piOutInPrFile;
  //  piOutInPrFileSet = true;
  //} 

  //if( paramLookUp( jobOptionsFile, (char*)"prOutInPiFile" ) ){
  //  prOutInPiFile = new char[200];
  //  *jobOptionsFile >> prOutInPiFile;
  //  prOutInPiFileSet = true;
  //} 
  

}


//=============================================================================
// readToyMcFlags( )
//=============================================================================
void UserInputs::readToyMcFlags( ifstream *jobOptionsFile ){

  //if( paramLookUp( jobOptionsFile, (char*)"reconUniSim" ) ){
  //  *jobOptionsFile >> reconUniSim;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"absorbUniSimStat" ) ){
  //  *jobOptionsFile >> absorbUniSimStat;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"absorbUniSimSyst" ) ){
  //  *jobOptionsFile >> absorbUniSimSyst;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tertiariesUniSimStat" ) ){
  //  *jobOptionsFile >> tertiariesUniSimStat;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"tertiariesUniSimSyst" ) ){
  //  *jobOptionsFile >> tertiariesUniSimSyst;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"momShiftUniSim" ) ){
  //  *jobOptionsFile >> momShiftUniSim;
  //}
  //
  //if( paramLookUp( jobOptionsFile, (char*)"absorptionSystError" ) ){
  //  absorptionSystErrorSet = true;
  //  *jobOptionsFile >> absorptionSystError;
  //} 
  //
  //if( paramLookUp( jobOptionsFile, (char*)"tertiariesSystError" ) ){
  //  tertiariesSystErrorSet = true;
  //  *jobOptionsFile >> tertiariesSystError;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"momShift" ) ){
  //  *jobOptionsFile >> momShift;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"migrationUniSimStat" ) ){
  //  *jobOptionsFile >> migrationUniSimStat;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"electronVetoUniSim" ) ){
  //  *jobOptionsFile >> electronVetoUniSim;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"pidMatrixUniSim" ) ){
  //  *jobOptionsFile >> pidMatrixUniSim;
  //}

}


//=============================================================================
// readErrorAnalysisParams( )
//=============================================================================
void UserInputs::readErrorAnalysisParams( ifstream *jobOptionsFile ){

  //if( paramLookUp( jobOptionsFile, (char*)"centralValueXsecFile" ) ){
  //  centralValueXsecFile = new char[1000];
  //  *jobOptionsFile >> centralValueXsecFile;
  //  centralValueXsecFileSet = true;
  //} 

  //bool found;
  //bool found_mig;
  //char words[200];
  //int i = 0, j = 0;
  //numErrMats = 0;
  //numMigUniSims = 0;

  //do{
  //  
  //  found = false;

  //  sprintf( words, "errmat%d", i+1 );
  //  if( paramLookUp( jobOptionsFile, words ) ){
  //    
  //    errmatName[i] = new char[200];
  //    uniSimFileName[i] = new char[1000];
  //    *jobOptionsFile >> errmatName[i];

  //    if( strcmp( errmatName[i], "migration_syst" ) == 0 ){
  //  numErrMats++;
  //  do{

  //    found_mig = false;

  //    sprintf( words, "migrationFile%d", j+1 );
  //    if( paramLookUp( jobOptionsFile, words ) ){
  //      migrationUniSimName[j] = new char[1000];
  //      *jobOptionsFile >> migrationUniSimName[j];
  //      found_mig = true;
  //      j++;
  //      numMigUniSims++;
  //    }
  //  }while( found_mig );
  //  numUniSimFiles[i] = numMigUniSims;
  //    }
  //    
  //    else if( strcmp( errmatName[i], "none" ) == 0 || strcmp( errmatName[i], "empty_subtraction" ) == 0 ){
  //  numUniSimFiles[i] = 0;
  //  numErrMats++;
  //    }

  //    else{
  //  sprintf( words, "uniSimFileName%d", i+1 );
  //  if( paramLookUp( jobOptionsFile, words ) )
  //    *jobOptionsFile >> uniSimFileName[i];
  //  else
  //    uniSimFileName[i] = (char*)"none";
  //  sprintf( words, "numUniSimFiles%d", i+1 );
  //  if( paramLookUp( jobOptionsFile, words ) )
  //    *jobOptionsFile >> numUniSimFiles[i]; 
  //  else
  //    numUniSimFiles[i] = 0;
  //  numErrMats++;
  //    }
 
  //    found = true;
  //    i++;
  //  }    
  //  
  //}while( found );

  // 

  //if( paramLookUp( jobOptionsFile, (char*)"analysis_p_min" ) ){
  //  *jobOptionsFile >> analysis_p_min;
  //}
  //if( paramLookUp( jobOptionsFile, (char*)"analysis_p_max" ) ){
  //  *jobOptionsFile >> analysis_p_max;
  //}
  //if( paramLookUp( jobOptionsFile, (char*)"analysis_th_min" ) ){
  //  *jobOptionsFile >> analysis_th_min;
  //}
  //if( paramLookUp( jobOptionsFile, (char*)"analysis_th_max" ) ){
  //  *jobOptionsFile >> analysis_th_max;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"errorsAboutMean" ) ){
  //  *jobOptionsFile >> errorsAboutMean;
  //}

  //if( paramLookUp( jobOptionsFile, (char*)"emptySubtractionError" ) ){
  //  *jobOptionsFile >> emptySubtractionError;
  //}

}


//=============================================================================
// printUserInputs( ) -- print out the final settings
//=============================================================================
void UserInputs::printUserInputs( ){

  cout << endl;
  cout << "Options input summary : " << endl << endl; 
  cout << " Process a maximum of " << numEventsToProcess;
  cout << " Centering with " << numCenteringEvents;
  if( isData ) cout << " data";
  else if( isMC ) { 
    cout << " monte carlo";
    if( modelSet ) cout << " using " << Model << " model ";
  }
  cout << " events." << endl << endl;

  //if( targetSet ){
  //  cout << "-------------- Target ----------------" << endl;
  //  cout << "  target : " << target << endl;
  //  cout << "  target thickness : " << targetThickness << " cm" << endl;
  //  cout << "  target density : " << targetDensity << endl;
  //  cout << "  target position : ( " << targetStartZ << ", " << targetStopZ << " ) mm " << endl;
  //  cout << "  target radius : " << targetRadius << " mm " << endl;
  //  cout << endl << endl;
  //}
  //else{
  //  cout << "-------------- Target -----------------" << endl;
  //  cout << " setting default target parameters: " << endl;
  //  cout << "  target thickness : " << targetThickness << " cm" << endl;
  //  cout << "  target density : " << targetDensity << endl;
  //  cout << "  target position : ( " << targetStartZ << ", " << targetStopZ << " ) mm " << endl;
  //  cout << "  target radius : " << targetRadius << " mm " << endl;
  //  cout << endl << endl;
  //}    
  if( inputFilesSet ){
    cout << "------------- Data Files -------------" << endl;    
    cout << "  Input file list = " << inputFiles << endl;
  }
  //if( potTargetSet )      cout << "      p.o.t = " << potTarget << endl;
  //if( potEmptySet )       cout << "      p.o.t = " << potEmpty << endl;
  cout << endl << endl;
  cout << "------------- Output Files -----------" << endl;
  if( rootOutputFileSet ) cout << "  Root output file = " << rootOutputFile << endl;
  if( SelEventListSet ) cout << "  Selected Event list file = " << SelEventList << endl;
  if( psOutputFileSet)    cout << "  PostScript output file = " << psOutputFile << endl;
  if( outputFileNameSet ) cout << "  Output files name = " << outputFileName << endl;  
  if( plotIndividualSet ) cout << "  Printing individual event plots to = " << plotIndividual << endl;  

  //cout << endl << endl;
  //cout << "------------- Analysis Methods -------" << endl;
  //cout << "  Pass over the data " << (int)multiplePass + 1 << " times" << endl;
  //if( correctionMethodSet ) cout << "  Corrections will be applied according to : " << correctionMethod << endl;
  //if( unsmearingMethodSet ) cout << "  Unsmearing will be applied according to : " << unsmearingMethod << endl;
  //else cout << "  Correction application method not set" << endl; 
  //cout << endl << endl;
  cout << "-------- Beam Particle Selection -----" << endl;
  cout << "  Entering track Z cut = " << zTPCCutoff << " cm" <<  endl;
  cout << "  Entering track circular cut = " << rCircleCut << " cm" <<  endl;
  cout << "  Mass cut range = (" << MassCutMin<< "," << MassCutMax<<") MeV" <<  endl;
  cout << "  Alpha Angle cut: " <<  alphaCut << " degrees" << endl;
  //cout << "  Beam particle radius cut = " << mwpcTargRadiusCut << " cm." << endl;
  //cout << "  Beam particle angle cut = " << mwpcTargAngleCut << " mrad" << endl;
  //cout << "  t0 cut = " << t0Cut << " ns" << endl;
  //cout << endl;

  cout << "------- XS Event Selection -------" << endl;
  cout << "  Slab Size Z = " << zSlabSize << " cm" << endl;
  cout << "  Branch Maximum distance  = " << branchMaxDist << " cm" << endl;
  cout << " Branch Cluster counting distance  = " << clusterMaxDist << " cm" << endl; 
  //if( trackTypeSet ) cout << "  Reconstructing tracks using : " << trackType << endl;
  //else cout << "  Track reconstruction method not set" << endl; 
  //cout << "  Perform momentum scale shift for data: " << momScaleShift << endl;
  //cout << "  Charge selection = " << chargeCut << endl;
  //cout << "  thx minumim fiducial cut = " << thxMinCut << " rad" << endl;
  //cout << "  thx maximum fiducial cut = " << thxMaxCut << " rad" << endl;
  //cout << "  thy minumim fiducial cut = " << thyMinCut << " rad" << endl;
  //cout << "  thy maximum fiducial cut = " << thyMaxCut << " rad" << endl;
  //cout << "  Hits in the road NDC2 cut = " << hitsRoadNdc2Cut << endl;
  //cout << "  Hits in the road NDC1 cut = " << hitsRoadNdc1Cut << endl;
  //cout << "  NDC1 track-hits chi-square cut = " << ndc1_c1dCut << endl;
  //cout << "  VERTEX1 radius cut = " << vertex1RadiusCut << endl;
  //cout << "  Beta cut = " << betaCut << endl;
  //cout << "  Ltof cut = " << ltofCut << " cm " << endl;
  //cout << "  Use careful tof hit selection : " << tofSelect << endl;
  //cout << "  Perform slab based tof correction on data : " << betaCorrect << endl;
  //cout << "  ToF chi2 cut = " << tofChi2Cut << endl;
  //cout << "  ToF pulse height cut = " << tofPulseHeightCut << endl;
  //cout << endl;
  //cout << "----- Track Efficiency Selection -----" << endl;
  //cout << "  Chi2 of VERTEX2 match cut = " << chi2Vertex2MatchCut << endl;
  //cout << "  Chi2 of VERTEX4 match cut = " << chi2Vertex4MatchCut << endl;
  //cout << "  VERTEX4 radius cut = " << vertex4RadiusCut << endl;
  //cout << "  Hits in the road NDC5 cut = " << hitsRoadNdc5Cut << endl;
  //cout << "  Chi2 of VERTEX2 match cut = " << chi2Vertex2MatchCut << endl;
  //cout << "  NDC track-hits chi-square cut = " << ndc_c1dCut << endl;
  //cout << endl;
  //cout << "----- Unsmearing Matrices Generation Options -----" << endl;
  //cout << "  Add non-Gaussian tails to MC = " << addNonGaussianTailsToMC << endl;
  //cout << "  Fraction of tracks to increase Prec-Ptrue diff = " << fractionToTails << endl;
  //cout << "  Prec-Ptrue increase scale factor = " << increasePrecPtrueDiffScaleFactor << endl;
  //cout << endl;
  //cout << "-------------- Particle ID -----------" << endl;
  //if( pidMethodSet )     cout << "  Particle ID Method : " << pidMethod << endl;
  //if( pionProbCutSet )   cout << "  Pion probability cut = " << pionProbCut << endl;
  //if( protonProbCutSet ) cout << "  Proton probability cut = " << protonProbCut << endl;
  //if( pidMethodSet && includeKaons ){
  //  cout << "  Include kaons in analysis from " << kaonMinP << " to " << kaonMaxP << " GeV/c" << endl;
  //  if( kaonProbCutSet ) cout << "  Kaon probability cut = " << kaonProbCut << endl;
  //}else
  //  cout << "  Do not include kaons" << endl;
  //if( pionPriorFileSet )   cout << "  Pion prior file = " << pionPriorFile << endl;
  //if( protonPriorFileSet ) cout << "  Proton prior file = " << protonPriorFile << endl;
  //if( kaonPriorFileSet )   cout << "  Kaon prior file = " << kaonPriorFile << endl;
  //cout << endl;
  //cout << "------------ Mirror Matching ---------" << endl;
  //cout << "  Fraction coincidences cut = " << fractionCoincidencesCut << endl;
  //cout << "  NDC1 number of true hits cut = " << mndc1nCut << endl;
  //cout << "  NDC2 number of true hits cut = " << mndc2nCut << endl;
  //cout << "  Max momentum of good mirror particle = " << mpCut << endl;
  //cout << "  Max thx of good mirror particle = " << mthxCut << endl;
  //cout << "  Max thy of good mirror particle = " << mthyCut << endl;
  //cout << "  Maximum momentum difference between rec and true track = " << pDiffCut << endl;
  //cout << "  Maximum angle difference between rec and true track = " << thDiffCut << endl;  
  //cout << endl;
  //cout << "--------- Detector PDF Files ---------" << endl;
  //if( tofPdfFileSet )              cout << "  Tof PDF file = " << tofPdfFile << endl;
  //if( tofPdfFileBinsSet )          cout << "  Tof PDF file bins = " << tofPdfFileBins << endl;
  //if( ckovPdfFileSet )             cout << "  Ckov PDF file = " << ckovPdfFile << endl;
  //cout << endl;
  //cout << "----- Analysis Correction Files -----" << endl;
  //if( reconEffFileSet )            cout << "  Reconstruction efficiency file = " << reconEffFile << endl;
  //if( absorptionFilePionsSet )     cout << "  Pion absorption file = " << absorptionFilePions << endl;
  //if( absorptionFileProtonsSet )   cout << "  Proton absorption file = " << absorptionFileProtons << endl;
  //if( absorptionFileKaonsSet )     cout << "  Kaon absorption file = " << absorptionFileKaons << endl;
  //if( tertiariesFilePionsSet )     cout << "  Tertiary pion file = " << tertiariesFilePions << endl;
  //if( tertiariesFileMuonsSet )     cout << "  Tertiary muon file = " << tertiariesFileMuons << endl;
  //if( tertiariesFileProtonsSet )   cout << "  Tertiary proton file = " << tertiariesFileProtons << endl;
  //if( unfoldingFile1Set )          cout << "  Momentum unfolding file 1 = " << unfoldingFile1 << endl;
  //if( unfoldingFile2Set )          cout << "  Momentum unfolding file 2 = " << unfoldingFile2 << endl;
  //if( unfoldingAngleFileSet )      cout << "  Angular unfolding file 1 = " << unfoldingAngleFile << endl;
  //if( pidMatrixFileSet )           cout << "  p/pi PID matrix file = " << pidMatrixFile << endl;
  //if( electronVetoFilePionsSet )   cout << "  Electron veto pion correction file = " << electronVetoFilePions << endl;
  //if( electronVetoFileProtonsSet ) cout << "  Electron veto proton correction file = " << electronVetoFileProtons << endl;
  //if( kaonMatrixFileSet )          cout << "  Kaon matrix file = " << kaonMatrixFile << endl;
  //cout << endl;
  //cout << "------- Unisim for Systematics ------" << endl;
  //cout << "  Perform reconstruction efficiency unisim : " << reconUniSim << endl;
  //cout << "  Perform absorbtion unisim (stat) : " << absorbUniSimStat << endl;
  //cout << "  Perform absorbtion unisim (syst) : " << absorbUniSimSyst << endl;
  //cout << "  Perform tertiaries unisim (stat) : " << tertiariesUniSimStat << endl;
  //cout << "  Perform tertiaries unisim (syst) : " << tertiariesUniSimSyst << endl;
  //cout << "  Perform migration  unisim (stat) : " << migrationUniSimStat << endl;
  //cout << "  Perform electron veto unisim : " << electronVetoUniSim << endl;
  //cout << "  Perform PID matrix unisim : " << pidMatrixUniSim << endl;
  //cout << "  Perform mom shift unisim : " << momShiftUniSim << endl;
  //if( absorbUniSimSyst ) cout << "  Systematic error on absorption correction : " << absorptionSystError << endl;
  //if( tertiariesUniSimSyst ) cout << "  Systematic error on tertiary correction : " << tertiariesSystError << endl;
  //if( momShiftUniSim ) cout << "  Mom shift scale : " << momShift << endl;
  //cout << endl;
  //if( centralValueXsecFileSet ){
  //  cout << "------- Error Analysis Parameters -------" << endl;
  //  cout << " Central value xsec file : " << centralValueXsecFile << endl << endl;
  //  cout << " number of error matrices : " << numErrMats << endl;
  //  for( int i = 0; i < numErrMats; i++ )
  //    cout << "   " << i+1 << " : " << numUniSimFiles[i] << " " << errmatName[i] << " : " 
  //     << uniSimFileName[i] << endl;
  //  cout << endl;
  //  cout << " number of migration unisims : " << numMigUniSims << endl;
  //  for( int i = 0; i < numMigUniSims; i++ )
  //    cout << "   " << i+1 << " : " << migrationUniSimName[i] << endl;
  //  cout << endl;
  //  cout << "Take error about the mean : " << errorsAboutMean << endl;
  //  cout << "Empty target subtraction error : " << 100*emptySubtractionError << "%" << endl;
  //  cout << endl;
  //  cout << "Analysis range : p = (" << analysis_p_min << "," << analysis_p_max 
  //   << ")  th = (" << analysis_th_min << "," << analysis_th_max << ")" << endl;
  //  cout << endl;
  //}

}

  
//=============================================================================
// paramLookUp( )
//=============================================================================
bool UserInputs::paramLookUp( ifstream *file, char* name ){

  char buffer[1000];
  
  if( file->is_open() )
  {

    file->seekg(0, ios::beg);
    
    do
    {
      *file >> buffer;
     
      if( strcmp( buffer, name ) == 0 )
	return true;
 
    }while( !file->eof() );
  
    file->seekg(0, ios::beg);
    file->clear( );
  }

  return false;

}



//#############################################################################
//
// END 
//
//#############################################################################
