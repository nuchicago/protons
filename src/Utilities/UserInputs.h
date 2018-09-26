
#ifndef USERINPUTS_H
#define USERINPUTS_H

#include "ROOTinclude.h"

class UserInputs {

 public:

  UserInputs( );
  UserInputs( char* jobOpts );


  //--------------------------------------------
  // member variables
  //--------------------------------------------

  //== input files
  char *inputFiles;
  bool inputFilesSet;
  char *correctionFile;
  bool correctionFileSet;
  bool isData;
  bool isMC;

  //== output files
  char *rootOutputFile;
  bool rootOutputFileSet;
  char *psOutputFile;
  bool psOutputFileSet;
  char *outputFileName;
  bool outputFileNameSet;
  char *SelEventList;
  bool SelEventListSet;
  char *plotIndividual;
  bool plotIndividualSet;

  int printMod;

  int verbose;
  int BendDirFilter;
  int numberIterations;
  bool multiplePass;
  char *correctionMethod;
  bool correctionMethodSet;
  char *unsmearingMethod;
  bool unsmearingMethodSet;

  //== data set description
  char *Model;
  bool modelSet;
  long numEventsToProcess;
  bool numEventsToProcessSet;
  int  potTarget;
  bool potTargetSet;
  int  potEmpty;
  bool potEmptySet;
  int RawWireVar;

  char *ParticleType;
  bool isProton;
  bool isPion;
  bool isKaon;
  

  //== BeamSelect Values
  bool zBeamCutoffSet;
  double zBeamCutoff;
  bool zTPCCutoffSet;
  double zTPCCutoff;
  bool zProjCutSet;
  double zProjCut;

  double rCircleCut;
  bool rCircleCutSet;
  double MassCutMin;
  bool MassCutMinSet;
  double MassCutMax;
  bool MassCutMaxSet;

  double pionMassCutMin;
  bool pionMassCutMinSet;
  double pionMassCutMax;
  bool pionMassCutMaxSet;

  double ThetaCut;
  bool ThetaCutSet;
  double PhiCut;
  bool PhiCutSet;

  int applyMassCut;
  int pickyTracksWC;

  // pion inclusive analysis cuts

  int pionCuts;
  double alphaCut;
  bool alphaCutSet;
  double pcPileupDist;
  bool pcPileupDistSet;
  double numPileupCut;
  bool numPileupCutSet;
  double numTracksShower;
  bool numTracksShowerSet;
  double lenTracksShower;
  bool lenTracksShowerSet;
  double tofMin;
  double tofMinSet;
  double tofMax;
  double tofMaxSet;



  //== EventSelector Values

  double zSlabSize;
  bool zSlabSizeSet;
  double branchMaxDist;
  bool branchMaxDistSet;
  double dedxNoBraggMax;
  bool dedxNoBraggMaxSet;
  double clusterMaxDist;
  bool clusterMaxDistSet;




    
  //== target parameters
  char   *target;
  bool   targetSet;
  double targetThickness;
  bool   targetThicknessSet;
  double targetDensity;
  bool   targetDensitySet;
  double targetStartZ;
  bool   targetStartZSet;
  double targetStopZ;
  bool   targetStopZSet;
  double targetRadius;
  bool   targetRadiusSet;
  
  //== Correction Files
  char *reconEffFile;
  bool reconEffFileSet;
  char *absorptionFilePions;
  bool absorptionFilePionsSet;
  char *absorptionFileProtons;
  bool absorptionFileProtonsSet;
  char *absorptionFileKaons;
  bool absorptionFileKaonsSet;
  char *tertiariesFilePions;
  bool tertiariesFilePionsSet;
  char *tertiariesFileMuons;
  bool tertiariesFileMuonsSet;
  char *tertiariesFileProtons;
  bool tertiariesFileProtonsSet;
  char *unfoldingFile1;
  bool unfoldingFile1Set;
  char *unfoldingFile2;
  bool unfoldingFile2Set;
  char *unfoldingAngleFile;
  bool unfoldingAngleFileSet;
  char *pidMatrixFile;
  bool pidMatrixFileSet;
  char *electronVetoFilePions;
  bool electronVetoFilePionsSet;
  char *electronVetoFileProtons;
  bool electronVetoFileProtonsSet;
  char *kaonMatrixFile;
  bool kaonMatrixFileSet;
  char *tofPdfFile;
  bool tofPdfFileSet;
  char *tofPdfFileBins;
  bool tofPdfFileBinsSet;
  char *ckovPdfFile;
  bool ckovPdfFileSet;
  char *piOutFile;
  bool piOutFileSet;
  char *prOutFile;
  bool prOutFileSet;
  char *piOutInPrFile;
  bool piOutInPrFileSet;
  char *prOutInPiFile;
  bool prOutInPiFileSet;

  //== Reconstruction Selection
  char *trackType;
  bool trackTypeSet;
  bool momScaleShift;

  //== Event Selection Cuts
  int    beamidCut;                 
  bool   beamidCutSet;
  double mwpcTargRadiusCut;         
  bool   mwpcTargRadiusCutSet;
  double mwpcTargAngleCut;          
  bool   mwpcTargAngleCutSet;
  double t0Cut;                     
  bool   t0CutSet;
 
  //== Track Selection Cuts
  bool   thxMinCutSet;
  double thxMinCut;
  bool   thxMaxCutSet;
  double thxMaxCut;
  bool   thyMinCutSet;
  double thyMinCut;
  bool   thyMaxCutSet;
  double thyMaxCut;
  int    hitsRoadNdc1Cut; 
  bool   hitsRoadNdc1CutSet;
  int    hitsRoadNdc2Cut;
  bool   hitsRoadNdc2CutSet; 
  int    ndc1_c1dCut;
  bool   ndc1_c1dCutSet;
  double ltofCut;
  bool   ltofCutSet;
  bool   tofSelect;
  bool   betaCorrect;
  double betaCut;
  bool   betaCutSet;
  double tofChi2Cut;
  bool   tofChi2CutSet;
  double tofPulseHeightCut;
  bool   tofPulseHeightCutSet;
  double chargeCut;
  bool   chargeCutSet;
  double vertex4RadiusCut;
  bool   vertex4RadiusCutSet;
  double vertex1RadiusCut;
  bool   vertex1RadiusCutSet;

  //== Efficiency Calculation Tracks Selection Cuts
  int    ndc_c1dCut;
  bool   ndc_c1dCutSet;
  double chi2Vertex2MatchCut;
  bool   chi2Vertex2MatchCutSet;
  double chi2Vertex4MatchCut;
  bool   chi2Vertex4MatchCutSet;
  int    hitsRoadNdc5Cut;
  bool   hitsRoadNdc5CutSet;

  //== Unsmearing matrix calculation
  bool   addNonGaussianTailsToMC;
  bool   addNonGaussianTailsToMCSet;
  double fractionToTails;
  bool   fractionToTailsSet;
  double increasePrecPtrueDiffScaleFactor;
  bool   increasePrecPtrueDiffScaleFactorSet;

  //== PID cuts
  char   *pidMethod;
  bool   pidMethodSet;
  double pionProbCut;
  bool   pionProbCutSet;
  double protonProbCut;
  bool   protonProbCutSet;
  bool   includeKaons;
  double kaonMinP;
  bool   kaonMinPSet;
  double kaonMaxP;
  bool   kaonMaxPSet;
  double kaonProbCut;
  bool   kaonProbCutSet;
  char   *pionPriorFile;
  bool   pionPriorFileSet;
  char   *protonPriorFile;
  bool   protonPriorFileSet;
  char   *kaonPriorFile;
  bool   kaonPriorFileSet;
  bool   useElectronVeto;
  double electronVetoThreshold;
  bool   electronVetoThresholdSet;
  double electronCkovCut;
  bool   electronCkovCutSet;
  double useCkovThreshold;
  bool   useCkovThresholdSet;
  double NpeCut;
  bool   NpeCutSet;
  char   *pdfReadFormat;
  bool   pdfReadFormatSet;
  int    pionNumGaus;
  bool   pionNumGausSet;
  int    protonNumGaus;
  bool   protonNumGausSet;
  int    kaonNumGaus;
  bool   kaonNumGausSet;

  //== Mirror matching cuts
  double fractionCoincidencesCut;
  bool   fractionCoincidencesCutSet;
  int    mndc1nCut;
  bool   mndc1nCutSet;
  int    mndc2nCut;
  bool   mndc2nCutSet;
  double mpCut;
  bool   mpCutSet;
  double mthxCut;
  bool   mthxCutSet;
  double mthyCut;
  bool   mthyCutSet;
  double pDiffCut;
  bool   pDiffCutSet;
  double thDiffCut;
  bool   thDiffCutSet;

  //== Unisim flags
  bool   reconUniSim;
  bool   absorbUniSimStat;
  bool   absorbUniSimSyst;
  bool   tertiariesUniSimStat;
  bool   tertiariesUniSimSyst;
  double absorptionSystError;
  bool   absorptionSystErrorSet;
  double tertiariesSystError;
  bool   tertiariesSystErrorSet;
  bool   momShiftUniSim;
  double momShift;
  bool   momShiftSet;
  bool   migrationUniSimStat;
  bool   electronVetoUniSim;
  bool   pidMatrixUniSim;

  //== Error Analysis Parameters
  char *centralValueXsecFile;
  bool centralValueXsecFileSet;

  int numErrMats, numMigUniSims;
  char* errmatName[20];
  int   numUniSimFiles[20];
  char* uniSimFileName[20];
  char* migrationUniSimName[20];

  double analysis_p_min;
  double analysis_p_max;
  double analysis_th_min;
  double analysis_th_max;

  bool errorsAboutMean;
  double emptySubtractionError;

 private:
  
  //--------------------------------------------
  // function signatures
  //--------------------------------------------

  void initialize( );

  void readJobOptions( ifstream *jobOptionsFile );

  void readIoFiles( ifstream *jobOptionsFile );
  void readDataSetParams( ifstream *jobOptionsFile );
  void readAnalysisCuts( ifstream *jobOptionsFile );
  void readCorrectionFiles( ifstream *jobOptionsFile );
  void readToyMcFlags( ifstream *jobOptionsFile );
  void readErrorAnalysisParams( ifstream *jobOptionsFile );

  void printUserInputs( );

  bool paramLookUp( ifstream *file, char* name );


};

#endif

 
