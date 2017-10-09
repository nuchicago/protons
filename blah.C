#define blah_cxx
#include "blah.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// ===================================================================================================================
// ====================================       PUT HISTOGRAMS HERE           ==========================================
// ===================================================================================================================

// ### mc ###

// ### all the angles (in the tpc) ###
TH1D *hMCEveryAngle = new TH1D("hMCEveryAngle", "Every angle (degrees)", 180, 0, 180);
// ### all the interacting angles ###
TH1D *hMCEveryIntAngle = new TH1D("hMCEveryIntAngle", "Every interacting angle", 90, 0, 180);
// ### all non interacting angles ###
TH1D *hMCEveryNonIntAngle = new TH1D("hMCEveryNonIntAngle", "Every non interacting angle", 90, 0, 180);
// ### inelastic angles -_- ###
TH1D *hMCInelasticAngle = new TH1D("hMCInelasticAngle", "Inelastic Angles", 180, 0, 180);
// ### elastic angles -_- ###
TH1D *hMCElasticAngle = new TH1D("hMCElasticAngle", "Elastic Angles", 180, 0, 180);
// ### elastic angles -_- ###
TH1D *hMCCoulombAngle = new TH1D("hMCCoulombAngle", "Coulombs Angles", 180, 0, 180);

// ### daughter info ###
TH1D *hMCNDaughters  = new TH1D("hMCNDaughters", "Number of daughters", 11, 0, 11);
TH1D *hMCDaughterpdg = new TH1D("hMCDaughterpdg", "daughter particle pdg codes", 3000, -500, 2500);
TH1D *hMCNPDaughters = new TH1D("hMCNPDaughters", "Number of proton daughters", 5, 0, 5);
    TH1D *hMC0Pntracks = new TH1D("hMC0Pntracks", "number of tracks reconstructed (0 proton daughter)", 5, 0, 5);
    TH1D *hMC1Pntracks = new TH1D("hMC1Pntracks", "number of tracks reconstructed (1 proton daughter)", 5, 0, 5);
    TH1D *hMC2Pntracks = new TH1D("hMC2Pntracks", "number of tracks reconstructed (2 proton daughter)", 5, 0, 5);
    TH1D *hMC3Pntracks = new TH1D("hMC3Pntracks", "number of tracks reconstructed (3 proton daughter)", 5, 0, 5);


// ### how many interactions per event ###
TH1D *hMCnInteractions = new TH1D("hMCnInteraction", "Number of interactions per event", 5, 0, 5);
// ### for 1 interaction events, how many tracks ###
TH1D *hMCOneIntNTracks = new TH1D("hMCnOneIntNTracks", "Number of tracks per event", 4, 1, 5);
// ### for (1 int, 1 track) events, event type ###
TH1D *hMCOneIntOneTrackInteractionType = new TH1D("hMCOneIntOneTrackInteractions", "Interaction types", 3, 0, 3);

// ### step size distance ###
TH1D *hMCGeant4StepSize = new TH1D("hMCG4StepSize", "Geant 4 Step Size",    45,  0, 45);
TH1D *hMCGeant4dEdX     = new TH1D("hMCG4dEdX",     "Geant 4 dE/dX (hack)", 100, 0, 10);
TH1D *hRecodEdX         = new TH1D("hRecodEdX",     "Reco dE/dX",           100, 0, 50);
TH1D *hRecodX           = new TH1D("hRecodX",       "Reco dX (pitch hit)",  100, 3, 7);

// ### STE ###

// ### all the angles -_- ###
TH1D *hRecoSTEEveryAngle = new TH1D("hRecoSTEEveryAngle", "Every spacepoint angle", 100, 0, 100);
// ### interaction angles ###
TH1D *hRecoSTEEveryAngleInts = new TH1D("hRecoSTEEveryAngleInts", "Every interacting spacepoint angle", 100, 0, 100);
  // ### inelastic ###
TH1D *hRecoSTEEveryAngleIntsInelast = new TH1D("hRecoSTEEveryAngleIntsInelast", "Every Inelast interacting spacepoint angle", 100, 0, 100);
  // ### coulomb  ###
TH1D *hRecoSTEEveryAngleIntsCoulomb = new TH1D("hRecoSTEEveryAngleIntsCoulomb", "Every Coulomb interacting spacepoint angle", 100, 0, 100);
  // ### elastic ###
TH1D *hRecoSTEEveryAngleIntsElastic = new TH1D("hRecoSTEEveryAngleIntsElastic", "Every Elastic interacting spacepoint angle", 100, 0, 100);
// ### non interaction angles ###
TH1D *hRecoSTEEveryAngleNoInts = new TH1D("hRecoSTEEveryAngleNoInts", "Every non interacting spacepoint angle", 100, 0, 100);

// ### DTE ###

// ### reco end to mc end ###
TH1D *hRecoDTETrackEndDelta = new TH1D("hRecoDTETrackEndDelta", "Primary track end distance to 2nd track end", 100, 0, 50);
// ### reco 1-2 angle if close ###
TH1D *hRecoDTETrackAngle = new TH1D("hRecoDTETrackAngle", "Angle between tracks if broken up", 90, 0, 180);
// ### reco 1-2 angle if close interacted ###
TH1D *hRecoDTETrackAngleInts = new TH1D("hRecoDTETrackAngleInts", "Angle between tracks if broken up (interactions)", 90, 0, 180);
// ### reco 1-2 angle if close not interacted ###
TH1D *hRecoDTETrackAngleNoInts = new TH1D("hRecoDTETrackAngleNoInts", "Angle between tracks if broken up (not interactions)", 90, 0, 180);
// ### inelastic ###
TH1D *hRecoDTETrackAngleIntsInelast = new TH1D("hRecoDTETrackAngleIntsInelast", "Angle between tracks if broken up (inelastics)", 90, 0, 180);
// ### elastic ###
TH1D *hRecoDTETrackAngleIntsCoulomb = new TH1D("hRecoDTETrackAngleIntsCoulomb", "Angle between tracks if broken up (coulombs)", 90, 0, 180);
// ### coulomb  ###
TH1D *hRecoDTETrackAngleIntsElastic = new TH1D("hRecoDTETrackAngleIntsElastic", "Angle between tracks if broken up (interactions)", 90, 0, 180);
// ### reco 1 position closest to reco 2 start if not near end ###
TH1D *hRecoDTEDistanceToMCIntFromTrack2 = new TH1D("hRecoDTEDistanceToMCIntFromTrack2", "Distance to mc int pt from track 2 start", 50, 0, 25);
// ### reco 23 avg distance to inelastic scatters ###
TH2D *hRecoTTE23vsg4 = new TH2D("hRecoTTE23g4", "Distance between tracks 2 and 3 vs (<start> - g4 int pt)", 50, 0, 50, 50, 0, 50);

// ### effish and danny devito ###
TH1D *hXSIntAnglePurity = new TH1D("hXSIntAnglePurity", "Purity of Selection by Angle", 9, 0, 9);
TH1D *hXSIntAngleEffish = new TH1D("hXSIntAngleEffish", "Effeciency of Selection by Angle", 9, 0, 9);


// ### some crazy shit ###
TH1D *hMCEnergyLossIntoTPC                  = new TH1D("hMCEnergyLossIntoTPC", "Kinetic Energy Loss", 100, 0, 1000);
TH1D *hMCEnergyLossIntoTPCRatio             = new TH1D("hMCEnergyLossIntoTPCRatio", "Kinetic Energy Loss", 50, 0, 1);
TH1D *hMCNormalizedAngleSignal              = new TH1D("hMCNormalizedAngleSignal",       "Truth Normalized Hadronic Angle Distributions",            180, 0, 180);
TH1D *hMCNormalizedAngleBackground          = new TH1D("hMCNormalizedAngleBackground",   "Truth Normalized non-Hadronic Angle Distributions",        180, 0, 180);
TH2D *hMCInitialInteractingKE               = new TH2D("hMCInitialInteractingKE",  "Initial vs Interacting KE", 500, 0, 1000, 500, 0, 1000);
TH2D *hMCElasticInitialInteractingKE        = new TH2D("hMCElasticInitialInteractingKE",  "Elastic Initial vs Interacting KE", 500, 0, 1000, 500, 0, 1000);
TH2D *hMCInelastInitialInteractingKE        = new TH2D("hMCInelastInitialInteractingKE",  "Inelast Initial vs Interacting KE", 500, 0, 1000, 500, 0, 1000);
TH1D *hRecoNormalizedAngleSignal            = new TH1D("hRecoNormalizedAngleSignal",     "Reconstucted Normalized Hadronic Angle Distributions",     180, 0, 180);
TH1D *hRecoNormalizedAngleBackground        = new TH1D("hRecoNormalizedAngleBackground", "Reconstucted Normalized non-Hadronic Angle Distributions", 180, 0, 180);
THStack *hMCNormalizedAngleStacked          = new THStack("hMCNormalizedAngleStacked",   "Truth Normalized Signal Combined");
THStack *hRecoNormalizedAngleStacked        = new THStack("hRecoNormalizedAngleStacked", "Reconstructed Normalized Signal Combined");
THStack *hRecoDeltaEnergyStacked            = new THStack("hRecoDeltaEnergyStacked",     "Initial KE Minus Total KE Loss (MeV)");
TH1D *hRecoMatchedInteractingKERatio        = new TH1D("hRecoMatchedInteractingKERatio", "Ratio of Interacting KE When Matched", 50, 0, 5);
TH2D *hRecoMatchedInteractingKECmpsn        = new TH2D("hRecoMatchedInteractingKECmpsn",         "Reconstructed vs True KE inelast", 100, 0, 1000, 100, 0, 1000);
TH2D *hRecoMatchedInteractingKECmpsnInelast = new TH2D("hRecoMatchedInteractingKECmpsnInelast",  "Reconstructed vs True KE"        , 100, 0, 1000, 100, 0, 1000);
TH2D *hRecoMatchedElasticAngleCmpsn         = new TH2D("hRecoMatchedElasticAngleCmpsn",   "Reconstructed vs True Elastic Angle", 45, 0, 45, 20, 0, 20);
TH1D *hRecoDeltaEnergy                      = new TH1D("hRecoDeltaEnergy", "Initial KE Minus Total KE Loss (MeV)"           , 50, 0, 1000);
TH1D *hRecoDeltaEnergyChannel0              = new TH1D("hRecoDeltaEnergyChannel0",         "True - Reco Loss Channel 0"     , 100, -200, 200);
TH2D *hRecoTrueVSRecoEnergyChannel0         = new TH2D("hRecoTrueVSRecoEnergyChannel0", "True vs Reco Loss"                 , 100, 0   , 800 , 100, 0, 800);
TH2D *hRecoDeltaEnergyvsIntKE               = new TH2D("hRecoDeltaEnergyvsIntKE", "Delta KE vs Int KE"                      , 100, 0   , 1000, 50, 0, 1000);
TH1D *hRecoDeltaEnergyElastic               = new TH1D("hRecoDeltaEnergyElastic", "Elastic Initial KE Minus Total KE Loss (MeV)"    , 50, 0, 1000);
TH1D *hRecoDeltaEnergyInelast               = new TH1D("hRecoDeltaEnergyInelast", "Inelastic Initial KE Minus Total KE Loss (MeV)"  , 50, 0, 1000);
TH2D *hRecoDeltaEnergyInelastvsIntKE        = new TH2D("hRecoDeltaEnergyInelastvsIntKE", "Delta KE vs Int KE"          , 50, 0, 1000, 50, 0, 1000);
TH1D *hRecoLowMeanDeltaEnergy               = new TH1D("hRecoLowMeanDeltaEnergy", "Initial KE Minus Total KE Loss (MeV)"            , 50, 0, 1000);
TH1D *hRecoLowMeanDeltaEnergySignal         = new TH1D("hRecoLowMeanDeltaEnergySignal", "Initial KE Minus Total KE Loss (MeV)"      , 50, 0, 1000);
TH1D *hRecoLowMeanDeltaEnergyBackground     = new TH1D("hRecoLowMeanDeltaEnergyBackground", "Initial KE Minus Total KE Loss (MeV)"  , 50, 0, 1000);
TH1D *hRecoEndMeandedx                      = new TH1D("hRecoEndMeandedx", "End Mean dE/dX for straight tracks"       , 45, 0, 45);
TH1D *hRecoEndMeandedxElastic               = new TH1D("hRecoEndMeandedxElastic", "End Mean dE/dX for straight tracks (elastics)", 45, 0, 45);
TH1D *hRecoEndMeandedxInelast               = new TH1D("hRecoEndMeandedxInelast", "End Mean dE/dX for straight tracks (inelasts)", 45, 0, 45);
TH1D *hRecoEndMeandedxNeither               = new TH1D("hRecoEndMeandedxNeither", "End Mean dE/dX for straight tracks (neithers)", 45, 0, 45);

// ### xs plots ###
//TH1D *hRecoIncidentKE           = new TH1D("hRecoIncidentKE",    "Incident Kinetic Energy"   , 22, -100, 1000);
//TH1D *hRecoInteractingKE        = new TH1D("hRecoInteractingKE", "Interacting Kinetic Energy", 22, -100, 1000);
TH1D *hRecoIncidentKE           = new TH1D("hRecoIncidentKE",    "Incident Kinetic Energy"   , 20, 0, 1000);
TH1D *hRecoInteractingKE        = new TH1D("hRecoInteractingKE", "Interacting Kinetic Energy", 20, 0, 1000);
    TH1D *hRecoInteractingKEChan1             = new TH1D("hRecoInteractingKEChan1",      "Interacting Kinetic Energy Channel 1"  , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan2             = new TH1D("hRecoInteractingKEChan2",      "Interacting Kinetic Energy Channel 2"  , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan3             = new TH1D("hRecoInteractingKEChan3",      "Interacting Kinetic Energy Channel 3"  , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan4             = new TH1D("hRecoInteractingKEChan4",      "Interacting Kinetic Energy Channel 4"  , 20, 0, 1000);
    TH1D *hRecoInteractingKESignal            = new TH1D("hRecoInteractingKESingal",     "Interacting Kinetic Energy Signal"     , 20, 0, 1000);
    TH1D *hRecoInteractingKESignalElastic     = new TH1D("hRecoInteractingKESingalElastic", "Interacting Kinetic Energy Signal Elastic"     , 20, 0, 1000);
    TH1D *hRecoInteractingKESignalElasticObs  = new TH1D("hRecoInteractingKESingalElasticObs", "Interacting Kinetic Energy Signal Elastic (>6)"     , 20, 0, 1000);
    TH1D *hRecoInteractingKESignalInelast     = new TH1D("hRecoInteractingKESingalInelast", "Interacting Kinetic Energy Signal Inelast"     , 20, 0, 1000);
    TH1D *hRecoInteractingKEBackground        = new TH1D("hRecoInteractingKEBackground", "Interacting Kinetic Energy Background" , 20, 0, 1000);
    TH1D *hRecoInteractingKEBackgroundElastic = new TH1D("hRecoInteractingKEBackgroundElastic", "Elastic Interacting Kinetic Energy Background" , 20, 0, 1000);
    TH1D *hRecoInteractingKEBackgroundInelast = new TH1D("hRecoInteractingKEBackgroundInelast", "Inelastic Interacting Kinetic Energy Background" , 20, 0, 1000);
    TH1D *hRecoInteractingKETrueSignal        = new TH1D("hRecoInteractingKETrueSignal", "Interacting Kinetic Energy True Signal" , 20, 0, 1000);
        TH1D *hRecoInteractingKETrueSignalInelast     = new TH1D("hRecoInteractingKETrueSignalInelast",   "True Signal"     , 20, 0, 1000);
        TH1D *hRecoInteractingKETrueSignalElastic     = new TH1D("hRecoInteractingKETrueSignalElastic",   "Interacting KE True Signal Elastic"       , 20, 0, 1000);
        TH1D *hRecoInteractingKETrueSignalElasticG6   = new TH1D("hRecoInteractingKETrueSignalElasticG6", "Interacting KE True Signal Elastic(>6)"   , 20, 0, 1000);
        TH1D *hRecoInteractingKETrueSignalInelastG6   = new TH1D("hRecoInteractingKETrueSignalInelastG6", "Interacting KE True Signal Inelastic(>6)" , 20, 0, 1000);
        TH1D *hRecoInteractingKETrueSignalG6          = new TH1D("hRecoInteractingKETrueSignalG6",        "Interacting KE True Signal (>6)"          , 20, 0, 1000);
        TH2D *hRecoInteractingKETrueSignalEA = new TH2D("hRecoInteractingKETrueSignalEA", "Interacting KE True Signal Elastic vs Angle", 500, 0, 1000, 20, 0, 20);
        TH2D *hRecoInteractingKETrueSignalIA = new TH2D("hRecoInteractingKETrueSignalIA", "Interacting KE True Signal Inelast vs Angle", 500, 0, 1000, 20, 0, 20);
    THStack *hRecoInteractingKEStacked   = new THStack("hRecoInteractingKEStacked",  "Interacting Kinetic Energy (combined)");
    THStack *hRecoInteractingKEStacked2  = new THStack("hRecoInteractingKEStacked2", "Interacting Kinetic Energy (combined)");
    THStack *hRecoInteractingKEStacked3  = new THStack("hRecoInteractingKEStacked3", "Interacting KE - All Topologies");
    THStack *hRecoInteractingKEStackedP  = new THStack("hRecoInteractingKEStackedP", "Interacting Kinetic Energy (signal and noise)");
    THStack *hRecoInteractingKEStackedE  = new THStack("hRecoInteractingKEStackedE", "Interacting Kinetic Energy (true and observed signal)");
        TH1D *hRecoInteractingKEChan1TrueSignal        = new TH1D("hRecoInteractingKEChan1TrueSignal","Interacting Kinetic Energy Channel 1 True Signal"          ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan1TrueSignalG6      = new TH1D("hRecoInteractingKEChan1TrueSignalG6","Interacting Kinetic Energy Channel 1 True Signal Visible",20, 0, 1000);
        TH1D *hRecoInteractingKEChan1Signal            = new TH1D("hRecoInteractingKEChan1Signal",        "Interacting Kinetic Energy Channel 1 Signal"           ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan1SignalElastic     = new TH1D("hRecoInteractingKEChan1SignalElastic", "Interacting Kinetic Energy Channel 1 Elastic Signal"   ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan1SignalInelast     = new TH1D("hRecoInteractingKEChan1SignalInelast","Interacting Kinetic Energy Channel 1 Inelast Signal"    ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan1Background        = new TH1D("hRecoInteractingKEChan1Background",    "Interacting Kinetic Energy Channel 1 Background"       ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan1BackgroundElastic = new TH1D("hRecoInteractingKEChan1BackgroundElastic","Interacting Kinetic Energy Channel 1 Background"    ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan1BackgroundInelast = new TH1D("hRecoInteractingKEChan1BackgroundInelast","Interacting Kinetic Energy Channel 1 Background"    ,20, 0, 1000);
        THStack *hRecoInteractingKEChan1StackedP       = new THStack("hRecoInteractingKEChan1StackedP",   "Interacting Kinetic Energy Channel 1(signal and noise)");
        THStack *hRecoInteractingKEChan1StackedE       = new THStack("hRecoInteractingKEChan1StackedE",   "Interacting Kinetic Energy Channel 1(true and observed signal)");
        TH1D *hRecoInteractingKEChan2Signal            = new TH1D("hRecoInteractingKEChan2Signal",        "Interacting Kinetic Energy Channel 2 Signal"           ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan2SignalElastic     = new TH1D("hRecoInteractingKEChan2SignalElastic", "Interacting Kinetic Energy Channel 2 Elastic Signal"   ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan2SignalInelast     = new TH1D("hRecoInteractingKEChan2SignalInelast","Interacting Kinetic Energy Channel  2 Inelast Signal"   ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan2Background        = new TH1D("hRecoInteractingKEChan2Background",    "Interacting Kinetic Energy Channel 2 Background"       ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan2BackgroundElastic = new TH1D("hRecoInteractingKEChan2BackgroundElastic","Interacting Kinetic Energy Channel 2 Background"    ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan2BackgroundInelast = new TH1D("hRecoInteractingKEChan2BackgroundInelast","Interacting Kinetic Energy Channel 2 Background"    ,20, 0, 1000);
        THStack *hRecoInteractingKEChan2StackedP       = new THStack("hRecoInteractingKEChan2StackedP",   "Interacting Kinetic Energy Channel 2(signal and noise)");
        THStack *hRecoInteractingKEChan2StackedE       = new THStack("hRecoInteractingKEChan2StackedE",   "Interacting Kinetic Energy Channel 2(true and observed signal)");
        TH1D *hRecoInteractingKEChan3TrueSignal        = new TH1D("hRecoInteractingKEChan3TrueSignal",    "Interacting Kinetic Energy Channel 3 True Signal"      ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan3TrueSignalG6      = new TH1D("hRecoInteractingKEChan3TrueSignalG6","Interacting Kinetic Energy Channel 3 True Signal Visible",20, 0, 1000);
        TH1D *hRecoInteractingKEChan3Signal            = new TH1D("hRecoInteractingKEChan3Signal",        "Interacting Kinetic Energy Channel 3 Signal"           ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan3SignalElastic     = new TH1D("hRecoInteractingKEChan3SignalElastic", "Interacting Kinetic Energy Channel 3 Signal Elastic"   ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan3SignalInelast     = new TH1D("hRecoInteractingKEChan3SignalInelast", "Interacting Kinetic Energy Channel 3 Signal Inelast"   ,20, 0, 1000);
        THStack *hRecoInteractingKEChan3StackedP       = new THStack("hRecoInteractingKEChan3StackedP",   "Interacting Kinetic Energy Channel 3(signal and noise)");
        THStack *hRecoInteractingKEChan3StackedE       = new THStack("hRecoInteractingKEChan3StackedE",   "Interacting Kinetic Energy Channel 3(true and observed signal)");
        TH1D *hRecoInteractingKEChan3Background        = new TH1D("hRecoInteractingKEChan3Background",    "Interacting Kinetic Energy Channel 3 Background"      ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan3BackgroundElastic = new TH1D("hRecoInteractingKEChan3BackgroundElastic","Interacting Kinetic Energy Channel 3 Background"   ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan3BackgroundInelast = new TH1D("hRecoInteractingKEChan3BackgroundInelast","Interacting Kinetic Energy Channel 3 Background"   ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan4TrueSignal    = new TH1D("hRecoInteractingKEChan4TrueSignal",    "Interacting Kinetic Energy Channel 4 True Signal"         ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan4TrueSignalG6  = new TH1D("hRecoInteractingKEChan4TrueSignalG6",  "Interacting Kinetic Energy Channel 4 True Signal Visible" ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan4Signal        = new TH1D("hRecoInteractingKEChan4Signal",        "Interacting Kinetic Energy Channel 4 Signal"              ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan4SignalElastic = new TH1D("hRecoInteractingKEChan4SignalElastic", "Interacting Kinetic Energy Channel 4 Signal Elastic"      ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan4SignalInelast = new TH1D("hRecoInteractingKEChan4SignalInelast", "Interacting Kinetic Energy Channel 4 Signal Inelast"      ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan4Background    = new TH1D("hRecoInteractingKEChan4Background",    "Interacting Kinetic Energy Channel 4 Background"          ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan4BackgroundElastic = new TH1D("hRecoInteractingKEChan4BackgroundElastic","Interacting Kinetic Energy Channel 4 Background"   ,20, 0, 1000);
        TH1D *hRecoInteractingKEChan4BackgroundInelast = new TH1D("hRecoInteractingKEChan4BackgroundInelast","Interacting Kinetic Energy Channel 4 Background"   ,20, 0, 1000);
        THStack *hRecoInteractingKEChan4StackedP   = new THStack("hRecoInteractingKEChan4StackedP",   "Interacting Kinetic Energy Channel 4(signal and noise)");
        THStack *hRecoInteractingKEChan4StackedE   = new THStack("hRecoInteractingKEChan4StackedE",   "Interacting Kinetic Energy Channel 4(true and observed signal)");
TH1D *hRecoInteractingKEPurity        = new TH1D("hRecoInteractingKEPurity", "Interacting Kinetic Energy (MeV) Purity"    , 20, 0, 1000);
TH1D *hRecoInteractingKEPurityElastic = new TH1D("hRecoInteractingKEPurityElastic", "Elastic Interacting Kinetic Energy (MeV) Purity"                 , 20, 0, 1000);
TH1D *hRecoInteractingKEPurityInelast = new TH1D("hRecoInteractingKEPurityInelast", "Inelastic Interacting Kinetic Energy (MeV) Purity"               , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan1Purity = new TH1D("hRecoInteractingKEChan1Purity", "Interacting Kinetic Energy (MeV) Channel 1 Purity"               , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan1PurityElastic = new TH1D("hRecoInteractingKEChan1PurityElastic", "Interacting Kinetic Energy (MeV) Channel 1 Purity" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan1PurityInelast = new TH1D("hRecoInteractingKEChan1PurityInelast", "Interacting Kinetic Energy (MeV) Channel 1 Purity" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan2Purity = new TH1D("hRecoInteractingKEChan2Purity", "Interacting Kinetic Energy (MeV) Channel 2 Purity"               , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan2PurityElastic = new TH1D("hRecoInteractingKEChan2PurityElastic", "Interacting Kinetic Energy (MeV) Channel 2 Purity" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan2PurityInelast = new TH1D("hRecoInteractingKEChan2PurityInelast", "Interacting Kinetic Energy (MeV) Channel 2 Purity" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan3Purity = new TH1D("hRecoInteractingKEChan3Purity", "Interacting Kinetic Energy (MeV) Channel 3 Purity"               , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan3PurityElastic = new TH1D("hRecoInteractingKEChan3PurityElastic", "Interacting Kinetic Energy (MeV) Channel 3 Purity" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan3PurityInelast = new TH1D("hRecoInteractingKEChan3PurityInelast", "Interacting Kinetic Energy (MeV) Channel 3 Purity" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan4Purity = new TH1D("hRecoInteractingKEChan4Purity", "Interacting Kinetic Energy (MeV) Channel 4 Purity"               , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan4PurityElastic = new TH1D("hRecoInteractingKEChan4PurityElastic", "Interacting Kinetic Energy (MeV) Channel 4 Purity" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan4PurityInelast = new TH1D("hRecoInteractingKEChan4PurityInelast", "Interacting Kinetic Energy (MeV) Channel 4 Purity" , 20, 0, 1000);
TH1D *hRecoInteractingKEEffish = new TH1D("hRecoInteractingKEEffish", "Interacting Kinetic Energy (MeV) Efficiency", 20, 0, 1000);
TH1D *hRecoInteractingKEElasticEffish    = new TH1D("hRecoInteractingKEElasticEffish",    "Interacting Kinetic Energy (MeV) Efficiency Elastic", 20, 0, 1000);
TH1D *hRecoInteractingKEElasticObsEffish = new TH1D("hRecoInteractingKEElasticObsEffish", "Interacting Kinetic Energy (MeV) Efficiency Elastic (>6)", 20, 0, 1000);
TH1D *hRecoInteractingKEInelastEffish    = new TH1D("hRecoInteractingKEInelastEffish",    "Interacting Kinetic Energy (MeV) Efficiency Inelastic", 20, 0, 1000);
    TH1D *hRecoInteractingKEChan1Effish  = new TH1D("hRecoInteractingKEChan1Effish",      "Interacting Kinetic Energy (MeV) Channel 1 Effish"    , 20, 0, 1000);
      TH1D *hRecoInteractingKEChan1ElasticEffish = new TH1D("hRecoInteractingKEChan1ElasticEffish", "Interacting Kinetic Energy (MeV) Channel 1 Elastic Effish" , 20, 0, 1000);
      TH1D *hRecoInteractingKEChan1InelastEffish = new TH1D("hRecoInteractingKEChan1InelastEffish", "Interacting Kinetic Energy (MeV) Channel 1 Inelast Effish" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan2Effish  = new TH1D("hRecoInteractingKEChan2Effish",      "Interacting Kinetic Energy (MeV) Channel 2 Effish"    , 20, 0, 1000);
      TH1D *hRecoInteractingKEChan2ElasticEffish = new TH1D("hRecoInteractingKEChan2ElasticEffish", "Interacting Kinetic Energy (MeV) Channel 2 Elastic Effish" , 20, 0, 1000);
      TH1D *hRecoInteractingKEChan2InelastEffish = new TH1D("hRecoInteractingKEChan2InelastEffish", "Interacting Kinetic Energy (MeV) Channel 2 Inelast Effish" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan3Effish = new TH1D("hRecoInteractingKEChan3Effish", "Interacting Kinetic Energy (MeV) Channel 3 Effish"    , 20, 0, 1000);
      TH1D *hRecoInteractingKEChan3ElasticEffish = new TH1D("hRecoInteractingKEChan3ElasticEffish", "Interacting Kinetic Energy (MeV) Channel 3 Elastic Effish" , 20, 0, 1000);
      TH1D *hRecoInteractingKEChan3InelastEffish = new TH1D("hRecoInteractingKEChan3InelastEffish", "Interacting Kinetic Energy (MeV) Channel 3 Inelast Effish" , 20, 0, 1000);
    TH1D *hRecoInteractingKEChan4Effish = new TH1D("hRecoInteractingKEChan4Effish", "Interacting Kinetic Energy (MeV) Channel 4 Effish"    , 20, 0, 1000);
      TH1D *hRecoInteractingKEChan4ElasticEffish = new TH1D("hRecoInteractingKEChan4ElasticEffish", "Interacting Kinetic Energy (MeV) Channel 4 Elastic Effish" , 20, 0, 1000);
      TH1D *hRecoInteractingKEChan4InelastEffish = new TH1D("hRecoInteractingKEChan4InelastEffish", "Interacting Kinetic Energy (MeV) Channel 4 Inelast Effish" , 20, 0, 1000);
TH1F *fCrossSection            = new TH1F("fCrossSection",      "Total Cross-Section"             , 20, 0, 1000);
TH1F *fCrossSectionElastic     = new TH1F("fElasticCrossSection",      "Elastic Cross-Section"     , 20, 0, 1000);
TH1F *fCrossSectionInelast     = new TH1F("fInelastCrossSection",      "Inelastic Cross-Section"   , 20, 0, 1000);

// ### add ###
TH1D *hRecoResolutionDiff     = new TH1D("hRecoResolutionDiff",    "(Reco - True) KE", 100, -500, 500 );
TH1D *hRecoResolutionDiff200  = new TH1D("hRecoResolutionDiff200", "(Reco - True) KE", 100, -500, 500 );
// ## by channel ##
TH1D *hRecoResolution0to50         = new TH1D("hRecoResolution0to50",         "Reco Interacting KE for True Slice (0, 50)"       ,    125 , -200 , 500);
TH1D *hRecoResolution0to50Chan1    = new TH1D("hRecoResolution0to50Chan1",    "Reco Interacting KE for True Slice (0, 50) Chan1" ,    125 , -200 , 500);
TH1D *hRecoResolution0to50Chan2    = new TH1D("hRecoResolution0to50Chan2",    "Reco Interacting KE for True Slice (0, 50) Chan2" ,    125 , -200 , 500);
TH1D *hRecoResolution0to50Chan3    = new TH1D("hRecoResolution0to50Chan3",    "Reco Interacting KE for True Slice (0, 50) Chan3" ,    125 , -200 , 500);
TH1D *hRecoResolution0to50Chan4    = new TH1D("hRecoResolution0to50Chan4",    "Reco Interacting KE for True Slice (0, 50) Chan4" ,    125 , -200 , 500);
TH1D *hRecoResolution50to100       = new TH1D("hRecoResolution50to100",       "Reco Interacting KE for True Slice (50, 100)"       ,  125 , -200 , 500);
TH1D *hRecoResolution50to100Chan1  = new TH1D("hRecoResolution50to100Chan1",  "Reco Interacting KE for True Slice (50, 100) Chan1" ,  125 , -200 , 500);
TH1D *hRecoResolution50to100Chan2  = new TH1D("hRecoResolution50to100Chan2",  "Reco Interacting KE for True Slice (50, 100) Chan2" ,  125 , -200 , 500);
TH1D *hRecoResolution50to100Chan3  = new TH1D("hRecoResolution50to100Chan3",  "Reco Interacting KE for True Slice (50, 100) Chan3" ,  125 , -200 , 500);
TH1D *hRecoResolution50to100Chan4  = new TH1D("hRecoResolution50to100Chan4",  "Reco Interacting KE for True Slice (50, 100) Chan4" ,  125 , -200 , 500);
TH1D *hRecoResolution100to150      = new TH1D("hRecoResolution100to150",      "Reco Interacting KE for True Slice (100, 150)"       , 125 , -200 , 500);
TH1D *hRecoResolution100to150Chan1 = new TH1D("hRecoResolution100to150Chan1", "Reco Interacting KE for True Slice (100, 150) Chan1" , 125 , -200 , 500);
TH1D *hRecoResolution100to150Chan2 = new TH1D("hRecoResolution100to150Chan2", "Reco Interacting KE for True Slice (100, 150) Chan2" , 125 , -200 , 500);
TH1D *hRecoResolution100to150Chan3 = new TH1D("hRecoResolution100to150Chan3", "Reco Interacting KE for True Slice (100, 150) Chan3" , 125 , -200 , 500);
TH1D *hRecoResolution100to150Chan4 = new TH1D("hRecoResolution100to150Chan4", "Reco Interacting KE for True Slice (100, 150) Chan4" , 125 , -200 , 500);
TH1D *hRecoResolution150to200      = new TH1D("hRecoResolution150to200",      "Reco Interacting KE for True Slice (150, 200)"      ,  125 , -200 , 500);
TH1D *hRecoResolution150to200Chan1 = new TH1D("hRecoResolution150to200Chan1", "Reco Interacting KE for True Slice (150, 200) Chan1",  125 , -200 , 500);
TH1D *hRecoResolution150to200Chan2 = new TH1D("hRecoResolution150to200Chan2", "Reco Interacting KE for True Slice (150, 200) Chan2",  125 , -200 , 500);
TH1D *hRecoResolution150to200Chan3 = new TH1D("hRecoResolution150to200Chan3", "Reco Interacting KE for True Slice (150, 200) Chan3",  125 , -200 , 500);
TH1D *hRecoResolution150to200Chan4 = new TH1D("hRecoResolution150to200Chan4", "Reco Interacting KE for True Slice (150, 200) Chan4",  125 , -200 , 500);
TH1D *hRecoResolution200to250      = new TH1D("hRecoResolution200to250",      "Reco Interacting KE for True Slice (200, 250)"       , 125 , 0    , 500);
TH1D *hRecoResolution200to250Chan1 = new TH1D("hRecoResolution200to250Chan1", "Reco Interacting KE for True Slice (200, 250) Chan1" , 125 , 0    , 500);
TH1D *hRecoResolution200to250Chan2 = new TH1D("hRecoResolution200to250Chan2", "Reco Interacting KE for True Slice (200, 250) Chan2" , 125 , 0    , 500);
TH1D *hRecoResolution200to250Chan3 = new TH1D("hRecoResolution200to250Chan3", "Reco Interacting KE for True Slice (200, 250) Chan3" , 125 , 0    , 500);
TH1D *hRecoResolution200to250Chan4 = new TH1D("hRecoResolution200to250Chan4", "Reco Interacting KE for True Slice (200, 250) Chan4" , 125 , 0    , 500);
TH1D *hRecoResolution250to300      = new TH1D("hRecoResolution250to300",      "Reco Interacting KE for True Slice (250, 300)"       , 125 , 0    , 500);
TH1D *hRecoResolution250to300Chan1 = new TH1D("hRecoResolution250to300Chan1", "Reco Interacting KE for True Slice (250, 300) Chan1" , 125 , 0    , 500);
TH1D *hRecoResolution250to300Chan2 = new TH1D("hRecoResolution250to300Chan2", "Reco Interacting KE for True Slice (250, 300) Chan2" , 125 , 0    , 500);
TH1D *hRecoResolution250to300Chan3 = new TH1D("hRecoResolution250to300Chan3", "Reco Interacting KE for True Slice (250, 300) Chan3" , 125 , 0    , 500);
TH1D *hRecoResolution250to300Chan4 = new TH1D("hRecoResolution250to300Chan4", "Reco Interacting KE for True Slice (250, 300) Chan4" , 125 , 0    , 500);
TH1D *hRecoResolution300to350      = new TH1D("hRecoResolution300to350",      "Reco Interacting KE for True Slice (300, 350)"       , 125 , 100  , 600);
TH1D *hRecoResolution300to350Chan1 = new TH1D("hRecoResolution300to350Chan1", "Reco Interacting KE for True Slice (300, 350) Chan1" , 125 , 100  , 600);
TH1D *hRecoResolution300to350Chan2 = new TH1D("hRecoResolution300to350Chan2", "Reco Interacting KE for True Slice (300, 350) Chan2" , 125 , 100  , 600);
TH1D *hRecoResolution300to350Chan3 = new TH1D("hRecoResolution300to350Chan3", "Reco Interacting KE for True Slice (300, 350) Chan3" , 125 , 100  , 600);
TH1D *hRecoResolution300to350Chan4 = new TH1D("hRecoResolution300to350Chan4", "Reco Interacting KE for True Slice (300, 350) Chan4" , 125 , 100  , 600);
TH1D *hRecoResolution350to400      = new TH1D("hRecoResolution350to400",      "Reco Interacting KE for True Slice (350, 400)"       , 125 , 100  , 600);
TH1D *hRecoResolution350to400Chan1 = new TH1D("hRecoResolution350to400Chan1", "Reco Interacting KE for True Slice (350, 400) Chan1" , 125 , 100  , 600);
TH1D *hRecoResolution350to400Chan2 = new TH1D("hRecoResolution350to400Chan2", "Reco Interacting KE for True Slice (350, 400) Chan2" , 125 , 100  , 600);
TH1D *hRecoResolution350to400Chan3 = new TH1D("hRecoResolution350to400Chan3", "Reco Interacting KE for True Slice (350, 400) Chan3" , 125 , 100  , 600);
TH1D *hRecoResolution350to400Chan4 = new TH1D("hRecoResolution350to400Chan4", "Reco Interacting KE for True Slice (350, 400) Chan4" , 125 , 100  , 600);
TH1D *hRecoResolution400to450      = new TH1D("hRecoResolution400to450",      "Reco Interacting KE for True Slice (400, 450)"       , 125 , 200  , 700);
TH1D *hRecoResolution400to450Chan1 = new TH1D("hRecoResolution400to450Chan1", "Reco Interacting KE for True Slice (400, 450) Chan1" , 125 , 200  , 700);
TH1D *hRecoResolution400to450Chan2 = new TH1D("hRecoResolution400to450Chan2", "Reco Interacting KE for True Slice (400, 450) Chan2" , 125 , 200  , 700);
TH1D *hRecoResolution400to450Chan3 = new TH1D("hRecoResolution400to450Chan3", "Reco Interacting KE for True Slice (400, 450) Chan3" , 125 , 200  , 700);
TH1D *hRecoResolution400to450Chan4 = new TH1D("hRecoResolution400to450Chan4", "Reco Interacting KE for True Slice (400, 450) Chan4" , 125 , 200  , 700);
TH1D *hRecoResolution450to500      = new TH1D("hRecoResolution450to500",      "Reco Interacting KE for True Slice (450, 500)"       , 125 , 200  , 700);
TH1D *hRecoResolution450to500Chan1 = new TH1D("hRecoResolution450to500Chan1", "Reco Interacting KE for True Slice (450, 500) Chan1" , 125 , 200  , 700);
TH1D *hRecoResolution450to500Chan2 = new TH1D("hRecoResolution450to500Chan2", "Reco Interacting KE for True Slice (450, 500) Chan2" , 125 , 200  , 700);
TH1D *hRecoResolution450to500Chan3 = new TH1D("hRecoResolution450to500Chan3", "Reco Interacting KE for True Slice (450, 500) Chan3" , 125 , 200  , 700);
TH1D *hRecoResolution450to500Chan4 = new TH1D("hRecoResolution450to500Chan4", "Reco Interacting KE for True Slice (450, 500) Chan4" , 125 , 200  , 700);
TH1D *hRecoResolution500to550      = new TH1D("hRecoResolution500to550",      "Reco Interacting KE for True Slice (500, 550)"       , 125 , 200  , 700);
TH1D *hRecoResolution500to550Chan1 = new TH1D("hRecoResolution500to550Chan1", "Reco Interacting KE for True Slice (500, 550) Chan1" , 125 , 200  , 700);
TH1D *hRecoResolution500to550Chan2 = new TH1D("hRecoResolution500to550Chan2", "Reco Interacting KE for True Slice (500, 550) Chan2" , 125 , 200  , 700);
TH1D *hRecoResolution500to550Chan3 = new TH1D("hRecoResolution500to550Chan3", "Reco Interacting KE for True Slice (500, 550) Chan3" , 125 , 200  , 700);
TH1D *hRecoResolution500to550Chan4 = new TH1D("hRecoResolution500to550Chan4", "Reco Interacting KE for True Slice (500, 550) Chan4" , 125 , 200  , 700);
TH1D *hRecoResolution550to600      = new TH1D("hRecoResolution550to600",      "Reco Interacting KE for True Slice (550, 600)"       , 125 , 200  , 700);
TH1D *hRecoResolution550to600Chan1 = new TH1D("hRecoResolution550to600Chan1", "Reco Interacting KE for True Slice (550, 600) Chan1" , 125 , 200  , 700);
TH1D *hRecoResolution550to600Chan2 = new TH1D("hRecoResolution550to600Chan2", "Reco Interacting KE for True Slice (550, 600) Chan2" , 125 , 200  , 700);
TH1D *hRecoResolution550to600Chan3 = new TH1D("hRecoResolution550to600Chan3", "Reco Interacting KE for True Slice (550, 600) Chan3" , 125 , 200  , 700);
TH1D *hRecoResolution550to600Chan4 = new TH1D("hRecoResolution550to600Chan4", "Reco Interacting KE for True Slice (550, 600) Chan4" , 125 , 200  , 700);

// ## everything I need to calculate a cross section ##
TH2D *hMCRecoIncidentMatrix              = new TH2D("hMCRecoIncidentMatrix",              "Incident Smearing Matrix",    20, 0, 1000, 20, 0, 1000);
TH2D *hMCRecoIncidentMatrixNormalized    = new TH2D("hMCRecoIncidentMatrixNormalized",    "Incident Smearing Matrix",    20, 0, 1000, 20, 0, 1000);
TH2D *hMCRecoInteractingMatrix           = new TH2D("hMCRecoInteractingMatrix",           "Interacting Smearing Matrix", 20, 0, 1000, 20, 0, 1000);
TH2D *hMCRecoInteractingMatrixNormalized = new TH2D("hMCRecoInteractingMatrixNormalized", "Interacting Smearing Matrix", 20, 0, 1000, 20, 0, 1000);
TH1D *hMCRecoIncidentKE        = new TH1D("hMCRecoIncidentKE",        "Truth Incident KE",           20, 0, 1000);
TH1D *hMCRecoInteractingKE     = new TH1D("hMCRecoInteractingKE",     "Truth Selected KE",         20, 0, 1000);
TH1D *hMCRecoInteractingKESignalInelast     = new TH1D("hMCRecoInteractingKESignalInealst",     "Observed Signal", 20, 0, 1000);
TH1D *hMCRecoInteractingKEBackgroundInelast = new TH1D("hMCRecoInteractingKEBackgroundInealst", "Observed Background", 20, 0, 1000);
TH1D *fN_selected                 = new TH1D("fN_selected", "N selected",                           20, 0, 1000);
TH1D *fN_incident                 = new TH1D("fN_incident", "N incident",                           20, 0, 1000);
TH1D *fInteractingKEPurityInelast = new TH1D("fInteractingKEPurityInelast", "Purity",               20, 0, 1000);
TH1D *fInteractingKEEffishInelast = new TH1D("fInteractingKEEffishInelast", "Efficiency",           20, 0, 1000);
TH1F *fRatio                      = new TH1F("fRatio",                 "Ratio",                     20, 0, 1000);
TH1F *fRatioCorrected             = new TH1F("fRatioCorrected",        "Ratio Energy Corrected",    20, 0, 1000);
TH1F *fInelasticCrossSection      = new TH1F("fInelasticCrossSection", "Inelastic XS",              20, 0, 1000);

void blah::Loop()
{

float particle_mass = 938.28; //<---Mass of Proton in MeV

// ########################################################################
// ### Fiducial Boundry Cuts (used to determine if a track is stopping) ###
// ########################################################################
double XLowerFid = 0;
double XUpperFid = 47;

double YLowerFid = -20;
double YUpperFid = 20;

double ZLowerFid = 0;
double ZUpperFid = 90;


// ### The assumed energy loss between the cryostat and the TPC ###
float entryTPCEnergyLoss = 36; //MeV


// ### The assumed mass of the incident particle (here we assume a proton) ###
float mass = 938.28;

//float RecoLength = 0;

double interaction_angle = 10;
//double tpc_interaction_angle = 5;
double tpc_poi_dedx_mean = 10;

if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes = 0, nb = 0;

// ### counters ###
int nTotalEvents=0;  int nEvntsTPC=0;   int nRecoMatched=0;
int nSimpleMCTopo=0; int nNoIntb4TPC=0; int nSimpleRecoTopo=0; int nTotalReco=0; 
int nMoreThan3=0; int nExactly4=0; 
int nOneReco=0; int nTwoReco=0;
int nInteractions=0; int nInteractions_afterAllCuts=0; int nTracks_oneInt=0;
int wInitialEvents=0; int wInteractions_afterAllCuts=0; int wInteractions=0; int wTracks_oneInt=0;
int nTotalTrueSignal=0;
int nOneTrackTrueSignal=0;
int nOneTrackBackground=0; nOneTrackSignal=0;
int nOneTrackSignalArray     [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
int nOneTrackBackgroundArray [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
int nTwoTrackTrueSignal=0;
int nTwoTrackTrueBackground=0;
int nFirstRegime=0;
int nTwoTrackSignal=0; 
int nTwoTrackBackground=0;
int nOtherRegime=0;
int nSanityCheck=0;
int nOtherTwoTrackTrueSignal=0;
int nOtherTwoTrackSignal=0;
int nOtherTwoTrackBackground=0;

// ### xs counters ###
int nChannel0=0; nChannel0True=0; int nChannel0Signal=0; int nChannel0Background=0;
int nChannel1=0; nChannel1True=0; int nChannel1Signal=0; int nChannel1Background=0;
int nChannel2=0; nChannel2True=0; int nChannel2Signal=0; int nChannel2Background=0;
int nChannel3=0; nChannel3True=0; int nChannel3Signal=0; int nChannel3Background=0;
int nChannel4=0; nChannel4True=0; int nChannel4Signal=0; int nChannel4Background=0;
int nChannel2TopTrue=0;
int nChannel0InelastTrue=0; int nChannel0InelastSignal=0; int nChannel0InelastBackground=0;
int nChannel1InelastTrue=0; int nChannel1InelastSignal=0; int nChannel1InelastBackground=0;
int nChannel2InelastTrue=0; int nChannel2InelastSignal=0; int nChannel2InelastBackground=0;
int nChannel3InelastTrue=0; int nChannel3InelastSignal=0; int nChannel3InelastBackground=0;
int nChannel4InelastTrue=0; int nChannel4InelastSignal=0; int nChannel4InelastBackground=0;
//int nTwoTrackSignalArray     [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//int nTwoTrackBackgroundArray [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// ### resolution study vars ###
std::vector< std::pair<   double, double > > matched_ke;
std::vector< std::vector< double> >          matched_ke_tuple;

float rho = 1396; //kg/m^3
//  float cm_per_m = 100;
float molar_mass = 39.95; //g/mol
float g_per_kg = 1000;
float avogadro = 6.022e+23; //number/mol
float number_density = rho*g_per_kg/molar_mass*avogadro;
float slab_width = 0.01;//in m
//float recip_number_density = 2.10 * 10^22;

// #### #### ####
//   event loop
// #### #### ####

for (Long64_t jentry=0; jentry<nentries;jentry++) {
//for (Long64_t jentry=0; jentry<5000;jentry++) {
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   nbytes += nb;

   // #############################
   // ### Counting Total Events ###
   // #############################
   nTotalEvents++;
   
   // === Outputting every nEvents to the screen ===
   if(nTotalEvents % 1000 == 0){std::cout<<"Event = "<<nTotalEvents<<std::endl;}   
   //std::cout << "--------------" << std::endl;
   //std::cout << "EVT NUM: " << nTotalEvents << std::endl;
   //std::cout << "--------------" << std::endl;
   
   bool ExitingTrack = false;
   
   // blah for protons reaching tpc
   if((*EndPointz)[0] < 0.0){continue;} 
   //if(nhits < 1){continue;}
   if(ntracks_reco == 0){continue;}
   nEvntsTPC++;
   
   // ============================================================
   //                           MC study area 
   // ============================================================
   int whyamilooping = 0;
   double primary_momentum = 1000.*primary_p;
   //std::cout << "primary momentum: " << primary_momentum << std::endl;

  //int mc_process = get_primary_process();
  for(size_t nG4 = 0; nG4 < geant_list_size; nG4++)
    {
    if((*process_primary)[nG4] == 10){continue;}//<--- 0-nonprimary, 1-primary
    //if(nTotalEvents == 218){std::cout << "  what is going on here: " << (*process_primary)[nG4] << std::endl;}
    //if(nTotalEvents == 320){std::cout << "  what is going on here (320): " << (*process_primary)[nG4] << std::endl;}
    if((*process_primary)[nG4] != 1){continue;}//<--- 0-nonprimary, 1-primary
    whyamilooping++;
    // ### some weighting stuff ###
    double mofit[] = {1.0098, 699.4, 108.1};
    double weight = 1.;
    //    weight = mofit[0]*exp(-1*(pow(1000.*(*StartPz)[nG4] - mofit[1], 2))/(2*pow(mofit[2], 2)));
    bool use_weight = true;
    int mc_end_process = 0;
    // ### energy ###
    //    double mc_primary_energy = 1000.*(*StartPz)[nG4]; 

    // ### table?? ###
    wInitialEvents += weight;

    // ### end process histogram ###
    if((*G4FinalProcess)[nG4] == "LArVoxelReadoutScoringProcess"){mc_end_process=1;}
    else if((*G4FinalProcess)[nG4] == "protonInelastic"){mc_end_process=2;}
    else{mc_end_process=3;}

    // ### the vars ###
    double mc_int_x = 999; double mc_int_y = 999; double mc_int_z = 999;
    double mc_ch2_int_x = 999; double mc_ch2_int_y = 999; double mc_ch2_int_z = 999;
    bool mc_ch2_flag = false;
    int interaction_vector_size = InteractionPoint->size();
    bool did_interact;
    double interaction_p = 0;
    double ep0 [3]; double ep2 [3]; double ep1 [3];
    double ev0 [3]; double ev1 [3];
    double eangle;
    double cp0 [3]; double cp2 [3]; double cp1 [3];
    double cv0 [3]; double cv1 [3];
    double cangle;
    double ip0 [3]; double ip2 [3]; double ip1 [3];
    double iv0 [3]; double iv1 [3];
    double iangle;
    int t13one = 0;
    int nd1_trpts = -999;
    int nd2_trpts = -999;

    // ## some logic ##
    if(interaction_vector_size > 1)
      {
      nInteractions++;
      wInteractions += weight;
      hMCnInteractions->Fill(interaction_vector_size);
      continue;
      }  
    if(interaction_vector_size == 1)
      {
      if((*InteractionPointType)[0] == 0)
        {
        did_interact = false;
        hMCnInteractions->Fill(0);
        }
      if((*InteractionPointType)[0] != 0)
        {
        did_interact = true;    
        nInteractions++;
        wInteractions += weight;
        hMCnInteractions->Fill(interaction_vector_size);
        hMCOneIntNTracks->Fill(ntracks_reco);
        if(ntracks_reco == 1)
          {
          nTracks_oneInt ++;
          wTracks_oneInt += weight;
          if( (*InteractionPointType)[0] == 13 ) 
            {hMCOneIntOneTrackInteractionType->Fill(0);}
          if( (*InteractionPointType)[0] == 8 ) 
            {hMCOneIntOneTrackInteractionType->Fill(1);}
          if( (*InteractionPointType)[0] == 3 ) 
            {hMCOneIntOneTrackInteractionType->Fill(2);}
          }
        }
      }
    else
      {did_interact = true;}
    nSimpleMCTopo++;

    // ### if it did interact then get the interaction point ###
    if(did_interact == true)
      {
      if(is_in_tpc((*InteractionPoint)[0]) == false)
        {continue;}
      // ### int pt ###
      mc_int_x = (*MidPosX)[0][(*InteractionPoint)[0]];
      mc_int_y = (*MidPosY)[0][(*InteractionPoint)[0]];
      mc_int_z = (*MidPosZ)[0][(*InteractionPoint)[0]];
      nInteractions_afterAllCuts ++;
      wInteractions_afterAllCuts += weight;
      // ### getting g4 angles by interaction type ###
      if( (*InteractionPointType)[0] == 8)
        {
        cp0[0] = (*MidPosX)[0][(*InteractionPoint)[0]-1]; 
        cp0[1] = (*MidPosY)[0][(*InteractionPoint)[0]-1]; 
        cp0[2] = (*MidPosZ)[0][(*InteractionPoint)[0]-1]; 
        cp1[0] = mc_int_x; cp1[1] = mc_int_y; cp1[2] = mc_int_z;
        cp2[0] = (*MidPosX)[0][(*InteractionPoint)[0]+1]; 
        cp2[1] = (*MidPosY)[0][(*InteractionPoint)[0]+1]; 
        cp2[2] = (*MidPosZ)[0][(*InteractionPoint)[0]+1]; 
        cv0[0] = cp1[0] - cp0[0]; cv0[1] = cp1[1] - cp0[1]; cv0[2] = cp1[2] - cp0[2];
        cv1[0] = cp2[0] - cp1[0]; cv1[1] = cp2[1] - cp1[1]; cv1[2] = cp2[2] - cp1[2];

        cangle = get_angle(cv1, cv0);
        hMCCoulombAngle->Fill(cangle);
        hMCNormalizedAngleBackground->Fill(cangle);
        }
      if( (*InteractionPointType)[0] == 3)
        {
        ep0[0] = (*MidPosX)[0][(*InteractionPoint)[0]-1]; 
        ep0[1] = (*MidPosY)[0][(*InteractionPoint)[0]-1]; 
        ep0[2] = (*MidPosZ)[0][(*InteractionPoint)[0]-1]; 
        ep1[0] = mc_int_x; ep1[1] = mc_int_y; ep1[2] = mc_int_z;
        ep2[0] = (*MidPosX)[0][(*InteractionPoint)[0]+1]; 
        ep2[1] = (*MidPosY)[0][(*InteractionPoint)[0]+1]; 
        ep2[2] = (*MidPosZ)[0][(*InteractionPoint)[0]+1]; 
        ev0[0] = ep1[0] - ep0[0]; ev0[1] = ep1[1] - ep0[1]; ev0[2] = ep1[2] - ep0[2];
        ev1[0] = ep2[0] - ep1[0]; ev1[1] = ep2[1] - ep1[1]; ev1[2] = ep2[2] - ep1[2];

        eangle = get_angle(ev1, ev0);
        hMCElasticAngle->Fill(eangle);
        hMCNormalizedAngleSignal->Fill(eangle);
        interaction_p = sqrt( pow( 1000*(*MidPx)[nG4][(*InteractionPoint)[0]-1], 2 ) 
                            + pow( 1000*(*MidPy)[nG4][(*InteractionPoint)[0]-1], 2 ) 
                            + pow( 1000*(*MidPz)[nG4][(*InteractionPoint)[0]-1], 2 ) );
        //std::cout << "      elastic mc interaction p: " << interaction_p << std::endl;
        //std::cout << "      ELASTIC INTERACTION" << std::endl;
        //std::cout << "          primary momentum            : " << primary_momentum << std::endl;
        //std::cout << "          number of trajectory points : " << (*NTrTrajPts)[nG4] << std::endl; 
        //std::cout << "          interaction point           : " << (*InteractionPoint)[0] << std::endl;
        //std::cout << "          interaction momentum        : " << interaction_p << std::endl;
        //std::cout << "          interaction px              : " << (*MidPx)[0][(*InteractionPoint)[0]] << std::endl;
        //std::cout << "          interaction py              : " << (*MidPy)[0][(*InteractionPoint)[0]] << std::endl;
        //std::cout << "          interaction pz              : " << (*MidPz)[0][(*InteractionPoint)[0]] << std::endl;
        //std::cout << "          interaction momentum -1     : " << sqrt( pow((*MidPx)[nG4][(*InteractionPoint)[0]-1], 2)
        //                                                               + pow((*MidPy)[nG4][(*InteractionPoint)[0]-1], 2) +
        //                                                               + pow((*MidPz)[nG4][(*InteractionPoint)[0]-1], 2) ) << std::endl; 
        //for(int i = 0; i < (*NTrTrajPts)[nG4]; i++)
        //  {
        //  if(i == (*InteractionPoint)[0]){std::cout << "                     interaction" << std::endl;}
        //  std::cout << "            momentum: " <<  sqrt( pow( 1000*(*MidPx)[0][i], 2 ) 
        //                                                + pow( 1000*(*MidPy)[0][i], 2 ) 
        //                                                + pow( 1000*(*MidPz)[0][i], 2 ) ) << std::endl;

        //  }
        }
      if( (*InteractionPointType)[0] == 13)
        {
        if(whyamilooping>1){continue;}
        interaction_p = sqrt( pow( 1000*(*MidPx)[nG4][(*InteractionPoint)[0]-1], 2 ) 
                            + pow( 1000*(*MidPy)[nG4][(*InteractionPoint)[0]-1], 2 ) 
                            + pow( 1000*(*MidPz)[nG4][(*InteractionPoint)[0]-1], 2 ) );
        //std::cout << "      inelastic mc interaction p: " << interaction_p << std::endl;
        //hMCNDaughters->Fill();
        hMCNPDaughters->Fill(NProtonDaughters);
        if(NProtonDaughters == 0)
          {
          t13one = -1;
          if( abs((*InteractionPoint)[0] - (*NTrTrajPts)[0]) < 2 )
            {
            mc_ch2_flag = true;
            if(ntracks_reco == 1){nChannel2True++;}
            }
          mc_ch2_flag = true;
          hMC0Pntracks->Fill(ntracks_reco);
          mc_ch2_int_x = (*MidPosX)[0][(*InteractionPoint)[0]];
          mc_ch2_int_y = (*MidPosY)[0][(*InteractionPoint)[0]];
          mc_ch2_int_z = (*MidPosZ)[0][(*InteractionPoint)[0]];
          //for(int i = 0; i < (*NumberDaughters)[nG4]; i++)
          //  {
          //  std::cout << "  inelastic no proton daughter pdg code: " << (*DPdgCode)[i] << std::endl;
          //  }
          }
        if(NProtonDaughters == 1)
          {
          t13one = 1;
          ip0[0] = (*MidPosX)[0][(*InteractionPoint)[0]-1]; 
          ip0[1] = (*MidPosY)[0][(*InteractionPoint)[0]-1]; 
          ip0[2] = (*MidPosZ)[0][(*InteractionPoint)[0]-1];
          ip1[0] = mc_int_x; ip1[1] = mc_int_y; ip1[2] = mc_int_z;
          ip2[0] = (*DPMidPosX)[0][1];
          ip2[1] = (*DPMidPosY)[0][1]; 
          ip2[2] = (*DPMidPosZ)[0][1]; 
          iv0[0] = ip1[0] - ip0[0]; iv0[1] = ip1[1] - ip0[1]; iv0[2] = ip1[2] - ip0[2]; 
          iv1[0] = ip2[0] - ip1[0]; iv1[1] = ip2[1] - ip1[1]; iv1[2] = ip2[2] - ip1[2]; 

          iangle = get_angle(iv1, iv0);
          hMCInelasticAngle->Fill(iangle);
          hMCNormalizedAngleSignal->Fill(iangle);
          hMC1Pntracks->Fill(ntracks_reco);
          }
        if(NProtonDaughters == 2)
          {
          //std::cout << "Inelastic collision. Two daughter protons. " << std::endl;
          for(int i = 0; i < (*NumberDaughters)[nG4]; i++)
            {
            if((*DPdgCode)[i] == 22 || (*DPdgCode)[i] == 2112){continue;}
            //double daughter_momentum = sqrt( pow( 1000*(*StartPz)[i], 2) + pow( 1000*(*StartPy)[i], 2) + pow( 1000*(*StartPx)[i], 2) );
            //std::cout << "  Particle: " << (*DPdgCode)[i] << "  momentum: " << daughter_momentum << std::endl;
            }
          nd1_trpts = (*NDTrTrajPts)[0];
          nd2_trpts = (*NDTrTrajPts)[1];
          hMC2Pntracks->Fill(ntracks_reco); 
          t13one = -2;
          }
        if(NProtonDaughters == 3)
          {hMC3Pntracks->Fill(ntracks_reco);}
        }
      }//<---End if did interact
    nNoIntb4TPC++;
    double mc_end_x = (*EndPointx)[nG4]; double mc_end_y = (*EndPointy)[nG4]; double mc_end_z = (*EndPointz)[nG4];
    // ### matching ###
    int tpc_spt_cnt = 0;
    double mc_energy_loss = 0;
    double mc_energy_loss_ratio = 0;
    double mc_match [3];
    double t = 0;
    double planar_x = 0; double planar_y = 0; double planar_z = 0;
    double temp_currentKE  = 9999;
    double temp_currentKE2 = 9999;
   
    if((*NTrTrajPts)[0] < 2){continue;} 
    //std::cout << "new particle!" << std::endl;
    //std::cout << "  total number of trajectory points: " << (*NTrTrajPts)[0] << std::endl;
    double dE = 0; double dEdX = 1.; double mc_dX = 0; double dumb_mc_dX = 0; 
    std::vector<double> g4_ke_at_slabs;
    for(size_t npt = 1; npt < (*NTrTrajPts)[0]-1; npt++)
      {
      double p0[] = {(*MidPosX)[0][npt-1], (*MidPosY)[0][npt-1], (*MidPosZ)[0][npt-1]};
      double p1[] = {(*MidPosX)[0][npt],   (*MidPosY)[0][npt],   (*MidPosZ)[0][npt]};
      double p2[] = {(*MidPosX)[0][npt+1], (*MidPosY)[0][npt+1], (*MidPosZ)[0][npt+1]};
      double v0[] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
      double v1[] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
      double angle = get_angle(v1, v0);

      if(is_in_tpc(npt) == false && (angle > interaction_angle)){continue;}
      if(is_in_tpc(npt))
        {
        if(did_interact == true && npt == (*InteractionPoint)[0])
          {hMCEveryIntAngle->Fill(angle);}
        if(did_interact == true && npt != (*InteractionPoint)[0])
          {hMCEveryNonIntAngle->Fill(angle); hMCNormalizedAngleBackground->Fill(angle);}
        if(did_interact == false)
          {hMCEveryNonIntAngle->Fill(angle); hMCNormalizedAngleBackground->Fill(angle);}
        hMCEveryAngle->Fill(angle);
        if(tpc_spt_cnt == 0)
          {
          tpc_spt_cnt++;
          for(int i = 0; i < 3; i++){mc_match[i] = v1[i];} 
          t = (0 - p1[2]) / (v1[2]);
          planar_x = p1[0] + t*v1[0];
          planar_y = p1[1] + t*v1[1];
          double temp_initialKE = sqrt(pow(primary_momentum, 2) + pow(938.28, 2)) - 938.28;
          temp_currentKE  = sqrt(pow(1000*(*MidPx)[0][npt], 2) + pow(1000*(*MidPy)[0][npt], 2) 
                               + pow(1000*(*MidPz)[0][npt], 2) + pow(938.28, 2)) - 938.28; 
          temp_currentKE2 = sqrt(pow(1000*(*MidPx)[0][npt-1], 2) + pow(1000*(*MidPy)[0][npt-1], 2) 
                               + pow(1000*(*MidPz)[0][npt-1], 2) + pow(938.28, 2)) - 938.28; 
          mc_energy_loss = temp_initialKE - temp_currentKE;
          mc_energy_loss_ratio = mc_energy_loss / temp_initialKE;
          dumb_mc_dX = sqrt(pow(p1[0]-(*EndPointx)[0],2) + pow(p1[1]-(*EndPointy)[0],2) + pow(p1[2]-(*EndPointz)[0],2));
          }//<--- End for spt in tpc
        // ### look at the distance between trajectory points ###
        double g4dp      = sqrt( pow(p1[0] - p0[0], 2) 
                               + pow(p1[1] - p0[1], 2) 
                               + pow(p1[2] - p0[2], 2) );
        double next_g4dp = sqrt( pow(p2[0] - p1[0], 2) 
                               + pow(p2[1] - p1[1], 2) 
                               + pow(p2[2] - p1[2], 2) );
        hMCGeant4StepSize->Fill(g4dp);

        // ### grabbing incident KE ###
        double midp       = sqrt( pow(1000*(*MidPx)[0][npt], 2) 
                                + pow(1000*(*MidPy)[0][npt], 2) 
                                + pow(1000*(*MidPz)[0][npt], 2) );
        double midke      = sqrt( pow(1000*(*MidPx)[0][npt], 2) + pow(1000*(*MidPy)[0][npt], 2) 
                                + pow(1000*(*MidPz)[0][npt], 2) + pow(938.28, 2)) - 938.28;
        double next_midke = sqrt( pow(1000*(*MidPx)[0][npt+1], 2) + pow(1000*(*MidPy)[0][npt+1], 2) 
                                + pow(1000*(*MidPz)[0][npt+1], 2) + pow(938.28, 2)) - 938.28;
        double preke      = sqrt( pow(1000*(*MidPx)[0][npt-1], 2) + pow(1000*(*MidPy)[0][npt-1], 2) 
                                + pow(1000*(*MidPz)[0][npt-1], 2) + pow(938.28, 2)) - 938.28;
        double delta_ke = preke - midke;
        mc_dX += g4dp;
        int g4_slab = (int)mc_dX + 1;
        double mc_next_dX = mc_dX + next_g4dp; 
        
        // ### g4 de/dx ###
        //std::cout << "  mc_dX: "   << mc_dX << std::endl;        
        //std::cout << "  next_dX: " << mc_next_dX << std::endl;        
        //if(mc_dX < g4_slab && g4_slab < mc_next_dX)
        //  {
        //  std::cout << "  surrounds slab!" << std::endl;
        //  double mc_ke_slope   = (next_midke - midke) / (mc_next_dX - mc_dX);
        //  double mc_ke_at_slab = mc_ke_slope*(g4_slab - mc_dX) + midke;
        //  g4_ke_at_slabs.push_back(1.);
        //  }

        }//<--- End if in tpc 
      }//<--- End mc spts loop
      hMCEnergyLossIntoTPC->Fill(mc_energy_loss);
      hMCEnergyLossIntoTPCRatio->Fill(mc_energy_loss_ratio);


// ### reco info ###
      
      // ### reco info ###
      if(ntracks_reco > 0){nTotalReco++;}
      if(ntracks_reco ==4){nExactly4++;}
      if(ntracks_reco > 3){nMoreThan3++; continue;} 
      nSimpleRecoTopo++;
    
      int reco_primary = get_primary(planar_x, planar_y, mc_match[0], mc_match[1], mc_match[2]);
      if(reco_primary == 999)
        {continue;}
      nRecoMatched++;
      int poi_space_point = 0;
      // ### primary kinetic energy and like things ###


      double momentum = primary_momentum;
      double mass = 938.28;
      //double ke = 1.;
      //double ke             = sqrt( pow(momentum, 2) + pow(mass, 2) ) - mass;
      //double ke             = temp_currentKE;
      //double ke             = temp_currentKE2;
      //double temp_init_ke   = temp_currentKE;
      double temp_init_ke   = temp_currentKE2;
      double init_ke   = sqrt( pow(momentum, 2) + pow(mass, 2) ) - mass;
      double ke = init_ke - 30;
      double true_ke = init_ke - 30;
      double interaction_ke = sqrt( pow(interaction_p, 2) + pow(mass, 2) ) - mass;
      //if(interaction_ke == 0){std::cout << "    zero interacting ke...?" << std::endl;}
      // ### stupid fix attempt for initial energy ###
      //double ke = 1.;
      //if(init_ke > 600){ke = init_ke - 30;}
      //if(init_ke < 600 && init_ke > 100){ke = init_ke - 45;}
      //if(init_ke < 100){continue;}
      //if(temp_init_ke < 10){ke = 0.89*init_ke;}
      //else{ke = temp_init_ke;}
      if(init_ke > 500){continue;}

      // ### reco track loop ###
      if(did_interact == true && (*InteractionPointType)[0] != 8)
        {
        //std::cout << "reco True interaction energy! " << interaction_p << std::endl;
        hMCInitialInteractingKE->Fill(ke, interaction_ke);
        hRecoInteractingKETrueSignal->Fill(interaction_ke);
        nTotalTrueSignal++;
        if((*InteractionPointType)[0] == 13)
          {
          hMCInelastInitialInteractingKE->Fill(ke, interaction_ke);
          hRecoInteractingKETrueSignalInelast->Fill(interaction_ke);
          if(t13one == 1){hRecoInteractingKETrueSignalIA->Fill(interaction_ke, iangle);}
          if(t13one == 1 && iangle > 6)
            {hRecoInteractingKETrueSignalInelastG6->Fill(interaction_ke); hRecoInteractingKETrueSignalG6->Fill(interaction_ke);}
          if(t13one == -1)
            {hRecoInteractingKETrueSignalG6->Fill(interaction_ke);}
          if(t13one == -2)
            {hRecoInteractingKETrueSignalG6->Fill(interaction_ke);}
          }
        if((*InteractionPointType)[0] == 3)     
          {
          hMCElasticInitialInteractingKE->Fill(ke, interaction_ke);
          hRecoInteractingKETrueSignalElastic->Fill(interaction_ke);
          hRecoInteractingKETrueSignalEA->Fill(interaction_ke, eangle);
          if(eangle > 6)
            {hRecoInteractingKETrueSignalElasticG6->Fill(interaction_ke); hRecoInteractingKETrueSignalG6->Fill(interaction_ke);}
          }
        }

      // ### reco track loop ###
      for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
        {
        if(nTPCtrk != reco_primary){continue;}
        //if(did_interact == true && (*InteractionPointType)[0] != 8)
        //  {
        //  std::cout << "reco True interaction energy! " << interaction_p << std::endl;
        //  nTotalTrueSignal++;
        //  }
        if(ntracks_reco == 1)
          {
          if(did_interact == true && (*InteractionPointType)[0] != 8){nOneTrackTrueSignal++;}
          nOneReco++;
          double tempx = 999; double tempy = 999; double tempz = 999; double tempr = 999;
          double min_r = 999; double min_r_angle = 0; double kink_close_to_int = -1.;
          int spts_left = 999; double dist_to_mc_end = 999; double XSAngle = 7.; double dummy_angle = -1;
          int dummycnt = 0;
          double XSAngleArray []      = {0,   1,  2,  3,  4,  5,  6,  7,  8,  9};  
          double dummy_angle_array [] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
          int    dummy_spt_array []   = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
          double dummy_rr_array []    = {0,   0,  0,  0,  0,  0,  0,  0,  0,  0};
          double dummy_ph_array []    = {0,   0,  0,  0,  0,  0,  0,  0,  0,  0};
          int    dummycnt_array []    = {0,   0,  0,  0,  0,  0,  0,  0,  0,  0};
          int    dummy_ch1_sn   []    = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
          int    dummy_inel_t13 []    = {0,   0,  0,  0,  0,  0,  0,  0,  0,  0};
          double dummy_elas_ang []    = {0,   0,  0,  0,  0,  0,  0,  0,  0,  0};
          //bool chan1sig=false; bool chan1back=false; 
          //bool chan2sig=false; bool chan2back=false; 
          //bool chan3sig=false; bool chan3back=false;
          for(size_t nspts = 1; nspts < (*ntrack_hits)[nTPCtrk]-1; nspts++)
            {
            if(did_interact == true)
              {
              tempx = (*track_xpos)[nTPCtrk][nspts] - mc_int_x; 
              tempy = (*track_ypos)[nTPCtrk][nspts] - mc_int_y; 
              tempz = (*track_zpos)[nTPCtrk][nspts] - mc_int_z; 
              tempr = sqrt( pow(tempx, 2) + pow(tempy, 2) + pow(tempz, 2) );          
              }
            double reco_p0[] = {(*track_xpos)[nTPCtrk][nspts-1], (*track_ypos)[nTPCtrk][nspts-1], (*track_zpos)[nTPCtrk][nspts-1]};   
            double reco_p1[] = {(*track_xpos)[nTPCtrk][nspts], (*track_ypos)[nTPCtrk][nspts], (*track_zpos)[nTPCtrk][nspts]};
            double reco_p2[] = {(*track_xpos)[nTPCtrk][nspts+1], (*track_ypos)[nTPCtrk][nspts+1], (*track_zpos)[nTPCtrk][nspts+1]};
            double reco_v0[] = {reco_p1[0] - reco_p0[0], reco_p1[1] - reco_p0[1], reco_p1[2] - reco_p0[2]};
            double reco_v1[] = {reco_p2[0] - reco_p1[0], reco_p2[1] - reco_p1[1], reco_p2[2] - reco_p1[2]};
            double reco_spt_angle = get_angle(reco_v1, reco_v0);
            if(TMath::IsNaN(reco_spt_angle)){reco_spt_angle = 0.;}
            hRecoSTEEveryAngle->Fill(reco_spt_angle);
            if(did_interact == true && tempr < 2 && reco_spt_angle > kink_close_to_int)
              {
              kink_close_to_int = reco_spt_angle; 
              spts_left = (*ntrack_hits)[nTPCtrk] - nspts;
              dist_to_mc_end = sqrt( pow( ((*track_xpos)[nTPCtrk][nspts] - mc_end_x), 2)
                                   + pow( ((*track_ypos)[nTPCtrk][nspts] - mc_end_y), 2)
                                   + pow( ((*track_zpos)[nTPCtrk][nspts] - mc_end_z), 2) );
              }
            if(did_interact == true && tempr > 2)
              {hRecoNormalizedAngleBackground->Fill(reco_spt_angle);}
            if(did_interact == false)
              {
              hRecoSTEEveryAngleNoInts->Fill(reco_spt_angle);
              hRecoNormalizedAngleBackground->Fill(reco_spt_angle);
              }
            if(reco_spt_angle > XSAngle && dummycnt == 0)
              {dummy_angle = reco_spt_angle; dummycnt++;}
            for(int i = 0; i < 10; i++)
              {
              if(reco_spt_angle > XSAngleArray[i] && dummycnt_array[i] == 0)
                {
                dummy_angle_array[i] = reco_spt_angle; 
                dummy_spt_array[i] = nspts; 
                dummycnt_array[i] = 1;
                // compute ph here to match up over there for(){}
                for(size_t ospts = 1; ospts < nspts; ospts++)
                  {
                  dummy_ph_array[i] += sqrt(    pow( (*track_xpos)[nTPCtrk][ospts-1] - (*track_xpos)[nTPCtrk][ospts], 2)
                                              + pow( (*track_ypos)[nTPCtrk][ospts-1] - (*track_ypos)[nTPCtrk][ospts], 2)
                                              + pow( (*track_zpos)[nTPCtrk][ospts-1] - (*track_zpos)[nTPCtrk][ospts], 2) );
                  
                  }
                // compute rr here to match up over there for(){}
                for(size_t mspts = nspts; mspts < (*ntrack_hits)[nTPCtrk]; mspts++)
                  {
                  dummy_rr_array[i] += sqrt( pow( (*track_xpos)[nTPCtrk][mspts-1] - (*track_xpos)[nTPCtrk][mspts], 2)
                                              + pow( (*track_ypos)[nTPCtrk][mspts-1] - (*track_ypos)[nTPCtrk][mspts], 2)
                                              + pow( (*track_zpos)[nTPCtrk][mspts-1] - (*track_zpos)[nTPCtrk][mspts], 2) );
                  }
                }
              }//<--- End binning loop for xs angle loop
            }//<--- End reco track spts loop     

          if(dummy_angle != -1)
            {
            if(did_interact == true && (*InteractionPointType)[0] != 8){nOneTrackSignal++;}
            else{nOneTrackBackground++;}
            }
          for(int i = 0; i < 10; i++)
            {
            if(dummy_angle_array[i] != -1)
              {
              // ## blahing g4 thinks there was an interaciton
              if(did_interact == true && (*InteractionPointType)[0] != 8)
                {
                double delta_g4recox = (*track_xpos)[nTPCtrk][dummy_spt_array[i]] - mc_int_x;
                double delta_g4recoy = (*track_ypos)[nTPCtrk][dummy_spt_array[i]] - mc_int_y;
                double delta_g4recoz = (*track_zpos)[nTPCtrk][dummy_spt_array[i]] - mc_int_z;
                double delta_g4recor = sqrt( pow(delta_g4recox, 2) + pow(delta_g4recoy, 2) + pow(delta_g4recoz, 2) );
                if(delta_g4recor < 2)
                  {
                  // ## blahing `visible` interaction
                  if((*InteractionPointType)[0] == 3 && eangle >  6)
                    {nOneTrackSignalArray[i]++; dummy_ch1_sn[i] = 1;}
                  if((*InteractionPointType)[0] == 3 && eangle <= 6)
                    {nOneTrackBackgroundArray[i]++; dummy_ch1_sn[i] = 0;}

                  if((*InteractionPointType)[0] == 13)
                    {
                    nOneTrackSignalArray[i]++;
                    dummy_ch1_sn[i] = 1;
                    //if(t13one == 1 && iangle > 6)
                    //  {nOneTrackSignalArray[i]++; dummy_ch1_sn[i] = 1;}
                    //if(t13one == 1 && iangle < 6)
                    //  {nOneTrackBackgroundArray[i]++; dummy_ch1_sn[i] = 0;}
                    //else
                    //  {nOneTrackSignalArray[i]++; dummy_ch1_sn[i] = 1;}
                    }
                  }
                else{nOneTrackBackgroundArray[i]++; dummy_ch1_sn[i] = 0;}// chan1back=true;}
                }
              if(did_interact == false)
                {nOneTrackBackgroundArray[i]++; dummy_ch1_sn[i] = 0;}
              }//<--- End if this bin got tripped
            }//<--- End looping over each angle bin 

          if(did_interact == true && kink_close_to_int != -1)
            {
            hRecoSTEEveryAngleInts->Fill(kink_close_to_int);
            // ### fill by type here ###
            if( (*InteractionPointType)[0] == 13 ) 
              {
              hRecoSTEEveryAngleIntsInelast->Fill(kink_close_to_int); 
              hRecoNormalizedAngleSignal->Fill(kink_close_to_int);
              }
            if( (*InteractionPointType)[0] == 8 ) 
              {
              hRecoSTEEveryAngleIntsCoulomb->Fill(kink_close_to_int);
              hRecoNormalizedAngleBackground->Fill(kink_close_to_int);
              }
            if( (*InteractionPointType)[0] == 3 ) 
              {
              hRecoSTEEveryAngleIntsElastic->Fill(kink_close_to_int);
              hRecoNormalizedAngleSignal->Fill(kink_close_to_int);
              }

            }//<--- End if did interact

          // ### probably where I need to collect 1 track topology blah info ###
          // ### need to loop over every spacepoint up until interaction candidate space point ###
          
          bool kick = false;

          // ### Channel 1 ###
          double calph = 0;
          double ph = 0;
          if(dummy_spt_array[5] != -1)
            {
            nChannel1++;
            if(did_interact == true && (*InteractionPointType)[0] != 8)
              {
              nChannel1True++;
              if((*InteractionPointType)[0] == 13){nChannel1InelastTrue++;}
              hRecoInteractingKEChan1TrueSignal->Fill(interaction_ke);
              if((*InteractionPointType)[0] == 3)
                {
                if(eangle > 6)
                  {hRecoInteractingKEChan1TrueSignalG6->Fill(interaction_ke);}
                }
              if((*InteractionPointType)[0] == 13)
                {
                if(t13one == 1 && iangle > 6)
                  {hRecoInteractingKEChan1TrueSignalG6->Fill(interaction_ke);}
                if(t13one != 1)
                  {hRecoInteractingKEChan1TrueSignalG6->Fill(interaction_ke);}
                }
              }//<--- End true interaction info
            ph = dummy_ph_array[5];
            bool ch1_had_signal = false;
            double reco_total_lost = 0;
            double reco_total_dX = 0;

            // ### filling channel 1 incident interacting ###
            for(int nspts = 0; nspts < ((*col_track_hits)[nTPCtrk] - 1); nspts++)
              {
              if(kick == true){continue;}
              // ### slab change ###
              double simChannel_E = 0;
              for(int simChannel = 0; simChannel < primary_num_simChannel; simChannel++)
                {
                if((*primary_simChannel)[simChannel] != (*col_track_wire)[0][nspts]){continue;}
                simChannel_E = (*primary_simChannel_E)[simChannel];
                }
              if( (ph - calph) > 1.)
                {
                double energyLoss = (*col_track_dedx)[nTPCtrk][nspts] * (*col_track_pitch_hit)[nTPCtrk][nspts];
                double reco_dX    = (*col_track_pitch_hit)[nTPCtrk][nspts];
                reco_total_lost += energyLoss; 
                reco_total_dX   += reco_dX; 
                hRecodX->Fill(10*reco_dX);
                int middle_slab = ((int) reco_total_dX) + 1;
                double reco_next_dX = reco_total_dX + (*col_track_pitch_hit)[nTPCtrk][nspts+1];
                double next_ke = ke - (*col_track_dedx)[nTPCtrk][nspts+1]*(*col_track_pitch_hit)[nTPCtrk][nspts+1];
                if(reco_total_dX < middle_slab && middle_slab < reco_next_dX)
                  {
                  double ke_slope   = (next_ke - ke) / (reco_next_dX - reco_total_dX);
                  double ke_at_slab = ke_slope*(middle_slab - reco_total_dX) + ke;
                  hRecoIncidentKE->Fill(ke);
                  hMCRecoIncidentKE->Fill(true_ke);
                  hMCRecoIncidentMatrix->Fill(ke, true_ke);
                  }
                calph += (*col_track_pitch_hit)[nTPCtrk][nspts];
                ke -= energyLoss;
                true_ke -= simChannel_E;
                }
              else
                { 
                kick = true;
                double last_ke = ke - (*col_track_dedx)[nTPCtrk][nspts+1]*(*col_track_pitch_hit)[nTPCtrk][nspts+1];
                int    last_slab = ((int) reco_total_dX) + 1; 
                double int_ke_slope   = (last_ke - ke) / (ph - calph);
                double int_ke_at_slab = int_ke_slope*(last_slab - reco_total_dX) + ke;
                ke = int_ke_at_slab;
                if(dummy_spt_array[5] != -1)
                  {
                  if(ke < 0){ke = 1.;}
                  //ke += 25;
                  hRecoInteractingKE->Fill(ke);
                  if(interaction_ke == 0){interaction_ke = true_ke;}
                  hMCRecoInteractingKE->Fill(interaction_ke);
                  hMCRecoInteractingMatrix->Fill(ke, interaction_ke);
                  hRecoInteractingKEChan1->Fill(ke);
                  if(dummy_ch1_sn[5] == 1)
                    {
                    hRecoMatchedInteractingKERatio->Fill(ke / interaction_ke);
                    hRecoMatchedInteractingKECmpsn->Fill(interaction_ke, ke);
                    hRecoInteractingKESignal->Fill(ke);
                    hRecoInteractingKEChan1Signal->Fill(ke);
                    if((*InteractionPointType)[0] == 3 || (*InteractionPointType)[0] == 13){ch1_had_signal = true;}
                    if((*InteractionPointType)[0] == 3) 
                      {
                      nChannel1InelastBackground++;
                      hRecoInteractingKEBackgroundInelast->Fill(ke);
                      hRecoInteractingKESignalElastic->Fill(ke);
                      hRecoInteractingKEChan1SignalElastic->Fill(ke);
                      hRecoInteractingKEChan1BackgroundInelast->Fill(ke);
                      hRecoMatchedElasticAngleCmpsn->Fill(dummy_angle_array[6], eangle);
                      if(eangle > 6)
                        {hRecoInteractingKESignalElasticObs->Fill(ke);}
                      hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                      }
                    if((*InteractionPointType)[0] == 13)
                      {
                      hRecoMatchedInteractingKECmpsnInelast->Fill(interaction_ke, ke);
                      // ### test ###
                      matched_ke.push_back(std::make_pair(ke, interaction_ke));
                      std::vector< double > ke_tuple ; 
                      ke_tuple.push_back(ke); ke_tuple.push_back(interaction_ke); ke_tuple.push_back(1.5);
                      matched_ke_tuple.push_back(ke_tuple);
                      nChannel1InelastSignal++;
                      hRecoInteractingKEBackgroundElastic->Fill(ke);
                      hRecoInteractingKESignalInelast->Fill(ke);
                      hRecoInteractingKEChan1SignalInelast->Fill(ke);
                      hRecoInteractingKEChan1BackgroundElastic->Fill(ke);
                      hMCRecoInteractingKESignalInelast->Fill(interaction_ke);
                      }
                    }
                  if(dummy_ch1_sn[5] == 0)
                    {
                    nChannel1InelastBackground++;
                    hRecoInteractingKEBackground->Fill(ke);
                    hRecoInteractingKEBackgroundElastic->Fill(ke);
                    hRecoInteractingKEBackgroundInelast->Fill(ke);
                    hRecoInteractingKEChan1Background->Fill(ke);
                    hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                    }
                  }
                }
              }//<--- End Channel 1 loop  

            double total_energy_loss = 0;
            for(int nspts = 0; nspts < ((*col_track_hits)[nTPCtrk] - 1); nspts++)
              {
              if(ch1_had_signal == false){continue;}
              total_energy_loss += (*col_track_dedx)[nTPCtrk][nspts] * (*col_track_pitch_hit)[nTPCtrk][nspts];
              }//<--- End energy delta loop test  
            double ch1_energy_delta = temp_init_ke - total_energy_loss;
            hRecoDeltaEnergy->Fill(ch1_energy_delta);
            hRecoDeltaEnergyvsIntKE->Fill(ch1_energy_delta, ke);
            if((*InteractionPointType)[0] == 3)
              {hRecoDeltaEnergyElastic->Fill(ch1_energy_delta);}
            if((*InteractionPointType)[0] == 13)
              {
              hRecoDeltaEnergyInelast->Fill(ch1_energy_delta);
              hRecoDeltaEnergyInelastvsIntKE->Fill(ch1_energy_delta, ke);
              }
            //std::cout << "      channel 1 delta ke: " << ch1_energy_delta << std::endl;

            }//<--- End if kinked track  
          // ### Straight tracks ###
          if(dummy_spt_array[5] == -1)
            {
            double ch20_endx = (*track_end_x)[nTPCtrk];
            double ch20_endy = (*track_end_y)[nTPCtrk];
            double ch20_endz = (*track_end_z)[nTPCtrk];
            bool straight_punchthru = false;
            if(  ch20_endx < 2   || ch20_endx > 45 
              || ch20_endy < -18 || ch20_endy > 18
              || ch20_endz < 2   || ch20_endz > 88 )
              {straight_punchthru = true;}
            // ## getting mean de/dx ##
            if(did_interact == true && (*InteractionPointType)[0] != 8){nChannel2TopTrue++;}
            double end_dedx_mean = 0; double end_dedx_temp = 0; int end_dedx_counter = 0;
            double end_track_temp = 0;
            for(int pspts = ((*col_track_hits)[nTPCtrk] - 2); pspts > 0; pspts--)
              {
              end_track_temp += (*col_track_pitch_hit)[nTPCtrk][pspts];
              if(end_track_temp < 2.5)
                {
                end_dedx_temp += (*col_track_dedx)[nTPCtrk][pspts];
                end_dedx_counter ++;
                }//<---End if end of track
              }//<--- End backwards blah loop
            end_dedx_mean = (end_dedx_temp) / (end_dedx_counter);
            // ### trying to optimize the end de/dx separation
            if(straight_punchthru == false)
              {
              hRecoEndMeandedx->Fill(end_dedx_mean);
              if((*InteractionPointType)[0] == 13)
                {hRecoEndMeandedxInelast->Fill(end_dedx_mean);}
              if((*InteractionPointType)[0] == 3)
                {hRecoEndMeandedxElastic->Fill(end_dedx_mean);}
              if((*InteractionPointType)[0] != 3 && (*InteractionPointType)[0] != 13)
                {hRecoEndMeandedxNeither->Fill(end_dedx_mean);}
              }

            // ### Channel 2 ###
            if(end_dedx_mean < 13 && straight_punchthru == false)
              {
              nChannel2++;
              double ch2_total_energy_loss = 0;
              bool ch2_signal = false;
              //for(int rspts = 0; rspts < ((*col_track_hits)[nTPCtrk] - 2); rspts++)
              //  {
              //  //hRecoIncidentKE->Fill(ke);
              //  //if(rspts == ((*col_track_hits)[nTPCtrk] - 3) ){hRecoInteractingKE->Fill(ke); hRecoInteractingKEChan2->Fill(ke);}
              //  double energyLoss = (*col_track_dedx)[nTPCtrk][rspts] * (*col_track_pitch_hit)[nTPCtrk][rspts];
              //  ch2_total_energy_loss += energyLoss;
              //  ke -= energyLoss;
              //  }  
              double ch2_endx = (*track_end_x)[nTPCtrk];
              double ch2_endy = (*track_end_y)[nTPCtrk];
              double ch2_endz = (*track_end_z)[nTPCtrk];
              double ch2_energy_delta = temp_init_ke - ch2_total_energy_loss;
              // ## just low mean events ##
              hRecoLowMeanDeltaEnergy->Fill(ch2_energy_delta);
              if(did_interact == true && (*InteractionPointType)[0] != 8)
                {
                nChannel2Signal++;
                double ch2_reco_mc_dR = sqrt( pow( (ch2_endx - mc_int_x), 2) 
                                            + pow( (ch2_endy - mc_int_y), 2) 
                                            + pow( (ch2_endz - mc_int_z), 2) );
                if((*InteractionPointType)[0] == 13)
                  {
                  nChannel2InelastTrue++;
                  if(ch2_reco_mc_dR < 2)
                    {
                    ch2_signal = true;
                    nChannel2InelastSignal++;
                    hRecoLowMeanDeltaEnergySignal->Fill(ch2_energy_delta);
                    }
                  if(ch2_reco_mc_dR >= 2)
                    {
                    ch2_signal = false;
                    nChannel2Background++;
                    nChannel2InelastBackground++;
                    hRecoLowMeanDeltaEnergyBackground->Fill(ch2_energy_delta);
                    } 
                  }
                if((*InteractionPointType)[0] == 3)
                  {
                  nChannel2InelastBackground++;
                  if(ch2_reco_mc_dR < 2)
                    {
                    ch2_signal = true;
                    hRecoLowMeanDeltaEnergySignal->Fill(ch2_energy_delta);
                    }
                  if(ch2_reco_mc_dR >= 2)
                    {
                    ch2_signal = false;
                    nChannel2Background++;
                    hRecoLowMeanDeltaEnergyBackground->Fill(ch2_energy_delta);
                    } 
                  }
                }
              else{ch2_signal = false; nChannel2Background++; nChannel2InelastBackground++;}
              // ## filling channel 2 ke ##
              //std::cout << "        channel 2 investigation.." << std::endl;
              //std::cout << "            number of blah obj: " << (*col_track_hits)[nTPCtrk] << std::endl;
              double before_ke = ke;
              double reco_total_lost = 0;
              double reco_total_dX = 0;
              double path_length = 0;
              for(int rspts = 0; rspts < (*col_track_hits)[nTPCtrk]; rspts++){path_length += (*col_track_pitch_hit)[nTPCtrk][rspts];}
              int num_slabs = (int) path_length;
              //std::cout << "    number of slabs: " << num_slabs << std::endl;
              for(int rspts = 0; rspts < ((*col_track_hits)[nTPCtrk]); rspts++)
                {
                // ### slab change ###
                double energyLoss = (*col_track_dedx)[nTPCtrk][rspts] * (*col_track_pitch_hit)[nTPCtrk][rspts];
                double reco_dX    = (*col_track_pitch_hit)[nTPCtrk][rspts];
                reco_total_lost += energyLoss; 
                reco_total_dX   += reco_dX; 
                hRecodX->Fill(10*reco_dX);
                double simChannel_E = 0;
                for(int simChannel = 0; simChannel < primary_num_simChannel; simChannel++)
                  {
                  if((*primary_simChannel)[simChannel] != (*col_track_wire)[0][rspts]){continue;}
                  simChannel_E = (*primary_simChannel_E)[simChannel];
                  }
                // ### incident! ###
                if(reco_total_dX < num_slabs)
                  { 
                  int middle_slab = ((int) reco_total_dX) + 1;
                  double previous_dX = reco_total_dX - (*col_track_pitch_hit)[nTPCtrk][rspts];
                  double reco_next_dX = reco_total_dX + (*col_track_pitch_hit)[nTPCtrk][rspts+1];
                  double next_ke = ke - (*col_track_dedx)[nTPCtrk][rspts+1]*(*col_track_pitch_hit)[nTPCtrk][rspts+1];
                  if(reco_total_dX < middle_slab && middle_slab < reco_next_dX)
                    {
                    double ke_slope   = (next_ke - ke) / (reco_next_dX - reco_total_dX);
                    double ke_at_slab = ke_slope*(middle_slab - reco_total_dX) + ke;
                    hRecoIncidentKE->Fill(ke_at_slab);
                    hMCRecoIncidentKE->Fill(true_ke);
                    hMCRecoIncidentMatrix->Fill(ke, true_ke);
                    before_ke = ke_at_slab;
                    //std::cout << "          filling incident... " << reco_total_dX << std::endl;
                    //std::cout << "                  ke: " << ke_at_slab << std::endl;
                    }
                  }
                ke -= energyLoss;
                true_ke -= simChannel_E;
                // ### interacting! ###
                if(rspts == (*col_track_hits)[nTPCtrk] -1)
                  {
                  // ### stupid I know ###
                  //ke = interaction_ke;
                  if(ke < 0){ke = before_ke - 5;}
                  //if(ke < 100){std::cout << "       ke: " << ke << std::endl; std::cout << " before ke: " << before_ke << std::endl;}
                  //ke += 25;
                  //ke += 25;
                  //std::cout << "                filling interacting... " << reco_total_dX << std::endl;
                  //std::cout << "                  ke: " << ke << std::endl;
                  hRecoInteractingKE->Fill(ke); 
                  if(interaction_ke == 0){interaction_ke = true_ke;}
                  hMCRecoInteractingKE->Fill(interaction_ke);
                  hMCRecoInteractingMatrix->Fill(ke, interaction_ke);
                  hRecoInteractingKEChan2->Fill(ke);
                  if(ch2_signal == true) 
                    {
                    hRecoMatchedInteractingKECmpsn->Fill(interaction_ke, ke);
                    hRecoInteractingKESignal->Fill(ke);
                    hRecoInteractingKEChan2Signal->Fill(ke);
                    if((*InteractionPointType)[0] == 13)
                      {
                      hRecoMatchedInteractingKECmpsnInelast->Fill(interaction_ke, ke);
                      // ### test ###
                      matched_ke.push_back(std::make_pair(ke, interaction_ke));
                      std::vector< double > ke_tuple ; 
                      ke_tuple.push_back(ke); ke_tuple.push_back(interaction_ke); ke_tuple.push_back(2.5);
                      matched_ke_tuple.push_back(ke_tuple);
                      hRecoInteractingKESignalInelast->Fill(ke);
                      hRecoInteractingKEChan2SignalInelast->Fill(ke);
                      hRecoInteractingKEBackgroundElastic->Fill(ke);
                      hRecoInteractingKEChan2BackgroundElastic->Fill(ke);
                      hMCRecoInteractingKESignalInelast->Fill(interaction_ke);
                      }
                    if((*InteractionPointType)[0] == 3)
                      {
                      hRecoInteractingKESignalElastic->Fill(ke);
                      hRecoInteractingKEChan2SignalElastic->Fill(ke);
                      hRecoInteractingKEBackgroundInelast->Fill(ke);
                      hRecoInteractingKEChan2BackgroundInelast->Fill(ke);
                      hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                      }
                    }
                  if(ch2_signal == false)
                    {
                    hRecoInteractingKEBackground->Fill(ke);
                    hRecoInteractingKEChan2Background->Fill(ke);
                    hRecoInteractingKEBackgroundInelast->Fill(ke);
                    hRecoInteractingKEChan2BackgroundInelast->Fill(ke);
                    hRecoInteractingKEBackgroundElastic->Fill(ke);
                    hRecoInteractingKEChan2BackgroundElastic->Fill(ke);
                    hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                    }
                  }
                }  
              }//<---End if Channel 2

            if(end_dedx_mean >= 13 || straight_punchthru == true)
              {
              nChannel0++;
              if((*InteractionPointType)[0] == 13){nChannel0InelastTrue++;}
              double reco_total_lost = 0; double reco_total_dX = 0;
              std::vector<double> stopping_ke_at_slabs;
              for(int qspts = 0; qspts < ((*col_track_hits)[nTPCtrk] - 2); qspts++)
                {
                double energyLoss = (*col_track_dedx)[nTPCtrk][qspts] * (*col_track_pitch_hit)[nTPCtrk][qspts];
                double reco_dX    = (*col_track_pitch_hit)[nTPCtrk][qspts];
                reco_total_lost += energyLoss; 
                reco_total_dX   += reco_dX; 
                hRecodX->Fill(10*reco_dX);
                int middle_slab = ((int) reco_total_dX) + 1;
                double reco_next_dX = reco_total_dX + (*col_track_pitch_hit)[nTPCtrk][qspts+1];
                double next_ke = ke - (*col_track_dedx)[nTPCtrk][qspts+1]*(*col_track_pitch_hit)[nTPCtrk][qspts+1];
                double simChannel_E = 0;
                for(int simChannel = 0; simChannel < primary_num_simChannel; simChannel++)
                  {
                  if((*primary_simChannel)[simChannel] != (*col_track_wire)[0][qspts]){continue;}
                  simChannel_E = (*primary_simChannel_E)[simChannel];
                  }
                if(reco_total_dX < middle_slab && middle_slab < reco_next_dX)
                  {
                  double ke_slope   = (next_ke - ke) / (reco_next_dX - reco_total_dX);
                  double ke_at_slab = ke_slope*(middle_slab - reco_total_dX) + ke;
                  stopping_ke_at_slabs.push_back(ke_at_slab);
                  hRecoIncidentKE->Fill(ke_at_slab);
                  hMCRecoIncidentKE->Fill(true_ke);
                  hMCRecoIncidentMatrix->Fill(ke, true_ke);
                  }
                ke -= energyLoss;
                true_ke -= simChannel_E;
                }  
              if((*InteractionPointType)[0] == 0 && straight_punchthru == false)
                {
                hRecoDeltaEnergyChannel0->Fill(temp_init_ke - reco_total_lost);
                hRecoTrueVSRecoEnergyChannel0->Fill(temp_init_ke, reco_total_lost);
                if(abs(mc_dX - reco_total_dX) < 2)
                  {
                  //std::cout << "                            for some reason this one worked..." << std::endl;
                  //std::cout << "    g4 ke at slabs! " << g4_ke_at_slabs.size() << std::endl;;
                  for(int i = 1; i <= stopping_ke_at_slabs.size(); i++)
                    {
                    double delta_ke_slabs = stopping_ke_at_slabs[i-1] - stopping_ke_at_slabs[i];
                    //std::cout << "                reco incident ke: " << stopping_ke_at_slabs[i] << std::endl;
                    //std::cout << "                reco incident delta: " << delta_ke_slabs << std::endl;
                    hRecodEdX->Fill(delta_ke_slabs);
                    } 
                  }
                }
              }//<---End if straight stopping track
              
            }//<--- End if straight track

          }//<--- End if only one track
        if(ntracks_reco == 2)
          {
          nTwoReco++;
          if(did_interact == true && (*InteractionPointType)[0] != 8){nTwoTrackTrueSignal++;}
          else nTwoTrackTrueBackground++;
          double other_track = 999.;
          int    ch3_sn = -1;
          if(nTPCtrk == 0){other_track = 1;}
          if(nTPCtrk == 1){other_track = 0;}
          double reco1_reco2_dx = (*track_end_x)[nTPCtrk] - (*track_xpos)[other_track][0];
          double reco1_reco2_dy = (*track_end_y)[nTPCtrk] - (*track_ypos)[other_track][0];
          double reco1_reco2_dz = (*track_end_z)[nTPCtrk] - (*track_zpos)[other_track][0];
          double reco1_reco2_dr = sqrt ( pow(reco1_reco2_dx, 2) + pow(reco1_reco2_dy, 2) + pow(reco1_reco2_dz, 2) ); 
          hRecoDTETrackEndDelta->Fill(reco1_reco2_dr);
          // ### Channel 3 ###
          if(reco1_reco2_dr < 3)
            {
            if(did_interact == true && (*InteractionPointType)[0] != 8)
              {
              nChannel3True++;
              hRecoInteractingKEChan3TrueSignal->Fill(interaction_ke);
              if((*InteractionPointType)[0] == 3)
                {
                if(eangle > 6)
                  {hRecoInteractingKEChan3TrueSignalG6->Fill(interaction_ke);}
                }
              if((*InteractionPointType)[0] == 13)
                {
                nChannel3InelastTrue++;
                if(t13one == 1 && iangle > 6)
                  {hRecoInteractingKEChan3TrueSignalG6->Fill(interaction_ke);}
                if(t13one != 1)
                  {hRecoInteractingKEChan3TrueSignalG6->Fill(interaction_ke);}
                }
              }
            nChannel3++;
            nFirstRegime++;
            double reco1_last[] = {(*track_xpos)[nTPCtrk][(*ntrack_hits)[nTPCtrk]-1], 
                                   (*track_ypos)[nTPCtrk][(*ntrack_hits)[nTPCtrk]-1], 
                                   (*track_zpos)[nTPCtrk][(*ntrack_hits)[nTPCtrk]-1]};
            double reco1_pnul[] = {(*track_xpos)[nTPCtrk][(*ntrack_hits)[nTPCtrk]-2], 
                                   (*track_ypos)[nTPCtrk][(*ntrack_hits)[nTPCtrk]-2], 
                                   (*track_zpos)[nTPCtrk][(*ntrack_hits)[nTPCtrk]-2]};
            double reco1_v[] = {reco1_last[0] - reco1_pnul[0], reco1_last[1] - reco1_pnul[1], reco1_pnul[2] - reco1_pnul[2]};
            double reco2_last[] = {(*track_xpos)[other_track][(*ntrack_hits)[other_track]-1], 
                                   (*track_ypos)[other_track][(*ntrack_hits)[other_track]-1], 
                                   (*track_zpos)[other_track][(*ntrack_hits)[other_track]-1]};
            double reco2_pnul[] = {(*track_xpos)[other_track][(*ntrack_hits)[other_track]-2], 
                                   (*track_ypos)[other_track][(*ntrack_hits)[other_track]-2], 
                                   (*track_zpos)[other_track][(*ntrack_hits)[other_track]-2]};
            double reco2_v[] = {reco2_last[0] - reco2_pnul[0], reco2_last[1] - reco2_pnul[1], reco2_pnul[2] - reco2_pnul[2]};
            double r12_angle = get_angle(reco1_v, reco2_v);
            hRecoDTETrackAngle->Fill(r12_angle);
            if(did_interact == true)
              {
              double delta_g4_reco12 = sqrt(pow(mc_int_x - reco1_last[0], 2) + 
                                            pow(mc_int_y - reco1_last[1], 2) + 
                                            pow(mc_int_z - reco1_last[2], 2));
              if((*InteractionPointType)[0] != 8) 
                { 
                if(delta_g4_reco12 < 3)
                  {
                  //nTwoTrackSignal++; ch3_sn = 1;}

                  // ## blahing `visible` interaction
                  if((*InteractionPointType)[0] == 3 && eangle > 6)
                    {nTwoTrackSignal++; ch3_sn = 1;}
                  if((*InteractionPointType)[0] == 3 && eangle < 6)
                    {nTwoTrackBackground++; ch3_sn = 0;}

                  if((*InteractionPointType)[0] == 13)
                    {
                    //if(t13one == 1 && iangle > 6)
                    //  {nTwoTrackSignal++; ch3_sn = 1;}
                    //if(t13one == 1 && iangle < 6)
                    //  {nTwoTrackBackground++; ch3_sn = 0;}
                    //else
                    //  {nTwoTrackSignal++; ch3_sn = 1;}
                    nTwoTrackSignal++;
                    ch3_sn = 1;
                    }
                  }

                else if(delta_g4_reco12 >= 3){nTwoTrackBackground++; ch3_sn = 0;}
                }
              // ### mb adding back missing background ###
              else{nTwoTrackBackground++; ch3_sn = 0;}

              hRecoDTETrackAngleInts->Fill(r12_angle);
              if( (*InteractionPointType)[0] == 13 ) 
                {hRecoDTETrackAngleIntsInelast->Fill(r12_angle);}
              if( (*InteractionPointType)[0] == 8 ) 
                {hRecoDTETrackAngleIntsCoulomb->Fill(r12_angle);}
              if( (*InteractionPointType)[0] == 3 ) 
                {hRecoDTETrackAngleIntsElastic->Fill(r12_angle);}
              }
            if(did_interact == false)
              {hRecoDTETrackAngleNoInts->Fill(r12_angle); nTwoTrackBackground++; ch3_sn = 0;}
            
            bool kick2 = false;
            double path_length = 0;
            double reco_total_lost = 0;
            double reco_total_dX = 0;
            double before_ke = 0;
            for(int rspts = 0; rspts < (*col_track_hits)[nTPCtrk]; rspts++){path_length += (*col_track_pitch_hit)[nTPCtrk][rspts];}
            int num_slabs = (int) path_length;
            // ### filling channel 3 ke ###
            for(int nspts = 0; nspts < (*col_track_hits)[nTPCtrk]; nspts++)
              {
              if(kick2 == true){continue;}
              double simChannel_E = 0;
              for(int simChannel = 0; simChannel < primary_num_simChannel; simChannel++)
                {
                if((*primary_simChannel)[simChannel] != (*col_track_wire)[0][nspts]){continue;}
                simChannel_E = (*primary_simChannel_E)[simChannel];
                }
              // ### slab change ###
              double energyLoss = (*col_track_dedx)[nTPCtrk][nspts] * (*col_track_pitch_hit)[nTPCtrk][nspts];
              double reco_dX    = (*col_track_pitch_hit)[nTPCtrk][nspts];
              reco_total_lost += energyLoss; 
              reco_total_dX   += reco_dX; 
              hRecodX->Fill(10*reco_dX);
              // ### incident! ###
              if(reco_total_dX < num_slabs)
                { 
                int middle_slab = ((int) reco_total_dX) + 1;
                double previous_dX = reco_total_dX - (*col_track_pitch_hit)[nTPCtrk][nspts];
                double reco_next_dX = reco_total_dX + (*col_track_pitch_hit)[nTPCtrk][nspts+1];
                double next_ke = ke - (*col_track_dedx)[nTPCtrk][nspts+1]*(*col_track_pitch_hit)[nTPCtrk][nspts+1];
                if(reco_total_dX < middle_slab && middle_slab < reco_next_dX)
                  {
                  double ke_slope   = (next_ke - ke) / (reco_next_dX - reco_total_dX);
                  double ke_at_slab = ke_slope*(middle_slab - reco_total_dX) + ke;
                  hRecoIncidentKE->Fill(ke_at_slab);
                  hMCRecoIncidentKE->Fill(true_ke);
                  hMCRecoIncidentMatrix->Fill(ke, true_ke);
                  before_ke = ke_at_slab;
                  }
                }
              ke -= energyLoss;
              true_ke -= simChannel_E;

              // ### interacting! ###
              if(nspts == ((*col_track_hits)[nTPCtrk] - 1))
                {
                // ### stupid I know ###
                if(ke < 0){ke = before_ke - 5;}
                //ke += 25;
                //ke = interaction_ke;
                //if(ke < 0){ke = 1.;}
                //ke += 25;
                kick2 = true; 
                hRecoInteractingKE->Fill(ke);
                if(interaction_ke == 0){interaction_ke = true_ke;}
                hMCRecoInteractingKE->Fill(interaction_ke);
                hMCRecoInteractingMatrix->Fill(ke, interaction_ke);
                hRecoInteractingKEChan3->Fill(ke);
                //if(did_interact == true && (*InteractionPointType)[0] != 8){hRecoInteractingKESignal->Fill(ke);}
                //else{hRecoInteractingKEBackground->Fill(ke);}
                if(ch3_sn == 1)
                  {
                  hRecoInteractingKESignal->Fill(ke);
                  hRecoMatchedInteractingKERatio->Fill(ke / interaction_ke);
                  hRecoMatchedInteractingKECmpsn->Fill(interaction_ke, ke);
                  hRecoInteractingKEChan3Signal->Fill(ke);
                  if((*InteractionPointType)[0] == 3) 
                    {
                    nChannel3InelastBackground++;
                    hRecoInteractingKEBackgroundInelast->Fill(ke);
                    hRecoInteractingKESignalElastic->Fill(ke); 
                    hRecoInteractingKEChan3SignalElastic->Fill(ke);
                    hRecoInteractingKEChan3BackgroundInelast->Fill(ke);
                    hRecoMatchedElasticAngleCmpsn->Fill(r12_angle, eangle);
                    if(eangle > 6)
                      {hRecoInteractingKESignalElasticObs->Fill(ke);}
                    hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                    }
                  if((*InteractionPointType)[0] == 13)
                    {
                    nChannel3InelastSignal++;
                    hRecoMatchedInteractingKECmpsnInelast->Fill(interaction_ke, ke);
                      // ### test ###
                      matched_ke.push_back(std::make_pair(ke, interaction_ke));
                      std::vector< double > ke_tuple ; 
                      ke_tuple.push_back(ke); ke_tuple.push_back(interaction_ke); ke_tuple.push_back(3.5);
                      matched_ke_tuple.push_back(ke_tuple);
                    hRecoInteractingKEBackgroundElastic->Fill(ke);
                    hRecoInteractingKESignalInelast->Fill(ke);
                    hRecoInteractingKEChan3SignalInelast->Fill(ke);
                    hRecoInteractingKEChan3BackgroundElastic->Fill(ke);
                    hMCRecoInteractingKESignalInelast->Fill(interaction_ke);
                    }
                  }
                else if(ch3_sn == 0)
                  { 
                  nChannel3InelastBackground++;
                  hRecoInteractingKEBackground->Fill(ke);
                  hRecoInteractingKEBackgroundElastic->Fill(ke);
                  hRecoInteractingKEBackgroundInelast->Fill(ke);
                  hRecoInteractingKEChan3Background->Fill(ke);
                  hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                  }
                // ### fill by interacting by signal/background ###
                }
              }//<---Channel 3  

            }//<--- End if elastic or inelastic candidate
          // ### Channel 4 ###
          if(reco1_reco2_dr >= 3)
            {
            if(did_interact == true && (*InteractionPointType)[0] != 8)
              {
              nChannel4True++;
              hRecoInteractingKEChan4TrueSignal->Fill(interaction_ke);
              if((*InteractionPointType)[0] == 3)
                {
                if(eangle > 6)
                  {hRecoInteractingKEChan4TrueSignalG6->Fill(interaction_ke);}
                }
              if((*InteractionPointType)[0] == 13)
                {
                nChannel4InelastTrue++;
                if(t13one == 1 && iangle > 6)
                  {hRecoInteractingKEChan4TrueSignalG6->Fill(interaction_ke);}
                if(t13one != 1)
                  {hRecoInteractingKEChan4TrueSignalG6->Fill(interaction_ke);}
                }
              }
            //nChannel4++;
            nOtherRegime++;
            int ch4_sn = -1;
            int min_to_reco2 = 999; double large_number = 999; 
            for(size_t nspts = 0; nspts < (*ntrack_hits)[nTPCtrk]; nspts++)
              {
              double mid_reco2_dx = (*track_xpos)[nTPCtrk][nspts] - (*track_xpos)[other_track][0];
              double mid_reco2_dy = (*track_ypos)[nTPCtrk][nspts] - (*track_ypos)[other_track][0];
              double mid_reco2_dz = (*track_zpos)[nTPCtrk][nspts] - (*track_zpos)[other_track][0];
              double mid_reco2_dr = sqrt ( pow(mid_reco2_dx, 2) + pow(mid_reco2_dy, 2) + pow(mid_reco2_dz, 2) ); 
              if(mid_reco2_dr < large_number){large_number = mid_reco2_dr; min_to_reco2 = nspts;}
              }
            if(did_interact == true)
              {
              double mc_int_reco_dx = mc_int_x - (*track_xpos)[nTPCtrk][min_to_reco2];
              double mc_int_reco_dy = mc_int_y - (*track_ypos)[nTPCtrk][min_to_reco2];
              double mc_int_reco_dz = mc_int_z - (*track_zpos)[nTPCtrk][min_to_reco2];
              double mc_int_reco_dr = sqrt ( pow(mc_int_reco_dx, 2) + pow(mc_int_reco_dy, 2) + pow(mc_int_reco_dz, 2) ); 
              hRecoDTEDistanceToMCIntFromTrack2->Fill(mc_int_reco_dr);
              double mc_int_reco_dxo = mc_int_x - (*track_xpos)[other_track][0];
              double mc_int_reco_dyo = mc_int_y - (*track_ypos)[other_track][0];
              double mc_int_reco_dzo = mc_int_z - (*track_zpos)[other_track][0];
              double mc_int_reco_dro = sqrt ( pow(mc_int_reco_dxo, 2) + pow(mc_int_reco_dyo, 2) + pow(mc_int_reco_dzo, 2) ); 
              if((*InteractionPointType)[0] != 8)
                {
                //nOtherTwoTrackTrueSignal++;
                nSanityCheck++;
                //std::cout << "      ch4 g4 reco dist: " << mc_int_reco_dr << std::endl;
                if(mc_int_reco_dr  < 2)
                  {
                  //nOtherTwoTrackSignal++; nChannel4++; ch4_sn = 1;}
                  nChannel4++;
                  // ## blahing `visible` interaction
                  if((*InteractionPointType)[0] == 3 && eangle > 6)
                    {nOtherTwoTrackSignal++; ch4_sn = 1;}
                  if((*InteractionPointType)[0] == 3 && eangle < 6)
                    {nOtherTwoTrackBackground++; ch4_sn = 0;}

                  if((*InteractionPointType)[0] == 13)
                    {
                    nOtherTwoTrackSignal++;
                    ch4_sn = 1;
                    }
                  }
                // ### tryna fix something ###
                if(mc_int_reco_dr >= 2)
                  {nOtherTwoTrackBackground++; nChannel4++; ch4_sn = 0;}
                }
              else if((*InteractionPointType)[0] == 8){nOtherTwoTrackBackground++; nChannel4++; ch4_sn = 0;}
              }
            if(did_interact == false){nOtherTwoTrackBackground++; nChannel4++; ch4_sn = 0;}
            bool kick3 = false;
            double midph = 0;
            double midcalph = 0; double calph = 0;
            for(size_t mspts = 1; mspts < min_to_reco2; mspts++)
              {
              midph += sqrt(  pow( (*track_xpos)[nTPCtrk][mspts-1] - (*track_xpos)[nTPCtrk][mspts], 2)
                            + pow( (*track_ypos)[nTPCtrk][mspts-1] - (*track_ypos)[nTPCtrk][mspts], 2)
                            + pow( (*track_zpos)[nTPCtrk][mspts-1] - (*track_zpos)[nTPCtrk][mspts], 2) );
              }  
            // ### filling channel 4 ke ###
            double reco_total_lost = 0; double reco_total_dX = 0;
            double before_ke = 0;
            for(int nspts = 0; nspts < (*col_track_hits)[nTPCtrk]; nspts++)
              {
              if(kick3 == true){continue;}
              double simChannel_E = 0;
              for(int simChannel = 0; simChannel < primary_num_simChannel; simChannel++)
                {
                if((*primary_simChannel)[simChannel] != (*col_track_wire)[0][nspts]){continue;}
                simChannel_E = (*primary_simChannel_E)[simChannel];
                }
              // ### slab change ###
              if( (midph - calph) > 1.)
                {
                double energyLoss = (*col_track_dedx)[nTPCtrk][nspts] * (*col_track_pitch_hit)[nTPCtrk][nspts];
                double reco_dX    = (*col_track_pitch_hit)[nTPCtrk][nspts];
                reco_total_lost += energyLoss; 
                reco_total_dX   += reco_dX; 
                hRecodX->Fill(10*reco_dX);
                int middle_slab = ((int) reco_total_dX) + 1;
                double reco_next_dX = reco_total_dX + (*col_track_pitch_hit)[nTPCtrk][nspts+1];
                double next_ke = ke - (*col_track_dedx)[nTPCtrk][nspts+1]*(*col_track_pitch_hit)[nTPCtrk][nspts+1];
                if(reco_total_dX < middle_slab && middle_slab < reco_next_dX)
                  {
                  double ke_slope   = (next_ke - ke) / (reco_next_dX - reco_total_dX);
                  double ke_at_slab = ke_slope*(middle_slab - reco_total_dX) + ke;
                  hRecoIncidentKE->Fill(ke_at_slab);
                  hMCRecoIncidentKE->Fill(true_ke);
                  hMCRecoIncidentMatrix->Fill(ke, true_ke);
                  before_ke = ke_at_slab;
                  //std::cout << "        filling incident: " << reco_total_dX << std::endl;
                  //std::cout << "            ke: "           << ke_at_slab    << std::endl;
                  }
                calph += (*col_track_pitch_hit)[nTPCtrk][nspts];
                ke -= energyLoss;
                true_ke -= simChannel_E;
                }
              else
                {
                // ### stupid I know ###
                if(ke < 0){ke = before_ke - 5;}
                  //std::cout << "            filling interacting: " << reco_total_dX << std::endl;
                  //std::cout << "                ke: "           << ke    << std::endl;
                //ke += 25;
                //ke = interaction_ke;
                kick3 = true;
                hRecoInteractingKE->Fill(ke);
                if(interaction_ke == 0){interaction_ke = true_ke;}
                hMCRecoInteractingKE->Fill(interaction_ke);
                hMCRecoInteractingMatrix->Fill(ke, interaction_ke);
                hRecoInteractingKEChan4->Fill(ke);
                if(ch4_sn == 1)
                  {
                  hRecoMatchedInteractingKECmpsn->Fill(interaction_ke, ke);
                  hRecoInteractingKESignal->Fill(ke);
                  hRecoMatchedInteractingKERatio->Fill(ke / interaction_ke);
                  hRecoInteractingKEChan4Signal->Fill(ke);
                  if((*InteractionPointType)[0] == 3) 
                    {
                    nChannel4InelastBackground++;
                    hRecoInteractingKEBackgroundInelast->Fill(ke);
                    hRecoInteractingKESignalElastic->Fill(ke);
                    hRecoInteractingKEChan4SignalElastic->Fill(ke);
                    hRecoInteractingKEChan4BackgroundInelast->Fill(ke);
                    if(eangle > 6)
                      {hRecoInteractingKESignalElasticObs->Fill(ke);}
                    hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                    }
                  if((*InteractionPointType)[0] == 13)
                    {
                    hRecoMatchedInteractingKECmpsnInelast->Fill(interaction_ke, ke);
                      // ### test ###
                      matched_ke.push_back(std::make_pair(ke, interaction_ke));
                      std::vector< double > ke_tuple ; 
                      ke_tuple.push_back(ke); ke_tuple.push_back(interaction_ke); ke_tuple.push_back(4.5);
                      matched_ke_tuple.push_back(ke_tuple);
                    nChannel4InelastSignal++;
                    hRecoInteractingKEBackgroundElastic->Fill(ke);
                    hRecoInteractingKESignalInelast->Fill(ke); 
                    hRecoInteractingKEChan4SignalInelast->Fill(ke);
                    hRecoInteractingKEChan4BackgroundElastic->Fill(ke);
                    hMCRecoInteractingKESignalInelast->Fill(interaction_ke);
                    }
                  }
                else if(ch4_sn == 0)
                  {
                  nChannel4InelastBackground++;
                  hRecoInteractingKEBackground->Fill(ke);
                  hRecoInteractingKEBackgroundElastic->Fill(ke);
                  hRecoInteractingKEBackgroundInelast->Fill(ke);
                  hRecoInteractingKEChan4Background->Fill(ke);
                  hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                  }
                }
              //double energyLoss = (*col_track_dedx)[nTPCtrk][nspts] * (*col_track_pitch_hit)[nTPCtrk][nspts];
              //ke -= energyLoss;
              //double reco_dX    = (*col_track_pitch_hit)[nTPCtrk][nspts];
              //hRecodX->Fill(10*reco_dX);
              }//<--- Channel 4  
            
            }//<--End if inelastic candidate
          }//<--- End if only two tracks


        if(ntracks_reco == 3)
          {
          int track2 = 999;
          int track3 = 999;
          int ch4_sn = -1;
          if(nTPCtrk == 0){track2 = 1; track3 = 2;}
          if(nTPCtrk == 1){track2 = 0; track3 = 2;}
          if(nTPCtrk == 2){track2 = 0; track3 = 1;}
          double track1_endx   = (*track_end_x)[nTPCtrk];  double track1_endy   = (*track_end_y)[nTPCtrk];  double track1_endz   = (*track_end_z)[nTPCtrk];
          double track2_startx = (*track_xpos)[track2][0]; double track2_starty = (*track_ypos)[track2][0]; double track2_startz = (*track_zpos)[track2][0];
          double track3_startx = (*track_xpos)[track3][0]; double track3_starty = (*track_ypos)[track3][0]; double track3_startz = (*track_zpos)[track3][0];
          //std::cout << "Reco event, 3 tracks... " << std::endl;
          double reco23_dr = sqrt(pow( track2_startx - track3_startx, 2) + 
                                  pow( track2_starty - track3_starty, 2) + 
                                  pow( track2_startz - track3_startz, 2)); 
          double reco12_dr = sqrt(pow( track2_startx - track1_endx, 2) + 
                                  pow( track2_starty - track1_endy, 2) + 
                                  pow( track2_startz - track1_endz, 2)); 
          double reco13_dr = sqrt(pow( track3_startx - track1_endx, 2) + 
                                  pow( track3_starty - track1_endy, 2) + 
                                  pow( track3_startz - track1_endz, 2)); 
          //std::cout << "    delta(start2,start3) :" << reco23_dr << std::endl;
          //std::cout << "        delta(start2,end1) :" << reco12_dr << std::endl;
          //std::cout << "        delta(start3,end1) :" << reco13_dr << std::endl;
          if((*InteractionPointType)[0] == 13)
            {
            //std::cout << "          there was an inelastic scatter " << std::endl;
            double track23_mean_startx = (track2_startx + track3_startx) / 2;
            double track23_mean_starty = (track2_starty + track3_starty) / 2;
            double track23_mean_startz = (track2_startz + track3_startz) / 2;
            double reco23_g4_dr = sqrt(pow( track23_mean_startx - mc_int_x, 2) + 
                                       pow( track23_mean_starty - mc_int_y, 2) + 
                                       pow( track23_mean_startz - mc_int_z, 2)); 
            hRecoTTE23vsg4->Fill(reco23_dr, reco23_g4_dr);
            //std::cout << "              mc interaction - (avg start of 2,3): " << reco23_g4_dr << std::endl;
            }
          // ### Channel 3 ###
          if(reco23_dr > 4)
            {
            // ### end to end tracks followed by end to end tracks ###
            nChannel3++;
            double closest_primary_mc_int = 999;
            double min_dist12 = 999;     double min_dist13 = 999;
            int    min_dist12_spt = 999; int    min_dist13_spt = 999;
            int ch33_sn = -1;
            int first_int_spt = 999;
            for(size_t nspts = 0; nspts < (*ntrack_hits)[nTPCtrk]; nspts++)
              {
              double track1_xpos = (*track_xpos)[nTPCtrk][nspts];
              double track1_ypos = (*track_ypos)[nTPCtrk][nspts];
              double track1_zpos = (*track_zpos)[nTPCtrk][nspts];
              double dist_temp12 = sqrt( pow(track1_xpos - track2_startx, 2) 
                                       + pow(track1_ypos - track2_starty, 2) 
                                       + pow(track1_zpos - track2_startz, 2) );
              double dist_temp13 = sqrt( pow(track1_xpos - track3_startx, 2) 
                                       + pow(track1_ypos - track3_starty, 2) 
                                       + pow(track1_zpos - track3_startz, 2) );

              if(dist_temp12 < min_dist12)
                {min_dist12 = dist_temp12; min_dist12_spt = nspts;}
              if(dist_temp13 < min_dist13)
                {min_dist13 = dist_temp13; min_dist13_spt = nspts;}
              }
            if( min_dist12_spt - min_dist13_spt > 0) {first_int_spt = min_dist13_spt;} 
            if( min_dist12_spt - min_dist13_spt <= 0){first_int_spt = min_dist12_spt;} 
            

            // ### truth info ###
            if(did_interact == true && (*InteractionPointType)[0] != 8)
              {
              nChannel3InelastTrue++;
              closest_primary_mc_int = sqrt(pow( (*track_xpos)[nTPCtrk][first_int_spt] - mc_int_x, 2) 
                                          + pow( (*track_ypos)[nTPCtrk][first_int_spt] - mc_int_y, 2) 
                                          + pow( (*track_zpos)[nTPCtrk][first_int_spt] - mc_int_z, 2));
              if(closest_primary_mc_int < 2){ch33_sn = 1;}
              else{ch33_sn = 0;}
              }
            else{ch33_sn = 0;}

            // ### reco info  ###
            bool kick33 = false;
            double midph = 0;
            double calph = 0;
            double reco_total_lost = 0; double reco_total_dX = 0;
            double before_ke = 0;
            for(size_t mspts = 1; mspts < first_int_spt; mspts++)
              {
              midph += sqrt(  pow( (*track_xpos)[nTPCtrk][mspts-1] - (*track_xpos)[nTPCtrk][mspts], 2)
                            + pow( (*track_ypos)[nTPCtrk][mspts-1] - (*track_ypos)[nTPCtrk][mspts], 2)
                            + pow( (*track_zpos)[nTPCtrk][mspts-1] - (*track_zpos)[nTPCtrk][mspts], 2) );
              }  
            for(int nspts = 0; nspts < (*col_track_hits)[nTPCtrk]; nspts++)
              {
              if(kick33 == true){continue;}
              double simChannel_E = 0;
              for(int simChannel = 0; simChannel < primary_num_simChannel; simChannel++)
                {
                if((*primary_simChannel)[simChannel] != (*col_track_wire)[0][nspts]){continue;}
                simChannel_E = (*primary_simChannel_E)[simChannel];
                }
              // ### slab change ###
              //hRecoIncidentKE->Fill(ke);
              if( (midph - calph) > 1.)
                {
                double energyLoss = (*col_track_dedx)[nTPCtrk][nspts] * (*col_track_pitch_hit)[nTPCtrk][nspts];
                double reco_dX    = (*col_track_pitch_hit)[nTPCtrk][nspts];
                reco_total_lost += energyLoss; 
                reco_total_dX   += reco_dX; 
                hRecodX->Fill(10*reco_dX);
                int middle_slab = ((int) reco_total_dX) + 1;
                double reco_next_dX = reco_total_dX + (*col_track_pitch_hit)[nTPCtrk][nspts+1];
                double next_ke = ke - (*col_track_dedx)[nTPCtrk][nspts+1]*(*col_track_pitch_hit)[nTPCtrk][nspts+1];
                if(reco_total_dX < middle_slab && middle_slab < reco_next_dX)
                  {
                  double ke_slope   = (next_ke - ke) / (reco_next_dX - reco_total_dX);
                  double ke_at_slab = ke_slope*(middle_slab - reco_total_dX) + ke;
                  hRecoIncidentKE->Fill(ke_at_slab);
                  hMCRecoIncidentKE->Fill(true_ke);
                  hMCRecoIncidentMatrix->Fill(ke, true_ke);
                  before_ke = ke_at_slab;
                  }
                calph += (*col_track_pitch_hit)[nTPCtrk][nspts];
                ke -= energyLoss;
                true_ke -= simChannel_E;
                }
              else
                {
                // ### stupid I know ###
                if(ke < 0){ke = before_ke - 5;}
                //ke += 25;
                //ke = interaction_ke;
                kick33 = true;
                hRecoInteractingKE->Fill(ke);
                if(interaction_ke == 0){interaction_ke = true_ke;}
                hMCRecoInteractingKE->Fill(interaction_ke);
                hMCRecoInteractingMatrix->Fill(ke, interaction_ke);
                hRecoInteractingKEChan3->Fill(ke);
                if(ch33_sn == 1)
                  {
                  hRecoMatchedInteractingKECmpsn->Fill(interaction_ke, ke);
                  hRecoInteractingKESignal->Fill(ke);
                  hRecoMatchedInteractingKERatio->Fill(ke / interaction_ke);
                  hRecoInteractingKEChan3Signal->Fill(ke);
                  if((*InteractionPointType)[0] == 3) 
                    {
                    nChannel3InelastBackground++;
                    hRecoInteractingKEBackgroundInelast->Fill(ke);
                    hRecoInteractingKESignalElastic->Fill(ke);
                    hRecoInteractingKEChan3SignalElastic->Fill(ke);
                    hRecoInteractingKEChan3BackgroundInelast->Fill(ke);
                    hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                    }
                  if((*InteractionPointType)[0] == 13)
                    {
                    hRecoMatchedInteractingKECmpsnInelast->Fill(interaction_ke, ke);
                    // ### test ###
                    matched_ke.push_back(std::make_pair(ke, interaction_ke));
                    std::vector< double > ke_tuple ; 
                    ke_tuple.push_back(ke); ke_tuple.push_back(interaction_ke); ke_tuple.push_back(4.5);
                    matched_ke_tuple.push_back(ke_tuple);
                    nChannel3InelastSignal++;
                    hRecoInteractingKEBackgroundElastic->Fill(ke);
                    hRecoInteractingKESignalInelast->Fill(ke); 
                    hRecoInteractingKEChan3SignalInelast->Fill(ke);
                    hRecoInteractingKEChan3BackgroundElastic->Fill(ke);
                    hMCRecoInteractingKESignalInelast->Fill(interaction_ke);
                    }
                  }
                else if(ch33_sn == 0)
                  {
                  nChannel3InelastBackground++;
                  hRecoInteractingKEBackground->Fill(ke);
                  hRecoInteractingKEBackgroundElastic->Fill(ke);
                  hRecoInteractingKEBackgroundInelast->Fill(ke);
                  hRecoInteractingKEChan3Background->Fill(ke);
                  hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                  }
                }
              }//<--- Channel 3  


            }//<--- End broken up track topology

          // ### more channel 4 ###
          if(reco23_dr <= 4)
            {
            nChannel4++;
            double end_primary_mc_int = 999;
            double min_dist23 = 999; int min_dist23_spt = 999;
            int ch34_sn = -1;
            for(size_t nspts = 0; nspts < (*ntrack_hits)[nTPCtrk]; nspts++)
              {
              double track1_xpos = (*track_xpos)[nTPCtrk][nspts];
              double track1_ypos = (*track_ypos)[nTPCtrk][nspts];
              double track1_zpos = (*track_zpos)[nTPCtrk][nspts];
              double dist_temp   = sqrt( pow(track1_xpos - track2_startx, 2) 
                                       + pow(track1_ypos - track2_starty, 2) 
                                       + pow(track1_zpos - track2_startz, 2) );

              if(dist_temp < min_dist23)
                {min_dist23 = dist_temp; min_dist23_spt = nspts;}
              }
            // ### truth info ###
            if(did_interact == true && (*InteractionPointType)[0] != 8)
              {
              nChannel4InelastTrue++;
              end_primary_mc_int = sqrt(pow( (*track_xpos)[nTPCtrk][min_dist23_spt] - mc_int_x, 2) 
                                      + pow( (*track_ypos)[nTPCtrk][min_dist23_spt] - mc_int_y, 2) 
                                      + pow( (*track_zpos)[nTPCtrk][min_dist23_spt] - mc_int_z, 2));
              if(end_primary_mc_int < 2){ch34_sn = 1;}
              else{ch34_sn = 0;}
              }
            else{ch34_sn = 0;}

            // ### reco info  ###
            bool kick4 = false;
            double midph = 0;
            double calph = 0;
            double reco_total_lost = 0; double reco_total_dX = 0;
            double before_ke = 0;
            for(size_t mspts = 1; mspts < min_dist23_spt; mspts++)
              {
              midph += sqrt(  pow( (*track_xpos)[nTPCtrk][mspts-1] - (*track_xpos)[nTPCtrk][mspts], 2)
                            + pow( (*track_ypos)[nTPCtrk][mspts-1] - (*track_ypos)[nTPCtrk][mspts], 2)
                            + pow( (*track_zpos)[nTPCtrk][mspts-1] - (*track_zpos)[nTPCtrk][mspts], 2) );
              }  
            for(int nspts = 0; nspts < (*col_track_hits)[nTPCtrk]; nspts++)
              {
              if(kick4 == true){continue;}
              double simChannel_E = 0;
              for(int simChannel = 0; simChannel < primary_num_simChannel; simChannel++)
                {
                if((*primary_simChannel)[simChannel] != (*col_track_wire)[0][nspts]){continue;}
                simChannel_E = (*primary_simChannel_E)[simChannel];
                }
              // ### slab change ###
              if( (midph - calph) > 1.)
                {
                double energyLoss = (*col_track_dedx)[nTPCtrk][nspts] * (*col_track_pitch_hit)[nTPCtrk][nspts];
                double reco_dX    = (*col_track_pitch_hit)[nTPCtrk][nspts];
                reco_total_lost += energyLoss; 
                reco_total_dX   += reco_dX; 
                hRecodX->Fill(10*reco_dX);
                int middle_slab = ((int) reco_total_dX) + 1;
                double reco_next_dX = reco_total_dX + (*col_track_pitch_hit)[nTPCtrk][nspts+1];
                double next_ke = ke - (*col_track_dedx)[nTPCtrk][nspts+1]*(*col_track_pitch_hit)[nTPCtrk][nspts+1];
                if(reco_total_dX < middle_slab && middle_slab < reco_next_dX)
                  {
                  double ke_slope   = (next_ke - ke) / (reco_next_dX - reco_total_dX);
                  double ke_at_slab = ke_slope*(middle_slab - reco_total_dX) + ke;
                  hRecoIncidentKE->Fill(ke_at_slab);
                  hMCRecoIncidentKE->Fill(true_ke);
                  hMCRecoIncidentMatrix->Fill(ke, true_ke);
                  before_ke = ke_at_slab;
                  }
                calph += (*col_track_pitch_hit)[nTPCtrk][nspts];
                ke -= energyLoss;
                true_ke -= simChannel_E;
                }
              else
                {
                // ### stupid I know ###
                if(ke < 0){ke = before_ke - 5;}
                //ke += 25;
                //ke = interaction_ke;
                kick4 = true;
                hRecoInteractingKE->Fill(ke);
                if(interaction_ke == 0){interaction_ke = true_ke;}
                hMCRecoInteractingKE->Fill(interaction_ke);
                hMCRecoInteractingMatrix->Fill(ke, interaction_ke);
                hRecoInteractingKEChan4->Fill(ke);
                if(ch34_sn == 1)
                  {
                  hRecoMatchedInteractingKECmpsn->Fill(interaction_ke, ke);
                  hRecoInteractingKESignal->Fill(ke);
                  hRecoMatchedInteractingKERatio->Fill(ke / interaction_ke);
                  hRecoInteractingKEChan4Signal->Fill(ke);
                  if((*InteractionPointType)[0] == 3) 
                    {
                    nChannel4InelastBackground++;
                    hRecoInteractingKEBackgroundInelast->Fill(ke);
                    hRecoInteractingKESignalElastic->Fill(ke);
                    hRecoInteractingKEChan4SignalElastic->Fill(ke);
                    hRecoInteractingKEChan4BackgroundInelast->Fill(ke);
                    hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                    }
                  if((*InteractionPointType)[0] == 13)
                    {
                    hRecoMatchedInteractingKECmpsnInelast->Fill(interaction_ke, ke);
                    // ### test ###
                    matched_ke.push_back(std::make_pair(ke, interaction_ke));
                    std::vector< double > ke_tuple ; 
                    ke_tuple.push_back(ke); ke_tuple.push_back(interaction_ke); ke_tuple.push_back(4.5);
                    matched_ke_tuple.push_back(ke_tuple);
                    nChannel4InelastSignal++;
                    hRecoInteractingKEBackgroundElastic->Fill(ke);
                    hRecoInteractingKESignalInelast->Fill(ke); 
                    hRecoInteractingKEChan4SignalInelast->Fill(ke);
                    hRecoInteractingKEChan4BackgroundElastic->Fill(ke);
                    hMCRecoInteractingKESignalInelast->Fill(interaction_ke);
                    }
                  }
                else if(ch34_sn == 0)
                  {
                  nChannel4InelastBackground++;
                  hRecoInteractingKEBackground->Fill(ke);
                  hRecoInteractingKEBackgroundElastic->Fill(ke);
                  hRecoInteractingKEBackgroundInelast->Fill(ke);
                  hRecoInteractingKEChan4Background->Fill(ke);
                  hMCRecoInteractingKEBackgroundInelast->Fill(interaction_ke);
                  }
                // ### fill by signal noise ###
                }
              }//<--- Channel 4  
            }//<--- Fork topology


          }//<--- End if three tracks
        }//<--- End reco track loop
    }//<---End MC primary loop 
}//<---End jentry loop  

//// ### end of loop printouts ###

std::cout << "Events... " << std::endl;
std::cout << "nTotal Events... " << nTotalEvents << std::endl;
std::cout << "n Events with tracks / end point is in tpc... " << nEvntsTPC << std::endl;
std::cout << "n Simple MC Topologies... " << nSimpleMCTopo << std::endl;
//std::cout << "n Interactions after all cuts... " << nInteractions_afterAllCuts << std::endl;
std::cout << std::endl;
std::cout << "n Total  Reco Topologies... " << nTotalReco << std::endl;
std::cout << "n Simple Reco Topologies... " << nSimpleRecoTopo << std::endl;
std::cout << "n More than 3 Tracsk Rec... " << nMoreThan3      << std::endl;
std::cout << "n Exactly 4 Tracsk Rec... " << nExactly4      << std::endl;
std::cout << "n Reco Events... " << nRecoMatched << std::endl;
std::cout << "n One Reco Events... " << nOneReco << std::endl;
std::cout << "n Two Reco Events... " << nTwoReco << std::endl;
//std::cout << "n Two Reco Interaction Events... " << nTwoTrackSignal << std::endl;
//std::cout << "  n One Reco Interaction Fakes.... " << nOneTrackBackground << std::endl;
//std::cout << "  n One Reco Interaction Events... " << nOneTrackSignal << std::endl;
//std::cout << "looping over bins for one track signal/noise... " << std::endl;
std::cout << "  true signal: " << nOneTrackTrueSignal << std::endl;
for(int i = 0; i < 10; i++)
  {
  double angle_selection_purity = 1.*(nOneTrackSignalArray[i]) / (nOneTrackSignalArray[i] + nOneTrackBackgroundArray[i] );
  double angle_selection_effish = 1.*(nOneTrackSignalArray[i]) / nOneTrackTrueSignal;
  double snratio = 1.*(nOneTrackSignalArray[i]) / nOneTrackBackgroundArray[i];
  if(i == 6)
    {
    std::cout << "    bin " << i << " ... signal: " << nOneTrackSignalArray[i]     << std::endl;
    std::cout << "    bin " << i << " ... noise : " << nOneTrackBackgroundArray[i] << std::endl;
    std::cout << "    bin " << i << " ... purity: " << angle_selection_purity      << std::endl;
    std::cout << "    bin " << i << " ... effish: " << angle_selection_effish      << std::endl;
    std::cout << "    bin " << i << " ... s/n:    " << snratio                     << std::endl;
    std::cout << std::endl;
    }
  hXSIntAnglePurity->Fill(i, angle_selection_purity);
  hXSIntAngleEffish->Fill(i, angle_selection_effish);
  }

double norm = 1;
hMCNormalizedAngleSignal->Scale(norm/hMCNormalizedAngleSignal->Integral(), "width");
hMCNormalizedAngleBackground->Scale(norm/hMCNormalizedAngleBackground->Integral(), "width");
hRecoNormalizedAngleSignal->Scale(norm/hRecoNormalizedAngleSignal->Integral(), "width");
hRecoNormalizedAngleBackground->Scale(norm/hRecoNormalizedAngleBackground->Integral(), "width");

std::cout << "Two tracksss" << std::endl;
std::cout << "True s: " << nTwoTrackTrueSignal << std::endl;
std::cout << "True b: " << nTwoTrackTrueBackground << std::endl;
std::cout << "  __/ regime " << std::endl;
std::cout << "  signal: " << nTwoTrackSignal << std::endl;
std::cout << "  background: " << nTwoTrackBackground << std::endl;
std::cout << "  __/__ regime " << std::endl;
std::cout << "  signal: " << nOtherTwoTrackSignal << std::endl;
std::cout << "  background: " << nOtherTwoTrackBackground << std::endl;

double two_track_selection_purity = 1.*(nTwoTrackSignal + nOtherTwoTrackSignal) / 
                                    (nTwoTrackSignal + nOtherTwoTrackSignal + nTwoTrackBackground + nOtherTwoTrackBackground);
double two_track_selection_effish = 1.*(nTwoTrackSignal + nOtherTwoTrackSignal) / (nTwoTrackTrueSignal);
std::cout << "two track selection purity: " << two_track_selection_purity << std::endl;
std::cout << "two track selection effish: " << two_track_selection_effish << std::endl;


double total_selection_purity = 1.*(nOneTrackSignalArray[6] + nTwoTrackSignal + nOtherTwoTrackSignal) 
                                /  (nOneTrackSignalArray[6] + nOneTrackBackgroundArray[7]
                                  + nTwoTrackSignal + nTwoTrackBackground
                                  + nOtherTwoTrackSignal + nOtherTwoTrackBackground );
double total_selection_effish = 1.*(nOneTrackSignalArray[6] + nTwoTrackSignal + nOtherTwoTrackSignal)
                                /  (nOneTrackTrueSignal + nTwoTrackTrueSignal);

std::cout << "TOTAL SELECTION PURITY: " << total_selection_purity << std::endl;
std::cout << "TOTAL SELECTION EFFISH: " << total_selection_effish << std::endl;

nChannel1Signal=nOneTrackSignalArray[6];
nChannel1Background=nOneTrackBackgroundArray[6];
nChannel3Signal=nTwoTrackSignal;
nChannel3Background=nTwoTrackBackground;
nChannel4Signal=nOtherTwoTrackSignal;
nChannel4Background=nOtherTwoTrackBackground;
std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;
std::cout << "double blahing area" << std::endl;
std::cout << "one track true signal (wholewise): " << nOneTrackTrueSignal << std::endl;
std::cout << "one track true signal (piecewise): " << (nChannel1True + nChannel2True) << std::endl;
std::cout << "total true signal... " << nTotalTrueSignal << std::endl;
std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;
std::cout << "### Channel ###   true     ###     obs      ###     signal      ###     background  ###" << std::endl;
std::cout << "      0           " << nChannel0True << "             "  << nChannel0 << "                " << nChannel0Signal << "               " << nChannel0Background << "   " << std::endl;

std::cout << "      1           " << nChannel1True << "             "  << nChannel1 << "                " << nChannel1Signal << "               " << nChannel1Background << "   " << std::endl;

std::cout << "      2           " << nChannel2True << "             "  << nChannel2 << "                " << nChannel2Signal << "               " << nChannel2Background << "   " << std::endl;

std::cout << "      3           " << nChannel3True << "             "  << nChannel3 << "                " << nChannel3Signal << "               " << nChannel3Background << "   " << std::endl;

std::cout << "      4           " << nChannel4True << "             "  << nChannel4 << "                " << nChannel4Signal << "               " << nChannel4Background << "   " << std::endl;


double ch1pur = 1.*(nChannel1Signal) / (nChannel1Signal + nChannel1Background);
double ch1eff = 1.*(nChannel1Signal) / (nChannel1True);
double ch2pur = 1.*(nChannel2Signal) / (nChannel2Signal + nChannel2Background);
double ch2eff = 1.*(nChannel2Signal) / (nChannel2True);
double ch3pur = 1.*(nChannel3Signal) / (nChannel3Signal + nChannel3Background);
double ch3eff = 1.*(nChannel3Signal) / (nChannel3True);
double ch4pur = 1.*(nChannel4Signal) / (nChannel4Signal + nChannel4Background);
double ch4eff = 1.*(nChannel4Signal) / (nChannel4True);

double chtotpur = 1.*(nChannel1Signal + nChannel3Signal + nChannel4Signal)
                    /(nChannel1Signal + nChannel3Signal + nChannel4Signal
                    + nChannel1Background + nChannel3Background + nChannel4Background);
double chtoteff = 1.*(nChannel1Signal + nChannel3Signal + nChannel4Signal)
                    /(nChannel1True   + nChannel3True   + nChannel4True);
double chtoteffall = 1.*(nChannel1Signal + nChannel3Signal + nChannel4Signal)
                       /(nChannel1True   + nChannel2TopTrue + nChannel3True   + nChannel4True);


std::cout << "ch1 p - e: " << ch1pur << " - " << ch1eff << std::endl;
//std::cout << "ch2 p - e: " << ch2pur << " - " << ch2eff << std::endl;
std::cout << "ch3 p - e: "    << ch3pur << " - " << ch3eff << std::endl;
std::cout << "ch4 p - e: "    << ch4pur << " - " << ch4eff << std::endl;
std::cout << "total p-e-e': " << chtotpur << " - " << chtoteff << " - " << chtoteffall << std::endl;


std::cout << std::endl;
std::cout << std::endl;
std::cout << "# Ch #  true  #  obs   #  sig  #  back  #" << std::endl;
std::cout << "  0     " << nChannel0InelastTrue << "      "  << nChannel0 << "      " << nChannel0InelastSignal << "      " << nChannel0InelastBackground << std::endl;
std::cout << "  1     " << nChannel1InelastTrue << "      "  << nChannel1 << "      " << nChannel1InelastSignal << "      " << nChannel1InelastBackground << std::endl;
std::cout << "  2     " << nChannel2InelastTrue << "      "  << nChannel2 << "      " << nChannel2InelastSignal << "      " << nChannel2InelastBackground << std::endl;
std::cout << "  3     " << nChannel3InelastTrue << "      "  << nChannel3 << "      " << nChannel3InelastSignal << "      " << nChannel3InelastBackground << std::endl;
std::cout << "  4     " << nChannel4InelastTrue << "      "  << nChannel4 << "      " << nChannel4InelastSignal << "      " << nChannel4InelastBackground << std::endl;
std::cout << std::endl;
std::cout << std::endl;

  // ### purity plots ###
  for(int jBin = 1; jBin < hRecoInteractingKESignal->GetNbinsX(); jBin++)
    {
    if(hRecoInteractingKESignal->GetBinContent(jBin) == 0 && hRecoInteractingKEBackground->GetBinContent(jBin) == 0){continue;}
    double s = 1.*hRecoInteractingKESignal->GetBinContent(jBin);
    double b = 1.*hRecoInteractingKEBackground->GetBinContent(jBin);
    double purity_temp = (s)/(s + b);
    double delta_purity_temp = (1. / (s + b)) * sqrt( s*(1 - s/(s + b)) );
    //std::cout << "purity: " << purity_temp << " +/- " << delta_purity_temp << std::endl;;
    hRecoInteractingKEPurity->SetBinContent(jBin, purity_temp);
    hRecoInteractingKEPurity->SetBinError(jBin, delta_purity_temp);
    }
  // ### elastic purity plots ###
  for(int jBin = 1; jBin < hRecoInteractingKESignalElastic->GetNbinsX(); jBin++)
    {
    if(hRecoInteractingKESignalElastic->GetBinContent(jBin) == 0 && hRecoInteractingKEBackgroundElastic->GetBinContent(jBin) == 0){continue;}
    double s = 1.*hRecoInteractingKESignalElastic->GetBinContent(jBin);
    double b = 1.*hRecoInteractingKEBackgroundElastic->GetBinContent(jBin);
    double purity_temp = (s)/(s + b);
    double delta_purity_temp = (1. / (s + b)) * sqrt( s*(1 - s/(s + b)) );
    //std::cout << "purity: " << purity_temp << " +/- " << delta_purity_temp << std::endl;;
    hRecoInteractingKEPurityElastic->SetBinContent(jBin, purity_temp);
    hRecoInteractingKEPurityElastic->SetBinError(jBin, delta_purity_temp);
    }
  // ### inelastic purity plots ###
  for(int jBin = 1; jBin < hRecoInteractingKESignalInelast->GetNbinsX(); jBin++)
    {
    if(hRecoInteractingKESignalInelast->GetBinContent(jBin) == 0 && hRecoInteractingKEBackgroundInelast->GetBinContent(jBin) == 0){continue;}
    double s = 1.*hRecoInteractingKESignalInelast->GetBinContent(jBin);
    double b = 1.*hRecoInteractingKEBackgroundInelast->GetBinContent(jBin);
    double purity_temp = (s)/(s + b);
    double delta_purity_temp = (1. / (s + b)) * sqrt( s*(1 - s/(s + b)) );
    //std::cout << "purity: " << purity_temp << " +/- " << delta_purity_temp << std::endl;;
    hRecoInteractingKEPurityInelast->SetBinContent(jBin, purity_temp);
    hRecoInteractingKEPurityInelast->SetBinError(jBin, delta_purity_temp);
    }

    // ## ch1 purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan1Signal->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan1Signal->GetBinContent(jBin) == 0 && hRecoInteractingKEChan1Background->GetBinContent(jBin) == 0){continue;}
      double ch1_s = 1.*hRecoInteractingKEChan1Signal->GetBinContent(jBin);
      double ch1_b = 1.*hRecoInteractingKEChan1Background->GetBinContent(jBin);
      double ch1_purity_temp = (ch1_s)/(ch1_s + ch1_b);
      double ch1_delta_purity_temp = (1. / (ch1_s + ch1_b)) * sqrt( ch1_s*(1 - ch1_s/(ch1_s + ch1_b)) );
      hRecoInteractingKEChan1Purity->SetBinContent(jBin, ch1_purity_temp);
      hRecoInteractingKEChan1Purity->SetBinError(jBin, ch1_delta_purity_temp);
      }
    // ## ch1 elastic purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan1SignalElastic->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan1SignalElastic->GetBinContent(jBin) == 0 && hRecoInteractingKEChan1BackgroundElastic->GetBinContent(jBin) == 0){continue;}
      double ch1_s = 1.*hRecoInteractingKEChan1SignalElastic->GetBinContent(jBin);
      double ch1_b = 1.*hRecoInteractingKEChan1BackgroundElastic->GetBinContent(jBin);
      double ch1_purity_temp = (ch1_s)/(ch1_s + ch1_b);
      double ch1_delta_purity_temp = (1. / (ch1_s + ch1_b)) * sqrt( ch1_s*(1 - ch1_s/(ch1_s + ch1_b)) );
      hRecoInteractingKEChan1PurityElastic->SetBinContent(jBin, ch1_purity_temp);
      hRecoInteractingKEChan1PurityElastic->SetBinError(jBin, ch1_delta_purity_temp);
      }
    // ## ch1 inelastic purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan1SignalInelast->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan1SignalInelast->GetBinContent(jBin) == 0 && hRecoInteractingKEChan1BackgroundInelast->GetBinContent(jBin) == 0){continue;}
      double ch1_s = 1.*hRecoInteractingKEChan1SignalInelast->GetBinContent(jBin);
      double ch1_b = 1.*hRecoInteractingKEChan1BackgroundInelast->GetBinContent(jBin);
      double ch1_purity_temp = (ch1_s)/(ch1_s + ch1_b);
      double ch1_delta_purity_temp = (1. / (ch1_s + ch1_b)) * sqrt( ch1_s*(1 - ch1_s/(ch1_s + ch1_b)) );
      hRecoInteractingKEChan1PurityInelast->SetBinContent(jBin, ch1_purity_temp);
      hRecoInteractingKEChan1PurityInelast->SetBinError(jBin, ch1_delta_purity_temp);
      }


    // ## ch2 purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan2Signal->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan2Signal->GetBinContent(jBin) == 0 && hRecoInteractingKEChan2Background->GetBinContent(jBin) == 0){continue;}
      double ch2_s = 1.*hRecoInteractingKEChan2Signal->GetBinContent(jBin);
      double ch2_b = 1.*hRecoInteractingKEChan2Background->GetBinContent(jBin);
      double ch2_purity_temp = (ch2_s)/(ch2_s + ch2_b);
      double ch2_delta_purity_temp = (1. / (ch2_s + ch2_b)) * sqrt( ch2_s*(1 - ch2_s/(ch2_s + ch2_b)) );
      hRecoInteractingKEChan2Purity->SetBinContent(jBin, ch2_purity_temp);
      hRecoInteractingKEChan2Purity->SetBinError(jBin, ch2_delta_purity_temp);
      }
    // ## ch2 elastic purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan2SignalElastic->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan2SignalElastic->GetBinContent(jBin) == 0 && hRecoInteractingKEChan2BackgroundElastic->GetBinContent(jBin) == 0){continue;}
      double ch2_s = 1.*hRecoInteractingKEChan2SignalElastic->GetBinContent(jBin);
      double ch2_b = 1.*hRecoInteractingKEChan2BackgroundElastic->GetBinContent(jBin);
      double ch2_purity_temp = (ch2_s)/(ch2_s + ch2_b);
      double ch2_delta_purity_temp = (1. / (ch2_s + ch2_b)) * sqrt( ch2_s*(1 - ch2_s/(ch2_s + ch2_b)) );
      hRecoInteractingKEChan2PurityElastic->SetBinContent(jBin, ch2_purity_temp);
      hRecoInteractingKEChan2PurityElastic->SetBinError(jBin, ch2_delta_purity_temp);
      }
    // ## ch1 inelastic purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan2SignalInelast->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan2SignalInelast->GetBinContent(jBin) == 0 && hRecoInteractingKEChan2BackgroundInelast->GetBinContent(jBin) == 0){continue;}
      double ch2_s = 1.*hRecoInteractingKEChan2SignalInelast->GetBinContent(jBin);
      double ch2_b = 1.*hRecoInteractingKEChan2BackgroundInelast->GetBinContent(jBin);
      double ch2_purity_temp = (ch2_s)/(ch2_s + ch2_b);
      double ch2_delta_purity_temp = (1. / (ch2_s + ch2_b)) * sqrt( ch2_s*(1 - ch2_s/(ch2_s + ch2_b)) );
      hRecoInteractingKEChan2PurityInelast->SetBinContent(jBin, ch2_purity_temp);
      hRecoInteractingKEChan2PurityInelast->SetBinError(jBin, ch2_delta_purity_temp);
      }


    // ## ch3 purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan3Signal->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan3Signal->GetBinContent(jBin) == 0 && hRecoInteractingKEChan3Background->GetBinContent(jBin) == 0){continue;}
      double ch3_s = 1.*hRecoInteractingKEChan3Signal->GetBinContent(jBin);
      double ch3_b = 1.*hRecoInteractingKEChan3Background->GetBinContent(jBin);
      double ch3_purity_temp = (ch3_s)/(ch3_s + ch3_b);
      double ch3_delta_purity_temp = (1. / (ch3_s + ch3_b)) * sqrt( ch3_s*(1 - ch3_s/(ch3_s + ch3_b)) );
      hRecoInteractingKEChan3Purity->SetBinContent(jBin, ch3_purity_temp);
      hRecoInteractingKEChan3Purity->SetBinError(jBin, ch3_delta_purity_temp);
      }
    // ## ch3 elastic purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan3SignalElastic->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan3SignalElastic->GetBinContent(jBin) == 0 && hRecoInteractingKEChan3BackgroundElastic->GetBinContent(jBin) == 0){continue;}
      double ch3_s = 1.*hRecoInteractingKEChan3SignalElastic->GetBinContent(jBin);
      double ch3_b = 1.*hRecoInteractingKEChan3BackgroundElastic->GetBinContent(jBin);
      double ch3_purity_temp = (ch3_s)/(ch3_s + ch3_b);
      double ch3_delta_purity_temp = (1. / (ch3_s + ch3_b)) * sqrt( ch3_s*(1 - ch3_s/(ch3_s + ch3_b)) );
      hRecoInteractingKEChan3PurityElastic->SetBinContent(jBin, ch3_purity_temp);
      hRecoInteractingKEChan3PurityElastic->SetBinError(jBin, ch3_delta_purity_temp);
      }
    // ## ch3 inelastic purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan3SignalInelast->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan3SignalInelast->GetBinContent(jBin) == 0 && hRecoInteractingKEChan3BackgroundInelast->GetBinContent(jBin) == 0){continue;}
      double ch3_s = 1.*hRecoInteractingKEChan3SignalInelast->GetBinContent(jBin);
      double ch3_b = 1.*hRecoInteractingKEChan3BackgroundInelast->GetBinContent(jBin);
      double ch3_purity_temp = (ch3_s)/(ch3_s + ch3_b);
      double ch3_delta_purity_temp = (1. / (ch3_s + ch3_b)) * sqrt( ch3_s*(1 - ch3_s/(ch3_s + ch3_b)) );
      hRecoInteractingKEChan3PurityInelast->SetBinContent(jBin, ch3_purity_temp);
      hRecoInteractingKEChan3PurityInelast->SetBinError(jBin, ch3_delta_purity_temp);
      }

    // ## ch4 purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan4Signal->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan4Signal->GetBinContent(jBin) == 0 && hRecoInteractingKEChan4Background->GetBinContent(jBin) == 0){continue;}
      double ch4_s = 1.*hRecoInteractingKEChan4Signal->GetBinContent(jBin);
      double ch4_b = 1.*hRecoInteractingKEChan4Background->GetBinContent(jBin);
      double ch4_purity_temp = (ch4_s)/(ch4_s + ch4_b);
      double ch4_delta_purity_temp = (1. / (ch4_s + ch4_b)) * sqrt( ch4_s*(1 - ch4_s/(ch4_s + ch4_b)) );
      hRecoInteractingKEChan4Purity->SetBinContent(jBin, ch4_purity_temp);
      hRecoInteractingKEChan4Purity->SetBinError(jBin, ch4_delta_purity_temp);
      }
    // ## ch4 elastic purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan4SignalElastic->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan4SignalElastic->GetBinContent(jBin) == 0 && hRecoInteractingKEChan4BackgroundElastic->GetBinContent(jBin) == 0){continue;}
      double ch4_s = 1.*hRecoInteractingKEChan4SignalElastic->GetBinContent(jBin);
      double ch4_b = 1.*hRecoInteractingKEChan4BackgroundElastic->GetBinContent(jBin);
      double ch4_purity_temp = (ch4_s)/(ch4_s + ch4_b);
      double ch4_delta_purity_temp = (1. / (ch4_s + ch4_b)) * sqrt( ch4_s*(1 - ch4_s/(ch4_s + ch4_b)) );
      hRecoInteractingKEChan4PurityElastic->SetBinContent(jBin, ch4_purity_temp);
      hRecoInteractingKEChan4PurityElastic->SetBinError(jBin, ch4_delta_purity_temp);
      }
    // ## ch4 inelastic purity ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan4SignalInelast->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKEChan4SignalInelast->GetBinContent(jBin) == 0 && hRecoInteractingKEChan4BackgroundInelast->GetBinContent(jBin) == 0){continue;}
      double ch4_s = 1.*hRecoInteractingKEChan4SignalInelast->GetBinContent(jBin);
      double ch4_b = 1.*hRecoInteractingKEChan3BackgroundInelast->GetBinContent(jBin);
      double ch4_purity_temp = (ch4_s)/(ch4_s + ch4_b);
      double ch4_delta_purity_temp = (1. / (ch4_s + ch4_b)) * sqrt( ch4_s*(1 - ch4_s/(ch4_s + ch4_b)) );
      hRecoInteractingKEChan4PurityInelast->SetBinContent(jBin, ch4_purity_temp);
      hRecoInteractingKEChan4PurityInelast->SetBinError(jBin, ch4_delta_purity_temp);
      }

    // ## truth level purity inelastic ##
    for(int jBin = 1; jBin < hMCRecoInteractingKESignalInelast->GetNbinsX(); jBin++)
      {
      if(hMCRecoInteractingKESignalInelast->GetBinContent(jBin) == 0 && hMCRecoInteractingKEBackgroundInelast->GetBinContent(jBin) == 0){continue;}
      double s = 1.*hMCRecoInteractingKESignalInelast->GetBinContent(jBin);
      double b = 1.*hMCRecoInteractingKEBackgroundInelast->GetBinContent(jBin);
      double purity_temp = (s)/(s + b);
      double delta_purity_temp = (1. / (s + b)) * sqrt( s*(1 - s/(s + b)) );
      //std::cout << "purity: " << purity_temp << " +/- " << delta_purity_temp << std::endl;;
      fInteractingKEPurityInelast->SetBinContent(jBin, purity_temp);
      fInteractingKEPurityInelast->SetBinError(jBin, delta_purity_temp);
      }




  
  // ### efficiency plots ###
  for(int jBin = 1; jBin < hRecoInteractingKESignal->GetNbinsX(); jBin++)
    {
    if(hRecoInteractingKETrueSignal->GetBinContent(jBin) == 0){continue;}
    double os = 1.*hRecoInteractingKESignal->GetBinContent(jBin);
    double ts = 1.*hRecoInteractingKETrueSignal->GetBinContent(jBin);
    double effish_temp = (os)/(ts);
    double delta_effish_temp = (1. / (ts)) * sqrt( os* abs(1 - os/ts) );
    hRecoInteractingKEEffish->SetBinContent(jBin, effish_temp);
    hRecoInteractingKEEffish->SetBinError(jBin, delta_effish_temp);
    }
    // ## inelastic effish ##
    for(int jBin = 1; jBin < hRecoInteractingKESignalInelast->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKETrueSignalInelast->GetBinContent(jBin) == 0){continue;}
      double osi = 1.*hRecoInteractingKESignalInelast->GetBinContent(jBin);
      double tsi = 1.*hRecoInteractingKETrueSignalInelast->GetBinContent(jBin);
      double effish_tempi = (osi)/(tsi);
      double delta_effish_tempi = (1. / (tsi)) * sqrt( osi* abs(1 - osi/tsi) );
      hRecoInteractingKEInelastEffish->SetBinContent(jBin, effish_tempi);
      hRecoInteractingKEInelastEffish->SetBinError(jBin, delta_effish_tempi);
      }
    // ## elastic effish    
    for(int jBin = 1; jBin < hRecoInteractingKESignalElastic->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKETrueSignalElastic->GetBinContent(jBin) == 0){continue;}
      double ose = 1.*hRecoInteractingKESignalElastic->GetBinContent(jBin);
      double tse = 1.*hRecoInteractingKETrueSignalElastic->GetBinContent(jBin);
      double effish_tempe = (ose)/(tse);
      double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
      hRecoInteractingKEElasticEffish->SetBinContent(jBin, effish_tempe);
      hRecoInteractingKEElasticEffish->SetBinError(jBin, delta_effish_tempe);
      }
    // ## elastic effish    
    for(int jBin = 1; jBin < hRecoInteractingKESignalElasticObs->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKETrueSignalElasticG6->GetBinContent(jBin) == 0){continue;}
      double ose = 1.*hRecoInteractingKESignalElasticObs->GetBinContent(jBin);
      double tse = 1.*hRecoInteractingKETrueSignalElasticG6->GetBinContent(jBin);
      double effish_tempe = (ose)/(tse);
      double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
      hRecoInteractingKEElasticObsEffish->SetBinContent(jBin, effish_tempe);
      hRecoInteractingKEElasticObsEffish->SetBinError(jBin, delta_effish_tempe);
      }
    // ## ch1 effish ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan1Signal->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKETrueSignal->GetBinContent(jBin) == 0){continue;}
      double ch1_os = 1.*hRecoInteractingKEChan1Signal->GetBinContent(jBin);
      double ch1_ts = 1.*hRecoInteractingKETrueSignal->GetBinContent(jBin);
      double ch1_effish_temp = (ch1_os)/(ch1_ts);
      double ch1_delta_effish_temp = (1. / (ch1_ts)) * sqrt( ch1_os* abs(1 - ch1_os/ch1_ts) );
      hRecoInteractingKEChan1Effish->SetBinContent(jBin, ch1_effish_temp);
      hRecoInteractingKEChan1Effish->SetBinError(jBin, ch1_delta_effish_temp);
      }
      // ## ch1 effish elastic ##
      for(int jBin = 1; jBin < hRecoInteractingKEChan1SignalElastic->GetNbinsX(); jBin++)
        {
        if(hRecoInteractingKETrueSignalElastic->GetBinContent(jBin) == 0){continue;}
        double ose = 1.*hRecoInteractingKEChan1SignalElastic->GetBinContent(jBin);
        double tse = 1.*hRecoInteractingKETrueSignalElastic->GetBinContent(jBin);
        double effish_tempe = (ose)/(tse);
        double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
        hRecoInteractingKEChan1ElasticEffish->SetBinContent(jBin, effish_tempe);
        hRecoInteractingKEChan1ElasticEffish->SetBinError(jBin, delta_effish_tempe);
        }
      // ## ch1 effish inelast ##
      for(int jBin = 1; jBin < hRecoInteractingKEChan1SignalInelast->GetNbinsX(); jBin++)
        {
        if(hRecoInteractingKETrueSignalInelast->GetBinContent(jBin) == 0){continue;}
        double ose = 1.*hRecoInteractingKEChan1SignalInelast->GetBinContent(jBin);
        double tse = 1.*hRecoInteractingKETrueSignalInelast->GetBinContent(jBin);
        double effish_tempe = (ose)/(tse);
        double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
        hRecoInteractingKEChan1InelastEffish->SetBinContent(jBin, effish_tempe);
        hRecoInteractingKEChan1InelastEffish->SetBinError(jBin, delta_effish_tempe);
        }
    // ## ch2 effish ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan2Signal->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKETrueSignal->GetBinContent(jBin) == 0){continue;}
      double ch2_os = 1.*hRecoInteractingKEChan2Signal->GetBinContent(jBin);
      double ch2_ts = 1.*hRecoInteractingKETrueSignal->GetBinContent(jBin);
      double ch2_effish_temp = (ch2_os)/(ch2_ts);
      double ch2_delta_effish_temp = (1. / (ch2_ts)) * sqrt( ch2_os* abs(1 - ch2_os/ch2_ts) );
      hRecoInteractingKEChan2Effish->SetBinContent(jBin, ch2_effish_temp);
      hRecoInteractingKEChan2Effish->SetBinError(jBin, ch2_delta_effish_temp);
      }
      // ## ch2 effish elastic ##
      for(int jBin = 1; jBin < hRecoInteractingKEChan2SignalElastic->GetNbinsX(); jBin++)
        {
        if(hRecoInteractingKETrueSignalElastic->GetBinContent(jBin) == 0){continue;}
        double ose = 1.*hRecoInteractingKEChan2SignalElastic->GetBinContent(jBin);
        double tse = 1.*hRecoInteractingKETrueSignalElastic->GetBinContent(jBin);
        double effish_tempe = (ose)/(tse);
        double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
        hRecoInteractingKEChan2ElasticEffish->SetBinContent(jBin, effish_tempe);
        hRecoInteractingKEChan2ElasticEffish->SetBinError(jBin, delta_effish_tempe);
        }
      // ## ch1 effish inelast ##
      for(int jBin = 1; jBin < hRecoInteractingKEChan2SignalInelast->GetNbinsX(); jBin++)
        {
        if(hRecoInteractingKETrueSignalInelast->GetBinContent(jBin) == 0){continue;}
        double ose = 1.*hRecoInteractingKEChan2SignalInelast->GetBinContent(jBin);
        double tse = 1.*hRecoInteractingKETrueSignalInelast->GetBinContent(jBin);
        double effish_tempe = (ose)/(tse);
        double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
        hRecoInteractingKEChan2InelastEffish->SetBinContent(jBin, effish_tempe);
        hRecoInteractingKEChan2InelastEffish->SetBinError(jBin, delta_effish_tempe);
        }
    // ## ch3 effish ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan3Signal->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKETrueSignal->GetBinContent(jBin) == 0){continue;}
      double ch3_os = 1.*hRecoInteractingKEChan3Signal->GetBinContent(jBin);
      double ch3_ts = 1.*hRecoInteractingKETrueSignal->GetBinContent(jBin);
      double ch3_effish_temp = (ch3_os)/(ch3_ts);
      double ch3_delta_effish_temp = (1. / (ch3_ts)) * sqrt( ch3_os* abs(1 - ch3_os/ch3_ts) );
      hRecoInteractingKEChan3Effish->SetBinContent(jBin, ch3_effish_temp);
      hRecoInteractingKEChan3Effish->SetBinError(jBin, ch3_delta_effish_temp);
      }
      // ## ch3 effish elastic ##
      for(int jBin = 1; jBin < hRecoInteractingKEChan3SignalElastic->GetNbinsX(); jBin++)
        {
        if(hRecoInteractingKETrueSignalElastic->GetBinContent(jBin) == 0){continue;}
        double ose = 1.*hRecoInteractingKEChan3SignalElastic->GetBinContent(jBin);
        double tse = 1.*hRecoInteractingKETrueSignalElastic->GetBinContent(jBin);
        double effish_tempe = (ose)/(tse);
        double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
        hRecoInteractingKEChan3ElasticEffish->SetBinContent(jBin, effish_tempe);
        hRecoInteractingKEChan3ElasticEffish->SetBinError(jBin, delta_effish_tempe);
        }
      // ## ch3 effish inelast ##
      for(int jBin = 1; jBin < hRecoInteractingKEChan3SignalInelast->GetNbinsX(); jBin++)
        {
        if(hRecoInteractingKETrueSignalInelast->GetBinContent(jBin) == 0){continue;}
        double ose = 1.*hRecoInteractingKEChan3SignalInelast->GetBinContent(jBin);
        double tse = 1.*hRecoInteractingKETrueSignalInelast->GetBinContent(jBin);
        double effish_tempe = (ose)/(tse);
        double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
        hRecoInteractingKEChan3InelastEffish->SetBinContent(jBin, effish_tempe);
        hRecoInteractingKEChan3InelastEffish->SetBinError(jBin, delta_effish_tempe);
        }
    // ## ch4 effish ##
    for(int jBin = 1; jBin < hRecoInteractingKEChan4Signal->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKETrueSignal->GetBinContent(jBin) == 0){continue;}
      double ch4_os = 1.*hRecoInteractingKEChan4Signal->GetBinContent(jBin);
      double ch4_ts = 1.*hRecoInteractingKETrueSignal->GetBinContent(jBin);
      double ch4_effish_temp = (ch4_os)/(ch4_ts);
      double ch4_delta_effish_temp = (1. / (ch4_ts)) * sqrt( ch4_os* abs(1 - ch4_os/ch4_ts) );
      hRecoInteractingKEChan4Effish->SetBinContent(jBin, ch4_effish_temp);
      hRecoInteractingKEChan4Effish->SetBinError(jBin, ch4_delta_effish_temp);
      }
      // ## ch4 effish elastic ##
      for(int jBin = 1; jBin < hRecoInteractingKEChan4SignalElastic->GetNbinsX(); jBin++)
        {
        if(hRecoInteractingKETrueSignalElastic->GetBinContent(jBin) == 0){continue;}
        double ose = 1.*hRecoInteractingKEChan4SignalElastic->GetBinContent(jBin);
        double tse = 1.*hRecoInteractingKETrueSignalElastic->GetBinContent(jBin);
        double effish_tempe = (ose)/(tse);
        double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
        hRecoInteractingKEChan4ElasticEffish->SetBinContent(jBin, effish_tempe);
        hRecoInteractingKEChan4ElasticEffish->SetBinError(jBin, delta_effish_tempe);
        }
      // ## ch4 effish inelast ##
      for(int jBin = 1; jBin < hRecoInteractingKEChan4SignalInelast->GetNbinsX(); jBin++)
        {
        if(hRecoInteractingKETrueSignalInelast->GetBinContent(jBin) == 0){continue;}
        double ose = 1.*hRecoInteractingKEChan4SignalInelast->GetBinContent(jBin);
        double tse = 1.*hRecoInteractingKETrueSignalInelast->GetBinContent(jBin);
        double effish_tempe = (ose)/(tse);
        double delta_effish_tempe = (1. / (tse)) * sqrt( ose* abs(1 - ose/tse) );
        hRecoInteractingKEChan4InelastEffish->SetBinContent(jBin, effish_tempe);
        hRecoInteractingKEChan4InelastEffish->SetBinError(jBin, delta_effish_tempe);
        }
      // ## ch1 effish elastic ##
      // ## ch1 effish inelast ##

    // ## true inelastic effish ##
    for(int jBin = 1; jBin < hMCRecoInteractingKESignalInelast->GetNbinsX(); jBin++)
      {
      if(hRecoInteractingKETrueSignalInelast->GetBinContent(jBin) == 0){continue;}
      double osi = 1.*hMCRecoInteractingKESignalInelast->GetBinContent(jBin);
      double tsi = 1.*hRecoInteractingKETrueSignalInelast->GetBinContent(jBin);
      double effish_tempi = (osi)/(tsi);
      double delta_effish_tempi = (1. / (tsi)) * sqrt( osi* abs(1 - osi/tsi) );
      fInteractingKEEffishInelast->SetBinContent(jBin, effish_tempi);
      fInteractingKEEffishInelast->SetBinError(jBin, delta_effish_tempi);
      }





// ### cross section plot making goes here ###
// ### incident and interacting kinetic energy plots need to go above this somewhere ###

// ### total ###
for( int iBin = 1; iBin <= hRecoInteractingKE->GetNbinsX(); ++iBin )
   {
   // ### If an incident bin is equal to zero then skip that bin ###
   if(hRecoIncidentKE->GetBinContent(iBin) == 0){continue;} //Temporary fix to ensure that no Infinities are propagated to pad
   if(hRecoInteractingKEEffish->GetBinContent(iBin) == 0){continue;}
   if(hRecoInteractingKEPurity->GetBinContent(iBin) == 0){continue;}

   // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
   float TempCrossSection = ( (hRecoInteractingKE->GetBinContent(iBin)*hRecoInteractingKEPurity->GetBinContent(iBin)*(1./hRecoInteractingKEEffish->GetBinContent(iBin)))
                            /(hRecoIncidentKE->GetBinContent(iBin)) ) * (1/number_density) * (1/slab_width);

   // ### Covert this into Barns ###
   float crossSection = TempCrossSection * (1/1e-28);

   // ### Putting the value on the plot
   fCrossSection->SetBinContent(iBin,crossSection);

   // ###########################################################
   // ### Calculating the error on the numerator of the ratio ###
   // ###########################################################
   float numError = pow(hRecoInteractingKE->GetBinContent(iBin),0.5);
   float num = hRecoInteractingKE->GetBinContent(iBin);

   // ### Putting in a protection against dividing by zero ###   
   if(num == 0){num = 1;}
   float term1 = numError/num;

   // #################################################
   // ### Calculating the error on the denomentator ###
   // #################################################
   float denomError = pow(hRecoIncidentKE->GetBinContent(iBin),0.5);
   float denom = hRecoIncidentKE->GetBinContent(iBin);
   if(denom == 0){denom = 1;}

   // ### Putting in a protection against dividing by zero ###
   float term2 = denomError/denom;

   float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);
   //std::cout<<"totalError = "<<totalError<<std::endl;
   fCrossSection->SetBinError(iBin,totalError);

   }//<---End iBin Loop
// ### elastic ###
for( int iBin = 1; iBin <= hRecoInteractingKE->GetNbinsX(); ++iBin )
   {
   // ### If an incident bin is equal to zero then skip that bin ###
   if(hRecoIncidentKE->GetBinContent(iBin) == 0){continue;} //Temporary fix to ensure that no Infinities are propagated to pad
   if(hRecoInteractingKEElasticEffish->GetBinContent(iBin) == 0){continue;}
   if(hRecoInteractingKEPurityElastic->GetBinContent(iBin) == 0){continue;}

   // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
   float TempCrossSection = ( (hRecoInteractingKE->GetBinContent(iBin)*hRecoInteractingKEPurityElastic->GetBinContent(iBin)
                            * (1./hRecoInteractingKEElasticEffish->GetBinContent(iBin)))
                            /(hRecoIncidentKE->GetBinContent(iBin)) ) * (1/number_density) * (1/slab_width);

   // ### Covert this into Barns ###
   float crossSection = TempCrossSection * (1/1e-28);

   // ### Putting the value on the plot
   fCrossSectionElastic->SetBinContent(iBin,crossSection);

   // ###########################################################
   // ### Calculating the error on the numerator of the ratio ###
   // ###########################################################
   float numError = pow(hRecoInteractingKE->GetBinContent(iBin),0.5);
   float num = hRecoInteractingKE->GetBinContent(iBin);

   // ### Putting in a protection against dividing by zero ###   
   if(num == 0){num = 1;}
   float term1 = numError/num;

   // #################################################
   // ### Calculating the error on the denomentator ###
   // #################################################
   float denomError = pow(hRecoIncidentKE->GetBinContent(iBin),0.5);
   float denom = hRecoIncidentKE->GetBinContent(iBin);
   if(denom == 0){denom = 1;}

   // ### Putting in a protection against dividing by zero ###
   float term2 = denomError/denom;

   float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);
   //std::cout<<"totalError = "<<totalError<<std::endl;
   fCrossSectionElastic->SetBinError(iBin,totalError);

   }//<---End iBin Loop
// ### inelastic ###
for( int iBin = 1; iBin <= hRecoInteractingKE->GetNbinsX(); ++iBin )
   {
   // ### If an incident bin is equal to zero then skip that bin ###
   if(hRecoIncidentKE->GetBinContent(iBin) == 0){continue;} //Temporary fix to ensure that no Infinities are propagated to pad
   if(hRecoInteractingKEInelastEffish->GetBinContent(iBin) == 0){continue;}
   if(hRecoInteractingKEPurityInelast->GetBinContent(iBin) == 0){continue;}

   // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
   float TempCrossSection = ( (hRecoInteractingKE->GetBinContent(iBin)*hRecoInteractingKEPurityInelast->GetBinContent(iBin)
                            * (1./hRecoInteractingKEInelastEffish->GetBinContent(iBin)))
                            /(hRecoIncidentKE->GetBinContent(iBin)) ) * (1/number_density) * (1/slab_width);

   float ratio = hRecoInteractingKE->GetBinContent(iBin) / (hRecoIncidentKE->GetBinContent(iBin) * number_density * slab_width * 1e-28);
   // ### Covert this into Barns ###
   float crossSection = TempCrossSection * (1/1e-28);

   // ### Putting the value on the plot
   fCrossSectionInelast->SetBinContent(iBin,crossSection);
   fRatio->SetBinContent(iBin, ratio);

   // ###########################################################
   // ### Calculating the error on the numerator of the ratio ###
   // ###########################################################
   float numError = pow(hRecoInteractingKE->GetBinContent(iBin),0.5);
   float num = hRecoInteractingKE->GetBinContent(iBin);

   // ### Putting in a protection against dividing by zero ###   
   if(num == 0){num = 1;}
   float term1 = numError/num;

   // #################################################
   // ### Calculating the error on the denomentator ###
   // #################################################
   float denomError = pow(hRecoIncidentKE->GetBinContent(iBin),0.5);
   float denom = hRecoIncidentKE->GetBinContent(iBin);
   if(denom == 0){denom = 1;}

   // ### Putting in a protection against dividing by zero ###
   float term2 = denomError/denom;

   float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);


   // ### playing with error ###
   //std::cout << "old calculation error: " << totalError <<std::endl;
   double n_int       = hRecoInteractingKE->GetBinContent(iBin);
   double delta_n_int = pow(hRecoInteractingKE->GetBinContent(iBin), .5);
   double eta         = hRecoInteractingKEPurityInelast->GetBinContent(iBin);
   double delta_eta   = hRecoInteractingKEPurityInelast->GetBinError(iBin);
   double n_inc       = hRecoIncidentKE->GetBinContent(iBin);
   double delta_n_inc = pow(hRecoIncidentKE->GetBinContent(iBin), .5);
   double eps         = hRecoInteractingKEInelastEffish->GetBinContent(iBin);
   double delta_eps   = hRecoInteractingKEInelastEffish->GetBinError(iBin);

   double term_a = pow( (delta_n_int) * ( (eta)         / (n_inc         * eps)         ), 2); 
   double term_b = pow( (delta_eta)   * ( (n_int)       / (n_inc         * eps)         ), 2); 
   double term_c = pow( (delta_n_inc) * ( (n_int * eta) / (pow(n_inc, 2) * eps)         ), 2);
   double term_d = pow( (delta_eps)   * ( (n_int * eta) / (n_inc         * pow(eps, 2)) ), 2);
   double quad_error = pow( pow( (delta_n_int) * ( (eta)         / (n_inc         * eps)         ), 2) 
                          + pow( (delta_eta)   * ( (n_int)       / (n_inc         * eps)         ), 2) 
                          + pow( (delta_n_inc) * ( (n_int * eta) / (pow(n_inc, 2) * eps)         ), 2)
                          + pow( (delta_eps)   * ( (n_int * eta) / (n_inc         * pow(eps, 2)) ), 2), .5);

   //std::cout << "       new error!: " << quad_error << std::endl;
   //std::cout << "               n_int = " << n_int << " +/- " << delta_n_int << std::endl;
   //std::cout << "               n_inc = " << n_inc << " +/- " << delta_n_inc << std::endl;
   //std::cout << "               eta   = " << eta   << " +/- " << delta_eta << std::endl;
   //std::cout << "               eps   = " << eps   << " +/- " << delta_eps << std::endl;
   //std::cout << "               terma = " << term_a << std::endl;
   //std::cout << "               termb = " << term_b << std::endl;
   //std::cout << "               termc = " << term_c << std::endl;
   //std::cout << "               termd = " << term_d << std::endl;

   fCrossSectionInelast->SetBinError(iBin,totalError);
   fRatio->SetBinError(iBin,totalError);

   }//<---End iBin Loop


  // ## smearing matrcies normalization ##
  for(int i = 0; i <= hMCRecoIncidentKE->GetNbinsX(); i++)
    {
    int column_total = 0;
    for(int j = 0; j < hMCRecoIncidentMatrix->GetNbinsY(); j++)
      {column_total += hMCRecoIncidentMatrix->GetBinContent(i,j);}
    if(column_total == 0){continue;}
    for(int j = 0; j < hMCRecoIncidentMatrixNormalized->GetNbinsY(); j++)
      {
      double normalized_value = hMCRecoIncidentMatrix->GetBinContent(i,j) / column_total;
      hMCRecoIncidentMatrixNormalized->SetBinContent(i,j,normalized_value); 
      }
    }
  for(int i = 0; i <= hMCRecoInteractingKE->GetNbinsX(); i++)
    {
    int column_total = 0;
    for(int j = 0; j < hMCRecoInteractingMatrix->GetNbinsY(); j++)
      {column_total += hMCRecoInteractingMatrix->GetBinContent(i,j);}
    if(column_total == 0){continue;}
    for(int j = 0; j < hMCRecoInteractingMatrixNormalized->GetNbinsY(); j++)
      {
      double normalized_value = hMCRecoInteractingMatrix->GetBinContent(i,j) / column_total;
      hMCRecoInteractingMatrixNormalized->SetBinContent(i,j,normalized_value); 
      }
    }
  // ## linear algebra ##
  for(int i = 0; i <= hMCRecoInteractingMatrixNormalized->GetNbinsX(); i++)
    {
    int reco_value = hRecoInteractingKE->GetBinContent(i);
    for(int j = 0; j < hMCRecoInteractingMatrixNormalized->GetNbinsY(); j++)
      {
      double weight = hMCRecoInteractingMatrixNormalized->GetBinContent(i,j);
      double weighted_value = reco_value * weight;
      fN_selected->AddBinContent(j, weighted_value);
      }
    }
  for(int i = 0; i <= hMCRecoIncidentMatrixNormalized->GetNbinsX(); i++)
    {
    int reco_value = hRecoIncidentKE->GetBinContent(i);
    for(int j = 0; j < hMCRecoIncidentMatrixNormalized->GetNbinsY(); j++)
      {
      double weight = hMCRecoIncidentMatrixNormalized->GetBinContent(i,j);
      double weighted_value = reco_value * weight;
      fN_incident->AddBinContent(j, weighted_value);
      }
    }
  // ## inelastic measurement with desmearing ##
  std::cout << "Inelastic Cross-Section Printouts..." << std::endl;
  for( int iBin = 1; iBin <= fN_selected->GetNbinsX(); ++iBin )
    {
    // ### If an incident bin is equal to zero then skip that bin ###
    if(fN_incident->GetBinContent(iBin) == 0){continue;} //Temporary fix to ensure that no Infinities are propagated to pad
    if(fInteractingKEEffishInelast->GetBinContent(iBin) == 0){continue;}
    if(fInteractingKEPurityInelast->GetBinContent(iBin) == 0){continue;}
  
    // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
    float n_sel = fN_selected->GetBinContent(iBin);
    float n_inc = fN_incident->GetBinContent(iBin);
    float pur   = fInteractingKEPurityInelast->GetBinContent(iBin); 
    float eff   = fInteractingKEEffishInelast->GetBinContent(iBin); 
    float TempCrossSection = ( n_sel*pur*(1./eff) )
                             /(n_inc) * (1/number_density) * (1/slab_width);
    float ratio = n_sel / (n_inc * number_density * slab_width * 1e-28);
    std::cout << "ratio corrected: " << ratio << std::endl;
    //float ratio_corrected
  
    // ### Covert this into Barns ###
    float crossSection = TempCrossSection * (1/1e-28);
    fInelasticCrossSection->SetBinContent(iBin, crossSection);
    fRatioCorrected->SetBinContent(iBin, ratio);
  
    // ### Calculating the error on the numerator of the ratio ###
    float numError = pow(n_sel,0.5);
    float num = n_sel;
    if(num == 0){num = 1;}
    float term1 = numError/num;
    float denomError = pow(n_inc, 0.5);
    float denom = n_inc;
    if(denom == 0){denom = 1;}
    float term2 = denomError/denom;
    float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);

    fInelasticCrossSection->SetBinError(iBin,totalError);
    fRatioCorrected->SetBinError(iBin, totalError);
    }//<---End iBin Loop



// ### resolution study ###
//for( i : pair )
for(auto const& tuple: matched_ke_tuple) 
    {
    double reco_int_ke = tuple[0]; 
    double true_int_ke = tuple[1]; 
    int topology_type  = tuple[2];
    if(true_int_ke < 50)
      {
      hRecoResolution0to50->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution0to50Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution0to50Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution0to50Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution0to50Chan4->Fill(reco_int_ke);}
      }
    if(true_int_ke < 100)
      {
      hRecoResolution50to100->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution50to100Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution50to100Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution50to100Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution50to100Chan4->Fill(reco_int_ke);}
      }
    if(true_int_ke < 150)
      {
      hRecoResolution100to150->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution100to150Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution100to150Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution100to150Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution100to150Chan4->Fill(reco_int_ke);}
      }
    if(true_int_ke < 200)
      {
      hRecoResolution150to200->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution150to200Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution150to200Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution150to200Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution150to200Chan4->Fill(reco_int_ke);}
      }
    else if(true_int_ke < 250)
      {
      hRecoResolution200to250->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution200to250Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution200to250Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution200to250Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution200to250Chan4->Fill(reco_int_ke);}
      }
    else if(true_int_ke < 300)
      {
      hRecoResolution250to300->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution250to300Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution250to300Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution250to300Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution250to300Chan4->Fill(reco_int_ke);}
      }
    else if(true_int_ke < 350)
      {
      hRecoResolution300to350->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution300to350Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution300to350Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution300to350Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution300to350Chan4->Fill(reco_int_ke);}
      }
    else if(true_int_ke < 400)
      {
      hRecoResolution350to400->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution350to400Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution350to400Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution350to400Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution350to400Chan4->Fill(reco_int_ke);}
      }
    else if(true_int_ke < 450)
      {
      hRecoResolution400to450->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution400to450Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution400to450Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution400to450Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution400to450Chan4->Fill(reco_int_ke);}
      }
    else if(true_int_ke < 500)
      {
      hRecoResolution450to500->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution450to500Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution450to500Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution450to500Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution450to500Chan4->Fill(reco_int_ke);}
      }
    else if(true_int_ke < 550)
      {
      hRecoResolution500to550->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution500to550Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution500to550Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution500to550Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution500to550Chan4->Fill(reco_int_ke);}
      }
    else if(true_int_ke < 600)
      {
      hRecoResolution550to600->Fill(reco_int_ke);
      if(topology_type == 1){hRecoResolution550to600Chan1->Fill(reco_int_ke);}
      if(topology_type == 2){hRecoResolution550to600Chan2->Fill(reco_int_ke);}
      if(topology_type == 3){hRecoResolution550to600Chan3->Fill(reco_int_ke);}
      if(topology_type == 4){hRecoResolution550to600Chan4->Fill(reco_int_ke);}
      }
    }

// ## resolution ###




for(auto const& pair: matched_ke) 
    {
    // ## diff ##
    hRecoResolutionDiff->Fill( pair.first - pair.second );
    if(pair.first > 200){hRecoResolutionDiff200->Fill( pair.first - pair.second);}
    }

//// ===========================================================================================
//// ============================  Writing out histograms to ROOT File =========================
//// ===========================================================================================
//
//// ###############################################
//// ### Creating a file to output my histograms ###
//// ###############################################
TFile myfile("mcblah_histos.root","RECREATE");

hMCEveryAngle->Write();
hMCEveryIntAngle->Write();
hMCEveryNonIntAngle->Write();
hMCInelasticAngle->Write();
hMCElasticAngle->Write();
hMCCoulombAngle->Write();
hMCNDaughters->Write();
hMCDaughterpdg->Write();
hMCNPDaughters->Write();
hMC0Pntracks->Write();
hMC1Pntracks->Write();
hMC2Pntracks->Write();
hMC3Pntracks->Write();
hMCnInteractions->Write();  
hMCOneIntNTracks->Write(); 
hMCOneIntOneTrackInteractionType->Write();
hMCGeant4StepSize->Write();
hMCGeant4dEdX->Write();
hRecodEdX->Write();
hRecodX->Write();
hMCEnergyLossIntoTPC->Write();
hMCEnergyLossIntoTPCRatio->Write();
hMCInitialInteractingKE->Write();
hMCElasticInitialInteractingKE->Write();
hMCInelastInitialInteractingKE->Write();
//// ### all angles ###
hRecoSTEEveryAngle->Write();
hRecoSTEEveryAngleInts->Write();
hRecoSTEEveryAngleIntsInelast->Write();
hRecoSTEEveryAngleIntsCoulomb->Write();
hRecoSTEEveryAngleIntsElastic->Write();
hRecoSTEEveryAngleNoInts->Write();
hRecoDTETrackEndDelta->Write();
hRecoDTETrackAngle->Write();
hRecoDTETrackAngleInts->Write();
hRecoDTETrackAngleIntsInelast->Write();
hRecoDTETrackAngleIntsCoulomb->Write();
hRecoDTETrackAngleIntsElastic->Write();
hRecoDTETrackAngleNoInts->Write();
hRecoDTEDistanceToMCIntFromTrack2->Write();
hRecoTTE23vsg4->Write();
hXSIntAnglePurity->Write();
hXSIntAngleEffish->Write();
hMCNormalizedAngleSignal->Write();
hMCNormalizedAngleBackground->Write();
hRecoNormalizedAngleSignal->Write();
hRecoNormalizedAngleBackground->Write();
  // ### signal ###
  hMCNormalizedAngleSignal->SetFillColor(kYellow+2);
  hMCNormalizedAngleSignal->SetLineColor(kYellow+2);
  hMCNormalizedAngleSignal->SetFillStyle(3006);
  hMCNormalizedAngleStacked->Add(hMCNormalizedAngleSignal);
  hRecoNormalizedAngleSignal->SetFillColor(kYellow+2);
  hRecoNormalizedAngleSignal->SetLineColor(kYellow+2);
  hRecoNormalizedAngleSignal->SetFillStyle(3006);
  hRecoNormalizedAngleStacked->Add(hRecoNormalizedAngleSignal);
  // ### no ints ###
  hMCNormalizedAngleBackground->SetFillColor(kGreen+2);
  hMCNormalizedAngleBackground->SetLineColor(kGreen+2);
  hMCNormalizedAngleBackground->SetFillStyle(3007);
  hMCNormalizedAngleStacked->Add(hMCNormalizedAngleBackground);
  hRecoNormalizedAngleBackground->SetFillColor(kGreen+2);
  hRecoNormalizedAngleBackground->SetLineColor(kGreen+2);
  hRecoNormalizedAngleBackground->SetFillStyle(3007);
  hRecoNormalizedAngleStacked->Add(hRecoNormalizedAngleBackground);
  // ### whole thing ###
hMCNormalizedAngleStacked->Write();
hRecoNormalizedAngleStacked->Write();
hRecoMatchedInteractingKERatio->Write();
hRecoMatchedInteractingKECmpsn->Write();
hRecoMatchedInteractingKECmpsnInelast->Write();
hRecoMatchedElasticAngleCmpsn->Write();
hRecoDeltaEnergy->Write();
hRecoDeltaEnergyChannel0->Write();
hRecoTrueVSRecoEnergyChannel0->Write();
hRecoDeltaEnergyvsIntKE->Write();
hRecoDeltaEnergyElastic->Write();
hRecoDeltaEnergyInelast->Write();
hRecoDeltaEnergyInelastvsIntKE->Write();
    // ### elastic ###
    hRecoDeltaEnergyElastic->SetFillColor(kYellow+2);
    hRecoDeltaEnergyElastic->SetLineColor(kYellow+2);
    hRecoDeltaEnergyElastic->SetFillStyle(3006);
    hRecoDeltaEnergyStacked->Add(hRecoDeltaEnergyElastic);
    // ### elastic ###
    hRecoDeltaEnergyInelast->SetFillColor(kGreen+2);
    hRecoDeltaEnergyInelast->SetLineColor(kGreen+2);
    hRecoDeltaEnergyInelast->SetFillStyle(3007);
    hRecoDeltaEnergyStacked->Add(hRecoDeltaEnergyInelast);
hRecoDeltaEnergyStacked->Write();
hRecoLowMeanDeltaEnergy->Write();
hRecoLowMeanDeltaEnergySignal->Write();
hRecoLowMeanDeltaEnergyBackground->Write();
hRecoEndMeandedx->Write();
hRecoEndMeandedxElastic->Write();
hRecoEndMeandedxInelast->Write();
hRecoEndMeandedxNeither->Write();

hRecoIncidentKE->Write();
hRecoInteractingKE->Write();
    hRecoInteractingKEChan1->Write();
    hRecoInteractingKEChan1TrueSignal->Write();
    hRecoInteractingKEChan1TrueSignalG6->Write();
    hRecoInteractingKEChan1Signal->Write();
    hRecoInteractingKEChan1SignalElastic->Write();
    hRecoInteractingKEChan1SignalInelast->Write();
    hRecoInteractingKEChan1Background->Write();
    hRecoInteractingKEChan1BackgroundElastic->Write();
    hRecoInteractingKEChan1BackgroundInelast->Write();
    hRecoInteractingKEChan2->Write();
    hRecoInteractingKEChan2Signal->Write();
    hRecoInteractingKEChan2SignalElastic->Write();
    hRecoInteractingKEChan2SignalInelast->Write();
    hRecoInteractingKEChan2Background->Write();
    hRecoInteractingKEChan2BackgroundElastic->Write();
    hRecoInteractingKEChan2BackgroundInelast->Write();
    hRecoInteractingKEChan3->Write();
    hRecoInteractingKEChan3TrueSignal->Write();
    hRecoInteractingKEChan3TrueSignalG6->Write();
    hRecoInteractingKEChan3Signal->Write();
    hRecoInteractingKEChan3SignalElastic->Write();
    hRecoInteractingKEChan3SignalInelast->Write();
    hRecoInteractingKEChan3Background->Write();
    hRecoInteractingKEChan3BackgroundElastic->Write();
    hRecoInteractingKEChan3BackgroundInelast->Write();
    hRecoInteractingKEChan4->Write();
    hRecoInteractingKEChan4TrueSignal->Write();
    hRecoInteractingKEChan4TrueSignalG6->Write();
    hRecoInteractingKEChan4Signal->Write();
    hRecoInteractingKEChan4SignalElastic->Write();
    hRecoInteractingKEChan4SignalInelast->Write();
    hRecoInteractingKEChan4Background->Write();
    hRecoInteractingKEChan4BackgroundElastic->Write();
    hRecoInteractingKEChan4BackgroundInelast->Write();
    hRecoInteractingKETrueSignal->Write();
    hRecoInteractingKETrueSignalInelast->Write();
    hRecoInteractingKETrueSignalInelastG6->Write();
    hRecoInteractingKETrueSignalElastic->Write();
    hRecoInteractingKETrueSignalElasticG6->Write();
    hRecoInteractingKETrueSignalEA->Write();
    hRecoInteractingKETrueSignalIA->Write();
    hRecoInteractingKETrueSignalG6->Write();
    hRecoInteractingKESignal->Write();
    hRecoInteractingKESignalElastic->Write();
    hRecoInteractingKESignalElasticObs->Write();
    hRecoInteractingKESignalInelast->Write();
    hRecoInteractingKEBackground->Write();
    hRecoInteractingKEBackgroundElastic->Write();
    hRecoInteractingKEBackgroundInelast->Write();
        // ## ch1 ##
        hRecoInteractingKEChan1->SetFillColor(kGreen+2);
        hRecoInteractingKEChan1->SetLineColor(kGreen+2);
        hRecoInteractingKEChan1->SetFillStyle(3007);
        hRecoInteractingKEStacked->Add(hRecoInteractingKEChan1);
        // ## ch2 ##
        hRecoInteractingKEChan2->SetFillColor(kBlue+2);
        hRecoInteractingKEChan2->SetLineColor(kBlue+2);
        hRecoInteractingKEChan2->SetFillStyle(3007);
        hRecoInteractingKEStacked->Add(hRecoInteractingKEChan2);
        // ## ch3 ##
        hRecoInteractingKEChan3->SetFillColor(kYellow+2);
        hRecoInteractingKEChan3->SetLineColor(kYellow+2);
        hRecoInteractingKEChan3->SetFillStyle(3006);
        hRecoInteractingKEStacked->Add(hRecoInteractingKEChan3);
        // ## ch4 ##
        hRecoInteractingKEChan4->SetFillColor(kRed+2);
        hRecoInteractingKEChan4->SetLineColor(kRed+2);
        hRecoInteractingKEChan4->SetFillStyle(3004);
        hRecoInteractingKEStacked->Add(hRecoInteractingKEChan4);
        // ### dumb test ###
        hRecoInteractingKEStacked2->Add(hRecoInteractingKEChan4);
        hRecoInteractingKEStacked2->Add(hRecoInteractingKEChan3);
        hRecoInteractingKEStacked2->Add(hRecoInteractingKEChan2);
        hRecoInteractingKEStacked2->Add(hRecoInteractingKEChan1);
        // ## signal ##
        hRecoInteractingKESignal->SetFillColor(kRed+2);
        hRecoInteractingKESignal->SetLineColor(kRed+3);
        hRecoInteractingKESignal->SetFillStyle(3007);
            hRecoInteractingKEChan1Signal->SetFillColor(kRed+2);
            hRecoInteractingKEChan1Signal->SetLineColor(kRed+3);
            hRecoInteractingKEChan1Signal->SetFillStyle(3007);
            hRecoInteractingKEChan3Signal->SetFillColor(kRed+2);
            hRecoInteractingKEChan3Signal->SetLineColor(kRed+3);
            hRecoInteractingKEChan3Signal->SetFillStyle(3007);
            hRecoInteractingKEChan4Signal->SetFillColor(kRed+2);
            hRecoInteractingKEChan4Signal->SetLineColor(kRed+3);
            hRecoInteractingKEChan4Signal->SetFillStyle(3007);
        // ## total ## 
        hRecoInteractingKEStackedP->Add(hRecoInteractingKESignal);
        hRecoInteractingKEStackedE->Add(hRecoInteractingKESignal);
        // ## by channels ##
          hRecoInteractingKEChan1StackedP->Add(hRecoInteractingKEChan1Signal);
          hRecoInteractingKEChan1StackedE->Add(hRecoInteractingKEChan1Signal);

          hRecoInteractingKEChan3StackedP->Add(hRecoInteractingKEChan3Signal);
          hRecoInteractingKEChan3StackedE->Add(hRecoInteractingKEChan3Signal);

          hRecoInteractingKEChan4StackedP->Add(hRecoInteractingKEChan4Signal);
          hRecoInteractingKEChan4StackedE->Add(hRecoInteractingKEChan4Signal);
        
        // ## true signal ##
        hRecoInteractingKETrueSignal->SetFillColor(kGreen+3);
        hRecoInteractingKETrueSignal->SetLineColor(kGreen+4);
        hRecoInteractingKETrueSignal->SetFillStyle(3006);
            hRecoInteractingKEChan1TrueSignal->SetFillColor(kGreen+2);
            hRecoInteractingKEChan1TrueSignal->SetLineColor(kGreen+3);
            hRecoInteractingKEChan1TrueSignal->SetFillStyle(3006);
            hRecoInteractingKEChan3TrueSignal->SetFillColor(kGreen+2);
            hRecoInteractingKEChan3TrueSignal->SetLineColor(kGreen+3);
            hRecoInteractingKEChan3TrueSignal->SetFillStyle(3006);
            hRecoInteractingKEChan4TrueSignal->SetFillColor(kGreen+2);
            hRecoInteractingKEChan4TrueSignal->SetLineColor(kGreen+3);
            hRecoInteractingKEChan4TrueSignal->SetFillStyle(3006);
        // ## total ##
        hRecoInteractingKEStackedE->Add(hRecoInteractingKETrueSignal);
        // ## by channels ##
        hRecoInteractingKEChan1StackedE->Add(hRecoInteractingKEChan1TrueSignal);
        hRecoInteractingKEChan3StackedE->Add(hRecoInteractingKEChan3TrueSignal);
        hRecoInteractingKEChan4StackedE->Add(hRecoInteractingKEChan4TrueSignal);


        // ## background ##
        hRecoInteractingKEBackground->SetFillColor(kCyan-1);
        hRecoInteractingKEBackground->SetLineColor(kCyan+3);
        hRecoInteractingKEBackground->SetFillStyle(3006);
            hRecoInteractingKEChan1Background->SetFillColor(kCyan-1);
            hRecoInteractingKEChan1Background->SetLineColor(kCyan+3);
            hRecoInteractingKEChan1Background->SetFillStyle(3006);
            hRecoInteractingKEChan3Background->SetFillColor(kCyan-1);
            hRecoInteractingKEChan3Background->SetLineColor(kCyan+3);
            hRecoInteractingKEChan3Background->SetFillStyle(3006);
            hRecoInteractingKEChan4Background->SetFillColor(kCyan-1);
            hRecoInteractingKEChan4Background->SetLineColor(kCyan+3);
            hRecoInteractingKEChan4Background->SetFillStyle(3006);
        // ## total ##
        hRecoInteractingKEStackedP->Add(hRecoInteractingKEBackground);
        // ## by channels ##
          hRecoInteractingKEChan1StackedP->Add(hRecoInteractingKEChan1Background);
          hRecoInteractingKEChan3StackedP->Add(hRecoInteractingKEChan3Background);
          hRecoInteractingKEChan4StackedP->Add(hRecoInteractingKEChan4Background);

    // ## stacked ##
    //hRecoInteractingKEStacked->GetXaxis()->SetTitle("Reconstructed KE (MeV)");
    //gPad->BuildLegend(0.75, 0.75, 0.95, 0.95, "");
    hRecoInteractingKEStacked->Write();
    hRecoInteractingKEStacked2->Write();
    hRecoInteractingKEStackedP->Write();
    hRecoInteractingKEStackedE->Write();
    hRecoInteractingKEChan1StackedP->Write();
    hRecoInteractingKEChan1StackedE->Write();
    hRecoInteractingKEChan3StackedP->Write();
    hRecoInteractingKEChan3StackedE->Write();
    hRecoInteractingKEChan4StackedP->Write();
    hRecoInteractingKEChan4StackedE->Write();
        
hRecoInteractingKEPurity->Write();
hRecoInteractingKEPurityElastic->Write();
hRecoInteractingKEPurityInelast->Write();
hRecoInteractingKEChan1Purity->Write();
hRecoInteractingKEChan1PurityElastic->Write();
hRecoInteractingKEChan1PurityInelast->Write();
hRecoInteractingKEChan2Purity->Write();
hRecoInteractingKEChan2PurityElastic->Write();
hRecoInteractingKEChan2PurityInelast->Write();
hRecoInteractingKEChan3Purity->Write();
hRecoInteractingKEChan3PurityElastic->Write();
hRecoInteractingKEChan3PurityInelast->Write();
hRecoInteractingKEChan4Purity->Write();
hRecoInteractingKEChan4PurityElastic->Write();
hRecoInteractingKEChan4PurityInelast->Write();
hRecoInteractingKEEffish->Write();
hRecoInteractingKEElasticEffish->Write();
hRecoInteractingKEElasticObsEffish->Write();
hRecoInteractingKEInelastEffish->Write();
hRecoInteractingKEChan1Effish->Write();
hRecoInteractingKEChan1ElasticEffish->Write();
hRecoInteractingKEChan1InelastEffish->Write();
hRecoInteractingKEChan2Effish->Write();
hRecoInteractingKEChan2ElasticEffish->Write();
hRecoInteractingKEChan2InelastEffish->Write();
hRecoInteractingKEChan3Effish->Write();
hRecoInteractingKEChan3ElasticEffish->Write();
hRecoInteractingKEChan3InelastEffish->Write();
hRecoInteractingKEChan4Effish->Write();
hRecoInteractingKEChan4ElasticEffish->Write();
hRecoInteractingKEChan4InelastEffish->Write();
fCrossSection->Write();
fCrossSectionElastic->Write();
fCrossSectionInelast->Write();
hRecoResolutionDiff->Write();
hRecoResolutionDiff200->Write();
hRecoResolution0to50->Write();
hRecoResolution0to50Chan1->Write();
hRecoResolution0to50Chan2->Write();
hRecoResolution0to50Chan3->Write();
hRecoResolution0to50Chan4->Write();
hRecoResolution50to100->Write();
hRecoResolution50to100Chan1->Write();
hRecoResolution50to100Chan2->Write();
hRecoResolution50to100Chan3->Write();
hRecoResolution50to100Chan4->Write();
hRecoResolution100to150->Write();
hRecoResolution100to150Chan1->Write();
hRecoResolution100to150Chan2->Write();
hRecoResolution100to150Chan3->Write();
hRecoResolution100to150Chan4->Write();
hRecoResolution150to200->Write();
hRecoResolution150to200Chan1->Write();
hRecoResolution150to200Chan2->Write();
hRecoResolution150to200Chan3->Write();
hRecoResolution150to200Chan4->Write();
hRecoResolution200to250->Write();
hRecoResolution200to250Chan1->Write();
hRecoResolution200to250Chan2->Write();
hRecoResolution200to250Chan3->Write();
hRecoResolution200to250Chan4->Write();
hRecoResolution250to300->Write();
hRecoResolution250to300Chan1->Write();
hRecoResolution250to300Chan2->Write();
hRecoResolution250to300Chan3->Write();
hRecoResolution250to300Chan4->Write();
hRecoResolution300to350->Write();
hRecoResolution300to350Chan1->Write();
hRecoResolution300to350Chan2->Write();
hRecoResolution300to350Chan3->Write();
hRecoResolution300to350Chan4->Write();
hRecoResolution350to400->Write();
hRecoResolution350to400Chan1->Write();
hRecoResolution350to400Chan2->Write();
hRecoResolution350to400Chan3->Write();
hRecoResolution350to400Chan4->Write();
hRecoResolution400to450->Write();
hRecoResolution400to450Chan1->Write();
hRecoResolution400to450Chan2->Write();
hRecoResolution400to450Chan3->Write();
hRecoResolution400to450Chan4->Write();
hRecoResolution450to500->Write();
hRecoResolution450to500Chan1->Write();
hRecoResolution450to500Chan2->Write();
hRecoResolution450to500Chan3->Write();
hRecoResolution450to500Chan4->Write();
hRecoResolution500to550->Write();
hRecoResolution500to550Chan1->Write();
hRecoResolution500to550Chan2->Write();
hRecoResolution500to550Chan3->Write();
hRecoResolution500to550Chan4->Write();
hRecoResolution550to600->Write();
hRecoResolution550to600Chan1->Write();
hRecoResolution550to600Chan2->Write();
hRecoResolution550to600Chan3->Write();
hRecoResolution550to600Chan4->Write();
hMCRecoIncidentKE->Write();
hMCRecoIncidentMatrix->Write();
hMCRecoIncidentMatrixNormalized->Write();
hMCRecoInteractingKE->Write();
hMCRecoInteractingMatrix->Write();
hMCRecoInteractingMatrixNormalized->Write();
hMCRecoInteractingKESignalInelast->Write();
hMCRecoInteractingKEBackgroundInelast->Write();
    // ## some stacked distributions ##
    THStack *hMCRecoStackedPInelast  = new THStack("hMCRecoStackedPInelast", "Purity");
    THStack *hMCRecoStackedEInelast  = new THStack("hMCRecoStackedEInelast", "Efficiency");
        // # sig #
        hMCRecoInteractingKESignalInelast->SetFillColor(kRed+2);
        hMCRecoInteractingKESignalInelast->SetLineColor(kRed+3);
        hMCRecoInteractingKESignalInelast->SetFillStyle(3007);
        hMCRecoStackedPInelast->Add(hMCRecoInteractingKESignalInelast);
        hMCRecoStackedEInelast->Add(hMCRecoInteractingKESignalInelast);
        // # bac #
        hMCRecoInteractingKEBackgroundInelast->SetFillColor(kCyan-1);
        hMCRecoInteractingKEBackgroundInelast->SetLineColor(kCyan+3);
        hMCRecoInteractingKEBackgroundInelast->SetFillStyle(3006);
        hMCRecoStackedPInelast->Add(hMCRecoInteractingKEBackgroundInelast);
        // # tru #
        hRecoInteractingKETrueSignalInelast->SetFillColor(kGreen+3);
        hRecoInteractingKETrueSignalInelast->SetLineColor(kGreen+4);
        hRecoInteractingKETrueSignalInelast->SetFillStyle(3006);
        hMCRecoStackedEInelast->Add(hRecoInteractingKETrueSignalInelast);
    hMCRecoStackedPInelast->Write();
    hMCRecoStackedEInelast->Write();
    
    // ## sanity check ##
    TH1D *hSanity1 = (TH1D*)fN_selected->Clone("hSanity1");
    hSanity1->Divide(hMCRecoInteractingKE);
    hSanity1->Write();
fN_selected->Write();
fN_incident->Write();
fInteractingKEPurityInelast->Write();
fInteractingKEEffishInelast->Write();
fInelasticCrossSection->Write();
fRatio->Write();
fRatioCorrected->Write();

}//<---End Loop() Function



