//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 22 12:12:17 2018 by ROOT version 6.12/04
// from TTree anatree/analysis tree
// found on file: ../files/positive_polarity_preproduction_test_reco_ntuple_20180821_1531.root
//////////////////////////////////////////////////////////

#ifndef PpolSummary_h
#define PpolSummary_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class PpolSummary {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Double_t        evttime;
   Double_t        efield[3];
   Int_t           t0;
   Int_t           ntracks_reco;
   vector<int>     *track_WC2TPC_match;
   vector<int>     *track_primary;
   vector<double>  *track_start_x;
   vector<double>  *track_start_y;
   vector<double>  *track_start_z;
   vector<double>  *track_end_x;
   vector<double>  *track_end_y;
   vector<double>  *track_end_z;
   vector<double>  *track_length;
   vector<int>     *ntrack_hits;
   vector<vector<double> > *track_xpos;
   vector<vector<double> > *track_ypos;
   vector<vector<double> > *track_zpos;
   vector<vector<int> > *nhit_ids;
   Double_t        track_spt_idarr[145][1000][100];   //[ntracks_reco]
   Double_t        track_spt_earr[145][1000][100];   //[ntracks_reco]
   vector<int>     *ind_track_hits;
   vector<double>  *ind_track_ke;
   vector<vector<double> > *ind_track_wire;
   vector<vector<double> > *ind_track_dedx;
   vector<vector<double> > *ind_track_dqdx;
   vector<vector<double> > *ind_track_rr;
   vector<vector<double> > *ind_track_pitch_hit;
   vector<int>     *col_track_hits;
   vector<double>  *col_track_ke;
   vector<vector<double> > *col_track_x;
   vector<vector<double> > *col_track_y;
   vector<vector<double> > *col_track_z;
   vector<vector<double> > *col_track_wire;
   vector<vector<double> > *col_track_dedx;
   vector<vector<double> > *col_track_dqdx;
   vector<vector<double> > *col_track_rr;
   vector<vector<double> > *col_track_pitch_hit;
   Int_t           nhits;
   vector<double>  *hit_time;
   vector<double>  *hit_amp;
   vector<double>  *hit_wire;
   vector<double>  *hit_view;
   Int_t           no_primaries;
   Int_t           geant_list_size;
   Double_t        primary_p;
   vector<int>     *PDG;
   vector<double>  *StartKE;
   vector<double>  *LastKE;
   vector<double>  *StartEnergy;
   vector<double>  *StartPx;
   vector<double>  *StartPy;
   vector<double>  *StartPz;
   vector<double>  *EndEnergy;
   vector<double>  *EndPx;
   vector<double>  *EndPy;
   vector<double>  *EndPz;
   vector<double>  *StartPointx;
   vector<double>  *StartPointy;
   vector<double>  *StartPointz;
   vector<double>  *EndPointx;
   vector<double>  *EndPointy;
   vector<double>  *EndPointz;
   vector<int>     *Process;
   vector<int>     *NumberDaughters;
   vector<int>     *Mother;
   vector<int>     *TrackId;
   vector<int>     *process_primary;
   vector<string>  *G4Process;
   vector<string>  *G4FinalProcess;
   vector<int>     *NTrTrajPts;
   Int_t           NProtonDaughters;
   Int_t           NNeutronDaughters;
   vector<double>  *NDTrTrajPts;
   vector<int>     *DTrackId;
   vector<int>     *DPdgCode;
   vector<double>  *DStartEnergy;
   vector<double>  *DStartP;
   vector<vector<double> > *MidPosX;
   vector<vector<double> > *MidPosY;
   vector<vector<double> > *MidPosZ;
   vector<vector<double> > *MidPx;
   vector<vector<double> > *MidPy;
   vector<vector<double> > *MidPz;
   vector<double>  *SlabX;
   vector<double>  *SlabY;
   vector<double>  *SlabZ;
   vector<double>  *SlapX;
   vector<double>  *SlapY;
   vector<double>  *SlapZ;
   vector<int>     *SlabN;
   vector<double>  *SlabE;
   Int_t           LastSlabInt;
   vector<vector<double> > *DMidPosX;
   vector<vector<double> > *DMidPosY;
   vector<vector<double> > *DMidPosZ;
   vector<int>     *InteractionPoint;
   vector<int>     *InteractionPointType;
   Int_t           num_tof_objects;
   Double_t        tofObject[1];   //[num_tof_objects]
   Int_t           num_wctracks;
   Double_t        wctrk_momentum[1];   //[num_wctracks]
   Double_t        wctrk_XFace[1];   //[num_wctracks]
   Double_t        wctrk_YFace[1];   //[num_wctracks]
   Double_t        wctrk_theta[1];   //[num_wctracks]
   Double_t        wctrk_phi[1];   //[num_wctracks]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_efield;   //!
   TBranch        *b_t0;   //!
   TBranch        *b_ntracks_reco;   //!
   TBranch        *b_track_WC2TPC_match;   //!
   TBranch        *b_track_primary;   //!
   TBranch        *b_track_start_x;   //!
   TBranch        *b_track_start_y;   //!
   TBranch        *b_track_start_z;   //!
   TBranch        *b_track_end_x;   //!
   TBranch        *b_track_end_y;   //!
   TBranch        *b_track_end_z;   //!
   TBranch        *b_track_length;   //!
   TBranch        *b_ntrack_hits;   //!
   TBranch        *b_track_xpos;   //!
   TBranch        *b_track_ypos;   //!
   TBranch        *b_track_zpos;   //!
   TBranch        *b_nhit_ids;   //!
   TBranch        *b_track_spt_idarr;   //!
   TBranch        *b_track_spt_earr;   //!
   TBranch        *b_ind_track_hits;   //!
   TBranch        *b_ind_track_ke;   //!
   TBranch        *b_ind_track_wire;   //!
   TBranch        *b_ind_track_dedx;   //!
   TBranch        *b_ind_track_dqdx;   //!
   TBranch        *b_ind_track_rr;   //!
   TBranch        *b_ind_track_pitch_hit;   //!
   TBranch        *b_col_track_hits;   //!
   TBranch        *b_col_track_ke;   //!
   TBranch        *b_col_track_x;   //!
   TBranch        *b_col_track_y;   //!
   TBranch        *b_col_track_z;   //!
   TBranch        *b_col_track_wire;   //!
   TBranch        *b_col_track_dedx;   //!
   TBranch        *b_col_track_dqdx;   //!
   TBranch        *b_col_track_rr;   //!
   TBranch        *b_col_track_pitch_hit;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_hit_time;   //!
   TBranch        *b_hit_amp;   //!
   TBranch        *b_hit_wire;   //!
   TBranch        *b_hit_view;   //!
   TBranch        *b_no_primaries;   //!
   TBranch        *b_geant_list_size;   //!
   TBranch        *b_primary;   //!
   TBranch        *b_PDG;   //!
   TBranch        *b_StartKE;   //!
   TBranch        *b_LastKE;   //!
   TBranch        *b_StartEnergy;   //!
   TBranch        *b_StartPx;   //!
   TBranch        *b_StartPy;   //!
   TBranch        *b_StartPz;   //!
   TBranch        *b_EndEnergy;   //!
   TBranch        *b_EndPx;   //!
   TBranch        *b_EndPy;   //!
   TBranch        *b_EndPz;   //!
   TBranch        *b_StartPointx;   //!
   TBranch        *b_StartPointy;   //!
   TBranch        *b_StartPointz;   //!
   TBranch        *b_EndPointx;   //!
   TBranch        *b_EndPointy;   //!
   TBranch        *b_EndPointz;   //!
   TBranch        *b_Process;   //!
   TBranch        *b_NumberDaughters;   //!
   TBranch        *b_Mother;   //!
   TBranch        *b_TrackId;   //!
   TBranch        *b_process_primary;   //!
   TBranch        *b_G4Process;   //!
   TBranch        *b_G4FinalProcess;   //!
   TBranch        *b_NTrTrajPts;   //!
   TBranch        *b_NProtonDaughters;   //!
   TBranch        *b_NNeutronDaughters;   //!
   TBranch        *b_NDTrTrajPts;   //!
   TBranch        *b_DTrackId;   //!
   TBranch        *b_DPdgCode;   //!
   TBranch        *b_DStartEnergy;   //!
   TBranch        *b_DStartP;   //!
   TBranch        *b_MidPosX;   //!
   TBranch        *b_MidPosY;   //!
   TBranch        *b_MidPosZ;   //!
   TBranch        *b_MidPx;   //!
   TBranch        *b_MidPy;   //!
   TBranch        *b_MidPz;   //!
   TBranch        *b_SlabX;   //!
   TBranch        *b_SlabY;   //!
   TBranch        *b_SlabZ;   //!
   TBranch        *b_SlapX;   //!
   TBranch        *b_SlapY;   //!
   TBranch        *b_SlapZ;   //!
   TBranch        *b_SlabN;   //!
   TBranch        *b_SlabE;   //!
   TBranch        *b_LastSlabInt;   //!
   TBranch        *b_DMidPosX;   //!
   TBranch        *b_DMidPosY;   //!
   TBranch        *b_DMidPosZ;   //!
   TBranch        *b_InteractionPoint;   //!
   TBranch        *b_InteractionPointType;   //!
   TBranch        *b_num_tof_objects;   //!
   TBranch        *b_tofObject;   //!
   TBranch        *b_num_wctracks;   //!
   TBranch        *b_wctrk_momentum;   //!
   TBranch        *b_wctrk_XFace;   //!
   TBranch        *b_wctrk_YFace;   //!
   TBranch        *b_wctrk_theta;   //!
   TBranch        *b_wctrk_phi;   //!

   PpolSummary(TTree *tree=0);
   virtual ~PpolSummary();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
      virtual double pointDistance(double x1, double y1, double z1, 
                        double x2, double y2, double z2);

};

#endif

#ifdef PpolSummary_cxx
PpolSummary::PpolSummary(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../files/positive_polarity_preproduction_test_reco_ntuple_20180821_1531.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../files/positive_polarity_preproduction_test_reco_ntuple_20180821_1531.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../files/positive_polarity_preproduction_test_reco_ntuple_20180821_1531.root:/anatree");
      dir->GetObject("anatree",tree);

   }
   Init(tree);
}

PpolSummary::~PpolSummary()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PpolSummary::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PpolSummary::LoadTree(Long64_t entry)
{
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

void PpolSummary::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   track_WC2TPC_match = 0;
   track_primary = 0;
   track_start_x = 0;
   track_start_y = 0;
   track_start_z = 0;
   track_end_x = 0;
   track_end_y = 0;
   track_end_z = 0;
   track_length = 0;
   ntrack_hits = 0;
   track_xpos = 0;
   track_ypos = 0;
   track_zpos = 0;
   nhit_ids = 0;
   ind_track_hits = 0;
   ind_track_ke = 0;
   ind_track_wire = 0;
   ind_track_dedx = 0;
   ind_track_dqdx = 0;
   ind_track_rr = 0;
   ind_track_pitch_hit = 0;
   col_track_hits = 0;
   col_track_ke = 0;
   col_track_x = 0;
   col_track_y = 0;
   col_track_z = 0;
   col_track_wire = 0;
   col_track_dedx = 0;
   col_track_dqdx = 0;
   col_track_rr = 0;
   col_track_pitch_hit = 0;
   hit_time = 0;
   hit_amp = 0;
   hit_wire = 0;
   hit_view = 0;
   PDG = 0;
   StartKE = 0;
   LastKE = 0;
   StartEnergy = 0;
   StartPx = 0;
   StartPy = 0;
   StartPz = 0;
   EndEnergy = 0;
   EndPx = 0;
   EndPy = 0;
   EndPz = 0;
   StartPointx = 0;
   StartPointy = 0;
   StartPointz = 0;
   EndPointx = 0;
   EndPointy = 0;
   EndPointz = 0;
   Process = 0;
   NumberDaughters = 0;
   Mother = 0;
   TrackId = 0;
   process_primary = 0;
   G4Process = 0;
   G4FinalProcess = 0;
   NTrTrajPts = 0;
   NDTrTrajPts = 0;
   DTrackId = 0;
   DPdgCode = 0;
   DStartEnergy = 0;
   DStartP = 0;
   MidPosX = 0;
   MidPosY = 0;
   MidPosZ = 0;
   MidPx = 0;
   MidPy = 0;
   MidPz = 0;
   SlabX = 0;
   SlabY = 0;
   SlabZ = 0;
   SlapX = 0;
   SlapY = 0;
   SlapZ = 0;
   SlabN = 0;
   SlabE = 0;
   DMidPosX = 0;
   DMidPosY = 0;
   DMidPosZ = 0;
   InteractionPoint = 0;
   InteractionPointType = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("efield", efield, &b_efield);
   fChain->SetBranchAddress("t0", &t0, &b_t0);
   fChain->SetBranchAddress("ntracks_reco", &ntracks_reco, &b_ntracks_reco);
   fChain->SetBranchAddress("track_WC2TPC_match", &track_WC2TPC_match, &b_track_WC2TPC_match);
   fChain->SetBranchAddress("track_primary", &track_primary, &b_track_primary);
   fChain->SetBranchAddress("track_start_x", &track_start_x, &b_track_start_x);
   fChain->SetBranchAddress("track_start_y", &track_start_y, &b_track_start_y);
   fChain->SetBranchAddress("track_start_z", &track_start_z, &b_track_start_z);
   fChain->SetBranchAddress("track_end_x", &track_end_x, &b_track_end_x);
   fChain->SetBranchAddress("track_end_y", &track_end_y, &b_track_end_y);
   fChain->SetBranchAddress("track_end_z", &track_end_z, &b_track_end_z);
   fChain->SetBranchAddress("track_length", &track_length, &b_track_length);
   fChain->SetBranchAddress("ntrack_hits", &ntrack_hits, &b_ntrack_hits);
   fChain->SetBranchAddress("track_xpos", &track_xpos, &b_track_xpos);
   fChain->SetBranchAddress("track_ypos", &track_ypos, &b_track_ypos);
   fChain->SetBranchAddress("track_zpos", &track_zpos, &b_track_zpos);
   fChain->SetBranchAddress("nhit_ids", &nhit_ids, &b_nhit_ids);
   fChain->SetBranchAddress("track_spt_idarr", track_spt_idarr, &b_track_spt_idarr);
   fChain->SetBranchAddress("track_spt_earr", track_spt_earr, &b_track_spt_earr);
   fChain->SetBranchAddress("ind_track_hits", &ind_track_hits, &b_ind_track_hits);
   fChain->SetBranchAddress("ind_track_ke", &ind_track_ke, &b_ind_track_ke);
   fChain->SetBranchAddress("ind_track_wire", &ind_track_wire, &b_ind_track_wire);
   fChain->SetBranchAddress("ind_track_dedx", &ind_track_dedx, &b_ind_track_dedx);
   fChain->SetBranchAddress("ind_track_dqdx", &ind_track_dqdx, &b_ind_track_dqdx);
   fChain->SetBranchAddress("ind_track_rr", &ind_track_rr, &b_ind_track_rr);
   fChain->SetBranchAddress("ind_track_pitch_hit", &ind_track_pitch_hit, &b_ind_track_pitch_hit);
   fChain->SetBranchAddress("col_track_hits", &col_track_hits, &b_col_track_hits);
   fChain->SetBranchAddress("col_track_ke", &col_track_ke, &b_col_track_ke);
   fChain->SetBranchAddress("col_track_x", &col_track_x, &b_col_track_x);
   fChain->SetBranchAddress("col_track_y", &col_track_y, &b_col_track_y);
   fChain->SetBranchAddress("col_track_z", &col_track_z, &b_col_track_z);
   fChain->SetBranchAddress("col_track_wire", &col_track_wire, &b_col_track_wire);
   fChain->SetBranchAddress("col_track_dedx", &col_track_dedx, &b_col_track_dedx);
   fChain->SetBranchAddress("col_track_dqdx", &col_track_dqdx, &b_col_track_dqdx);
   fChain->SetBranchAddress("col_track_rr", &col_track_rr, &b_col_track_rr);
   fChain->SetBranchAddress("col_track_pitch_hit", &col_track_pitch_hit, &b_col_track_pitch_hit);
   fChain->SetBranchAddress("nhits", &nhits, &b_nhits);
   fChain->SetBranchAddress("hit_time", &hit_time, &b_hit_time);
   fChain->SetBranchAddress("hit_amp", &hit_amp, &b_hit_amp);
   fChain->SetBranchAddress("hit_wire", &hit_wire, &b_hit_wire);
   fChain->SetBranchAddress("hit_view", &hit_view, &b_hit_view);
   fChain->SetBranchAddress("no_primaries", &no_primaries, &b_no_primaries);
   fChain->SetBranchAddress("geant_list_size", &geant_list_size, &b_geant_list_size);
   fChain->SetBranchAddress("primary_p", &primary_p, &b_primary);
   fChain->SetBranchAddress("PDG", &PDG, &b_PDG);
   fChain->SetBranchAddress("StartKE", &StartKE, &b_StartKE);
   fChain->SetBranchAddress("LastKE", &LastKE, &b_LastKE);
   fChain->SetBranchAddress("StartEnergy", &StartEnergy, &b_StartEnergy);
   fChain->SetBranchAddress("StartPx", &StartPx, &b_StartPx);
   fChain->SetBranchAddress("StartPy", &StartPy, &b_StartPy);
   fChain->SetBranchAddress("StartPz", &StartPz, &b_StartPz);
   fChain->SetBranchAddress("EndEnergy", &EndEnergy, &b_EndEnergy);
   fChain->SetBranchAddress("EndPx", &EndPx, &b_EndPx);
   fChain->SetBranchAddress("EndPy", &EndPy, &b_EndPy);
   fChain->SetBranchAddress("EndPz", &EndPz, &b_EndPz);
   fChain->SetBranchAddress("StartPointx", &StartPointx, &b_StartPointx);
   fChain->SetBranchAddress("StartPointy", &StartPointy, &b_StartPointy);
   fChain->SetBranchAddress("StartPointz", &StartPointz, &b_StartPointz);
   fChain->SetBranchAddress("EndPointx", &EndPointx, &b_EndPointx);
   fChain->SetBranchAddress("EndPointy", &EndPointy, &b_EndPointy);
   fChain->SetBranchAddress("EndPointz", &EndPointz, &b_EndPointz);
   fChain->SetBranchAddress("Process", &Process, &b_Process);
   fChain->SetBranchAddress("NumberDaughters", &NumberDaughters, &b_NumberDaughters);
   fChain->SetBranchAddress("Mother", &Mother, &b_Mother);
   fChain->SetBranchAddress("TrackId", &TrackId, &b_TrackId);
   fChain->SetBranchAddress("process_primary", &process_primary, &b_process_primary);
   fChain->SetBranchAddress("G4Process", &G4Process, &b_G4Process);
   fChain->SetBranchAddress("G4FinalProcess", &G4FinalProcess, &b_G4FinalProcess);
   fChain->SetBranchAddress("NTrTrajPts", &NTrTrajPts, &b_NTrTrajPts);
   fChain->SetBranchAddress("NProtonDaughters", &NProtonDaughters, &b_NProtonDaughters);
   fChain->SetBranchAddress("NNeutronDaughters", &NNeutronDaughters, &b_NNeutronDaughters);
   fChain->SetBranchAddress("NDTrTrajPts", &NDTrTrajPts, &b_NDTrTrajPts);
   fChain->SetBranchAddress("DTrackId", &DTrackId, &b_DTrackId);
   fChain->SetBranchAddress("DPdgCode", &DPdgCode, &b_DPdgCode);
   fChain->SetBranchAddress("DStartEnergy", &DStartEnergy, &b_DStartEnergy);
   fChain->SetBranchAddress("DStartP", &DStartP, &b_DStartP);
   fChain->SetBranchAddress("MidPosX", &MidPosX, &b_MidPosX);
   fChain->SetBranchAddress("MidPosY", &MidPosY, &b_MidPosY);
   fChain->SetBranchAddress("MidPosZ", &MidPosZ, &b_MidPosZ);
   fChain->SetBranchAddress("MidPx", &MidPx, &b_MidPx);
   fChain->SetBranchAddress("MidPy", &MidPy, &b_MidPy);
   fChain->SetBranchAddress("MidPz", &MidPz, &b_MidPz);
   fChain->SetBranchAddress("SlabX", &SlabX, &b_SlabX);
   fChain->SetBranchAddress("SlabY", &SlabY, &b_SlabY);
   fChain->SetBranchAddress("SlabZ", &SlabZ, &b_SlabZ);
   fChain->SetBranchAddress("SlapX", &SlapX, &b_SlapX);
   fChain->SetBranchAddress("SlapY", &SlapY, &b_SlapY);
   fChain->SetBranchAddress("SlapZ", &SlapZ, &b_SlapZ);
   fChain->SetBranchAddress("SlabN", &SlabN, &b_SlabN);
   fChain->SetBranchAddress("SlabE", &SlabE, &b_SlabE);
   fChain->SetBranchAddress("LastSlabInt", &LastSlabInt, &b_LastSlabInt);
   fChain->SetBranchAddress("DMidPosX", &DMidPosX, &b_DMidPosX);
   fChain->SetBranchAddress("DMidPosY", &DMidPosY, &b_DMidPosY);
   fChain->SetBranchAddress("DMidPosZ", &DMidPosZ, &b_DMidPosZ);
   fChain->SetBranchAddress("InteractionPoint", &InteractionPoint, &b_InteractionPoint);
   fChain->SetBranchAddress("InteractionPointType", &InteractionPointType, &b_InteractionPointType);
   fChain->SetBranchAddress("num_tof_objects", &num_tof_objects, &b_num_tof_objects);
   fChain->SetBranchAddress("tofObject", tofObject, &b_tofObject);
   fChain->SetBranchAddress("num_wctracks", &num_wctracks, &b_num_wctracks);
   fChain->SetBranchAddress("wctrk_momentum", wctrk_momentum, &b_wctrk_momentum);
   fChain->SetBranchAddress("wctrk_XFace", wctrk_XFace, &b_wctrk_XFace);
   fChain->SetBranchAddress("wctrk_YFace", wctrk_YFace, &b_wctrk_YFace);
   fChain->SetBranchAddress("wctrk_theta", wctrk_theta, &b_wctrk_theta);
   fChain->SetBranchAddress("wctrk_phi", wctrk_phi, &b_wctrk_phi);
   Notify();
}

Bool_t PpolSummary::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PpolSummary::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PpolSummary::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
double PpolSummary::pointDistance(double x1, double y1, double z1, double x2, double y2, double z2){

    double dist = sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
    return dist;}
#endif // #ifdef PpolSummary_cxx
