#ifndef PLOTMODULE_H
#define PLOTMODULE_H
#include <vector>
#include <iostream>
#include <stdio.h>
#include "../Utilities/UtilityFunctions.h"
#include "Riostream.h"
#include "math.h"
#include "../Utilities/ROOTinclude.h"

class PlotModule{

  public:
    PlotModule();
    PlotModule(char *pdfOutputFile);


    TCanvas *c1;
    TPad *pad1;
    TPad *pad2;
    TPad *pad3;
    TPad *pad4;
    TPad *pad5;
    TPad *pad6;
    TPad *pad7;

    TGraph2D *inductionHits;
    TGraph2D *collectionHits;
    TGraph2D *Event3dPrimary;
    TGraph *EventYZprimary;
    TGraph *EventXZprimary;

    void CloseSummary(char* pdfOutputFile);


    void EventSummary( char *pdfOutputFile,int run, int subrun, int event,
                    std::vector<double> beamMatchInfo, std::vector<double> interactionInfo, std::vector<int> *ClusterIDvect,
                     double InitialKE, double beamMass, std::vector<int> branchTracks,
                    std::vector< std::vector<double> > *track_xpos,
                    std::vector< std::vector<double> > *track_ypos,
                    std::vector< std::vector<double> > *track_zpos,
                    std::vector<double> *ntrack_hits,
                    std::vector< std::vector<double> > *col_track_x,
                    std::vector< std::vector<double> > *col_track_y,
                    std::vector< std::vector<double> > *col_track_z,
                    std::vector<double> *col_track_hits,
                    std::vector<std::vector<double>> *col_track_dedx,
                    std::vector<std::vector<double>> *col_track_pitch_hit,
                    double *wctrk_XFace, double *wctrk_YFace,
                    std::vector<int> *nhits,
                    std::vector<float> *hit_time,
                    std::vector<float> *hit_amp,
                    std::vector<int> *hit_wire);



};

#endif

