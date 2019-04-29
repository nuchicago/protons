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
    PlotModule(TPostScript *psOutputFile);


    void CloseSummary(TPostScript * psOutputFile);


    void EventSummary( TPostScript *psOutputFile,int run, int subrun, int event, int ntracks_reco,
                    std::vector<double> beamMatchInfo, double *interactionInfo, std::vector<int> ClusterIDvect,
                     double InitialKE, double IntKE , double beamMass, 
                    std::vector< std::vector<double> > *track_xpos,
                    std::vector< std::vector<double> > *track_ypos,
                    std::vector< std::vector<double> > *track_zpos,
                    std::vector<int> *ntrack_hits,
                    std::vector< std::vector<double> > *col_track_x,
                    std::vector< std::vector<double> > *col_track_y,
                    std::vector< std::vector<double> > *col_track_z,
                    std::vector<int> *col_track_hits,
                    std::vector<std::vector<double>> *col_track_dedx,
                    std::vector<std::vector<double>> *col_track_rr,
                    double wctrk_XFace, double wctrk_YFace,
                    int nhits,
                    std::vector<double> *hit_time,
                    std::vector<double> *hit_amp,
                    std::vector<double> *hit_wire);



};

#endif

