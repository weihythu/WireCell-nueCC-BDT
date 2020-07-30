#ifndef BDT_vis_2_h
#define BDT_vis_2_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_vis_2{
    public:
        BDT_vis_2();
        BDT_vis_2(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float vis_2_n_vtx_segs;
        float vis_2_min_angle;
        float vis_2_min_weak_track;
        float vis_2_angle_beam;
        float vis_2_min_angle1;
        float vis_2_iso_angle1;
        float vis_2_min_medium_dQ_dx;
        float vis_2_min_length;
        float vis_2_sg_length;
        float vis_2_max_angle;
        float vis_2_max_weak_track;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
