#ifndef BDT_vis_1_h
#define BDT_vis_1_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_vis_1{
    public:
        BDT_vis_1();
        BDT_vis_1(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float vis_1_n_vtx_segs;
        float vis_1_energy;
        float vis_1_num_good_tracks;
        float vis_1_max_angle;
        float vis_1_max_shower_angle;
        float vis_1_tmp_length1;
        float vis_1_tmp_length2;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
