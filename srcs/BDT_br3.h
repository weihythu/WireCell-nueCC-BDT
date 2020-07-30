#ifndef BDT_br3_h
#define BDT_br3_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_br3{
    public:
        BDT_br3();
        BDT_br3(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float br3_1_energy;
        float br3_1_n_shower_segments;
        float br3_1_sg_flag_trajectory;
        float br3_1_sg_direct_length;
        float br3_1_sg_length;
        float br3_1_total_main_length;
        float br3_1_total_length;
        float br3_1_iso_angle;
        float br3_1_sg_flag_topology;
        float br3_2_n_ele;
        float br3_2_n_other;
        float br3_2_other_fid;
        float br3_4_acc_length;
        float br3_4_total_length;
        float br3_7_min_angle;
        float br3_8_max_dQ_dx;
        float br3_8_n_main_segs;
   

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
