#ifndef BDT_br3_5_h
#define BDT_br3_5_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_br3_5{
    public:
        BDT_br3_5();
        BDT_br3_5(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float br3_5_v_dir_length;
        float br3_5_v_total_length;
        float br3_5_v_flag_avoid_muon_check;
        float br3_5_v_n_seg;
        float br3_5_v_angle;
        float br3_5_v_sg_length;
        float br3_5_v_energy;
        //float br3_5_v_n_main_segs;
        float br3_5_v_n_segs;
        float br3_5_v_shower_main_length;
        float br3_5_v_shower_total_length;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
