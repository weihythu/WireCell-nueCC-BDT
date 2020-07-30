#ifndef BDT_br4_tro_h
#define BDT_br4_tro_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_br4_tro{
    public:
        BDT_br4_tro();
        BDT_br4_tro(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float br4_1_shower_main_length;
        float br4_1_shower_total_length;
        float br4_1_min_dis;
        float br4_1_energy;
        float br4_1_flag_avoid_muon_check;
        float br4_1_n_vtx_segs;
        float br4_1_n_main_segs;
        float br4_2_ratio_45;
        float br4_2_ratio_35;
        float br4_2_ratio_25;
        float br4_2_ratio_15;
        float br4_2_ratio1_45;
        float br4_2_ratio1_35;
        float br4_2_ratio1_25;
        float br4_2_ratio1_15;
        float br4_2_iso_angle;
        float br4_2_iso_angle1;
        float br4_2_angle;
        float tro_3_stem_length;
        float tro_3_n_muon_segs;
        float tro_3_energy;


    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
