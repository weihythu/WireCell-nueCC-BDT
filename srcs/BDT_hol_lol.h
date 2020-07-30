#ifndef BDT_hol_lol_h
#define BDT_hol_lol_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_hol_lol{
    public:
        BDT_hol_lol();
        BDT_hol_lol(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float hol_1_n_valid_tracks;
        float hol_1_min_angle;
        float hol_1_energy;
        float hol_1_flag_all_shower;
        float hol_1_min_length;
        float hol_2_min_angle;
        float hol_2_medium_dQ_dx;
        float hol_2_ncount;
        float lol_3_angle_beam;
        float lol_3_n_valid_tracks;
        float lol_3_min_angle;
        float lol_3_vtx_n_segs;
        float lol_3_shower_main_length;
        float lol_3_n_out;
        float lol_3_n_sum;    

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
