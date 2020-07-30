#ifndef BDT_lol_2_h
#define BDT_lol_2_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_lol_2{
    public:
        BDT_lol_2();
        BDT_lol_2(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float lol_2_v_length;
        float lol_2_v_angle;
        float lol_2_v_type;
        float lol_2_v_vtx_n_segs;
        float lol_2_v_energy;
        float lol_2_v_shower_main_length;
        float lol_2_v_flag_dir_weak;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
