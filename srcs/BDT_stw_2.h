#ifndef BDT_stw_2_h
#define BDT_stw_2_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_stw_2{
    public:
        BDT_stw_2();
        BDT_stw_2(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float stw_2_v_medium_dQ_dx;
        float stw_2_v_energy;
        float stw_2_v_angle;
        float stw_2_v_dir_length;
        float stw_2_v_max_dQ_dx;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
