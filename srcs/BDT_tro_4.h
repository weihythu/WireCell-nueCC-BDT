#ifndef BDT_tro_4_h
#define BDT_tro_4_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_tro_4{
    public:
        BDT_tro_4();
        BDT_tro_4(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float tro_4_v_dir2_mag;
        float tro_4_v_angle;
        float tro_4_v_angle1;
        float tro_4_v_angle2;
        float tro_4_v_length;
        float tro_4_v_length1;
        float tro_4_v_medium_dQ_dx;
        float tro_4_v_end_dQ_dx;
        float tro_4_v_energy;
        float tro_4_v_shower_main_length;
        float tro_4_v_flag_shower_trajectory;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
