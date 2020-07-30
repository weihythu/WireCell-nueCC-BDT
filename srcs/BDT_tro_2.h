#ifndef BDT_tro_2_h
#define BDT_tro_2_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_tro_2{
    public:
        BDT_tro_2();
        BDT_tro_2(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float tro_2_v_energy;
        float tro_2_v_stem_length;
        float tro_2_v_iso_angle;
        float tro_2_v_max_length;
        float tro_2_v_angle;
    
    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
