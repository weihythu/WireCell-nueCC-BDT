#ifndef BDT_stw_4_h
#define BDT_stw_4_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_stw_4{
    public:
        BDT_stw_4();
        BDT_stw_4(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float stw_4_v_angle;
        float stw_4_v_dis;
        float stw_4_v_energy;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
