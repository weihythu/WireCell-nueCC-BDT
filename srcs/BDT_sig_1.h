#ifndef BDT_sig_1_h
#define BDT_sig_1_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_sig_1{
    public:
        BDT_sig_1();
        BDT_sig_1(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float sig_1_v_angle;
        float sig_1_v_flag_single_shower;
        float sig_1_v_energy;
        float sig_1_v_energy_1;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
