#ifndef BDT_sig_2_h
#define BDT_sig_2_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_sig_2{
    public:
        BDT_sig_2();
        BDT_sig_2(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float sig_2_v_energy;
        float sig_2_v_shower_angle;
        float sig_2_v_flag_single_shower;
        float sig_2_v_medium_dQ_dx;
        float sig_2_v_start_dQ_dx;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
