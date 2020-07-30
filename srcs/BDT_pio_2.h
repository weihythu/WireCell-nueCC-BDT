#ifndef BDT_pio_2_h
#define BDT_pio_2_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_pio_2{
    public:
        BDT_pio_2();
        BDT_pio_2(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float pio_2_v_dis2;
        float pio_2_v_angle2;
        float pio_2_v_acc_length;
        float pio_mip_id;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
