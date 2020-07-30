#ifndef BDT_pio_1_h
#define BDT_pio_1_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_pio_1{
    public:
        BDT_pio_1();
        BDT_pio_1(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float pio_1_mass;
        float pio_1_pio_type;
        float pio_1_energy_1;
        float pio_1_energy_2;
        float pio_1_dis_1;
        float pio_1_dis_2;
        float pio_mip_id;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
