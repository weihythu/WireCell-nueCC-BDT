#ifndef BDT_mgo_mgt_h
#define BDT_mgo_mgt_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_mgo_mgt{
    public:
        BDT_mgo_mgt();
        BDT_mgo_mgt(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float mgo_energy;
        float mgo_max_energy;
        float mgo_total_energy;
        float mgo_n_showers;
        float mgo_max_energy_1;
        float mgo_max_energy_2;
        float mgo_total_other_energy;
        float mgo_n_total_showers;
        float mgo_total_other_energy_1;
        float mgt_flag_single_shower;
        float mgt_max_energy;
        float mgt_total_other_energy;
        float mgt_max_energy_1;
        float mgt_e_indirect_max_energy;
        float mgt_e_direct_max_energy;
        float mgt_n_direct_showers;
        float mgt_e_direct_total_energy;
        float mgt_flag_indirect_max_pio;
        float mgt_e_indirect_total_energy;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
