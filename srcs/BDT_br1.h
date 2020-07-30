#ifndef BDT_br1_h
#define BDT_br1_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_br1{
    public:
        BDT_br1();
        BDT_br1(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float br1_1_shower_type;
        float br1_1_vtx_n_segs;
        float br1_1_energy;
        float br1_1_n_segs;
        float br1_1_flag_sg_topology;
        float br1_1_flag_sg_trajectory;
        float br1_1_sg_length;
        float br1_2_n_connected;
        float br1_2_max_length;
        float br1_2_n_connected_1;
        float br1_2_n_shower_segs;
        float br1_2_max_length_ratio;
        float br1_2_shower_length;
        float br1_3_n_connected_p;
        float br1_3_max_length_p;
        float br1_3_n_shower_main_segs; 
    
    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
