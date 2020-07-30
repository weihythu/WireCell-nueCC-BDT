#ifndef BDT_mipquality_h
#define BDT_mipquality_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_mipquality{
    public:
        BDT_mipquality();
        BDT_mipquality(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float mip_quality_energy;
        float mip_quality_overlap;
        float mip_quality_n_showers;
        float mip_quality_n_tracks;
        float mip_quality_flag_inside_pi0;
        float mip_quality_n_pi0_showers;
        float mip_quality_shortest_length;
        float mip_quality_acc_length;
        float mip_quality_shortest_angle;
        float mip_quality_flag_proton;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
