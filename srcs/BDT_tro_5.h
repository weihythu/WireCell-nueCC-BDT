#ifndef BDT_tro_5_h
#define BDT_tro_5_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_tro_5{
    public:
        BDT_tro_5();
        BDT_tro_5(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float tro_5_v_max_angle;
        float tro_5_v_min_angle;
        float tro_5_v_max_length;
        float tro_5_v_iso_angle;
        float tro_5_v_n_vtx_segs;
        float tro_5_v_min_count;
        float tro_5_v_max_count;
        float tro_5_v_energy;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
