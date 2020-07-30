#ifndef BDT_br3_6_h
#define BDT_br3_6_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_br3_6{
    public:
        BDT_br3_6();
        BDT_br3_6(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float br3_6_v_angle;
        float br3_6_v_angle1;
        float br3_6_v_flag_shower_trajectory;
        float br3_6_v_direct_length;
        float br3_6_v_length;
        float br3_6_v_n_other_vtx_segs;
        float br3_6_v_energy;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
