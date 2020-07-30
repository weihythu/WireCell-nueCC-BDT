#ifndef BDT_lol_1_h
#define BDT_lol_1_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_lol_1{
    public:
        BDT_lol_1();
        BDT_lol_1(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float lol_1_v_energy;
        float lol_1_v_vtx_n_segs;
        float lol_1_v_nseg;
        float lol_1_v_angle;


    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
