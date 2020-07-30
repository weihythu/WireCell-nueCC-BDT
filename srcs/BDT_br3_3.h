#ifndef BDT_br3_3_h
#define BDT_br3_3_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_br3_3{
    public:
        BDT_br3_3();
        BDT_br3_3(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float br3_3_v_energy;
        float br3_3_v_angle;
        float br3_3_v_dir_length;
        float br3_3_v_length;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
