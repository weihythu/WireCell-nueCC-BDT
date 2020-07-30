#ifndef BDT_gap_h
#define BDT_gap_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_gap{
    public:
        BDT_gap();
        BDT_gap(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float gap_flag_prolong_u;
        float gap_flag_prolong_v;
        float gap_flag_prolong_w;
        float gap_flag_parallel;
        float gap_n_points;
        float gap_n_bad;
        float gap_energy;
        float gap_num_valid_tracks;
        float gap_flag_single_shower;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
