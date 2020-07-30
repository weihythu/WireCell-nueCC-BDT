#ifndef BDT_stw_spt_h
#define BDT_stw_spt_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_stw_spt{
    public:
        BDT_stw_spt();
        BDT_stw_spt(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float stw_1_energy;
        float stw_1_dis;
        float stw_1_dQ_dx;
        float stw_1_flag_single_shower;
        float stw_1_n_pi0;
        float stw_1_num_valid_tracks;
        float spt_shower_main_length;
        float spt_shower_total_length;
        float spt_angle_beam;
        float spt_angle_vertical;
        float spt_max_dQ_dx;
        float spt_angle_beam_1;
        float spt_angle_drift;
        float spt_angle_drift_1;
        float spt_num_valid_tracks;
        float spt_n_vtx_segs;
        float spt_max_length;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
