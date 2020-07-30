#ifndef BDT_trimuon_h
#define BDT_trimuon_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_trimuon{
    public:
        BDT_trimuon();
        BDT_trimuon(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float stem_len_energy;
        float stem_len_length;
        float stem_len_flag_avoid_muon_check;
        float stem_len_num_daughters;
        float stem_len_daughter_length;
        float brm_n_mu_segs;
        float brm_Ep;
        float brm_acc_length;
        float brm_shower_total_length;
        float brm_connected_length;
        float brm_n_size;
        float brm_acc_direct_length;
        float brm_n_shower_main_segs;
        float brm_n_mu_main;
        float lem_shower_main_length;
        float lem_n_3seg;
        float lem_e_charge;
        float lem_e_dQdx;
        float lem_shower_num_main_segs;


    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
