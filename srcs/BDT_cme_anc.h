#ifndef BDT_cme_anc_h
#define BDT_cme_anc_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_cme_anc{
    public:
        BDT_cme_anc();
        BDT_cme_anc(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float cme_mu_energy;
        float cme_energy;
        float cme_mu_length;
        float cme_length;
        float cme_angle_beam;
        float anc_angle;
        float anc_max_angle;
        float anc_max_length;
        float anc_acc_forward_length;
        float anc_acc_backward_length;
        float anc_acc_forward_length1;
        float anc_shower_main_length;
        float anc_shower_total_length;
        float anc_flag_main_outside;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
