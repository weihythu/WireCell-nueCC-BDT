#include "BDT_cme_anc.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_cme_anc::BDT_cme_anc(){
    reset();
}

BDT_cme_anc::BDT_cme_anc(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_cme_anc::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("cme_mu_energy",&cme_mu_energy);
    reader->AddVariable("cme_energy",&cme_energy);
    reader->AddVariable("cme_mu_length",&cme_mu_length);
    reader->AddVariable("cme_length",&cme_length);
    reader->AddVariable("cme_angle_beam",&cme_angle_beam);
    reader->AddVariable("anc_angle",&anc_angle);
    reader->AddVariable("anc_max_angle",&anc_max_angle);
    reader->AddVariable("anc_max_length",&anc_max_length);
    reader->AddVariable("anc_acc_forward_length",&anc_acc_forward_length);
    reader->AddVariable("anc_acc_backward_length",&anc_acc_backward_length);
    reader->AddVariable("anc_acc_forward_length1",&anc_acc_forward_length1);
    reader->AddVariable("anc_shower_main_length",&anc_shower_main_length);
    reader->AddVariable("anc_shower_total_length",&anc_shower_total_length);
    reader->AddVariable("anc_flag_main_outside",&anc_flag_main_outside);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_cme_anc::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_cme_anc::reset(){
    cme_mu_energy=0;
    cme_energy=0;
    cme_mu_length=0;
    cme_length=0;
    cme_angle_beam=0;
    anc_angle=0;
    anc_max_angle=0;
    anc_max_length=0;
    anc_acc_forward_length=0;
    anc_acc_backward_length=0;
    anc_acc_forward_length1=0;
    anc_shower_main_length=0;
    anc_shower_total_length=0;
    anc_flag_main_outside=0;
}





