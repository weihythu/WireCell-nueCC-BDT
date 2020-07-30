#include "BDT_sig_2.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_sig_2::BDT_sig_2(){
    reset();
}

BDT_sig_2::BDT_sig_2(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_sig_2::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("sig_2_v_energy",&sig_2_v_energy);
    reader->AddVariable("sig_2_v_shower_angle",&sig_2_v_shower_angle);
    reader->AddVariable("sig_2_v_flag_single_shower",&sig_2_v_flag_single_shower);
    reader->AddVariable("sig_2_v_medium_dQ_dx",&sig_2_v_medium_dQ_dx);  
    reader->AddVariable("sig_2_v_start_dQ_dx",&sig_2_v_start_dQ_dx);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_sig_2::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_sig_2::reset(){
    sig_2_v_energy=0;
    sig_2_v_shower_angle=0;
    sig_2_v_flag_single_shower=0;
    sig_2_v_medium_dQ_dx=0;
    sig_2_v_start_dQ_dx=0;
}





