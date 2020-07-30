#include "BDT_sig_1.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_sig_1::BDT_sig_1(){
    reset();
}

BDT_sig_1::BDT_sig_1(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_sig_1::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("sig_1_v_angle",&sig_1_v_angle);
    reader->AddVariable("sig_1_v_flag_single_shower",&sig_1_v_flag_single_shower);
    reader->AddVariable("sig_1_v_energy",&sig_1_v_energy);
    reader->AddVariable("sig_1_v_energy_1",&sig_1_v_energy_1);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_sig_1::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_sig_1::reset(){
    sig_1_v_angle=0;
    sig_1_v_flag_single_shower=0;
    sig_1_v_energy=0;
    sig_1_v_energy_1=0;
}





