#include "BDT_stw_2.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_stw_2::BDT_stw_2(){
    reset();
}

BDT_stw_2::BDT_stw_2(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_stw_2::init(){

    reader = new TMVA::Reader();
    
    reader->AddVariable("stw_2_v_medium_dQ_dx",&stw_2_v_medium_dQ_dx);
    reader->AddVariable("stw_2_v_energy",&stw_2_v_energy);
    reader->AddVariable("stw_2_v_angle",&stw_2_v_angle);
    reader->AddVariable("stw_2_v_dir_length",&stw_2_v_dir_length);
    reader->AddVariable("stw_2_v_max_dQ_dx",&stw_2_v_max_dQ_dx);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_stw_2::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_stw_2::reset(){
    stw_2_v_medium_dQ_dx=0;
    stw_2_v_energy=0;
    stw_2_v_angle=0;
    stw_2_v_dir_length=0;
    stw_2_v_max_dQ_dx=0;
}





