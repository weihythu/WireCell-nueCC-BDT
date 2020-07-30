#include "BDT_stw_3.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_stw_3::BDT_stw_3(){
    reset();
}

BDT_stw_3::BDT_stw_3(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_stw_3::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("stw_3_v_angle",&stw_3_v_angle);
    reader->AddVariable("stw_3_v_dir_length",&stw_3_v_dir_length);
    reader->AddVariable("stw_3_v_energy",&stw_3_v_energy);
    reader->AddVariable("stw_3_v_medium_dQ_dx",&stw_3_v_medium_dQ_dx);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_stw_3::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_stw_3::reset(){
    stw_3_v_angle=0;
    stw_3_v_dir_length=0;
    stw_3_v_energy=0;
    stw_3_v_medium_dQ_dx=0;
}





