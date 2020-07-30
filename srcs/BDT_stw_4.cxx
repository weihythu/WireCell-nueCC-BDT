#include "BDT_stw_4.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_stw_4::BDT_stw_4(){
    reset();
}

BDT_stw_4::BDT_stw_4(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_stw_4::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("stw_4_v_angle",&stw_4_v_angle);
    reader->AddVariable("stw_4_v_dis",&stw_4_v_dis);
    reader->AddVariable("stw_4_v_energy",&stw_4_v_energy);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_stw_4::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_stw_4::reset(){
    stw_4_v_angle=0;
    stw_4_v_dis=0;
    stw_4_v_energy=0;
}





