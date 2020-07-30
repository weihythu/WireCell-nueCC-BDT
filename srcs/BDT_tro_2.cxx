#include "BDT_tro_2.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_tro_2::BDT_tro_2(){
    reset();
}

BDT_tro_2::BDT_tro_2(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_tro_2::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("tro_2_v_energy",&tro_2_v_energy);
    reader->AddVariable("tro_2_v_stem_length",&tro_2_v_stem_length);
    reader->AddVariable("tro_2_v_iso_angle",&tro_2_v_iso_angle);
    reader->AddVariable("tro_2_v_max_length",&tro_2_v_max_length);
    reader->AddVariable("tro_2_v_angle",&tro_2_v_angle);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_tro_2::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_tro_2::reset(){
    tro_2_v_energy=0;
    tro_2_v_stem_length=0;
    tro_2_v_iso_angle=0;
    tro_2_v_max_length=0;
    tro_2_v_angle=0;
}





