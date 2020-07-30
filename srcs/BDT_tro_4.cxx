#include "BDT_tro_4.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_tro_4::BDT_tro_4(){
    reset();
}

BDT_tro_4::BDT_tro_4(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_tro_4::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("tro_4_v_dir2_mag",&tro_4_v_dir2_mag);
    reader->AddVariable("tro_4_v_angle",&tro_4_v_angle);
    reader->AddVariable("tro_4_v_angle1",&tro_4_v_angle1);
    reader->AddVariable("tro_4_v_angle2",&tro_4_v_angle2);
    reader->AddVariable("tro_4_v_length",&tro_4_v_length);
    reader->AddVariable("tro_4_v_length1",&tro_4_v_length1);
    reader->AddVariable("tro_4_v_medium_dQ_dx",&tro_4_v_medium_dQ_dx);
    reader->AddVariable("tro_4_v_end_dQ_dx",&tro_4_v_end_dQ_dx);
    reader->AddVariable("tro_4_v_energy",&tro_4_v_energy);
    reader->AddVariable("tro_4_v_shower_main_length",&tro_4_v_shower_main_length);
    reader->AddVariable("tro_4_v_flag_shower_trajectory",&tro_4_v_flag_shower_trajectory);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_tro_4::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_tro_4::reset(){
    tro_4_v_dir2_mag=0;
    tro_4_v_angle=0;
    tro_4_v_angle1=0;
    tro_4_v_angle2=0;
    tro_4_v_length=0;
    tro_4_v_length1=0;
    tro_4_v_medium_dQ_dx=0;
    tro_4_v_end_dQ_dx=0;
    tro_4_v_energy=0;
    tro_4_v_shower_main_length=0;
    tro_4_v_flag_shower_trajectory=0;
}





