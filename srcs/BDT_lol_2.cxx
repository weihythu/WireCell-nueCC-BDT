#include "BDT_lol_2.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_lol_2::BDT_lol_2(){
    reset();
}

BDT_lol_2::BDT_lol_2(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_lol_2::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("lol_2_v_length",&lol_2_v_length);
    reader->AddVariable("lol_2_v_angle",&lol_2_v_angle);
    reader->AddVariable("lol_2_v_type",&lol_2_v_type);
    reader->AddVariable("lol_2_v_vtx_n_segs",&lol_2_v_vtx_n_segs);
    reader->AddVariable("lol_2_v_energy",&lol_2_v_energy);
    reader->AddVariable("lol_2_v_shower_main_length",&lol_2_v_shower_main_length);
    reader->AddVariable("lol_2_v_flag_dir_weak",&lol_2_v_flag_dir_weak);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_lol_2::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_lol_2::reset(){
    lol_2_v_length=0;
    lol_2_v_angle=0;
    lol_2_v_type=0;
    lol_2_v_vtx_n_segs=0;
    lol_2_v_energy=0;
    lol_2_v_shower_main_length=0;
    lol_2_v_flag_dir_weak=0;
}





