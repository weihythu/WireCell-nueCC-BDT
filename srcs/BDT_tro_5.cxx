#include "BDT_tro_5.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_tro_5::BDT_tro_5(){
    reset();
}

BDT_tro_5::BDT_tro_5(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_tro_5::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("tro_5_v_max_angle",&tro_5_v_max_angle);
    reader->AddVariable("tro_5_v_min_angle",&tro_5_v_min_angle);
    reader->AddVariable("tro_5_v_max_length",&tro_5_v_max_length);
    reader->AddVariable("tro_5_v_iso_angle",&tro_5_v_iso_angle);
    reader->AddVariable("tro_5_v_n_vtx_segs",&tro_5_v_n_vtx_segs);
    reader->AddVariable("tro_5_v_min_count",&tro_5_v_min_count);
    reader->AddVariable("tro_5_v_max_count",&tro_5_v_max_count);
    reader->AddVariable("tro_5_v_energy",&tro_5_v_energy);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_tro_5::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_tro_5::reset(){
    tro_5_v_max_angle=0;
    tro_5_v_min_angle=0;
    tro_5_v_max_length=0;
    tro_5_v_iso_angle=0;
    tro_5_v_n_vtx_segs=0;
    tro_5_v_min_count=0;
    tro_5_v_max_count=0;
    tro_5_v_energy=0;
}





