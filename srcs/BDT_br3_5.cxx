#include "BDT_br3_5.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_br3_5::BDT_br3_5(){
    reset();
}

BDT_br3_5::BDT_br3_5(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_br3_5::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("br3_5_v_dir_length",&br3_5_v_dir_length);
    reader->AddVariable("br3_5_v_total_length",&br3_5_v_total_length);
    reader->AddVariable("br3_5_v_flag_avoid_muon_check",&br3_5_v_flag_avoid_muon_check);
    reader->AddVariable("br3_5_v_n_seg",&br3_5_v_n_seg);
    reader->AddVariable("br3_5_v_angle",&br3_5_v_angle);
    reader->AddVariable("br3_5_v_sg_length",&br3_5_v_sg_length);
    reader->AddVariable("br3_5_v_energy",&br3_5_v_energy);
    //reader->AddVariable("br3_5_v_n_main_segs",&br3_5_v_n_main_segs);
    reader->AddVariable("br3_5_v_n_segs",&br3_5_v_n_segs);
    reader->AddVariable("br3_5_v_shower_main_length",&br3_5_v_shower_main_length);
    reader->AddVariable("br3_5_v_shower_total_length",&br3_5_v_shower_total_length);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_br3_5::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_br3_5::reset(){
        br3_5_v_dir_length=0;
        br3_5_v_total_length=0;
        br3_5_v_flag_avoid_muon_check=0;
        br3_5_v_n_seg=0;
        br3_5_v_angle=0;
        br3_5_v_sg_length=0;
        br3_5_v_energy=0;
        //br3_5_v_n_main_segs=0;
        br3_5_v_n_segs=0;
        br3_5_v_shower_main_length=0;
        br3_5_v_shower_total_length=0;
}





