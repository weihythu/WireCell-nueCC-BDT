#include "BDT_br4_tro.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_br4_tro::BDT_br4_tro(){
    reset();
}

BDT_br4_tro::BDT_br4_tro(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_br4_tro::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("br4_1_shower_main_length",&br4_1_shower_main_length);
    reader->AddVariable("br4_1_shower_total_length",&br4_1_shower_total_length);
    reader->AddVariable("br4_1_min_dis",&br4_1_min_dis);
    reader->AddVariable("br4_1_energy",&br4_1_energy);
    reader->AddVariable("br4_1_flag_avoid_muon_check",&br4_1_flag_avoid_muon_check);
    reader->AddVariable("br4_1_n_vtx_segs",&br4_1_n_vtx_segs);
    reader->AddVariable("br4_1_n_main_segs",&br4_1_n_main_segs);
    reader->AddVariable("br4_2_ratio_45",&br4_2_ratio_45);
    reader->AddVariable("br4_2_ratio_35",&br4_2_ratio_35);
    reader->AddVariable("br4_2_ratio_25",&br4_2_ratio_25);
    reader->AddVariable("br4_2_ratio_15",&br4_2_ratio_15);
    reader->AddVariable("br4_2_ratio1_45",&br4_2_ratio1_45);
    reader->AddVariable("br4_2_ratio1_35",&br4_2_ratio1_35);
    reader->AddVariable("br4_2_ratio1_25",&br4_2_ratio1_25);
    reader->AddVariable("br4_2_ratio1_15",&br4_2_ratio1_15);
    reader->AddVariable("br4_2_iso_angle",&br4_2_iso_angle);
    reader->AddVariable("br4_2_iso_angle1",&br4_2_iso_angle1);
    reader->AddVariable("br4_2_angle",&br4_2_angle);
    reader->AddVariable("tro_3_stem_length",&tro_3_stem_length);
    reader->AddVariable("tro_3_n_muon_segs",&tro_3_n_muon_segs);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_br4_tro::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_br4_tro::reset(){
    br4_1_shower_main_length=0;
    br4_1_shower_total_length=0;
    br4_1_min_dis=0;
    br4_1_energy=0;
    br4_1_flag_avoid_muon_check=0;
    br4_1_n_vtx_segs=0;
    br4_1_n_main_segs=0;
    br4_2_ratio_45=0;
    br4_2_ratio_35=0;
    br4_2_ratio_25=0;
    br4_2_ratio_15=0;
    br4_2_ratio1_45=0;
    br4_2_ratio1_35=0;
    br4_2_ratio1_25=0;
    br4_2_ratio1_15=0;
    br4_2_iso_angle=0;
    br4_2_iso_angle1=0;
    br4_2_angle=0;
    tro_3_stem_length=0;
    tro_3_n_muon_segs=0;
}





