#include "BDT_stemdir_br2.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_stemdir_br2::BDT_stemdir_br2(){
    reset();
}

BDT_stemdir_br2::BDT_stemdir_br2(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_stemdir_br2::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
    reader->AddVariable("stem_dir_angle",&stem_dir_angle);
    reader->AddVariable("stem_dir_energy",&stem_dir_energy);
    reader->AddVariable("stem_dir_angle1",&stem_dir_angle1);
    reader->AddVariable("stem_dir_angle2",&stem_dir_angle2);
    reader->AddVariable("stem_dir_angle3",&stem_dir_angle3);
    reader->AddVariable("stem_dir_ratio",&stem_dir_ratio);
    reader->AddVariable("br2_num_valid_tracks",&br2_num_valid_tracks);
    reader->AddVariable("br2_n_shower_main_segs",&br2_n_shower_main_segs);
    reader->AddVariable("br2_max_angle",&br2_max_angle);
    reader->AddVariable("br2_sg_length",&br2_sg_length);
    reader->AddVariable("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_stemdir_br2::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_stemdir_br2::reset(){
    stem_dir_flag_single_shower=0;
    stem_dir_angle=0;
    stem_dir_energy=0;
    stem_dir_angle1=0;
    stem_dir_angle2=0;
    stem_dir_angle3=0;
    stem_dir_ratio=0;
    br2_num_valid_tracks=0;
    br2_n_shower_main_segs=0;
    br2_max_angle=0;
    br2_sg_length=0;
    br2_flag_sg_trajectory=0;
}





