#include "BDT_mipid.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_mipid::BDT_mipid(){
    reset();
}

BDT_mipid::BDT_mipid(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_mipid::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("mip_energy", &mip_energy);
    reader->AddVariable("mip_n_end_reduction", &mip_n_end_reduction);
    reader->AddVariable("mip_n_first_mip", &mip_n_first_mip);
    reader->AddVariable("mip_n_first_non_mip", &mip_n_first_non_mip);
    reader->AddVariable("mip_n_first_non_mip_1", &mip_n_first_non_mip_1);
    reader->AddVariable("mip_n_first_non_mip_2", &mip_n_first_non_mip_2);
    reader->AddVariable("mip_vec_dQ_dx_0", &mip_vec_dQ_dx_0);
    reader->AddVariable("mip_vec_dQ_dx_1", &mip_vec_dQ_dx_1);
    reader->AddVariable("mip_max_dQ_dx_sample", &mip_max_dQ_dx_sample);
    reader->AddVariable("mip_n_below_threshold", &mip_n_below_threshold);
    reader->AddVariable("mip_n_below_zero", &mip_n_below_zero);
    reader->AddVariable("mip_n_lowest", &mip_n_lowest);
    reader->AddVariable("mip_n_highest", &mip_n_highest);
    reader->AddVariable("mip_lowest_dQ_dx", &mip_lowest_dQ_dx);
    reader->AddVariable("mip_highest_dQ_dx", &mip_highest_dQ_dx);
    reader->AddVariable("mip_medium_dQ_dx", &mip_medium_dQ_dx);
    reader->AddVariable("mip_stem_length", &mip_stem_length);
    reader->AddVariable("mip_length_main", &mip_length_main);
    reader->AddVariable("mip_length_total", &mip_length_total);
    reader->AddVariable("mip_angle_beam", &mip_angle_beam);
    reader->AddVariable("mip_iso_angle", &mip_iso_angle);
    reader->AddVariable("mip_n_vertex", &mip_n_vertex);
    reader->AddVariable("mip_n_good_tracks", &mip_n_good_tracks);
    reader->AddVariable("mip_E_indirect_max_energy", &mip_E_indirect_max_energy);
    reader->AddVariable("mip_flag_all_above", &mip_flag_all_above);
    reader->AddVariable("mip_min_dQ_dx_5", &mip_min_dQ_dx_5);
    reader->AddVariable("mip_n_other_vertex", &mip_n_other_vertex);
    reader->AddVariable("mip_n_stem_size", &mip_n_stem_size);
    reader->AddVariable("mip_flag_stem_trajectory", &mip_flag_stem_trajectory);
    reader->AddVariable("mip_min_dis", &mip_min_dis);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_mipid::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_mipid::reset(){
    mip_energy=0;
    mip_n_end_reduction=0;    
    mip_n_first_mip=0;
    mip_n_first_non_mip=0;
    mip_n_first_non_mip_1=0;
    mip_n_first_non_mip_2=0;
    mip_vec_dQ_dx_0=0;
    mip_vec_dQ_dx_1=0;
    mip_max_dQ_dx_sample=0;
    mip_n_below_threshold=0;
    mip_n_below_zero=0;
    mip_n_lowest=0;
    mip_n_highest=0;
    mip_lowest_dQ_dx=0;
    mip_highest_dQ_dx=0;
    mip_medium_dQ_dx=0;
    mip_stem_length=0;
    mip_length_main=0;
    mip_length_total=0;
    mip_angle_beam=0;
    mip_iso_angle=0;
    mip_n_vertex=0;
    mip_n_good_tracks=0;
    mip_E_indirect_max_energy=0;
    mip_flag_all_above=0;
    mip_min_dQ_dx_5=0;
    mip_n_other_vertex=0; 
    mip_n_stem_size=0;
    mip_flag_stem_trajectory=0;
    mip_min_dis=0;
}





