#include "BDT_br3.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_br3::BDT_br3(){
    reset();
}

BDT_br3::BDT_br3(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_br3::init(){

    reader = new TMVA::Reader();
    
    reader->AddVariable("br3_1_energy",&br3_1_energy);
    reader->AddVariable("br3_1_n_shower_segments",&br3_1_n_shower_segments);
    reader->AddVariable("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
    reader->AddVariable("br3_1_sg_direct_length",&br3_1_sg_direct_length);
    reader->AddVariable("br3_1_sg_length",&br3_1_sg_length);
    reader->AddVariable("br3_1_total_main_length",&br3_1_total_main_length);
    reader->AddVariable("br3_1_total_length",&br3_1_total_length);
    reader->AddVariable("br3_1_iso_angle",&br3_1_iso_angle);
    reader->AddVariable("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
    reader->AddVariable("br3_2_n_ele",&br3_2_n_ele);
    reader->AddVariable("br3_2_n_other",&br3_2_n_other);
    reader->AddVariable("br3_2_other_fid",&br3_2_other_fid);
    reader->AddVariable("br3_4_acc_length",&br3_4_acc_length);
    reader->AddVariable("br3_4_total_length",&br3_4_total_length);
    reader->AddVariable("br3_7_min_angle",&br3_7_min_angle);
    reader->AddVariable("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
    reader->AddVariable("br3_8_n_main_segs",&br3_8_n_main_segs);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_br3::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_br3::reset(){
    br3_1_energy=0;
    br3_1_n_shower_segments=0;
    br3_1_sg_flag_trajectory=0;
    br3_1_sg_direct_length=0;
    br3_1_sg_length=0;
    br3_1_total_main_length=0;
    br3_1_total_length=0;
    br3_1_iso_angle=0;
    br3_1_sg_flag_topology=0;
    br3_2_n_ele=0;
    br3_2_n_other=0;
    br3_2_other_fid=0;
    br3_4_acc_length=0;
    br3_4_total_length=0;
    br3_7_min_angle=0;
    br3_8_max_dQ_dx=0;
    br3_8_n_main_segs=0;
}





