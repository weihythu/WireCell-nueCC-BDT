#include "BDT_br1.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_br1::BDT_br1(){
    reset();
}

BDT_br1::BDT_br1(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_br1::init(){

    reader = new TMVA::Reader();
    
    reader->AddVariable("br1_1_shower_type",&br1_1_shower_type);
    reader->AddVariable("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
    reader->AddVariable("br1_1_energy",&br1_1_energy);
    reader->AddVariable("br1_1_n_segs",&br1_1_n_segs);
    reader->AddVariable("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
    reader->AddVariable("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
    reader->AddVariable("br1_1_sg_length",&br1_1_sg_length);
    reader->AddVariable("br1_2_n_connected",&br1_2_n_connected);
    reader->AddVariable("br1_2_max_length",&br1_2_max_length);
    reader->AddVariable("br1_2_n_connected_1",&br1_2_n_connected_1);
    reader->AddVariable("br1_2_n_shower_segs",&br1_2_n_shower_segs);
    reader->AddVariable("br1_2_max_length_ratio",&br1_2_max_length_ratio);
    reader->AddVariable("br1_2_shower_length",&br1_2_shower_length);
    reader->AddVariable("br1_3_n_connected_p",&br1_3_n_connected_p);
    reader->AddVariable("br1_3_max_length_p",&br1_3_max_length_p);
    reader->AddVariable("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_br1::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_br1::reset(){
    br1_1_shower_type=0;
    br1_1_vtx_n_segs=0;
    br1_1_energy=0;
    br1_1_n_segs=0;
    br1_1_flag_sg_topology=0;
    br1_1_flag_sg_trajectory=0;
    br1_1_sg_length=0;
    br1_2_n_connected=0;
    br1_2_max_length=0;
    br1_2_n_connected_1=0;
    br1_2_n_shower_segs=0;
    br1_2_max_length_ratio=0;
    br1_2_shower_length=0;
    br1_3_n_connected_p=0;
    br1_3_max_length_p=0;
    br1_3_n_shower_main_segs=0; 
}





