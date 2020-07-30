#include "BDT_tro_1.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_tro_1::BDT_tro_1(){
    reset();
}

BDT_tro_1::BDT_tro_1(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_tro_1::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("tro_1_v_particle_type",&tro_1_v_particle_type);
    reader->AddVariable("tro_1_v_flag_dir_weak",&tro_1_v_flag_dir_weak);
    reader->AddVariable("tro_1_v_min_dis",&tro_1_v_min_dis);
    reader->AddVariable("tro_1_v_sg1_length",&tro_1_v_sg1_length);
    reader->AddVariable("tro_1_v_shower_main_length",&tro_1_v_shower_main_length);
    reader->AddVariable("tro_1_v_max_n_vtx_segs",&tro_1_v_max_n_vtx_segs);
    reader->AddVariable("tro_1_v_tmp_length",&tro_1_v_tmp_length);
    reader->AddVariable("tro_1_v_medium_dQ_dx",&tro_1_v_medium_dQ_dx);
    reader->AddVariable("tro_1_v_dQ_dx_cut", &tro_1_v_dQ_dx_cut);
    reader->AddVariable("tro_1_v_flag_shower_topology", &tro_1_v_flag_shower_topology);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_tro_1::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_tro_1::reset(){
    tro_1_v_particle_type=0;
    tro_1_v_flag_dir_weak=0;
    tro_1_v_min_dis=0;
    tro_1_v_sg1_length=0;
    tro_1_v_shower_main_length=0;
    tro_1_v_max_n_vtx_segs=0;
    tro_1_v_tmp_length=0;
    tro_1_v_medium_dQ_dx=0;
    tro_1_v_dQ_dx_cut=0;
    tro_1_v_flag_shower_topology=0;
}





