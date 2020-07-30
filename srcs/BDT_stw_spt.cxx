#include "BDT_stw_spt.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_stw_spt::BDT_stw_spt(){
    reset();
}

BDT_stw_spt::BDT_stw_spt(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_stw_spt::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("stw_1_energy",&stw_1_energy);
    reader->AddVariable("stw_1_dis",&stw_1_dis);
    reader->AddVariable("stw_1_dQ_dx",&stw_1_dQ_dx);
    reader->AddVariable("stw_1_flag_single_shower",&stw_1_flag_single_shower);
    reader->AddVariable("stw_1_n_pi0",&stw_1_n_pi0);
    reader->AddVariable("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
    reader->AddVariable("spt_shower_main_length",&spt_shower_main_length);
    reader->AddVariable("spt_shower_total_length",&spt_shower_total_length);
    reader->AddVariable("spt_angle_beam",&spt_angle_beam);
    reader->AddVariable("spt_angle_vertical",&spt_angle_vertical);
    reader->AddVariable("spt_max_dQ_dx",&spt_max_dQ_dx);
    reader->AddVariable("spt_angle_beam_1",&spt_angle_beam_1);
    reader->AddVariable("spt_angle_drift",&spt_angle_drift);
    reader->AddVariable("spt_angle_drift_1",&spt_angle_drift_1);
    reader->AddVariable("spt_num_valid_tracks",&spt_num_valid_tracks);
    reader->AddVariable("spt_n_vtx_segs",&spt_n_vtx_segs);
    reader->AddVariable("spt_max_length",&spt_max_length);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_stw_spt::evaluate(){
   
    /* double stw_spt_score = reader->EvaluateMVA("MyBDT"); */

    /* if(stw_spt_score<-0.1) */ 
    /* std::cout */
    /* <<stw_1_energy<<" " */
    /* <<stw_1_dis<<" " */
    /* <<stw_1_dQ_dx<<" " */
    /* <<stw_1_flag_single_shower<<" " */
    /* <<stw_1_n_pi0<<" " */
    /* <<stw_1_num_valid_tracks<<" " */
    /* <<spt_shower_main_length<<" " */
    /* <<spt_shower_total_length<<" " */
    /* <<spt_angle_beam<<" " */
    /* <<spt_angle_vertical<<" " */
    /* <<spt_max_dQ_dx<<" " */
    /* <<spt_angle_beam_1<<" " */
    /* <<spt_angle_drift<<" " */
    /* <<spt_angle_drift_1<<" " */
    /* <<spt_num_valid_tracks<<" " */
    /* <<spt_n_vtx_segs<<" " */
    /* <<spt_max_length<<" " */
    /* <<std::endl; */
   
    //return stw_spt_score;
    return reader->EvaluateMVA("MyBDT");
}

void BDT_stw_spt::reset(){
    stw_1_energy=0;
    stw_1_dis=0;
    stw_1_dQ_dx=0;
    stw_1_flag_single_shower=0;
    stw_1_n_pi0=0;
    stw_1_num_valid_tracks=0;
    spt_shower_main_length=0;
    spt_shower_total_length=0;
    spt_angle_beam=0;
    spt_angle_vertical=0;
    spt_max_dQ_dx=0;
    spt_angle_beam_1=0;
    spt_angle_drift=0;
    spt_angle_drift_1=0;
    spt_num_valid_tracks=0;
    spt_n_vtx_segs=0;
    spt_max_length=0;
}





