#include "BDT_vis_1.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_vis_1::BDT_vis_1(){
    reset();
}

BDT_vis_1::BDT_vis_1(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_vis_1::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
    reader->AddVariable("vis_1_energy",&vis_1_energy);
    reader->AddVariable("vis_1_num_good_tracks",&vis_1_num_good_tracks);
    reader->AddVariable("vis_1_max_angle",&vis_1_max_angle);
    reader->AddVariable("vis_1_max_shower_angle",&vis_1_max_shower_angle);
    reader->AddVariable("vis_1_tmp_length1",&vis_1_tmp_length1);
    reader->AddVariable("vis_1_tmp_length2",&vis_1_tmp_length2);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_vis_1::evaluate(){
   
    /* double vis_1_score = reader->EvaluateMVA("MyBDT"); */

    /* if(vis_1_score<-0.1) */ 
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
   
    //return vis_1_score;
    return reader->EvaluateMVA("MyBDT");
}

void BDT_vis_1::reset(){
    vis_1_n_vtx_segs=0;
    vis_1_energy=0;
    vis_1_num_good_tracks=0;
    vis_1_max_angle=0;
    vis_1_max_shower_angle=0;
    vis_1_tmp_length1=0;
    vis_1_tmp_length2=0;
}





