#include "BDT_hol_lol.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_hol_lol::BDT_hol_lol(){
    reset();
}

BDT_hol_lol::BDT_hol_lol(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_hol_lol::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("hol_1_n_valid_tracks",&hol_1_n_valid_tracks);
    reader->AddVariable("hol_1_min_angle",&hol_1_min_angle);
    reader->AddVariable("hol_1_energy",&hol_1_energy);
    reader->AddVariable("hol_1_flag_all_shower",&hol_1_flag_all_shower);
    reader->AddVariable("hol_1_min_length",&hol_1_min_length);
    reader->AddVariable("hol_2_min_angle",&hol_2_min_angle);
    reader->AddVariable("hol_2_medium_dQ_dx",&hol_2_medium_dQ_dx);
    reader->AddVariable("hol_2_ncount",&hol_2_ncount);
    reader->AddVariable("lol_3_angle_beam",&lol_3_angle_beam);
    reader->AddVariable("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
    reader->AddVariable("lol_3_min_angle",&lol_3_min_angle);
    reader->AddVariable("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
    reader->AddVariable("lol_3_shower_main_length",&lol_3_shower_main_length);
    reader->AddVariable("lol_3_n_out",&lol_3_n_out);
    reader->AddVariable("lol_3_n_sum",&lol_3_n_sum);
    
    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_hol_lol::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_hol_lol::reset(){
    hol_1_n_valid_tracks=0;
    hol_1_min_angle=0;
    hol_1_energy=0;
    hol_1_flag_all_shower=0;
    hol_1_min_length=0;
    hol_2_min_angle=0;
    hol_2_medium_dQ_dx=0;
    hol_2_ncount=0;
    lol_3_angle_beam=0;
    lol_3_n_valid_tracks=0;
    lol_3_min_angle=0;
    lol_3_vtx_n_segs=0;
    lol_3_shower_main_length=0;
    lol_3_n_out=0;
    lol_3_n_sum=0;    
}





