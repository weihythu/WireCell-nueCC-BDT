#include "BDT_gap.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_gap::BDT_gap(){
    reset();
}

BDT_gap::BDT_gap(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_gap::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("gap_flag_prolong_u",&gap_flag_prolong_u);
    reader->AddVariable("gap_flag_prolong_v",&gap_flag_prolong_v);
    reader->AddVariable("gap_flag_prolong_w",&gap_flag_prolong_w);
    reader->AddVariable("gap_flag_parallel",&gap_flag_parallel);
    reader->AddVariable("gap_n_points",&gap_n_points);
    reader->AddVariable("gap_n_bad",&gap_n_bad);
    reader->AddVariable("gap_energy",&gap_energy);
    reader->AddVariable("gap_num_valid_tracks",&gap_num_valid_tracks);
    reader->AddVariable("gap_flag_single_shower",&gap_flag_single_shower); 

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_gap::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_gap::reset(){
    gap_flag_prolong_u=0;
    gap_flag_prolong_v=0;
    gap_flag_prolong_w=0;
    gap_flag_parallel=0;
    gap_n_points=0;
    gap_n_bad=0;
    gap_energy=0;
    gap_num_valid_tracks=0;
    gap_flag_single_shower=0;
}





