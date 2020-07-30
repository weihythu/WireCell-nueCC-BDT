#include "BDT_trimuon.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_trimuon::BDT_trimuon(){
    reset();
}

BDT_trimuon::BDT_trimuon(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_trimuon::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("stem_len_energy",&stem_len_energy);
    reader->AddVariable("stem_len_length",&stem_len_length);
    reader->AddVariable("stem_len_flag_avoid_muon_check",&stem_len_flag_avoid_muon_check);
    reader->AddVariable("stem_len_num_daughters",&stem_len_num_daughters);
    reader->AddVariable("stem_len_daughter_length",&stem_len_daughter_length);
    reader->AddVariable("brm_n_mu_segs",&brm_n_mu_segs);
    reader->AddVariable("brm_Ep",&brm_Ep);
    reader->AddVariable("brm_acc_length",&brm_acc_length);
    reader->AddVariable("brm_shower_total_length",&brm_shower_total_length);
    reader->AddVariable("brm_connected_length",&brm_connected_length);
    reader->AddVariable("brm_n_size",&brm_n_size);
    reader->AddVariable("brm_acc_direct_length",&brm_acc_direct_length);
    reader->AddVariable("brm_n_shower_main_segs",&brm_n_shower_main_segs);
    reader->AddVariable("brm_n_mu_main",&brm_n_mu_main);
    reader->AddVariable("lem_shower_main_length",&lem_shower_main_length);
    reader->AddVariable("lem_n_3seg",&lem_n_3seg);
    reader->AddVariable("lem_e_charge",&lem_e_charge);
    reader->AddVariable("lem_e_dQdx",&lem_e_dQdx);
    reader->AddVariable("lem_shower_num_main_segs",&lem_shower_num_main_segs);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_trimuon::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_trimuon::reset(){
    stem_len_energy=0;
    stem_len_length=0;
    stem_len_flag_avoid_muon_check=0;
    stem_len_num_daughters=0;
    stem_len_daughter_length=0;
    brm_n_mu_segs=0;
    brm_Ep=0;
    brm_acc_length=0;
    brm_shower_total_length=0;
    brm_connected_length=0;
    brm_n_size=0;
    brm_acc_direct_length=0;
    brm_n_shower_main_segs=0;
    brm_n_mu_main=0;
    lem_shower_main_length=0;
    lem_n_3seg=0;
    lem_e_charge=0;
    lem_e_dQdx=0;
    lem_shower_num_main_segs=0;

}
