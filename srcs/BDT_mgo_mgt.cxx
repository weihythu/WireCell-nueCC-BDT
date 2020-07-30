#include "BDT_mgo_mgt.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_mgo_mgt::BDT_mgo_mgt(){
    reset();
}

BDT_mgo_mgt::BDT_mgo_mgt(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_mgo_mgt::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("mgo_energy",&mgo_energy);
    reader->AddVariable("mgo_max_energy",&mgo_max_energy);
    reader->AddVariable("mgo_total_energy",&mgo_total_energy);
    reader->AddVariable("mgo_n_showers",&mgo_n_showers);
    reader->AddVariable("mgo_max_energy_1",&mgo_max_energy_1);
    reader->AddVariable("mgo_max_energy_2",&mgo_max_energy_2);
    reader->AddVariable("mgo_total_other_energy",&mgo_total_other_energy);
    reader->AddVariable("mgo_n_total_showers",&mgo_n_total_showers);
    reader->AddVariable("mgo_total_other_energy_1",&mgo_total_other_energy_1);
    reader->AddVariable("mgt_flag_single_shower",&mgt_flag_single_shower);
    reader->AddVariable("mgt_max_energy",&mgt_max_energy);
    reader->AddVariable("mgt_total_other_energy",&mgt_total_other_energy);
    reader->AddVariable("mgt_max_energy_1",&mgt_max_energy_1);
    reader->AddVariable("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
    reader->AddVariable("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
    reader->AddVariable("mgt_n_direct_showers",&mgt_n_direct_showers);
    reader->AddVariable("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
    reader->AddVariable("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
    reader->AddVariable("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_mgo_mgt::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_mgo_mgt::reset(){
    mgo_energy=0;
    mgo_max_energy=0;
    mgo_total_energy=0;
    mgo_n_showers=0;
    mgo_max_energy_1=0;
    mgo_max_energy_2=0;
    mgo_total_other_energy=0;
    mgo_n_total_showers=0;
    mgo_total_other_energy_1=0;
    mgt_flag_single_shower=0;
    mgt_max_energy=0;
    mgt_total_other_energy=0;
    mgt_max_energy_1=0;
    mgt_e_indirect_max_energy=0;
    mgt_e_direct_max_energy=0;
    mgt_n_direct_showers=0;
    mgt_e_direct_total_energy=0;
    mgt_flag_indirect_max_pio=0;
    mgt_e_indirect_total_energy=0;
}





