#include "BDT_mipquality.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_mipquality::BDT_mipquality(){
    reset();
}

BDT_mipquality::BDT_mipquality(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_mipquality::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("mip_quality_energy",&mip_quality_energy);
    reader->AddVariable("mip_quality_overlap",&mip_quality_overlap);
    reader->AddVariable("mip_quality_n_showers",&mip_quality_n_showers);
    reader->AddVariable("mip_quality_n_tracks",&mip_quality_n_tracks);
    reader->AddVariable("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
    reader->AddVariable("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
    reader->AddVariable("mip_quality_shortest_length",&mip_quality_shortest_length);
    reader->AddVariable("mip_quality_acc_length",&mip_quality_acc_length);
    reader->AddVariable("mip_quality_shortest_angle",&mip_quality_shortest_angle);
    reader->AddVariable("mip_quality_flag_proton",&mip_quality_flag_proton);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_mipquality::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_mipquality::reset(){
    mip_quality_energy=0;
    mip_quality_overlap=0;
    mip_quality_n_showers=0;
    mip_quality_n_tracks=0;
    mip_quality_flag_inside_pi0=0;
    mip_quality_n_pi0_showers=0;
    mip_quality_shortest_length=0;
    mip_quality_acc_length=0;
    mip_quality_shortest_angle=0;
    mip_quality_flag_proton=0;
}





