#include "BDT_pio_1.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_pio_1::BDT_pio_1(){
    reset();
}

BDT_pio_1::BDT_pio_1(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_pio_1::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("pio_1_mass",&pio_1_mass);
    reader->AddVariable("pio_1_pio_type",&pio_1_pio_type);
    reader->AddVariable("pio_1_energy_1",&pio_1_energy_1);
    reader->AddVariable("pio_1_energy_2",&pio_1_energy_2);
    reader->AddVariable("pio_1_dis_1",&pio_1_dis_1);
    reader->AddVariable("pio_1_dis_2",&pio_1_dis_2);
    reader->AddVariable("pio_mip_id",&pio_mip_id);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_pio_1::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_pio_1::reset(){
    pio_1_mass=0;
    pio_1_pio_type=0;
    pio_1_energy_1=0;
    pio_1_energy_2=0;
    pio_1_dis_1=0;
    pio_1_dis_2=0;
    pio_mip_id=0;
}





