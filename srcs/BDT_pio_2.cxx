#include "BDT_pio_2.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_pio_2::BDT_pio_2(){
    reset();
}

BDT_pio_2::BDT_pio_2(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_pio_2::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("pio_2_v_dis2",&pio_2_v_dis2);
    reader->AddVariable("pio_2_v_angle2",&pio_2_v_angle2);
    reader->AddVariable("pio_2_v_acc_length",&pio_2_v_acc_length);
    reader->AddVariable("pio_mip_id",&pio_mip_id);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_pio_2::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_pio_2::reset(){
    pio_2_v_dis2=0;
    pio_2_v_angle2=0;
    pio_2_v_acc_length=0;
    pio_mip_id=0;
}





