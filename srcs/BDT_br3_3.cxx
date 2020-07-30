#include "BDT_br3_3.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_br3_3::BDT_br3_3(){
    reset();
}

BDT_br3_3::BDT_br3_3(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_br3_3::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("br3_3_v_energy",&br3_3_v_energy);
    reader->AddVariable("br3_3_v_angle",&br3_3_v_angle);
    reader->AddVariable("br3_3_v_dir_length",&br3_3_v_dir_length);
    reader->AddVariable("br3_3_v_length",&br3_3_v_length);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_br3_3::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_br3_3::reset(){
    br3_3_v_energy=0;
    br3_3_v_angle=0;
    br3_3_v_dir_length=0;
    br3_3_v_length=0;
}





