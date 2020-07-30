#include "BDT_br3_6.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_br3_6::BDT_br3_6(){
    reset();
}

BDT_br3_6::BDT_br3_6(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_br3_6::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("br3_6_v_angle",&br3_6_v_angle);
    reader->AddVariable("br3_6_v_angle1",&br3_6_v_angle1);
    reader->AddVariable("br3_6_v_flag_shower_trajectory",&br3_6_v_flag_shower_trajectory);
    reader->AddVariable("br3_6_v_direct_length",&br3_6_v_direct_length);
    reader->AddVariable("br3_6_v_length",&br3_6_v_length);
    reader->AddVariable("br3_6_v_n_other_vtx_segs",&br3_6_v_n_other_vtx_segs);
    reader->AddVariable("br3_6_v_energy",&br3_6_v_energy);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_br3_6::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_br3_6::reset(){
    br3_6_v_angle=0;
    br3_6_v_angle1=0;
    br3_6_v_flag_shower_trajectory=0;
    br3_6_v_direct_length=0;
    br3_6_v_length=0;
    br3_6_v_n_other_vtx_segs=0;
    br3_6_v_energy=0;
}





