#include "BDT_lol_1.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_lol_1::BDT_lol_1(){
    reset();
}

BDT_lol_1::BDT_lol_1(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_lol_1::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("lol_1_v_energy",&lol_1_v_energy);
    reader->AddVariable("lol_1_v_vtx_n_segs",&lol_1_v_vtx_n_segs);
    reader->AddVariable("lol_1_v_nseg",&lol_1_v_nseg);
    reader->AddVariable("lol_1_v_angle",&lol_1_v_angle);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_lol_1::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_lol_1::reset(){
    lol_1_v_energy=0;
    lol_1_v_vtx_n_segs=0;
    lol_1_v_nseg=0;
    lol_1_v_angle=0;
}





