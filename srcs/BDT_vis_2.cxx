#include "BDT_vis_2.h"

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

BDT_vis_2::BDT_vis_2(){
    reset();
}

BDT_vis_2::BDT_vis_2(TString filename){
        /// BDT initilization
        _filename_ = filename;
        init();
        reset();
}

void BDT_vis_2::init(){

    reader = new TMVA::Reader();

    reader->AddVariable("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
    reader->AddVariable("vis_2_min_angle",&vis_2_min_angle);
    reader->AddVariable("vis_2_min_weak_track",&vis_2_min_weak_track);
    reader->AddVariable("vis_2_angle_beam", &vis_2_angle_beam);
    reader->AddVariable("vis_2_min_angle1",&vis_2_min_angle1);
    reader->AddVariable("vis_2_iso_angle1",&vis_2_iso_angle1);
    reader->AddVariable("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
    reader->AddVariable("vis_2_min_length",&vis_2_min_length);
    reader->AddVariable("vis_2_sg_length",&vis_2_sg_length);
    reader->AddVariable("vis_2_max_angle",&vis_2_max_angle);
    reader->AddVariable("vis_2_max_weak_track",&vis_2_max_weak_track);

    reader->BookMVA( "MyBDT", _filename_);
}

double BDT_vis_2::evaluate(){
    return reader->EvaluateMVA("MyBDT");
}

void BDT_vis_2::reset(){
    vis_2_n_vtx_segs=0;
    vis_2_min_angle=0;
    vis_2_min_weak_track=0;
    vis_2_angle_beam=0;
    vis_2_min_angle1=0;
    vis_2_iso_angle1=0;
    vis_2_min_medium_dQ_dx=0;
    vis_2_min_length=0;
    vis_2_sg_length=0;
    vis_2_max_angle=0;
    vis_2_max_weak_track=0;
}





