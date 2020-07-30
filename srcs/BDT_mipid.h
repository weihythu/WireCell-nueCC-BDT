#ifndef BDT_mipid_h
#define BDT_mipid_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_mipid{
    public:
        BDT_mipid();
        BDT_mipid(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float mip_energy;
        float mip_n_end_reduction;    
        float mip_n_first_mip;
        float mip_n_first_non_mip;
        float mip_n_first_non_mip_1;
        float mip_n_first_non_mip_2;
        float mip_vec_dQ_dx_0;
        float mip_vec_dQ_dx_1;
        float mip_max_dQ_dx_sample;
        float mip_n_below_threshold;
        float mip_n_below_zero;
        float mip_n_lowest;
        float mip_n_highest;
        float mip_lowest_dQ_dx;
        float mip_highest_dQ_dx;
        float mip_medium_dQ_dx;
        float mip_stem_length;
        float mip_length_main;
        float mip_length_total;
        float mip_angle_beam;
        float mip_iso_angle;
        float mip_n_vertex;
        float mip_n_good_tracks;
        float mip_E_indirect_max_energy;
        float mip_flag_all_above;
        float mip_min_dQ_dx_5;
        float mip_n_other_vertex; 
        float mip_n_stem_size;
        float mip_flag_stem_trajectory;
        float mip_min_dis;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
