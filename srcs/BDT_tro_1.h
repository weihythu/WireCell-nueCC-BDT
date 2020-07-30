#ifndef BDT_tro_1_h
#define BDT_tro_1_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_tro_1{
    public:
        BDT_tro_1();
        BDT_tro_1(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float tro_1_v_particle_type;
        float tro_1_v_flag_dir_weak;
        float tro_1_v_min_dis;
        float tro_1_v_sg1_length;
        float tro_1_v_shower_main_length;
        float tro_1_v_max_n_vtx_segs;
        float tro_1_v_tmp_length;
        float tro_1_v_medium_dQ_dx;
        float tro_1_v_dQ_dx_cut;
        float tro_1_v_flag_shower_topology;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
