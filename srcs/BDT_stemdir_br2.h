#ifndef BDT_stemdir_br2_h
#define BDT_stemdir_br2_h

#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

class BDT_stemdir_br2{
    public:
        BDT_stemdir_br2();
        BDT_stemdir_br2(TString filename);
        void init();
        void reset();
        double evaluate();

        //BDT input variables
        float stem_dir_flag_single_shower;
        float stem_dir_angle;
        float stem_dir_energy;
        float stem_dir_angle1;
        float stem_dir_angle2;
        float stem_dir_angle3;
        float stem_dir_ratio;
        float br2_num_valid_tracks;
        float br2_n_shower_main_segs;
        float br2_max_angle;
        float br2_sg_length;
        float br2_flag_sg_trajectory;

    private:
        TString _filename_;
        TMVA::Reader *reader;
};

#endif
