#include <vector>
#include <iostream>

#include "TString.h" 
#include "TTree.h"
#include "TFile.h"

#include "BDT_mipid.h"
#include "BDT_mipquality.h"
#include "BDT_gap.h"
#include "BDT_hol_lol.h"
#include "BDT_cme_anc.h"
#include "BDT_mgo_mgt.h"
#include "BDT_stw_spt.h"
#include "BDT_trimuon.h"
#include "BDT_br1.h"
#include "BDT_stemdir_br2.h"
#include "BDT_br3.h"
#include "BDT_br3_3.h"
#include "BDT_br3_5.h"
#include "BDT_br3_6.h"
#include "BDT_br4_tro.h"
#include "BDT_pio_1.h"
#include "BDT_pio_2.h"
#include "BDT_vis_1.h"
#include "BDT_vis_2.h"
#include "BDT_stw_2.h"
#include "BDT_stw_3.h"
#include "BDT_stw_4.h"
#include "BDT_sig_1.h"
#include "BDT_sig_2.h"
#include "BDT_lol_1.h"
#include "BDT_lol_2.h"
#include "BDT_tro_1.h"
#include "BDT_tro_2.h"
#include "BDT_tro_4.h"
#include "BDT_tro_5.h"





using namespace std;

int main(int argc, char* argv[]){

    cout << "BDT preprocessing!" << endl;
    if(argc!=4) {cout << "Incorrect input arguments!" << endl; return 1;} 
   
    bool _MC_ = atoi(argv[3]);

    TString filename = argv[1];    
    TFile *inputfile = new TFile(filename);
    TString type = argv[2];
    TTree *T = (TTree*)inputfile->Get(type);

    int run, subrun, event;
    int nueTag;
    int recoFC;
    int truth_nue;
    int truth_CC;
    int truth_inFV;
    int truth_cosmic;
    float trueEnu;
    float weight, lowEweight;
    float trueEdep;
    float nuvtx_diff;
    float showervtx_diff;

    T->SetBranchAddress("run",&run);
    T->SetBranchAddress("subrun",&subrun);
    T->SetBranchAddress("event",&event);
    T->SetBranchAddress("nueTag",&nueTag);
    T->SetBranchAddress("recoFC",&recoFC);
    if(_MC_){
    T->SetBranchAddress("truth_nue",&truth_nue);
    T->SetBranchAddress("truth_CC",&truth_CC);
    T->SetBranchAddress("truth_inFV",&truth_inFV);
    T->SetBranchAddress("truth_cosmic",&truth_cosmic);
    T->SetBranchAddress("trueEnu",&trueEnu);
    T->SetBranchAddress("weight",&weight);
    T->SetBranchAddress("lowEweight",&lowEweight);
    T->SetBranchAddress("trueEdep",&trueEdep);
    T->SetBranchAddress("nuvtx_diff",&nuvtx_diff);
    T->SetBranchAddress("showervtx_diff",&showervtx_diff);
    }
    /// BDT tagger variables
    // mipid
    int mip_flag;
    double mip_energy;
    int mip_n_end_reduction;    
    int mip_n_first_mip;
    int mip_n_first_non_mip;
    int mip_n_first_non_mip_1;
    int mip_n_first_non_mip_2;
    double mip_vec_dQ_dx_0;
    double mip_vec_dQ_dx_1;
    double mip_vec_dQ_dx_2;
    double mip_vec_dQ_dx_3;
    double mip_vec_dQ_dx_4;
    double mip_vec_dQ_dx_5;
    double mip_vec_dQ_dx_6;
    double mip_max_dQ_dx_sample;
    int mip_n_below_threshold;
    int mip_n_below_zero;
    int mip_n_lowest;
    int mip_n_highest;
    double mip_lowest_dQ_dx;
    double mip_highest_dQ_dx;
    double mip_medium_dQ_dx;
    double mip_stem_length;
    double mip_length_main;
    double mip_length_total;
    double mip_angle_beam;
    double mip_iso_angle;
    int mip_n_vertex;
    int mip_n_good_tracks;
    double mip_E_indirect_max_energy;
    int mip_flag_all_above;
    double mip_min_dQ_dx_5;
    int mip_n_other_vertex; 
    int mip_n_stem_size;
    int mip_flag_stem_trajectory;
    double mip_min_dis;
    int mip_filled;


    T->SetBranchAddress("mip_flag",&mip_flag);
    T->SetBranchAddress("mip_filled",&mip_filled);
    T->SetBranchAddress("mip_energy",&mip_energy);
    T->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
    T->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
    T->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
    T->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
    T->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
    T->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
    T->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
    T->SetBranchAddress("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2);
    T->SetBranchAddress("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3);
    T->SetBranchAddress("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4);
    T->SetBranchAddress("mip_vec_dQ_dx_5",&mip_vec_dQ_dx_5);
    T->SetBranchAddress("mip_vec_dQ_dx_6",&mip_vec_dQ_dx_6);
    T->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
    T->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
    T->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
    T->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
    T->SetBranchAddress("mip_n_highest",&mip_n_highest);
    T->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
    T->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
    T->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
    T->SetBranchAddress("mip_stem_length",&mip_stem_length);
    T->SetBranchAddress("mip_length_main",&mip_length_main);
    T->SetBranchAddress("mip_length_total",&mip_length_total);
    T->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
    T->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
    T->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
    T->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
    T->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
    T->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
    T->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
    T->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
    T->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
    T->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
    T->SetBranchAddress("mip_min_dis",&mip_min_dis);

    /// gap
    int gap_flag;
    int gap_filled;
    int gap_flag_prolong_u;
    int gap_flag_prolong_v;
    int gap_flag_prolong_w;
    int gap_flag_parallel;
    int gap_n_points;
    int gap_n_bad;
    double gap_energy;
    int gap_num_valid_tracks;
    int gap_flag_single_shower;

    T->SetBranchAddress("gap_flag",&gap_flag);
    T->SetBranchAddress("gap_filled",&gap_filled);
    T->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
    T->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
    T->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
    T->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
    T->SetBranchAddress("gap_n_points",&gap_n_points);
    T->SetBranchAddress("gap_n_bad",&gap_n_bad);
    T->SetBranchAddress("gap_energy",&gap_energy);
    T->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
    T->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);

    // muon energy and shower angle: cme_anc
    int cme_flag;
    double cme_mu_energy;
    double cme_energy;
    double cme_mu_length;
    double cme_length;
    double cme_angle_beam;

    int anc_flag;
    double anc_angle;
    double anc_max_angle;
    double anc_max_length;
    double anc_acc_forward_length;
    double anc_acc_backward_length;
    double anc_acc_forward_length1;
    double anc_shower_main_length;
    double anc_shower_total_length;
    int anc_flag_main_outside;

    T->SetBranchAddress("cme_flag",&cme_flag);
    T->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
    T->SetBranchAddress("cme_energy",&cme_energy);
    T->SetBranchAddress("cme_mu_length",&cme_mu_length);
    T->SetBranchAddress("cme_length",&cme_length);
    T->SetBranchAddress("cme_angle_beam",&cme_angle_beam);

    T->SetBranchAddress("anc_flag",&anc_flag);
    T->SetBranchAddress("anc_angle",&anc_angle);
    T->SetBranchAddress("anc_max_angle",&anc_max_angle);
    T->SetBranchAddress("anc_max_length",&anc_max_length);
    T->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
    T->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
    T->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
    T->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
    T->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
    T->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);

    // multiple gammas taggers
    // mgo_mgt
    int mgo_flag;
    double mgo_energy;
    double mgo_max_energy;
    double mgo_total_energy;
    int mgo_n_showers;
    double mgo_max_energy_1;
    double mgo_max_energy_2;
    double mgo_total_other_energy;
    int mgo_n_total_showers;
    double mgo_total_other_energy_1;

    int mgt_flag;
    int mgt_flag_single_shower;
    double mgt_max_energy;
    double mgt_total_other_energy;
    double mgt_max_energy_1;
    double mgt_e_indirect_max_energy;
    double mgt_e_direct_max_energy;
    int mgt_n_direct_showers;
    double mgt_e_direct_total_energy;
    int mgt_flag_indirect_max_pio;
    double mgt_e_indirect_total_energy;

    T->SetBranchAddress("mgo_flag",&mgo_flag);
    T->SetBranchAddress("mgo_energy",&mgo_energy);
    T->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
    T->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
    T->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
    T->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
    T->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
    T->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
    T->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
    T->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);

    T->SetBranchAddress("mgt_flag",&mgt_flag);
    T->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
    T->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
    T->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
    T->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
    T->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
    T->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
    T->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
    T->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
    T->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
    T->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);

    // br1
    int br1_flag;

    int br1_1_shower_type;
    int br1_1_vtx_n_segs;
    double br1_1_energy;
    int br1_1_n_segs;
    int br1_1_flag_sg_topology;
    int br1_1_flag_sg_trajectory;
    double br1_1_sg_length;

    int br1_2_n_connected;
    double br1_2_max_length;
    int br1_2_n_connected_1;
    int br1_2_n_shower_segs;
    double br1_2_max_length_ratio;
    double br1_2_shower_length;

    int br1_3_n_connected_p;
    double br1_3_max_length_p;
    int br1_3_n_shower_main_segs;

    T->SetBranchAddress("br1_flag",&br1_flag);

    T->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
    T->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
    T->SetBranchAddress("br1_1_energy",&br1_1_energy);
    T->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
    T->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
    T->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
    T->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);

    T->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
    T->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
    T->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
    T->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
    T->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
    T->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);

    T->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
    T->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
    T->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);

    // stem_dir_br2 merge
    int stem_dir_flag;
    int stem_dir_flag_single_shower;
    double stem_dir_angle;
    double stem_dir_energy;
    double stem_dir_angle1;
    double stem_dir_angle2;
    double stem_dir_angle3;
    double stem_dir_ratio;

    int br2_flag;
    int br2_num_valid_tracks;
    int br2_n_shower_main_segs;
    double br2_max_angle;
    double br2_sg_length;
    int br2_flag_sg_trajectory;

    T->SetBranchAddress("stem_dir_flag",&stem_dir_flag);
    T->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
    T->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
    T->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
    T->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
    T->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
    T->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
    T->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);

    T->SetBranchAddress("br2_flag",&br2_flag);
    T->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
    T->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
    T->SetBranchAddress("br2_max_angle",&br2_max_angle);
    T->SetBranchAddress("br2_sg_length",&br2_sg_length);
    T->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);

    // stem length, low energy michel, broken muon
    int stem_len_flag;
    double stem_len_energy;
    double stem_len_length;
    int stem_len_flag_avoid_muon_check;
    int stem_len_num_daughters;
    double stem_len_daughter_length;

    int brm_flag;
    int brm_n_mu_segs;
    double brm_Ep;
    double brm_acc_length;
    double brm_shower_total_length;
    double brm_connected_length;
    int brm_n_size;
    double brm_acc_direct_length;
    int brm_n_shower_main_segs;
    int brm_n_mu_main;

    int lem_flag;
    double lem_shower_main_length;
    int lem_n_3seg;
    double lem_e_charge;
    double lem_e_dQdx;
    int lem_shower_num_main_segs;

    T->SetBranchAddress("stem_len_flag", &stem_len_flag);
    T->SetBranchAddress("stem_len_energy", &stem_len_energy);
    T->SetBranchAddress("stem_len_length", &stem_len_length);
    T->SetBranchAddress("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check);
    T->SetBranchAddress("stem_len_num_daughters", &stem_len_num_daughters);
    T->SetBranchAddress("stem_len_daughter_length", &stem_len_daughter_length);

    T->SetBranchAddress("brm_n_flag",&brm_flag);
    T->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
    T->SetBranchAddress("brm_Ep",&brm_Ep);
    T->SetBranchAddress("brm_acc_length",&brm_acc_length);
    T->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
    T->SetBranchAddress("brm_connected_length",&brm_connected_length);
    T->SetBranchAddress("brm_n_size",&brm_n_size);
    T->SetBranchAddress("brm_nacc_direct_length",&brm_acc_direct_length);
    T->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
    T->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);

    T->SetBranchAddress("lem_flag",&lem_flag);
    T->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
    T->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
    T->SetBranchAddress("lem_e_charge",&lem_e_charge);
    T->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
    T->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);

    // br4 + tro single input variables
    int br4_flag;
    
    double br4_1_shower_main_length;
    double br4_1_shower_total_length;
    double br4_1_min_dis;
    double br4_1_energy;
    int br4_1_flag_avoid_muon_check;
    int br4_1_n_vtx_segs;
    int br4_1_n_main_segs;

    double br4_2_ratio_45;
    double br4_2_ratio_35;
    double br4_2_ratio_25;
    double br4_2_ratio_15;
    double br4_2_ratio1_45;
    double br4_2_ratio1_35;
    double br4_2_ratio1_25;
    double br4_2_ratio1_15;
    double br4_2_iso_angle;
    double br4_2_iso_angle1;
    double br4_2_angle;

    int tro_3_flag;
    double tro_3_stem_length;
    int tro_3_n_muon_segs;
    double tro_3_energy;

    T->SetBranchAddress("tro_3_flag",&tro_3_flag);
    T->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
    T->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
    T->SetBranchAddress("tro_3_energy",&tro_3_energy);

    T->SetBranchAddress("br4_flag", &br4_flag);
    T->SetBranchAddress("br4_1_shower_main_length", &br4_1_shower_main_length);
    T->SetBranchAddress("br4_1_shower_total_length", &br4_1_shower_total_length);
    T->SetBranchAddress("br4_1_min_dis", &br4_1_min_dis);
    T->SetBranchAddress("br4_1_energy", &br4_1_energy);
    T->SetBranchAddress("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check);
    T->SetBranchAddress("br4_1_n_vtx_segs", &br4_1_n_vtx_segs);
    T->SetBranchAddress("br4_1_br4_1_n_main_segs", &br4_1_n_main_segs);

    T->SetBranchAddress("br4_2_ratio_45", &br4_2_ratio_45);
    T->SetBranchAddress("br4_2_ratio_35", &br4_2_ratio_35);
    T->SetBranchAddress("br4_2_ratio_25", &br4_2_ratio_25);
    T->SetBranchAddress("br4_2_ratio_15", &br4_2_ratio_15);
    T->SetBranchAddress("br4_2_ratio1_45", &br4_2_ratio1_45);
    T->SetBranchAddress("br4_2_ratio1_35", &br4_2_ratio1_35);
    T->SetBranchAddress("br4_2_ratio1_25", &br4_2_ratio1_25);
    T->SetBranchAddress("br4_2_ratio1_15", &br4_2_ratio1_15);
    T->SetBranchAddress("br4_2_iso_angle", &br4_2_iso_angle);
    T->SetBranchAddress("br4_2_iso_angle1", &br4_2_iso_angle1);
    T->SetBranchAddress("br4_2_angle", &br4_2_angle);

    //mip quality
    int mip_quality_flag;
    int mip_quality_filled;
    double mip_quality_energy;
    int mip_quality_overlap;
    int mip_quality_n_showers;
    int mip_quality_n_tracks;
    int mip_quality_flag_inside_pi0;
    int mip_quality_n_pi0_showers;
    double mip_quality_shortest_length;
    double mip_quality_acc_length;
    double mip_quality_shortest_angle;
    int mip_quality_flag_proton;

    T->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
    T->SetBranchAddress("mip_quality_filled",&mip_quality_filled);
    T->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
    T->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
    T->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
    T->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
    T->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
    T->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
    T->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
    T->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
    T->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
    T->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);

    // pio tagger1
    int pio_flag;
    int pio_filled;
    int pio_mip_id;
    int pio_flag_pio;

    int pio_1_flag;
    double pio_1_mass;
    int pio_1_pio_type;
    double pio_1_energy_1;
    double pio_1_energy_2;
    double pio_1_dis_1;
    double pio_1_dis_2;

    T->SetBranchAddress("pio_flag",&pio_flag);
    T->SetBranchAddress("pio_filled",&pio_filled);
    T->SetBranchAddress("pio_mip_id",&pio_mip_id);
    T->SetBranchAddress("pio_flag_pio",&pio_flag_pio);

    T->SetBranchAddress("pio_1_flag",&pio_1_flag);
    T->SetBranchAddress("pio_1_mass",&pio_1_mass);
    T->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
    T->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
    T->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
    T->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
    T->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);

    // pio tagger 2
    // share flags with pio tagger 1
    std::vector<double> *pio_2_v_dis2 = new std::vector<double>;
    std::vector<double> *pio_2_v_angle2 = new std::vector<double>;
    std::vector<double> *pio_2_v_acc_length = new std::vector<double>;
    std::vector<int> *pio_2_v_flag = new std::vector<int>;

    T->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
    T->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
    T->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);
    T->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);

    // stw + spt for single gamma
    int stw_1_flag;
    double stw_1_energy;
    double stw_1_dis;
    double stw_1_dQ_dx;
    int stw_1_flag_single_shower;
    int stw_1_n_pi0;
    int stw_1_num_valid_tracks;

    int spt_flag;
    int spt_flag_single_shower;
    double spt_shower_main_length;
    double spt_shower_total_length;
    double spt_angle_beam;
    double spt_angle_vertical;
    double spt_max_dQ_dx;
    double spt_angle_beam_1;
    double spt_angle_drift;
    double spt_angle_drift_1;
    int spt_num_valid_tracks;
    double spt_n_vtx_segs;
    double spt_max_length;

    T->SetBranchAddress("stw_1_flag",&stw_1_flag); 
    T->SetBranchAddress("stw_1_energy",&stw_1_energy);
    T->SetBranchAddress("stw_1_dis",&stw_1_dis);
    T->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
    T->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
    T->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
    T->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
    
    T->SetBranchAddress("spt_flag", &spt_flag);
    T->SetBranchAddress("spt_flag_single_shower", &spt_flag_single_shower);
    T->SetBranchAddress("spt_shower_main_length", &spt_shower_main_length);
    T->SetBranchAddress("spt_shower_total_length", &spt_shower_total_length);
    T->SetBranchAddress("spt_angle_beam", &spt_angle_beam);
    T->SetBranchAddress("spt_angle_vertical", &spt_angle_vertical);
    T->SetBranchAddress("spt_max_dQ_dx", &spt_max_dQ_dx);
    T->SetBranchAddress("spt_angle_beam_1", &spt_angle_beam_1);
    T->SetBranchAddress("spt_angle_drift", &spt_angle_drift);
    T->SetBranchAddress("spt_angle_drift_1", &spt_angle_drift_1);
    T->SetBranchAddress("spt_num_valid_tracks", &spt_num_valid_tracks);
    T->SetBranchAddress("spt_n_vtx_segs", &spt_n_vtx_segs);
    T->SetBranchAddress("spt_max_length", &spt_max_length);
    
    // vis tagger 1
    int vis_1_flag;
    int vis_1_filled;
    int vis_1_n_vtx_segs;
    double vis_1_energy;
    int vis_1_num_good_tracks;
    double vis_1_max_angle;
    double vis_1_max_shower_angle;
    double vis_1_tmp_length1;
    double vis_1_tmp_length2;
    double vis_1_particle_type;

    T->SetBranchAddress("vis_1_flag",&vis_1_flag);
    T->SetBranchAddress("vis_1_filled",&vis_1_filled);
    T->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
    T->SetBranchAddress("vis_1_energy",&vis_1_energy);
    T->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
    T->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
    T->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
    T->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
    T->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
    T->SetBranchAddress("vis_1_particle_type",&vis_1_particle_type);

    // vis tagger 2
    int vis_2_flag;
    int vis_2_filled;
    int vis_2_n_vtx_segs;
    double vis_2_min_angle;
    int vis_2_min_weak_track;
    double vis_2_angle_beam;
    double vis_2_min_angle1;
    double vis_2_iso_angle1;
    double vis_2_min_medium_dQ_dx;
    double vis_2_min_length;
    double vis_2_sg_length;
    double vis_2_max_angle;
    int vis_2_max_weak_track;

    T->SetBranchAddress("vis_2_flab",&vis_2_flag);
    T->SetBranchAddress("vis_2_filled",&vis_2_filled);
    T->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
    T->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
    T->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
    T->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
    T->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
    T->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
    T->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
    T->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
    T->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
    T->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
    T->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);

    // stw tagger 2
    std::vector<int> *stw_2_v_flag = new std::vector<int>;
    std::vector<double> *stw_2_v_medium_dQ_dx = new std::vector<double>;
    std::vector<double> *stw_2_v_energy = new std::vector<double>;
    std::vector<double> *stw_2_v_angle = new std::vector<double>;
    std::vector<double> *stw_2_v_dir_length = new std::vector<double>;
    std::vector<double> *stw_2_v_max_dQ_dx = new std::vector<double>;

    T->SetBranchAddress("stw_2_v_flag", &stw_2_v_flag);
    T->SetBranchAddress("stw_2_v_medium_dQ_dx", &stw_2_v_medium_dQ_dx);
    T->SetBranchAddress("stw_2_v_energy", &stw_2_v_energy);
    T->SetBranchAddress("stw_2_v_angle", &stw_2_v_angle);
    T->SetBranchAddress("stw_2_v_dir_length", &stw_2_v_dir_length);
    T->SetBranchAddress("stw_2_v_max_dQ_dx", &stw_2_v_max_dQ_dx);

    // stw tagger 3
    std::vector<int> *stw_3_v_flag = new std::vector<int>;
    std::vector<double> *stw_3_v_angle = new std::vector<double>;
    std::vector<double> *stw_3_v_dir_length = new std::vector<double>;
    std::vector<double> *stw_3_v_energy = new std::vector<double>;
    std::vector<double> *stw_3_v_medium_dQ_dx = new std::vector<double>;
    
    T->SetBranchAddress("stw_3_v_flag",&stw_3_v_flag);
    T->SetBranchAddress("stw_3_v_angle",&stw_3_v_angle);
    T->SetBranchAddress("stw_3_v_dir_length",&stw_3_v_dir_length);
    T->SetBranchAddress("stw_3_v_energy",&stw_3_v_energy);
    T->SetBranchAddress("stw_3_v_medium_dQ_dx",&stw_3_v_medium_dQ_dx);
    
    // stw tagger 4
    std::vector<int> *stw_4_v_flag = new std::vector<int>;
    std::vector<double> *stw_4_v_angle = new std::vector<double>;
    std::vector<double> *stw_4_v_dis = new std::vector<double>;
    std::vector<double> *stw_4_v_energy = new std::vector<double>;

    T->SetBranchAddress("stw_4_v_flag",&stw_4_v_flag);
    T->SetBranchAddress("stw_4_v_angle",&stw_4_v_angle);
    T->SetBranchAddress("stw_4_v_dis",&stw_4_v_dis);
    T->SetBranchAddress("stw_4_v_energy",&stw_4_v_energy);

    // sig tagger 1
    std::vector<int> *sig_1_v_flag= new std::vector<int>;
    std::vector<double> *sig_1_v_angle = new std::vector<double>;
    std::vector<int> *sig_1_v_flag_single_shower= new std::vector<int>;
    std::vector<double> *sig_1_v_energy= new std::vector<double>;
    std::vector<double> *sig_1_v_energy_1= new std::vector<double>;

    T->SetBranchAddress("sig_1_v_flag",&sig_1_v_flag);
    T->SetBranchAddress("sig_1_v_angle",&sig_1_v_angle);
    T->SetBranchAddress("sig_1_v_flag_single_shower",&sig_1_v_flag_single_shower);
    T->SetBranchAddress("sig_1_v_energy",&sig_1_v_energy);
    T->SetBranchAddress("sig_1_v_energy_1",&sig_1_v_energy_1);

    // sig tagger 2
    std::vector<int> *sig_2_v_flag= new std::vector<int>;
    std::vector<double> *sig_2_v_energy= new std::vector<double>;
    std::vector<double> *sig_2_v_shower_angle= new std::vector<double>;
    std::vector<int> *sig_2_v_flag_single_shower= new std::vector<int>;
    std::vector<double> *sig_2_v_medium_dQ_dx= new std::vector<double>;
    std::vector<double> *sig_2_v_start_dQ_dx= new std::vector<double>;

    T->SetBranchAddress("sig_2_v_flag",&sig_2_v_flag);
    T->SetBranchAddress("sig_2_v_energy",&sig_2_v_energy);
    T->SetBranchAddress("sig_2_v_shower_angle",&sig_2_v_shower_angle);
    T->SetBranchAddress("sig_2_v_flag_single_shower",&sig_2_v_flag_single_shower);
    T->SetBranchAddress("sig_2_v_medium_dQ_dx",&sig_2_v_medium_dQ_dx);
    T->SetBranchAddress("sig_2_v_start_dQ_dx",&sig_2_v_start_dQ_dx);

    // br3 scalar
    int br3_1_flag;
    double br3_1_energy;
    int br3_1_n_shower_segments;
    int br3_1_sg_flag_trajectory;
    double br3_1_sg_direct_length;
    double br3_1_sg_length;
    double br3_1_total_main_length;
    double br3_1_total_length;
    double br3_1_iso_angle;
    int br3_1_sg_flag_topology;

    int br3_2_flag;
    int br3_2_n_ele;
    int br3_2_n_other;
    int br3_2_other_fid;

    int br3_4_flag;
    double br3_4_acc_length;
    double br3_4_total_length;

    int br3_7_flag;
    double br3_7_min_angle;

    int br3_8_flag;
    double br3_8_max_dQ_dx;
    int br3_8_n_main_segs;

    T->SetBranchAddress("br3_1_flag",&br3_1_flag);
    T->SetBranchAddress("br3_1_energy",&br3_1_energy);
    T->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
    T->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
    T->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
    T->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
    T->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
    T->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
    T->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
    T->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);

    T->SetBranchAddress("br3_2_flag",&br3_2_flag);
    T->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
    T->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
    T->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);

    T->SetBranchAddress("br3_4_flag", &br3_4_flag);
    T->SetBranchAddress("br3_4_acc_length", &br3_4_acc_length);
    T->SetBranchAddress("br3_4_total_length", &br3_4_total_length);

    T->SetBranchAddress("br3_7_flag",&br3_7_flag);
    T->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);

    T->SetBranchAddress("br3_8_flag",&br3_8_flag);
    T->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
    T->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);

    // br3 tagger 3
    std::vector<int> *br3_3_v_flag = new std::vector<int>;
    std::vector<double> *br3_3_v_energy = new std::vector<double>;
    std::vector<double> *br3_3_v_angle = new std::vector<double>;
    std::vector<double> *br3_3_v_dir_length = new std::vector<double>;
    std::vector<double> *br3_3_v_length = new std::vector<double>;
   
    T->SetBranchAddress("br3_3_v_flag",&br3_3_v_flag);
    T->SetBranchAddress("br3_3_v_energy",&br3_3_v_energy);
    T->SetBranchAddress("br3_3_v_angle",&br3_3_v_angle);
    T->SetBranchAddress("br3_3_v_dir_length",&br3_3_v_dir_length);
    T->SetBranchAddress("br3_3_v_length",&br3_3_v_length);

    // br3 tagger 5
    std::vector<int> *br3_5_v_flag = new std::vector<int>;
    std::vector<double> *br3_5_v_dir_length = new std::vector<double>;
    std::vector<double> *br3_5_v_total_length = new std::vector<double>;
    std::vector<int> *br3_5_v_flag_avoid_muon_check = new std::vector<int>;
    std::vector<int> *br3_5_v_n_seg = new std::vector<int>;
    std::vector<double> *br3_5_v_angle = new std::vector<double>;
    std::vector<double> *br3_5_v_sg_length = new std::vector<double>;
    std::vector<double> *br3_5_v_energy = new std::vector<double>;
    std::vector<int> *br3_5_v_n_main_segs = new std::vector<int>;
    std::vector<int> *br3_5_v_n_segs = new std::vector<int>;
    std::vector<double> *br3_5_v_shower_main_length = new std::vector<double>;
    std::vector<double> *br3_5_v_shower_total_length = new std::vector<double>;

    T->SetBranchAddress("br3_5_v_flag", &br3_5_v_flag); 
    T->SetBranchAddress("br3_5_v_dir_length", &br3_5_v_dir_length);
    T->SetBranchAddress("br3_5_v_total_length", &br3_5_v_total_length);
    T->SetBranchAddress("br3_5_v_flag_avoid_muon_check", &br3_5_v_flag_avoid_muon_check);
    T->SetBranchAddress("br3_5_v_n_seg", &br3_5_v_n_seg);
    T->SetBranchAddress("br3_5_v_angle", &br3_5_v_angle);
    T->SetBranchAddress("br3_5_v_sg_length", &br3_5_v_sg_length);
    T->SetBranchAddress("br3_5_v_energy", &br3_5_v_energy);
    T->SetBranchAddress("br3_5_v_n_main_segs", &br3_5_v_n_main_segs);
    T->SetBranchAddress("br3_5_v_n_segs", &br3_5_v_n_segs);
    T->SetBranchAddress("br3_5_v_shower_main_length", &br3_5_v_shower_main_length);
    T->SetBranchAddress("br3_5_v_shower_total_length", &br3_5_v_shower_total_length);

    //br3 tagger 6
    std::vector<int> *br3_6_v_flag = new std::vector<int>;
    std::vector<double> *br3_6_v_angle = new std::vector<double>;
    std::vector<double> *br3_6_v_angle1 = new std::vector<double>;
    std::vector<int> *br3_6_v_flag_shower_trajectory = new std::vector<int>;
    std::vector<double> *br3_6_v_direct_length = new std::vector<double>;
    std::vector<double> *br3_6_v_length = new std::vector<double>;
    std::vector<int> *br3_6_v_n_other_vtx_segs = new std::vector<int>;
    std::vector<double> *br3_6_v_energy = new std::vector<double>;

    T->SetBranchAddress("br3_6_v_flag",&br3_6_v_flag);
    T->SetBranchAddress("br3_6_v_angle",&br3_6_v_angle);
    T->SetBranchAddress("br3_6_v_angle1",&br3_6_v_angle1);
    T->SetBranchAddress("br3_6_v_flag_shower_trajectory",&br3_6_v_flag_shower_trajectory);
    T->SetBranchAddress("br3_6_v_direct_length",&br3_6_v_direct_length);
    T->SetBranchAddress("br3_6_v_length",&br3_6_v_length);
    T->SetBranchAddress("br3_6_v_n_other_vtx_segs",&br3_6_v_n_other_vtx_segs);
    T->SetBranchAddress("br3_6_v_energy",&br3_6_v_energy);

    // hol + lol merge
    int hol_1_n_valid_tracks;
    double hol_1_min_angle;
    double hol_1_energy;
    int hol_1_flag_all_shower;
    double hol_1_min_length;
    int hol_1_flag;

    double hol_2_min_angle;
    double hol_2_medium_dQ_dx;
    int hol_2_ncount;
    int hol_2_flag;

    double lol_3_angle_beam;
    int lol_3_n_valid_tracks;
    double lol_3_min_angle;
    int lol_3_vtx_n_segs;
    double lol_3_shower_main_length;
    int lol_3_n_out;
    int lol_3_n_sum;    
    int lol_3_flag;

    T->SetBranchAddress("hol_1_flag", &hol_1_flag);
    T->SetBranchAddress("hol_1_n_valid_tracks", &hol_1_n_valid_tracks);
    T->SetBranchAddress("hol_1_min_angle", &hol_1_min_angle);
    T->SetBranchAddress("hol_1_energy", &hol_1_energy);
    T->SetBranchAddress("hol_1_all_shower", &hol_1_flag_all_shower);
    T->SetBranchAddress("hol_1_min_length", &hol_1_min_length);

    T->SetBranchAddress("hol_2_flag", &hol_2_flag);
    T->SetBranchAddress("hol_2_min_angle", &hol_2_min_angle);
    T->SetBranchAddress("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx);
    T->SetBranchAddress("hol_2_ncount", &hol_2_ncount);
    
    T->SetBranchAddress("lol_3_flag",&lol_3_flag);
    T->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
    T->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
    T->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
    T->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
    T->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
    T->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
    T->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);

    // lol tagger 1
    std::vector<int> *lol_1_v_flag= new std::vector<int>;
    std::vector<double> *lol_1_v_energy = new std::vector<double>;
    std::vector<int> *lol_1_v_vtx_n_segs= new std::vector<int>;
    std::vector<int> *lol_1_v_nseg= new std::vector<int>;
    std::vector<double> *lol_1_v_angle= new std::vector<double>;

    T->SetBranchAddress("lol_1_v_energy",&lol_1_v_energy);
    T->SetBranchAddress("lol_1_v_vtx_n_segs",&lol_1_v_vtx_n_segs);
    T->SetBranchAddress("lol_1_v_nseg",&lol_1_v_nseg);
    T->SetBranchAddress("lol_1_v_angle",&lol_1_v_angle);
    T->SetBranchAddress("lol_1_v_flag",&lol_1_v_flag);

    // lol tagger 2
    std::vector<double> *lol_2_v_length= new std::vector<double>;
    std::vector<double> *lol_2_v_angle= new std::vector<double>;
    std::vector<int> *lol_2_v_type= new std::vector<int>;
    std::vector<int> *lol_2_v_vtx_n_segs= new std::vector<int>;
    std::vector<double> *lol_2_v_energy= new std::vector<double>;
    std::vector<double> *lol_2_v_shower_main_length= new std::vector<double>;
    std::vector<int> *lol_2_v_flag_dir_weak= new std::vector<int>;
    std::vector<int> *lol_2_v_flag = new std::vector<int>;

    T->SetBranchAddress("lol_2_v_flag",&lol_2_v_flag);
    T->SetBranchAddress("lol_2_v_length",&lol_2_v_length);
    T->SetBranchAddress("lol_2_v_angle",&lol_2_v_angle);
    T->SetBranchAddress("lol_2_v_type",&lol_2_v_type);
    T->SetBranchAddress("lol_2_v_vtx_n_segs",&lol_2_v_vtx_n_segs);
    T->SetBranchAddress("lol_2_v_energy",&lol_2_v_energy);
    T->SetBranchAddress("lol_2_v_shower_main_length",&lol_2_v_shower_main_length);
    T->SetBranchAddress("lol_2_v_flag_dir_weak",&lol_2_v_flag_dir_weak);

    // tro tagger 1
    std::vector<int> *tro_1_v_particle_type= new std::vector<int>;
    std::vector<int> *tro_1_v_flag_dir_weak= new std::vector<int>;
    std::vector<double> *tro_1_v_min_dis = new std::vector<double>;
    std::vector<double> *tro_1_v_sg1_length = new std::vector<double>;
    std::vector<double> *tro_1_v_shower_main_length = new std::vector<double>;
    std::vector<int> *tro_1_v_max_n_vtx_segs= new std::vector<int>;
    std::vector<double> *tro_1_v_tmp_length = new std::vector<double>;
    std::vector<double> *tro_1_v_medium_dQ_dx = new std::vector<double>;
    std::vector<double> *tro_1_v_dQ_dx_cut = new std::vector<double>;
    std::vector<int> *tro_1_v_flag_shower_topology= new std::vector<int>;
    std::vector<int> *tro_1_v_flag= new std::vector<int>;

    T->SetBranchAddress("tro_1_v_flag",&tro_1_v_flag);
    T->SetBranchAddress("tro_1_v_particle_type",&tro_1_v_particle_type);
    T->SetBranchAddress("tro_1_v_flag_dir_weak",&tro_1_v_flag_dir_weak);
    T->SetBranchAddress("tro_1_v_min_dis",&tro_1_v_min_dis);
    T->SetBranchAddress("tro_1_v_sg1_length",&tro_1_v_sg1_length);
    T->SetBranchAddress("tro_1_v_shower_main_length",&tro_1_v_shower_main_length);
    T->SetBranchAddress("tro_1_v_max_n_vtx_segs",&tro_1_v_max_n_vtx_segs);
    T->SetBranchAddress("tro_1_v_tmp_length",&tro_1_v_tmp_length);
    T->SetBranchAddress("tro_1_v_medium_dQ_dx",&tro_1_v_medium_dQ_dx);
    T->SetBranchAddress("tro_1_v_dQ_dx_cut",&tro_1_v_dQ_dx_cut);
    T->SetBranchAddress("tro_1_v_flag_shower_topology",&tro_1_v_flag_shower_topology);

    // tro tagger 2
    std::vector<int> *tro_2_v_flag= new std::vector<int>;
    std::vector<double> *tro_2_v_energy = new std::vector<double>;
    std::vector<double> *tro_2_v_stem_length = new std::vector<double>;
    std::vector<double> *tro_2_v_iso_angle = new std::vector<double>;
    std::vector<double> *tro_2_v_max_length = new std::vector<double>;
    std::vector<double> *tro_2_v_angle = new std::vector<double>;

    T->SetBranchAddress("tro_2_v_flag",&tro_2_v_flag);
    T->SetBranchAddress("tro_2_v_energy",&tro_2_v_energy);
    T->SetBranchAddress("tro_2_v_stem_length",&tro_2_v_stem_length);
    T->SetBranchAddress("tro_2_v_iso_angle",&tro_2_v_iso_angle);
    T->SetBranchAddress("tro_2_v_max_length",&tro_2_v_max_length);
    T->SetBranchAddress("tro_2_v_angle",&tro_2_v_angle);

    // tro tagger 4
    std::vector<int> *tro_4_v_flag= new std::vector<int>;
    std::vector<double> *tro_4_v_dir2_mag = new std::vector<double>;
    std::vector<double> *tro_4_v_angle = new std::vector<double>;
    std::vector<double> *tro_4_v_angle1 = new std::vector<double>;
    std::vector<double> *tro_4_v_angle2 = new std::vector<double>;
    std::vector<double> *tro_4_v_length = new std::vector<double>;
    std::vector<double> *tro_4_v_length1 = new std::vector<double>;
    std::vector<double> *tro_4_v_medium_dQ_dx = new std::vector<double>;
    std::vector<double> *tro_4_v_end_dQ_dx = new std::vector<double>;
    std::vector<double> *tro_4_v_energy = new std::vector<double>;
    std::vector<double> *tro_4_v_shower_main_length = new std::vector<double>;
    std::vector<int> *tro_4_v_flag_shower_trajectory= new std::vector<int>;

    T->SetBranchAddress("tro_4_v_flag",&tro_4_v_flag);
    T->SetBranchAddress("tro_4_v_dir2_mag",&tro_4_v_dir2_mag);
    T->SetBranchAddress("tro_4_v_angle",&tro_4_v_angle);
    T->SetBranchAddress("tro_4_v_angle1",&tro_4_v_angle1);
    T->SetBranchAddress("tro_4_v_angle2",&tro_4_v_angle2);
    T->SetBranchAddress("tro_4_v_length",&tro_4_v_length);
    T->SetBranchAddress("tro_4_v_length1",&tro_4_v_length1);
    T->SetBranchAddress("tro_4_v_medium_dQ_dx",&tro_4_v_medium_dQ_dx);
    T->SetBranchAddress("tro_4_v_end_dQ_dx",&tro_4_v_end_dQ_dx);
    T->SetBranchAddress("tro_4_v_energy",&tro_4_v_energy);
    T->SetBranchAddress("tro_4_v_shower_main_length",&tro_4_v_shower_main_length);
    T->SetBranchAddress("tro_4_v_flag_shower_trajectory",&tro_4_v_flag_shower_trajectory);

    // tro tagger 5
    std::vector<int> *tro_5_v_flag = new std::vector<int>;
    std::vector<double> *tro_5_v_max_angle = new std::vector<double>;
    std::vector<double> *tro_5_v_min_angle = new std::vector<double>;
    std::vector<double> *tro_5_v_max_length = new std::vector<double>;
    std::vector<double> *tro_5_v_iso_angle = new std::vector<double>;
    std::vector<int> *tro_5_v_n_vtx_segs= new std::vector<int>;
    std::vector<int> *tro_5_v_min_count= new std::vector<int>;
    std::vector<int> *tro_5_v_max_count= new std::vector<int>;
    std::vector<double> *tro_5_v_energy = new std::vector<double>;

    T->SetBranchAddress("tro_5_v_flag",&tro_5_v_flag);
    T->SetBranchAddress("tro_5_v_max_angle",&tro_5_v_max_angle);
    T->SetBranchAddress("tro_5_v_min_angle",&tro_5_v_min_angle);
    T->SetBranchAddress("tro_5_v_max_length",&tro_5_v_max_length);
    T->SetBranchAddress("tro_5_v_iso_angle",&tro_5_v_iso_angle);
    T->SetBranchAddress("tro_5_v_n_vtx_segs",&tro_5_v_n_vtx_segs);
    T->SetBranchAddress("tro_5_v_min_count",&tro_5_v_min_count);
    T->SetBranchAddress("tro_5_v_max_count",&tro_5_v_max_count);
    T->SetBranchAddress("tro_5_v_energy",&tro_5_v_energy);





    /// BDT initilization
    //BDT_mipid class_mipid("dataset/mipid_type3_BDT800_5.weights.xml");
    BDT_mipid class_mipid("dataset/mipid_BDT.weights.xml");
    BDT_gap class_gap("dataset/gap_BDT.weights.xml");
    BDT_cme_anc class_cme_anc("dataset/cme_anc_BDT.weights.xml");
    BDT_mgo_mgt class_mgo_mgt("dataset/mgo_mgt_BDT.weights.xml");
    BDT_br1 class_br1("dataset/br1_BDT.weights.xml");
    BDT_br3 class_br3("dataset/br3_BDT.weights.xml");
    BDT_br3_3 class_br3_3("dataset/br3_3_BDT.weights.xml");
    BDT_br3_5 class_br3_5("dataset/br3_5_BDT.weights.xml");
    BDT_br3_6 class_br3_6("dataset/br3_6_BDT.weights.xml");
    BDT_hol_lol class_hol_lol("dataset/hol_lol_BDT.weights.xml");
    BDT_stemdir_br2 class_stemdir_br2("dataset/stem_dir_br2_BDT.weights.xml");
    BDT_trimuon class_trimuon("dataset/stl_lem_brm_BDT.weights.xml");
    BDT_br4_tro class_br4_tro("dataset/br4_tro_BDT.weights.xml");
    BDT_mipquality class_mipquality("dataset/mipquality_BDT.weights.xml");
    BDT_pio_1 class_pio_1("dataset/pio_1_BDT.weights.xml");
    BDT_pio_2 class_pio_2("dataset/pio_2_BDT.weights.xml");
    BDT_stw_spt class_stw_spt("dataset/stw_spt_BDT.weights.xml");
    BDT_vis_1 class_vis_1("dataset/vis_1_BDT.weights.xml");
    BDT_vis_2 class_vis_2("dataset/vis_2_BDT.weights.xml");
    BDT_stw_2 class_stw_2("dataset/stw_2_BDT.weights.xml");
    BDT_stw_3 class_stw_3("dataset/stw_3_BDT.weights.xml");
    BDT_stw_4 class_stw_4("dataset/stw_4_BDT.weights.xml");
    BDT_sig_1 class_sig_1("dataset/sig_1_BDT.weights.xml");
    BDT_sig_2 class_sig_2("dataset/sig_2_BDT.weights.xml");
    BDT_lol_1 class_lol_1("dataset/lol_1_BDT.weights.xml");
    BDT_lol_2 class_lol_2("dataset/lol_2_BDT.weights.xml");
    BDT_tro_1 class_tro_1("dataset/tro_1_BDT.weights.xml");
    BDT_tro_2 class_tro_2("dataset/tro_2_BDT.weights.xml");
    BDT_tro_4 class_tro_4("dataset/tro_4_BDT.weights.xml");
    BDT_tro_5 class_tro_5("dataset/tro_5_BDT.weights.xml");




    /// outputfile
    float mipid_score = 0;
    float gap_score = 0;
    float hol_lol_score = 0;
    float cme_anc_score = 0;
    float mgo_mgt_score = 0;
    float br1_score = 0;
    float br3_score = 0;
    float br3_3_score = 0;
    float br3_5_score = 0;
    float br3_6_score = 0;
    float stemdir_br2_score = 0;
    float trimuon_score = 0;
    float br4_tro_score = 0;
    float mipquality_score = 0;
    float pio_1_score = 0;
    float pio_2_score = 0;
    float stw_spt_score = 0;
    float vis_1_score = 0;
    float vis_2_score = 0;
    float stw_2_score = 0;
    float stw_3_score = 0;
    float stw_4_score = 0;
    float sig_1_score = 0;
    float sig_2_score = 0;
    float lol_1_score = 0;
    float lol_2_score = 0;
    float tro_1_score = 0;
    float tro_2_score = 0;
    float tro_4_score = 0;
    float tro_5_score = 0;
    int temp_flag = 1; // current integrated tagger cut-based result
    int rest_flag = 1; // the rest tagger cut-based result

    // additional variable considered in the training
    float mip_energy_f = 0;
    /* float gap_energy_f = 0; */ 
    /* float gap_flag_single_shower_f = 0; */ 
    /* float gap_num_valid_tracks_f = 0; */ 
    /* float gap_n_bad_f = 0; */ 
    /* float mgt_total_other_energy_f = 0; */ 
    /* float mgt_e_indirect_total_energy_f = 0; */ 
    /* float mgt_e_direct_total_energy_f = 0; */  
    float mip_vec_dQ_dx_0_f = 0; 
    float mip_vec_dQ_dx_1_f = 0; 
    float mip_vec_dQ_dx_2_f = 0; 
    float mip_vec_dQ_dx_3_f = 0; 
    float mip_vec_dQ_dx_4_f = 0; 
    /* float mip_vec_dQ_dx_5_f = 0; */ 
    /* float mip_vec_dQ_dx_6_f = 0; */ 
    float spt_angle_beam_f = 0; 
    float spt_angle_drift_f = 0; 
    float cme_mu_energy_f = 0; 
    float cme_mu_length_f = 0; 
    float cme_length_f = 0; 

    TFile* outputfile = new TFile("bdtscore_"+type+".root","RECREATE");
    TTree* sig = new TTree(type, type);
    sig->Branch("run",&run, "run/I");
    sig->Branch("subrun",&subrun, "subrun/I");
    sig->Branch("event",&event, "event/I");
    sig->Branch("nueTag",&nueTag, "nueTag/I");
    sig->Branch("recoFC",&recoFC, "recoFC/I");
    if(_MC_){
    sig->Branch("truth_nue",&truth_nue, "truth_nue/I");
    sig->Branch("truth_CC",&truth_CC, "truth_CC/I");
    sig->Branch("truth_inFV",&truth_inFV, "truth_inFV/I");
    sig->Branch("truth_cosmic",&truth_cosmic, "truth_cosmic/I");
    sig->Branch("trueEnu",&trueEnu, "trueEnu/F");
    sig->Branch("weight",&weight, "weight/F");
    sig->Branch("lowEweight",&lowEweight, "lowEweight/F");
    sig->Branch("trueEdep",&trueEdep, "trueEdep/F");
    sig->Branch("nuvtx_diff",&nuvtx_diff, "nuvtx_diff/F");
    sig->Branch("showervtx_diff",&showervtx_diff, "showervtx_diff/F");
    }
    sig->Branch("mip_energy", &mip_energy_f, "mip_energy/F");
    /* sig->Branch("gap_energy", &gap_energy_f, "gap_energy/F"); */ 
    /* sig->Branch("gap_flag_single_shower", &gap_flag_single_shower_f, "gap_flag_single_shower/F"); */ 
    /* sig->Branch("gap_num_valid_tracks", &gap_num_valid_tracks_f, "gap_num_valid_tracks/F"); */ 
    /* sig->Branch("gap_n_bad", &gap_n_bad_f, "gap_n_bad/F"); */ 
    /* sig->Branch("mgt_total_other_energy", &mgt_total_other_energy_f, "mgt_total_other_energy/F"); */ 
    /* sig->Branch("mgt_e_indirect_total_energy", &mgt_e_indirect_total_energy_f, "mgt_e_indirect_total_energy/F"); */ 
    /* sig->Branch("mgt_e_direct_total_energy", &mgt_e_direct_total_energy_f, "mgt_e_direct_total_energy/F"); */  
    sig->Branch("mip_vec_dQ_dx_0", &mip_vec_dQ_dx_0_f, "mip_vec_dQ_dx_0/F"); 
    sig->Branch("mip_vec_dQ_dx_1", &mip_vec_dQ_dx_1_f, "mip_vec_dQ_dx_1/F"); 
    sig->Branch("mip_vec_dQ_dx_2", &mip_vec_dQ_dx_2_f, "mip_vec_dQ_dx_2/F"); 
    sig->Branch("mip_vec_dQ_dx_3", &mip_vec_dQ_dx_3_f, "mip_vec_dQ_dx_3/F"); 
    sig->Branch("mip_vec_dQ_dx_4", &mip_vec_dQ_dx_4_f, "mip_vec_dQ_dx_4/F"); 
    /* sig->Branch("mip_vec_dQ_dx_5", &mip_vec_dQ_dx_5_f, "mip_vec_dQ_dx_5/F"); */ 
    /* sig->Branch("mip_vec_dQ_dx_6", &mip_vec_dQ_dx_6_f, "mip_vec_dQ_dx_6/F"); */ 
    sig->Branch("spt_angle_beam", &spt_angle_beam_f, "spt_angle_beam/F"); 
    sig->Branch("spt_angle_drift", &spt_angle_drift_f, "spt_angle_drift/F"); 
    sig->Branch("cme_mu_energy", &cme_mu_energy_f, "cme_mu_energy/F"); 
    sig->Branch("cme_mu_length", &cme_mu_length_f, "cme_mu_length/F"); 
    sig->Branch("cme_length", &cme_length_f, "cme_length/F"); 


    sig->Branch("mipid_score",&mipid_score, "mipid_score/F");
    sig->Branch("gap_score",&gap_score, "gap_score/F");
    sig->Branch("hol_lol_score",&hol_lol_score, "hol_lol_score/F");
    sig->Branch("cme_anc_score", &cme_anc_score, "cme_anc_score/F");
    sig->Branch("mgo_mgt_score", &mgo_mgt_score, "mgo_mgt_score/F");
    sig->Branch("br1_score",&br1_score, "br1_score/F");
    sig->Branch("br3_score",&br3_score, "br3_score/F");
    sig->Branch("br3_3_score",&br3_3_score, "br3_3_score/F");
    sig->Branch("br3_5_score",&br3_5_score, "br3_5_score/F");
    sig->Branch("br3_6_score",&br3_6_score, "br3_6_score/F");
    sig->Branch("stemdir_br2_score",&stemdir_br2_score, "stemdir_br2_score/F");
    sig->Branch("trimuon_score",&trimuon_score, "trimuon_score/F");
    sig->Branch("br4_tro_score",&br4_tro_score, "br4_tro_score/F");
    sig->Branch("mipquality_score",&mipquality_score, "mipquality_score/F");
    sig->Branch("pio_1_score",&pio_1_score, "pio_1_score/F");
    sig->Branch("pio_2_score",&pio_2_score, "pio_2_score/F");
    sig->Branch("stw_spt_score",&stw_spt_score, "stw_spt_score/F");
    sig->Branch("vis_1_score",&vis_1_score, "vis_1_score/F");
    sig->Branch("vis_2_score",&vis_2_score, "vis_2_score/F");
    sig->Branch("stw_2_score",&stw_2_score, "stw_2_score/F");
    sig->Branch("stw_3_score",&stw_3_score, "stw_3_score/F");
    sig->Branch("stw_4_score",&stw_4_score, "stw_4_score/F");
    sig->Branch("sig_1_score",&sig_1_score, "sig_1_score/F");
    sig->Branch("sig_2_score",&sig_2_score, "sig_2_score/F");
    sig->Branch("lol_1_score",&lol_1_score, "lol_1_score/F");
    sig->Branch("lol_2_score",&lol_2_score, "lol_2_score/F");
    sig->Branch("tro_1_score",&tro_1_score, "tro_1_score/F");
    sig->Branch("tro_2_score",&tro_2_score, "tro_2_score/F");
    sig->Branch("tro_4_score",&tro_4_score, "tro_4_score/F");
    sig->Branch("tro_5_score",&tro_5_score, "tro_5_score/F");

    sig->Branch("temp_flag",&temp_flag, "temp_flag/I");
    sig->Branch("rest_flag",&rest_flag, "rest_flag/I");

    float default_scale = 1.0; //BDT default value scaling factor

    bool mipid_on = true;
    bool gap_on = true;
    bool cme_anc_on = true;
    bool mgo_mgt_on = true;
    bool br1_on = true;
    bool br3_on = true;
    bool br3_3_on = true; 
    bool br3_5_on = true; 
    bool br3_6_on = true; 
    bool hol_lol_on = true; 
    bool stemdir_br2_on = true;
    bool trimuon_on = true;
    bool br4_tro_on = true;
    bool mipquality_on = true;
    bool pio_1_on = true;
    bool pio_2_on = true;
    bool stw_spt_on = true;
    bool vis_1_on = true;
    bool vis_2_on = true;
    bool stw_2_on = true; 
    bool stw_3_on = true; 
    bool stw_4_on = true; 
    bool sig_1_on = true; 
    bool sig_2_on = true; 
    bool lol_1_on = true; 
    bool lol_2_on = true; 
    bool tro_1_on = true; 
    bool tro_2_on = true; 
    bool tro_4_on = true; 
    bool tro_5_on = true; 


    //// Loop over each event
    for(int i=0; i<T->GetEntries(); i++){

    T->GetEntry(i);

    /// read variables and calculate BDT score
    /// data type conversion within class if necessary

    // mipid
    if(mipid_on){
    class_mipid.reset();

    class_mipid.mip_energy=mip_energy;
    class_mipid.mip_n_end_reduction=mip_n_end_reduction;    
    class_mipid.mip_n_first_mip=mip_n_first_mip;
    class_mipid.mip_n_first_non_mip=mip_n_first_non_mip;
    class_mipid.mip_n_first_non_mip_1=mip_n_first_non_mip_1;
    class_mipid.mip_n_first_non_mip_2=mip_n_first_non_mip_2;
    class_mipid.mip_vec_dQ_dx_0=mip_vec_dQ_dx_0;
    class_mipid.mip_vec_dQ_dx_1=mip_vec_dQ_dx_1;
    class_mipid.mip_max_dQ_dx_sample=mip_max_dQ_dx_sample;
    class_mipid.mip_n_below_threshold=mip_n_below_threshold;
    class_mipid.mip_n_below_zero=mip_n_below_zero;
    class_mipid.mip_n_lowest=mip_n_lowest;
    class_mipid.mip_n_highest=mip_n_highest;
    class_mipid.mip_lowest_dQ_dx=mip_lowest_dQ_dx;
    class_mipid.mip_highest_dQ_dx=mip_highest_dQ_dx;
    class_mipid.mip_medium_dQ_dx=mip_medium_dQ_dx;
    class_mipid.mip_stem_length=mip_stem_length;
    class_mipid.mip_length_main=mip_length_main;
    class_mipid.mip_length_total=mip_length_total;
    class_mipid.mip_angle_beam=mip_angle_beam;
    class_mipid.mip_iso_angle=mip_iso_angle;
    class_mipid.mip_n_vertex=mip_n_vertex;
    class_mipid.mip_n_good_tracks=mip_n_good_tracks;
    class_mipid.mip_E_indirect_max_energy=mip_E_indirect_max_energy;
    class_mipid.mip_flag_all_above=mip_flag_all_above;
    class_mipid.mip_min_dQ_dx_5=mip_min_dQ_dx_5;
    class_mipid.mip_n_other_vertex=mip_n_other_vertex; 
    class_mipid.mip_n_stem_size=mip_n_stem_size;
    class_mipid.mip_flag_stem_trajectory=mip_flag_stem_trajectory;
    if(mip_min_dis>1000) mip_min_dis = 1000.0;
    class_mipid.mip_min_dis=mip_min_dis;
    
    mipid_score = class_mipid.evaluate();
    if(mip_filled==0) mipid_score = 0.26/default_scale; // average signal score
    //if(i<20) cout << "BDT score: "<<mipid_score<<endl;
    }

    //gap
    if(gap_on){
    class_gap.reset();

    class_gap.gap_flag_prolong_u=gap_flag_prolong_u;
    class_gap.gap_flag_prolong_v=gap_flag_prolong_v;
    class_gap.gap_flag_prolong_w=gap_flag_prolong_w;
    class_gap.gap_flag_parallel=gap_flag_parallel;
    class_gap.gap_n_points=gap_n_points;
    class_gap.gap_n_bad=gap_n_bad;
    class_gap.gap_energy=gap_energy;
    class_gap.gap_num_valid_tracks=gap_num_valid_tracks;
    class_gap.gap_flag_single_shower=gap_flag_single_shower;

    gap_score = class_gap.evaluate();
    if(gap_filled==0) gap_score = 0.32/default_scale; // average signal score
    }

    // cme_anc
    if(cme_anc_on){
    class_cme_anc.reset();

    class_cme_anc.cme_mu_energy=cme_mu_energy;
    class_cme_anc.cme_energy=cme_energy;
    class_cme_anc.cme_mu_length=cme_mu_length;
    class_cme_anc.cme_length=cme_length;
    class_cme_anc.cme_angle_beam=cme_angle_beam;
    class_cme_anc.anc_angle=anc_angle;
    class_cme_anc.anc_max_angle=anc_max_angle;
    class_cme_anc.anc_max_length=anc_max_length;
    class_cme_anc.anc_acc_forward_length=anc_acc_forward_length;
    class_cme_anc.anc_acc_backward_length=anc_acc_backward_length;
    class_cme_anc.anc_acc_forward_length1=anc_acc_forward_length1;
    class_cme_anc.anc_shower_main_length=anc_shower_main_length;
    class_cme_anc.anc_shower_total_length=anc_shower_total_length;
    class_cme_anc.anc_flag_main_outside=anc_flag_main_outside;
   
    cme_anc_score = class_cme_anc.evaluate();
    }

    // mgo_mgt
    if(mgo_mgt_on){
    class_mgo_mgt.reset();

    class_mgo_mgt.mgo_energy=mgo_energy;
    class_mgo_mgt.mgo_max_energy=mgo_max_energy;
    class_mgo_mgt.mgo_total_energy=mgo_total_energy;
    class_mgo_mgt.mgo_n_showers=mgo_n_showers;
    class_mgo_mgt.mgo_max_energy_1=mgo_max_energy_1;
    class_mgo_mgt.mgo_max_energy_2=mgo_max_energy_2;
    class_mgo_mgt.mgo_total_other_energy=mgo_total_other_energy;
    class_mgo_mgt.mgo_n_total_showers=mgo_n_total_showers;
    class_mgo_mgt.mgo_total_other_energy_1=mgo_total_other_energy_1;
    class_mgo_mgt.mgt_flag_single_shower=mgt_flag_single_shower;
    class_mgo_mgt.mgt_max_energy=mgt_max_energy;
    class_mgo_mgt.mgt_total_other_energy=mgt_total_other_energy;
    class_mgo_mgt.mgt_max_energy_1=mgt_max_energy_1;
    class_mgo_mgt.mgt_e_indirect_max_energy=mgt_e_indirect_max_energy;
    class_mgo_mgt.mgt_e_direct_max_energy=mgt_e_direct_max_energy;
    class_mgo_mgt.mgt_n_direct_showers=mgt_n_direct_showers;
    class_mgo_mgt.mgt_e_direct_total_energy=mgt_e_direct_total_energy;
    class_mgo_mgt.mgt_flag_indirect_max_pio=mgt_flag_indirect_max_pio;
    class_mgo_mgt.mgt_e_indirect_total_energy=mgt_e_indirect_total_energy;

    mgo_mgt_score = class_mgo_mgt.evaluate();
    }

    // br1
    if(br1_on){
    class_br1.reset();

    class_br1.br1_1_shower_type=br1_1_shower_type;
    class_br1.br1_1_vtx_n_segs=br1_1_vtx_n_segs;
    class_br1.br1_1_energy=br1_1_energy;
    class_br1.br1_1_n_segs=br1_1_n_segs;
    class_br1.br1_1_flag_sg_topology=br1_1_flag_sg_topology;
    class_br1.br1_1_flag_sg_trajectory=br1_1_flag_sg_trajectory;
    class_br1.br1_1_sg_length=br1_1_sg_length;
    class_br1.br1_2_n_connected=br1_2_n_connected;
    class_br1.br1_2_max_length=br1_2_max_length;
    class_br1.br1_2_n_connected_1=br1_2_n_connected_1;
    class_br1.br1_2_n_shower_segs=br1_2_n_shower_segs;
    class_br1.br1_2_max_length_ratio=br1_2_max_length_ratio;
    class_br1.br1_2_shower_length=br1_2_shower_length;
    class_br1.br1_3_n_connected_p=br1_3_n_connected_p;
    class_br1.br1_3_max_length_p=br1_3_max_length_p;
    class_br1.br1_3_n_shower_main_segs=br1_3_n_shower_main_segs;

    br1_score = class_br1.evaluate();
    }

    // br3
    if(br3_on){
    class_br3.reset();

    class_br3.br3_1_energy=br3_1_energy;
    class_br3.br3_1_n_shower_segments=br3_1_n_shower_segments;
    class_br3.br3_1_sg_flag_trajectory=br3_1_sg_flag_trajectory;
    class_br3.br3_1_sg_direct_length=br3_1_sg_direct_length;
    class_br3.br3_1_sg_length=br3_1_sg_length;
    class_br3.br3_1_total_main_length=br3_1_total_main_length;
    class_br3.br3_1_total_length=br3_1_total_length;
    class_br3.br3_1_iso_angle=br3_1_iso_angle;
    class_br3.br3_1_sg_flag_topology=br3_1_sg_flag_topology;
    class_br3.br3_2_n_ele=br3_2_n_ele;
    class_br3.br3_2_n_other=br3_2_n_other;
    class_br3.br3_2_other_fid=br3_2_other_fid;
    class_br3.br3_4_acc_length=br3_4_acc_length;
    class_br3.br3_4_total_length=br3_4_total_length;
    class_br3.br3_7_min_angle=br3_7_min_angle;
    class_br3.br3_8_max_dQ_dx=br3_8_max_dQ_dx;
    class_br3.br3_8_n_main_segs=br3_8_n_main_segs;

    br3_score = class_br3.evaluate();
    }





    // stemdir_br2
    if(stemdir_br2_on){
    class_stemdir_br2.reset();

    class_stemdir_br2.stem_dir_flag_single_shower=stem_dir_flag_single_shower;
    class_stemdir_br2.stem_dir_angle=stem_dir_angle;
    class_stemdir_br2.stem_dir_energy=stem_dir_energy;
    class_stemdir_br2.stem_dir_angle1=stem_dir_angle1;
    class_stemdir_br2.stem_dir_angle2=stem_dir_angle2;
    class_stemdir_br2.stem_dir_angle3=stem_dir_angle3;
    class_stemdir_br2.stem_dir_ratio=stem_dir_ratio;
    class_stemdir_br2.br2_num_valid_tracks=br2_num_valid_tracks;
    class_stemdir_br2.br2_n_shower_main_segs=br2_n_shower_main_segs;
    class_stemdir_br2.br2_max_angle=br2_max_angle;
    class_stemdir_br2.br2_sg_length=br2_sg_length;
    class_stemdir_br2.br2_flag_sg_trajectory=br2_flag_sg_trajectory;

    stemdir_br2_score = class_stemdir_br2.evaluate();
    }


    // three muon related tagger stem length, lem, brm
    if(trimuon_on){
    class_trimuon.reset();
    
    class_trimuon.stem_len_energy=stem_len_energy;
    class_trimuon.stem_len_length=stem_len_length;
    class_trimuon.stem_len_flag_avoid_muon_check=stem_len_flag_avoid_muon_check;
    class_trimuon.stem_len_num_daughters=stem_len_num_daughters;
    class_trimuon.stem_len_daughter_length=stem_len_daughter_length;
    class_trimuon.brm_n_mu_segs=brm_n_mu_segs;
    class_trimuon.brm_Ep=brm_Ep;
    class_trimuon.brm_acc_length=brm_acc_length;
    class_trimuon.brm_shower_total_length=brm_shower_total_length;
    class_trimuon.brm_connected_length=brm_connected_length;
    class_trimuon.brm_n_size=brm_n_size;
    class_trimuon.brm_acc_direct_length=brm_acc_direct_length;
    class_trimuon.brm_n_shower_main_segs=brm_n_shower_main_segs;
    class_trimuon.brm_n_mu_main=brm_n_mu_main;
    class_trimuon.lem_shower_main_length=lem_shower_main_length;
    class_trimuon.lem_n_3seg=lem_n_3seg;
    class_trimuon.lem_e_charge=lem_e_charge;
    class_trimuon.lem_e_dQdx=lem_e_dQdx;
    class_trimuon.lem_shower_num_main_segs=lem_shower_num_main_segs;

    trimuon_score = class_trimuon.evaluate();
    }

    // br4 + tro
    if(br4_tro_on){
    class_br4_tro.reset();

    class_br4_tro.br4_1_shower_main_length=br4_1_shower_main_length;
    class_br4_tro.br4_1_shower_total_length=br4_1_shower_total_length;
    class_br4_tro.br4_1_min_dis=br4_1_min_dis;
    class_br4_tro.br4_1_energy=br4_1_energy;
    class_br4_tro.br4_1_flag_avoid_muon_check=br4_1_flag_avoid_muon_check;
    class_br4_tro.br4_1_n_vtx_segs=br4_1_n_vtx_segs;
    class_br4_tro.br4_1_n_main_segs=br4_1_n_main_segs;
    class_br4_tro.br4_2_ratio_45=br4_2_ratio_45;
    class_br4_tro.br4_2_ratio_35=br4_2_ratio_35;
    class_br4_tro.br4_2_ratio_25=br4_2_ratio_25;
    class_br4_tro.br4_2_ratio_15=br4_2_ratio_15;
    class_br4_tro.br4_2_ratio1_45=br4_2_ratio1_45;
    class_br4_tro.br4_2_ratio1_35=br4_2_ratio1_35;
    class_br4_tro.br4_2_ratio1_25=br4_2_ratio1_25;
    class_br4_tro.br4_2_ratio1_15=br4_2_ratio1_15;
    class_br4_tro.br4_2_iso_angle=br4_2_iso_angle;
    class_br4_tro.br4_2_iso_angle1=br4_2_iso_angle1;
    class_br4_tro.br4_2_angle=br4_2_angle;
    class_br4_tro.tro_3_stem_length=tro_3_stem_length;
    class_br4_tro.tro_3_n_muon_segs=tro_3_n_muon_segs;

    br4_tro_score = class_br4_tro.evaluate();
    }

    // mip quality
    if(mipquality_on){
    class_mipquality.reset();
    
    class_mipquality.mip_quality_energy=mip_quality_energy;
    class_mipquality.mip_quality_overlap=mip_quality_overlap;
    class_mipquality.mip_quality_n_showers=mip_quality_n_showers;
    class_mipquality.mip_quality_n_tracks=mip_quality_n_tracks;
    class_mipquality.mip_quality_flag_inside_pi0=mip_quality_flag_inside_pi0;
    class_mipquality.mip_quality_n_pi0_showers=mip_quality_n_pi0_showers;
    if(mip_quality_shortest_length>1000) mip_quality_shortest_length = 1000;
    class_mipquality.mip_quality_shortest_length=mip_quality_shortest_length;
    class_mipquality.mip_quality_acc_length=mip_quality_acc_length;
    if(std::isnan(mip_quality_shortest_angle)) mip_quality_shortest_angle = 0; 
    class_mipquality.mip_quality_shortest_angle=mip_quality_shortest_angle;
    class_mipquality.mip_quality_flag_proton=mip_quality_flag_proton;
    
    mipquality_score = class_mipquality.evaluate();
    }

    // pio tagger1
    if(pio_1_on){
    class_pio_1.reset();
    
    class_pio_1.pio_mip_id=pio_mip_id;
    class_pio_1.pio_1_mass=pio_1_mass;
    class_pio_1.pio_1_pio_type=pio_1_pio_type;
    class_pio_1.pio_1_energy_1=pio_1_energy_1;
    class_pio_1.pio_1_energy_2=pio_1_energy_2;
    class_pio_1.pio_1_dis_1=pio_1_dis_1;
    class_pio_1.pio_1_dis_2=pio_1_dis_2; 

    pio_1_score = class_pio_1.evaluate();
    if(pio_flag_pio==0 || pio_filled==0) pio_1_score = 0.1/default_scale; //0.1 
    }

    // stw + spt
    if(stw_spt_on){
    class_stw_spt.reset();

    class_stw_spt.stw_1_energy=stw_1_energy;
    class_stw_spt.stw_1_dis=stw_1_dis;
    class_stw_spt.stw_1_dQ_dx=stw_1_dQ_dx;
    class_stw_spt.stw_1_flag_single_shower=stw_1_flag_single_shower;
    class_stw_spt.stw_1_n_pi0=stw_1_n_pi0;
    class_stw_spt.stw_1_num_valid_tracks=stw_1_num_valid_tracks;
    class_stw_spt.spt_shower_main_length=spt_shower_main_length;
    class_stw_spt.spt_shower_total_length=spt_shower_total_length;
    class_stw_spt.spt_angle_beam=spt_angle_beam;
    class_stw_spt.spt_angle_vertical=spt_angle_vertical;
    class_stw_spt.spt_max_dQ_dx=spt_max_dQ_dx;
    class_stw_spt.spt_angle_beam_1=spt_angle_beam_1;
    class_stw_spt.spt_angle_drift=spt_angle_drift;
    class_stw_spt.spt_angle_drift_1=spt_angle_drift_1;
    class_stw_spt.spt_num_valid_tracks=spt_num_valid_tracks;
    class_stw_spt.spt_n_vtx_segs=spt_n_vtx_segs;
    class_stw_spt.spt_max_length=spt_max_length;

    stw_spt_score = class_stw_spt.evaluate();

    }

    // vis tagger 1
    if(vis_1_on){
    class_vis_1.reset();

    class_vis_1.vis_1_n_vtx_segs=vis_1_n_vtx_segs;
    class_vis_1.vis_1_energy=vis_1_energy;
    class_vis_1.vis_1_num_good_tracks=vis_1_num_good_tracks;
    class_vis_1.vis_1_max_angle=vis_1_max_angle;
    class_vis_1.vis_1_max_shower_angle=vis_1_max_shower_angle;
    class_vis_1.vis_1_tmp_length1=vis_1_tmp_length1;
    class_vis_1.vis_1_tmp_length2=vis_1_tmp_length2;

    if(vis_1_filled == 0) vis_1_score = 0.6/default_scale; //0.6
    else vis_1_score = class_vis_1.evaluate();

    }

    // vis tagger 2
    if(vis_2_on){
        class_vis_2.reset();

        class_vis_2.vis_2_n_vtx_segs=vis_2_n_vtx_segs;
        class_vis_2.vis_2_min_angle=vis_2_min_angle;
        class_vis_2.vis_2_min_weak_track=vis_2_min_weak_track;
        class_vis_2.vis_2_angle_beam=vis_2_angle_beam;
        class_vis_2.vis_2_min_angle1=vis_2_min_angle1;
        class_vis_2.vis_2_iso_angle1=vis_2_iso_angle1;
        class_vis_2.vis_2_min_medium_dQ_dx=vis_2_min_medium_dQ_dx;
        class_vis_2.vis_2_min_length=vis_2_min_length;
        class_vis_2.vis_2_sg_length=vis_2_sg_length;
        class_vis_2.vis_2_max_angle=vis_2_max_angle;
        class_vis_2.vis_2_max_weak_track=vis_2_max_weak_track;
    
        if(vis_2_filled == 0) vis_2_score = 0.52/default_scale; //0.52
        else vis_2_score = class_vis_2.evaluate();

    }

    // hol + lol 
    if(hol_lol_on){
        class_hol_lol.reset();

        class_hol_lol.hol_1_n_valid_tracks=hol_1_n_valid_tracks;
        class_hol_lol.hol_1_min_angle=hol_1_min_angle;
        class_hol_lol.hol_1_energy=hol_1_energy;
        class_hol_lol.hol_1_flag_all_shower=hol_1_flag_all_shower;
        class_hol_lol.hol_1_min_length=hol_1_min_length;
        class_hol_lol.hol_2_min_angle=hol_2_min_angle;
        class_hol_lol.hol_2_medium_dQ_dx=hol_2_medium_dQ_dx;
        class_hol_lol.hol_2_ncount=hol_2_ncount;
        class_hol_lol.lol_3_angle_beam=lol_3_angle_beam;
        class_hol_lol.lol_3_n_valid_tracks=lol_3_n_valid_tracks;
        class_hol_lol.lol_3_min_angle=lol_3_min_angle;
        class_hol_lol.lol_3_vtx_n_segs=lol_3_vtx_n_segs;
        class_hol_lol.lol_3_shower_main_length=lol_3_shower_main_length;
        class_hol_lol.lol_3_n_out=lol_3_n_out;
        class_hol_lol.lol_3_n_sum=lol_3_n_sum;    

        hol_lol_score = class_hol_lol.evaluate();
        /* if(hol_lol_score<-0.1){ */
        /*     std::cout */
        /*     <<run<<" " */
        /*     <<subrun<<" " */
        /*     <<event<<" " */
        /*     <<hol_1_n_valid_tracks<<" " */
        /*     <<hol_1_min_angle<<" " */
        /*     <<hol_1_energy<<" " */
        /*     <<hol_1_flag_all_shower<<" " */
        /*     <<hol_1_min_length<<" " */
        /*     <<hol_2_min_angle<<" " */
        /*     <<hol_2_medium_dQ_dx<<" " */
        /*     <<hol_2_ncount<<" " */
        /*     <<lol_3_angle_beam<<" " */
        /*     <<lol_3_n_valid_tracks<<" " */
        /*     <<lol_3_min_angle<<" " */
        /*     <<lol_3_vtx_n_segs<<" " */
        /*     <<lol_3_shower_main_length<<" " */
        /*     <<lol_3_n_out<<" " */
        /*     <<lol_3_n_sum<<" " */   
        /*     <<std::endl; */
        /* } */
    }
    






    ///////// FOR vector input 
    /// vector: lowest BDT score and multiplication of flags
    std::vector<float> *bbddtt = new std::vector<float>;  
    
    // pio tagger2
    int pio_2_flag_temp = 1;
    if(pio_2_on){
    bbddtt->clear();
    for(size_t j=0; j<pio_2_v_flag->size(); j++)
    {
        class_pio_2.reset();

        class_pio_2.pio_2_v_dis2=pio_2_v_dis2->at(j);
        class_pio_2.pio_2_v_angle2=pio_2_v_angle2->at(j);
        class_pio_2.pio_2_v_acc_length=pio_2_v_acc_length->at(j);
        class_pio_2.pio_mip_id=pio_mip_id;

        float temp_score = class_pio_2.evaluate();
        bbddtt->push_back(temp_score);
        pio_2_flag_temp *= pio_2_v_flag->at(j);
    }

    if(pio_flag_pio==1 || pio_filled==0 || pio_2_v_flag->size()==0) pio_2_score = 0.2/default_scale; //0.2 
    else pio_2_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // stw tagger 2
    int stw_2_flag_temp = 1;  
    if(stw_2_on){
    bbddtt->clear();
    for(size_t j=0; j<stw_2_v_flag->size(); j++)
    {
        class_stw_2.reset();

        class_stw_2.stw_2_v_medium_dQ_dx=stw_2_v_medium_dQ_dx->at(j);
        class_stw_2.stw_2_v_energy=stw_2_v_energy->at(j);
        class_stw_2.stw_2_v_angle=stw_2_v_angle->at(j);
        class_stw_2.stw_2_v_dir_length=stw_2_v_dir_length->at(j);
        class_stw_2.stw_2_v_max_dQ_dx=stw_2_v_max_dQ_dx->at(j);
    
        float temp_score = class_stw_2.evaluate();
        bbddtt->push_back(temp_score);
        stw_2_flag_temp *= stw_2_v_flag->at(j);
    }

    if(stw_2_v_flag->size()==0) stw_2_score = 0.7/default_scale; //0.7
    else stw_2_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // stw tagger 3
    int stw_3_flag_temp = 1;  
    if(stw_3_on){
    bbddtt->clear();
    for(size_t j=0; j<stw_3_v_flag->size(); j++)
    {
        class_stw_3.reset();

        class_stw_3.stw_3_v_angle=stw_3_v_angle->at(j);
        class_stw_3.stw_3_v_dir_length=stw_3_v_dir_length->at(j);
        class_stw_3.stw_3_v_energy=stw_3_v_energy->at(j);
        class_stw_3.stw_3_v_medium_dQ_dx=stw_3_v_medium_dQ_dx->at(j);
    
        float temp_score = class_stw_3.evaluate();
        bbddtt->push_back(temp_score);
        stw_3_flag_temp *= stw_3_v_flag->at(j);
    }

    if(stw_3_v_flag->size()==0) stw_3_score = 0.5/default_scale; //0.5
    else stw_3_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // stw tagger 4 
    int stw_4_flag_temp = 1;  
    if(stw_4_on){
    bbddtt->clear();
    for(size_t j=0; j<stw_4_v_flag->size(); j++)
    {
        class_stw_4.reset();

        class_stw_4.stw_4_v_angle=stw_4_v_angle->at(j);
        class_stw_4.stw_4_v_dis=stw_4_v_dis->at(j);
        class_stw_4.stw_4_v_energy=stw_4_v_energy->at(j);
    
        float temp_score = class_stw_4.evaluate();
        bbddtt->push_back(temp_score);
        stw_4_flag_temp *= stw_4_v_flag->at(j);
    }

    if(stw_4_v_flag->size()==0) stw_4_score = 0.7/default_scale; //0.7
    else stw_4_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }
    
    // sig tagger 1 
    int sig_1_flag_temp = 1;  
    if(sig_1_on){
    bbddtt->clear();
    for(size_t j=0; j<sig_1_v_flag->size(); j++)
    {
        class_sig_1.reset();

        class_sig_1.sig_1_v_angle=sig_1_v_angle->at(j);
        class_sig_1.sig_1_v_flag_single_shower=sig_1_v_flag_single_shower->at(j);
        class_sig_1.sig_1_v_energy=sig_1_v_energy->at(j);
        class_sig_1.sig_1_v_energy_1=sig_1_v_energy_1->at(j);
    
        float temp_score = class_sig_1.evaluate();
        bbddtt->push_back(temp_score);
        sig_1_flag_temp *= sig_1_v_flag->at(j);
    }

    if(sig_1_v_flag->size()==0) sig_1_score = 0.59/default_scale; //0.59
    else sig_1_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // sig tagger 2 
    int sig_2_flag_temp = 1;  
    if(sig_2_on){
    bbddtt->clear();
    for(size_t j=0; j<sig_2_v_flag->size(); j++)
    {
        class_sig_2.reset();
        
        class_sig_2.sig_2_v_energy = sig_2_v_energy->at(j);
        class_sig_2.sig_2_v_shower_angle = sig_2_v_shower_angle->at(j);
        class_sig_2.sig_2_v_flag_single_shower = sig_2_v_flag_single_shower->at(j);
        class_sig_2.sig_2_v_medium_dQ_dx = sig_2_v_medium_dQ_dx->at(j);
        class_sig_2.sig_2_v_start_dQ_dx = sig_2_v_start_dQ_dx->at(j);
    
        float temp_score = class_sig_2.evaluate();
        bbddtt->push_back(temp_score);
        sig_2_flag_temp *= sig_2_v_flag->at(j);
    }

    if(sig_2_v_flag->size()==0) sig_2_score = 0.55/default_scale;
    else sig_2_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }
    
    // br3 tagger 3 
    int br3_3_flag_temp = 1;  
    if(br3_3_on){
    bbddtt->clear();
    for(size_t j=0; j<br3_3_v_flag->size(); j++)
    {
        class_br3_3.reset();
        
        class_br3_3.br3_3_v_energy=br3_3_v_energy->at(j);
        class_br3_3.br3_3_v_angle=br3_3_v_angle->at(j);
        class_br3_3.br3_3_v_dir_length=br3_3_v_dir_length->at(j);
        class_br3_3.br3_3_v_length=br3_3_v_length->at(j);

        float temp_score = class_br3_3.evaluate();
        bbddtt->push_back(temp_score);
        br3_3_flag_temp *= br3_3_v_flag->at(j);
    }

    if(br3_3_v_flag->size()==0) br3_3_score = 0.3/default_scale;
    else br3_3_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // br3 tagger 5 
    int br3_5_flag_temp = 1;  
    if(br3_5_on){
    bbddtt->clear();
    for(size_t j=0; j<br3_5_v_flag->size(); j++)
    {
        class_br3_5.reset();
        
        class_br3_5.br3_5_v_dir_length=br3_5_v_dir_length->at(j);
        class_br3_5.br3_5_v_total_length=br3_5_v_total_length->at(j);
        class_br3_5.br3_5_v_flag_avoid_muon_check=br3_5_v_flag_avoid_muon_check->at(j);
        class_br3_5.br3_5_v_n_seg=br3_5_v_n_seg->at(j);
        class_br3_5.br3_5_v_angle=br3_5_v_angle->at(j);
        class_br3_5.br3_5_v_sg_length=br3_5_v_sg_length->at(j);
        class_br3_5.br3_5_v_energy=br3_5_v_energy->at(j);
        //class_br3_5.br3_5_v_n_main_segs=br3_5_v_n_main_segs->at(j);
        class_br3_5.br3_5_v_n_segs=br3_5_v_n_segs->at(j);
        class_br3_5.br3_5_v_shower_main_length=br3_5_v_shower_main_length->at(j);
        class_br3_5.br3_5_v_shower_total_length=br3_5_v_shower_total_length->at(j);

        float temp_score = class_br3_5.evaluate();
        bbddtt->push_back(temp_score);
        br3_5_flag_temp *= br3_5_v_flag->at(j);
    }

    if(br3_5_v_flag->size()==0) br3_5_score = 0.42/default_scale; //0.42
    else br3_5_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // br3 tagger 6 
    int br3_6_flag_temp = 1;  
    if(br3_6_on){
    bbddtt->clear();
    for(size_t j=0; j<br3_6_v_flag->size(); j++)
    {
        class_br3_6.reset();
        
        class_br3_6.br3_6_v_angle=br3_6_v_angle->at(j);
        class_br3_6.br3_6_v_angle1=br3_6_v_angle1->at(j);
        class_br3_6.br3_6_v_flag_shower_trajectory=br3_6_v_flag_shower_trajectory->at(j);
        class_br3_6.br3_6_v_direct_length=br3_6_v_direct_length->at(j);
        class_br3_6.br3_6_v_length=br3_6_v_length->at(j);
        class_br3_6.br3_6_v_n_other_vtx_segs=br3_6_v_n_other_vtx_segs->at(j);
        class_br3_6.br3_6_v_energy=br3_6_v_energy->at(j);

        float temp_score = class_br3_6.evaluate();
        bbddtt->push_back(temp_score);
        br3_6_flag_temp *= br3_6_v_flag->at(j);
    }

    if(br3_6_v_flag->size()==0) br3_6_score = 0.75/default_scale;
    else br3_6_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // lol tagger 1 
    int lol_1_flag_temp = 1;  
    if(lol_1_on){
    bbddtt->clear();
    for(size_t j=0; j<lol_1_v_flag->size(); j++)
    {
        class_lol_1.reset();

        class_lol_1.lol_1_v_energy=lol_1_v_energy->at(j);
        class_lol_1.lol_1_v_vtx_n_segs=lol_1_v_vtx_n_segs->at(j);
        class_lol_1.lol_1_v_nseg=lol_1_v_nseg->at(j);
        class_lol_1.lol_1_v_angle=lol_1_v_angle->at(j);

        float temp_score = class_lol_1.evaluate();
        bbddtt->push_back(temp_score);
        lol_1_flag_temp *= lol_1_v_flag->at(j);
    }

    if(lol_1_v_flag->size()==0) lol_1_score = 0.85/default_scale;
    else lol_1_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // lol tagger 2 
    int lol_2_flag_temp = 1;  
    if(lol_2_on){
    bbddtt->clear();
    for(size_t j=0; j<lol_2_v_flag->size(); j++)
    {
        class_lol_2.reset();

        class_lol_2.lol_2_v_length=lol_2_v_length->at(j);
        class_lol_2.lol_2_v_angle=lol_2_v_angle->at(j);
        class_lol_2.lol_2_v_type=lol_2_v_type->at(j);
        class_lol_2.lol_2_v_vtx_n_segs=lol_2_v_vtx_n_segs->at(j);
        class_lol_2.lol_2_v_energy=lol_2_v_energy->at(j);
        class_lol_2.lol_2_v_shower_main_length=lol_2_v_shower_main_length->at(j);
        class_lol_2.lol_2_v_flag_dir_weak=lol_2_v_flag_dir_weak->at(j);

        float temp_score = class_lol_2.evaluate();
        bbddtt->push_back(temp_score);
        lol_2_flag_temp *= lol_2_v_flag->at(j);
    }

    if(lol_2_v_flag->size()==0) lol_2_score = 0.7/default_scale;
    else lol_2_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // tro tagger 1 
    int tro_1_flag_temp = 1;  
    if(tro_1_on){
    bbddtt->clear();
    for(size_t j=0; j<tro_1_v_flag->size(); j++)
    {
        class_tro_1.reset();

        class_tro_1.tro_1_v_particle_type=tro_1_v_particle_type->at(j);
        class_tro_1.tro_1_v_flag_dir_weak=tro_1_v_flag_dir_weak->at(j);
        class_tro_1.tro_1_v_min_dis=tro_1_v_min_dis->at(j);
        class_tro_1.tro_1_v_sg1_length=tro_1_v_sg1_length->at(j);
        class_tro_1.tro_1_v_shower_main_length=tro_1_v_shower_main_length->at(j);
        class_tro_1.tro_1_v_max_n_vtx_segs=tro_1_v_max_n_vtx_segs->at(j);
        class_tro_1.tro_1_v_tmp_length=tro_1_v_tmp_length->at(j);
        class_tro_1.tro_1_v_medium_dQ_dx=tro_1_v_medium_dQ_dx->at(j);
        float tro_1_v_dQ_dx_cut_temp = tro_1_v_dQ_dx_cut->at(j);
        if(std::isinf(tro_1_v_dQ_dx_cut_temp) || std::isnan(tro_1_v_dQ_dx_cut_temp)) tro_1_v_dQ_dx_cut_temp = 10;
        class_tro_1.tro_1_v_dQ_dx_cut=tro_1_v_dQ_dx_cut_temp;
        class_tro_1.tro_1_v_flag_shower_topology=tro_1_v_flag_shower_topology->at(j);

        float temp_score = class_tro_1.evaluate();
        bbddtt->push_back(temp_score);
        tro_1_flag_temp *= tro_1_v_flag->at(j);
    }

    if(tro_1_v_flag->size()==0) tro_1_score = 0.28/default_scale;
    else tro_1_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // tro tagger 2 
    int tro_2_flag_temp = 1;  
    if(tro_2_on){
    bbddtt->clear();
    for(size_t j=0; j<tro_2_v_flag->size(); j++)
    {
        class_tro_2.reset();

        class_tro_2.tro_2_v_energy=tro_2_v_energy->at(j);
        class_tro_2.tro_2_v_stem_length=tro_2_v_stem_length->at(j);
        class_tro_2.tro_2_v_iso_angle=tro_2_v_iso_angle->at(j);
        class_tro_2.tro_2_v_max_length=tro_2_v_max_length->at(j);
        class_tro_2.tro_2_v_angle=tro_2_v_angle->at(j);

        float temp_score = class_tro_2.evaluate();
        bbddtt->push_back(temp_score);
        tro_2_flag_temp *= tro_2_v_flag->at(j);
    }

    if(tro_2_v_flag->size()==0) tro_2_score = 0.35/default_scale;
    else tro_2_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // tro tagger 4 
    int tro_4_flag_temp = 1;  
    if(tro_4_on){
    bbddtt->clear();
    for(size_t j=0; j<tro_4_v_flag->size(); j++)
    {
        class_tro_4.reset();
    
        class_tro_4.tro_4_v_dir2_mag=tro_4_v_dir2_mag->at(j);
        class_tro_4.tro_4_v_angle=tro_4_v_angle->at(j);
        class_tro_4.tro_4_v_angle1=tro_4_v_angle1->at(j);
        class_tro_4.tro_4_v_angle2=tro_4_v_angle2->at(j);
        class_tro_4.tro_4_v_length=tro_4_v_length->at(j);
        class_tro_4.tro_4_v_length1=tro_4_v_length1->at(j);
        class_tro_4.tro_4_v_medium_dQ_dx=tro_4_v_medium_dQ_dx->at(j);
        class_tro_4.tro_4_v_end_dQ_dx=tro_4_v_end_dQ_dx->at(j);
        class_tro_4.tro_4_v_energy=tro_4_v_energy->at(j);
        class_tro_4.tro_4_v_shower_main_length=tro_4_v_shower_main_length->at(j);
        class_tro_4.tro_4_v_flag_shower_trajectory=tro_4_v_flag_shower_trajectory->at(j);

        float temp_score = class_tro_4.evaluate();
        bbddtt->push_back(temp_score);
        tro_4_flag_temp *= tro_4_v_flag->at(j);
    }

    if(tro_4_v_flag->size()==0) tro_4_score = 0.33/default_scale;
    else tro_4_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }

    // tro tagger 5 
    int tro_5_flag_temp = 1;  
    if(tro_5_on){
    bbddtt->clear();
    for(size_t j=0; j<tro_5_v_flag->size(); j++)
    {
        class_tro_5.reset();
    
        class_tro_5.tro_5_v_max_angle=tro_5_v_max_angle->at(j);
        class_tro_5.tro_5_v_min_angle=tro_5_v_min_angle->at(j);
        class_tro_5.tro_5_v_max_length=tro_5_v_max_length->at(j);
        class_tro_5.tro_5_v_iso_angle=tro_5_v_iso_angle->at(j);
        class_tro_5.tro_5_v_n_vtx_segs=tro_5_v_n_vtx_segs->at(j);
        class_tro_5.tro_5_v_min_count=tro_5_v_min_count->at(j);
        class_tro_5.tro_5_v_max_count=tro_5_v_max_count->at(j);
        class_tro_5.tro_5_v_energy=tro_5_v_energy->at(j);

        float temp_score = class_tro_5.evaluate();
        bbddtt->push_back(temp_score);
        tro_5_flag_temp *= tro_5_v_flag->at(j);
    }

    if(tro_5_v_flag->size()==0) tro_5_score = 0.5/default_scale;
    else tro_5_score = *std::min_element(bbddtt->begin(), bbddtt->end());
    }



    // TTree fill
    temp_flag = mip_flag*gap_flag*cme_flag*anc_flag*mgo_flag*mgt_flag*br1_flag*br3_1_flag*br3_2_flag*br3_4_flag*br3_7_flag*br3_8_flag
        *(stem_dir_flag||!stem_dir_flag_single_shower)*br2_flag
        *stem_len_flag*lem_flag*brm_flag*tro_3_flag*br4_flag*mip_quality_flag
        *(pio_mip_id||!pio_filled||!pio_flag_pio||pio_1_flag)*(pio_mip_id||!pio_filled||pio_flag_pio||pio_2_flag_temp)
        *stw_1_flag*spt_flag*vis_1_flag*vis_2_flag*stw_2_flag_temp*stw_3_flag_temp*stw_4_flag_temp*sig_1_flag_temp*sig_2_flag_temp
        *br3_3_flag_temp*br3_5_flag_temp*br3_6_flag_temp*hol_1_flag*hol_2_flag*lol_3_flag*lol_1_flag_temp*lol_2_flag_temp
        *tro_1_flag_temp*tro_2_flag_temp*tro_4_flag_temp*tro_5_flag_temp;
    rest_flag = 1; // for any specific test 

    mip_energy_f = mip_energy;
    /* gap_energy_f = gap_energy; */ 
    /* gap_flag_single_shower_f = gap_flag_single_shower; */ 
    /* gap_num_valid_tracks_f = gap_num_valid_tracks; */ 
    /* gap_n_bad_f = gap_n_bad; */ 
    /* mgt_total_other_energy_f = mgt_total_other_energy; */ 
    /* mgt_e_indirect_total_energy_f = mgt_e_indirect_total_energy; */ 
    /* mgt_e_direct_total_energy_f = mgt_e_direct_total_energy; */  
    mip_vec_dQ_dx_0_f = mip_vec_dQ_dx_0; 
    mip_vec_dQ_dx_1_f = mip_vec_dQ_dx_1; 
    mip_vec_dQ_dx_2_f = mip_vec_dQ_dx_2; 
    mip_vec_dQ_dx_3_f = mip_vec_dQ_dx_3; 
    mip_vec_dQ_dx_4_f = mip_vec_dQ_dx_4; 
    /* mip_vec_dQ_dx_5_f = mip_vec_dQ_dx_5; */ 
    /* mip_vec_dQ_dx_6_f = mip_vec_dQ_dx_6; */ 
    spt_angle_beam_f = spt_angle_beam; 
    spt_angle_drift_f = spt_angle_drift; 
    cme_mu_energy_f = cme_mu_energy; 
    cme_mu_length_f = cme_mu_length; 
    cme_length_f = cme_length; 

    sig->Fill();

    }

    sig->Write();
    outputfile->Close();
    inputfile->Close();
    return 0;
}
