// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

using namespace std;

TFile *input = 0;
TFile *output = 0;
TMVA::Factory *factory = 0;
TMVA::DataLoader *dataloader = 0;


void TestEvaluate(double min, double max, int nbins, float BDTcut);
void GetROC(TH1F* hs, TH1F* hb, TGraph* roc);
void GetEff(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TGraph* eff);
void GetPurity(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TH1F* hb, double bfrac, TGraph* purity);
void GetEffPurity(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TH1F* hb, double bfrac, TGraph* ep);

void TestEvaluate(TString filename, double min, double max, int nbins, float BDTcut = 0)
{
    gROOT->ProcessLine(".x DrawOption.cc");
    TH1F* hsig = new TH1F("hsig", "BDT score for sig", nbins/10, min, max);
    TH1F* hbkg1 = new TH1F("hbkg1", "BDT score for bkg1", nbins/10, min, max);
    TH1F* hbkg2 = new TH1F("hbkg2", "BDT score for bkg2", nbins/10, min, max);
    TH1F* hbkg3 = new TH1F("hbkg3", "BDT score for bkg3", nbins/10, min, max);
    
    TH1F* hsigROC = new TH1F("hsigROC", "BDT score for sigROC", nbins, min, max);
    TH1F* hbkgROC = new TH1F("hbkgROC", "BDT score for bkgROC", nbins, min, max);

    ///  two additional plots
    ///  signal efficiency vs true Edep
    TH1F* hs_cut = new TH1F("hs_cut", "Signal cut-based", 30, 0, 3000);
    TH1F* hs_allcut = new TH1F("hs_allcut", "Signal cut-based", 30, 0, 3000);
    TH1F* hs_bdt = new TH1F("hs_bdt", "Signal BDT", 30, 0, 3000);
    TH1F* hstotal = new TH1F("hstotal", "All signal events", 30, 0, 3000);
    TH1F* hseff_cut = new TH1F("hseff_cut", "Signal efficiency cut-based", 30, 0, 3000);
    TH1F* hseff_allcut = new TH1F("hseff_allcut", "Signal efficiency cut-based", 30, 0, 3000);
    TH1F* hseff_bdt = new TH1F("hseff_bdt", "Signal efficiency BDT", 30, 0, 3000);
    TH1F* hb_cut = new TH1F("hb_cut", "Bkg cut-based", 30, 0, 3000);
    TH1F* hb_allcut = new TH1F("hb_allcut", "Bkg cut-based", 30, 0, 3000);
    TH1F* hb_bdt = new TH1F("hb_bdt", "Bkg BDT", 30, 0, 3000);
    TH1F* hbtotal = new TH1F("hbtotal", "All bkg events", 30, 0, 3000);  
    TH1F* hbeff_cut = new TH1F("hbeff_cut", "Bkg passing rate cut-based", 30, 0, 3000);
    TH1F* hbeff_allcut = new TH1F("hbeff_allcut", "Bkg passing rate cut-based", 30, 0, 3000);
    TH1F* hbeff_bdt = new TH1F("hbeff_bdt", "Bkg passing rate BDT", 30, 0, 3000);
    TH1F* hpurity_cut = new TH1F("hpurity_cut", "Purity cut-based", 30, 0, 3000); 
    TH1F* hpurity_allcut = new TH1F("hpurity_allcut", "Purity cut-based", 30, 0, 3000); 
    TH1F* hpurity_bdt = new TH1F("hpurity_bdt", "Purity BDT", 30, 0, 3000); 
    ///

/// read variables from file

  Int_t run, subrun, event;
  Int_t nueTag;
  Int_t truth_nue;
  Int_t truth_CC;
  Int_t truth_inFV;
  Int_t truth_cosmic;
  Float_t trueEdep;  
  Float_t trueEnu;  
  Float_t nuvtx_diff;  
  Float_t showervtx_diff;  
  Float_t weight;
  Float_t lowEweight;

  double mip_energy;
  int mip_n_end_reduction;
  int mip_n_first_mip;
  int mip_n_first_non_mip;
  int mip_n_first_non_mip_1;
  int mip_n_first_non_mip_2;
  double mip_vec_dQ_dx_0;
  double mip_vec_dQ_dx_1;
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
  int gap_flag_prolong_u;
  int gap_flag_prolong_v;
  int gap_flag_prolong_w;
  int gap_flag_parallel;
  int gap_n_points;
  int gap_n_bad;
  double gap_energy;
  int gap_num_valid_tracks;
  int gap_flag_single_shower;
  double cme_mu_energy;
  double cme_energy;
  double cme_mu_length;
  double cme_length;
  double cme_angle_beam;
  double anc_angle;
  double anc_max_angle;
  double anc_max_length;
  double anc_acc_forward_length;
  double anc_acc_backward_length;
  double anc_acc_forward_length1;
  double anc_shower_main_length;
  double anc_shower_total_length;
  int anc_flag_main_outside;
  double mgo_energy;
  double mgo_max_energy;
  double mgo_total_energy;
  int mgo_n_showers;
  double mgo_max_energy_1;
  double mgo_max_energy_2;
  double mgo_total_other_energy;
  int mgo_n_total_showers;
  double mgo_total_other_energy_1;
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
  double br3_1_energy;
  int br3_1_n_shower_segments;
  int br3_1_sg_flag_trajectory;
  double br3_1_sg_direct_length;
  double br3_1_sg_length;
  double br3_1_total_main_length;
  double br3_1_total_length;
  double br3_1_iso_angle;
  int br3_1_sg_flag_topology;
  int br3_2_n_ele;
  int br3_2_n_other;
  int br3_2_other_fid;
  double br3_4_acc_length;
  double br3_4_total_length;
  double br3_7_min_angle;
  double br3_8_max_dQ_dx;
  int br3_8_n_main_segs;
  int stem_dir_flag_single_shower;
  double stem_dir_angle;
  double stem_dir_energy;
  double stem_dir_angle1;
  double stem_dir_angle2;
  double stem_dir_angle3;
  double stem_dir_ratio;
  int br2_num_valid_tracks;
  int br2_n_shower_main_segs;
  double br2_max_angle;
  double br2_sg_length;
  int br2_flag_sg_trajectory;
  double stem_len_energy;
  double stem_len_length;
  int stem_len_flag_avoid_muon_check;
  int stem_len_num_daughters;
  double stem_len_daughter_length;
  int brm_n_mu_segs;
  double brm_Ep;
  double brm_acc_length;
  double brm_shower_total_length;
  double brm_connected_length;
  int brm_n_size;
  double brm_acc_direct_length;
  int brm_n_shower_main_segs;
  int brm_n_mu_main;
  double lem_shower_main_length;
  int lem_n_3seg;
  double lem_e_charge;
  double lem_e_dQdx;
  int lem_shower_num_main_segs;
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
  double tro_3_stem_length;
  int tro_3_n_muon_segs;
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
  int pio_mip_id;
  double pio_1_mass;
  int pio_1_pio_type;
  double pio_1_energy_1;
  double pio_1_energy_2;
  double pio_1_dis_1;
  double pio_1_dis_2;
  double stw_1_energy;
  double stw_1_dis;
  double stw_1_dQ_dx;
  int stw_1_flag_single_shower;
  int stw_1_n_pi0;
  int stw_1_num_valid_tracks;
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
  int vis_1_n_vtx_segs;
  double vis_1_energy;
  int vis_1_num_good_tracks;
  double vis_1_max_angle;
  double vis_1_max_shower_angle;
  double vis_1_tmp_length1;
  double vis_1_tmp_length2;
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
  int hol_1_n_valid_tracks;
  double hol_1_min_angle;
  double hol_1_energy;
  int hol_1_flag_all_shower;
  double hol_1_min_length;
  double hol_2_min_angle;
  double hol_2_medium_dQ_dx;
  int hol_2_ncount;
  double lol_3_angle_beam;
  int lol_3_n_valid_tracks;
  double lol_3_min_angle;
  int lol_3_vtx_n_segs;
  double lol_3_shower_main_length;
  int lol_3_n_out;
  int lol_3_n_sum;
 
  double mip_vec_dQ_dx_2 = 0;
  double mip_vec_dQ_dx_3 = 0;
  double mip_vec_dQ_dx_4 = 0;
  double mip_vec_dQ_dx_5 = 0;
  double mip_vec_dQ_dx_6 = 0;
  double mip_vec_dQ_dx_7 = 0;
  double mip_vec_dQ_dx_8 = 0;
  double mip_vec_dQ_dx_9 = 0;
  double mip_vec_dQ_dx_10 = 0;
  double mip_vec_dQ_dx_11 = 0;
  double mip_vec_dQ_dx_12 = 0;
  double mip_vec_dQ_dx_13 = 0;
  double mip_vec_dQ_dx_14 = 0;
  double mip_vec_dQ_dx_15 = 0;
  double mip_vec_dQ_dx_16 = 0;
  double mip_vec_dQ_dx_17 = 0;
  double mip_vec_dQ_dx_18 = 0;
  double mip_vec_dQ_dx_19 = 0;

  float br3_3_score = 0;
  float br3_5_score = 0;
  float br3_6_score = 0;
  float pio_2_score = 0;
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

///

    TString fname = filename;
    if (!gSystem->AccessPathName( fname )) {
        input = TFile::Open( fname ); // check if file in local directory exists
    }
    TTree *sigTree = (TTree*)input->Get("sig");
    TTree *bkgTree = (TTree*)input->Get("bkg");

    sigTree->SetBranchAddress("run",&run);
    sigTree->SetBranchAddress("subrun",&subrun);
    sigTree->SetBranchAddress("event",&event);
    sigTree->SetBranchAddress("nueTag",&nueTag);
    sigTree->SetBranchAddress("truth_nue",&truth_nue);
    sigTree->SetBranchAddress("truth_CC",&truth_CC);
    sigTree->SetBranchAddress("truth_inFV",&truth_inFV);
    sigTree->SetBranchAddress("truth_cosmic",&truth_cosmic);
    sigTree->SetBranchAddress("weight",&weight);
    sigTree->SetBranchAddress("lowEweight",&lowEweight);
    sigTree->SetBranchAddress("trueEdep",&trueEdep);
    sigTree->SetBranchAddress("trueEnu",&trueEnu);
    sigTree->SetBranchAddress("nuvtx_diff",&nuvtx_diff);
    sigTree->SetBranchAddress("showervtx_diff",&showervtx_diff);
    sigTree->SetBranchAddress("nueTag",&nueTag);

    sigTree->SetBranchAddress("mip_energy",&mip_energy);
    sigTree->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
    sigTree->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
    sigTree->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
    sigTree->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
    sigTree->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
    sigTree->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
    sigTree->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
    sigTree->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
    sigTree->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
    sigTree->SetBranchAddress("mip_n_highest",&mip_n_highest);
    sigTree->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
    sigTree->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
    sigTree->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
    sigTree->SetBranchAddress("mip_stem_length",&mip_stem_length);
    sigTree->SetBranchAddress("mip_length_main",&mip_length_main);
    sigTree->SetBranchAddress("mip_length_total",&mip_length_total);
    sigTree->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
    sigTree->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
    sigTree->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
    sigTree->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
    sigTree->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
    sigTree->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
    sigTree->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
    sigTree->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
    sigTree->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
    sigTree->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
    sigTree->SetBranchAddress("mip_min_dis",&mip_min_dis);
    sigTree->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
    sigTree->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
    sigTree->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
    sigTree->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
    sigTree->SetBranchAddress("gap_n_points",&gap_n_points);
    sigTree->SetBranchAddress("gap_n_bad",&gap_n_bad);
    sigTree->SetBranchAddress("gap_energy",&gap_energy);
    sigTree->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
    sigTree->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);
    sigTree->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
    sigTree->SetBranchAddress("cme_energy",&cme_energy);
    sigTree->SetBranchAddress("cme_mu_length",&cme_mu_length);
    sigTree->SetBranchAddress("cme_length",&cme_length);
    sigTree->SetBranchAddress("cme_angle_beam",&cme_angle_beam);
    sigTree->SetBranchAddress("anc_angle",&anc_angle);
    sigTree->SetBranchAddress("anc_max_angle",&anc_max_angle);
    sigTree->SetBranchAddress("anc_max_length",&anc_max_length);
    sigTree->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
    sigTree->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
    sigTree->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
    sigTree->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
    sigTree->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
    sigTree->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);
    sigTree->SetBranchAddress("mgo_energy",&mgo_energy);
    sigTree->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
    sigTree->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
    sigTree->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
    sigTree->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
    sigTree->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
    sigTree->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
    sigTree->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
    sigTree->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);
    sigTree->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
    sigTree->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
    sigTree->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
    sigTree->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
    sigTree->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
    sigTree->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
    sigTree->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
    sigTree->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
    sigTree->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
    sigTree->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);
    sigTree->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
    sigTree->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
    sigTree->SetBranchAddress("br1_1_energy",&br1_1_energy);
    sigTree->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
    sigTree->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
    sigTree->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
    sigTree->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);
    sigTree->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
    sigTree->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
    sigTree->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
    sigTree->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
    sigTree->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
    sigTree->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);
    sigTree->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
    sigTree->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
    sigTree->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);
    sigTree->SetBranchAddress("br3_1_energy",&br3_1_energy);
    sigTree->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
    sigTree->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
    sigTree->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
    sigTree->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
    sigTree->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
    sigTree->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
    sigTree->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
    sigTree->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
    sigTree->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
    sigTree->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
    sigTree->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);
    sigTree->SetBranchAddress("br3_4_acc_length",&br3_4_acc_length);
    sigTree->SetBranchAddress("br3_4_total_length",&br3_4_total_length);
    sigTree->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);
    sigTree->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
    sigTree->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);
    sigTree->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
    sigTree->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
    sigTree->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
    sigTree->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
    sigTree->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
    sigTree->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
    sigTree->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);
    sigTree->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
    sigTree->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
    sigTree->SetBranchAddress("br2_max_angle",&br2_max_angle);
    sigTree->SetBranchAddress("br2_sg_length",&br2_sg_length);
    sigTree->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);
    sigTree->SetBranchAddress("stem_len_energy",&stem_len_energy);
    sigTree->SetBranchAddress("stem_len_length",&stem_len_length);
    sigTree->SetBranchAddress("stem_len_flag_avoid_muon_check",&stem_len_flag_avoid_muon_check);
    sigTree->SetBranchAddress("stem_len_num_daughters",&stem_len_num_daughters);
    sigTree->SetBranchAddress("stem_len_daughter_length",&stem_len_daughter_length);
    sigTree->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
    sigTree->SetBranchAddress("brm_Ep",&brm_Ep);
    sigTree->SetBranchAddress("brm_acc_length",&brm_acc_length);
    sigTree->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
    sigTree->SetBranchAddress("brm_connected_length",&brm_connected_length);
    sigTree->SetBranchAddress("brm_n_size",&brm_n_size);
    sigTree->SetBranchAddress("brm_nacc_direct_length",&brm_acc_direct_length); // naming issue
    sigTree->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
    sigTree->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);
    sigTree->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
    sigTree->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
    sigTree->SetBranchAddress("lem_e_charge",&lem_e_charge);
    sigTree->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
    sigTree->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);
    sigTree->SetBranchAddress("br4_1_shower_main_length",&br4_1_shower_main_length);
    sigTree->SetBranchAddress("br4_1_shower_total_length",&br4_1_shower_total_length);
    sigTree->SetBranchAddress("br4_1_min_dis",&br4_1_min_dis);
    sigTree->SetBranchAddress("br4_1_energy",&br4_1_energy);
    sigTree->SetBranchAddress("br4_1_flag_avoid_muon_check",&br4_1_flag_avoid_muon_check);
    sigTree->SetBranchAddress("br4_1_n_vtx_segs",&br4_1_n_vtx_segs);
    sigTree->SetBranchAddress("br4_1_br4_1_n_main_segs",&br4_1_n_main_segs); // naming issue
    sigTree->SetBranchAddress("br4_2_ratio_45",&br4_2_ratio_45);
    sigTree->SetBranchAddress("br4_2_ratio_35",&br4_2_ratio_35);
    sigTree->SetBranchAddress("br4_2_ratio_25",&br4_2_ratio_25);
    sigTree->SetBranchAddress("br4_2_ratio_15",&br4_2_ratio_15);
    sigTree->SetBranchAddress("br4_2_ratio1_45",&br4_2_ratio1_45);
    sigTree->SetBranchAddress("br4_2_ratio1_35",&br4_2_ratio1_35);
    sigTree->SetBranchAddress("br4_2_ratio1_25",&br4_2_ratio1_25);
    sigTree->SetBranchAddress("br4_2_ratio1_15",&br4_2_ratio1_15);
    sigTree->SetBranchAddress("br4_2_iso_angle",&br4_2_iso_angle);
    sigTree->SetBranchAddress("br4_2_iso_angle1",&br4_2_iso_angle1);
    sigTree->SetBranchAddress("br4_2_angle",&br4_2_angle);
    sigTree->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
    sigTree->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
    sigTree->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
    sigTree->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
    sigTree->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
    sigTree->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
    sigTree->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
    sigTree->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
    sigTree->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
    sigTree->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
    sigTree->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
    sigTree->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
    sigTree->SetBranchAddress("pio_mip_id",&pio_mip_id);
    sigTree->SetBranchAddress("pio_1_mass",&pio_1_mass);
    sigTree->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
    sigTree->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
    sigTree->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
    sigTree->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
    sigTree->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);
    sigTree->SetBranchAddress("stw_1_energy",&stw_1_energy);
    sigTree->SetBranchAddress("stw_1_dis",&stw_1_dis);
    sigTree->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
    sigTree->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
    sigTree->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
    sigTree->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
    sigTree->SetBranchAddress("spt_shower_main_length",&spt_shower_main_length);
    sigTree->SetBranchAddress("spt_shower_total_length",&spt_shower_total_length);
    sigTree->SetBranchAddress("spt_angle_beam",&spt_angle_beam);
    sigTree->SetBranchAddress("spt_angle_vertical",&spt_angle_vertical);
    sigTree->SetBranchAddress("spt_max_dQ_dx",&spt_max_dQ_dx);
    sigTree->SetBranchAddress("spt_angle_beam_1",&spt_angle_beam_1);
    sigTree->SetBranchAddress("spt_angle_drift",&spt_angle_drift);
    sigTree->SetBranchAddress("spt_angle_drift_1",&spt_angle_drift_1);
    sigTree->SetBranchAddress("spt_num_valid_tracks",&spt_num_valid_tracks);
    sigTree->SetBranchAddress("spt_n_vtx_segs",&spt_n_vtx_segs);
    sigTree->SetBranchAddress("spt_max_length",&spt_max_length);
    sigTree->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
    sigTree->SetBranchAddress("vis_1_energy",&vis_1_energy);
    sigTree->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
    sigTree->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
    sigTree->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
    sigTree->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
    sigTree->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
    sigTree->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
    sigTree->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
    sigTree->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
    sigTree->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
    sigTree->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
    sigTree->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
    sigTree->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
    sigTree->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
    sigTree->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
    sigTree->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
    sigTree->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);
    sigTree->SetBranchAddress("hol_1_n_valid_tracks",&hol_1_n_valid_tracks);
    sigTree->SetBranchAddress("hol_1_min_angle",&hol_1_min_angle);
    sigTree->SetBranchAddress("hol_1_energy",&hol_1_energy);
    sigTree->SetBranchAddress("hol_1_all_shower",&hol_1_flag_all_shower); // naming issue
    sigTree->SetBranchAddress("hol_1_min_length",&hol_1_min_length);
    sigTree->SetBranchAddress("hol_2_min_angle",&hol_2_min_angle);
    sigTree->SetBranchAddress("hol_2_medium_dQ_dx",&hol_2_medium_dQ_dx);
    sigTree->SetBranchAddress("hol_2_ncount",&hol_2_ncount);
    sigTree->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
    sigTree->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
    sigTree->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
    sigTree->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
    sigTree->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
    sigTree->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
    sigTree->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_5",&mip_vec_dQ_dx_5);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_6",&mip_vec_dQ_dx_6);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_7",&mip_vec_dQ_dx_7);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_8",&mip_vec_dQ_dx_8);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_9",&mip_vec_dQ_dx_9);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_10",&mip_vec_dQ_dx_10);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_11",&mip_vec_dQ_dx_11);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_12",&mip_vec_dQ_dx_12);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_13",&mip_vec_dQ_dx_13);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_14",&mip_vec_dQ_dx_14);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_15",&mip_vec_dQ_dx_15);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_16",&mip_vec_dQ_dx_16);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_17",&mip_vec_dQ_dx_17);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_18",&mip_vec_dQ_dx_18);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_19",&mip_vec_dQ_dx_19);
    sigTree->SetBranchAddress("br3_3_score",&br3_3_score);
    sigTree->SetBranchAddress("br3_5_score",&br3_5_score);
    sigTree->SetBranchAddress("br3_6_score",&br3_6_score);
    sigTree->SetBranchAddress("pio_2_score",&pio_2_score);
    sigTree->SetBranchAddress("stw_2_score",&stw_2_score);
    sigTree->SetBranchAddress("stw_3_score",&stw_3_score);
    sigTree->SetBranchAddress("stw_4_score",&stw_4_score);
    sigTree->SetBranchAddress("sig_1_score",&sig_1_score);
    sigTree->SetBranchAddress("sig_2_score",&sig_2_score);
    sigTree->SetBranchAddress("lol_1_score",&lol_1_score);
    sigTree->SetBranchAddress("lol_2_score",&lol_2_score);
    sigTree->SetBranchAddress("tro_1_score",&tro_1_score);
    sigTree->SetBranchAddress("tro_2_score",&tro_2_score);
    sigTree->SetBranchAddress("tro_4_score",&tro_4_score);
    sigTree->SetBranchAddress("tro_5_score",&tro_5_score);

    bkgTree->SetBranchAddress("run",&run);
    bkgTree->SetBranchAddress("subrun",&subrun);
    bkgTree->SetBranchAddress("event",&event);
    bkgTree->SetBranchAddress("nueTag",&nueTag);
    bkgTree->SetBranchAddress("truth_nue",&truth_nue);
    bkgTree->SetBranchAddress("truth_CC",&truth_CC);
    bkgTree->SetBranchAddress("truth_inFV",&truth_inFV);
    bkgTree->SetBranchAddress("truth_cosmic",&truth_cosmic);
    bkgTree->SetBranchAddress("weight",&weight);
    bkgTree->SetBranchAddress("lowEweight",&lowEweight);
    bkgTree->SetBranchAddress("trueEdep",&trueEdep);
    bkgTree->SetBranchAddress("trueEnu",&trueEnu);
    bkgTree->SetBranchAddress("nuvtx_diff",&nuvtx_diff);
    bkgTree->SetBranchAddress("showervtx_diff",&showervtx_diff);
    bkgTree->SetBranchAddress("nueTag",&nueTag);
   
    bkgTree->SetBranchAddress("mip_energy",&mip_energy);
    bkgTree->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
    bkgTree->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
    bkgTree->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
    bkgTree->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
    bkgTree->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
    bkgTree->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
    bkgTree->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
    bkgTree->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
    bkgTree->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
    bkgTree->SetBranchAddress("mip_n_highest",&mip_n_highest);
    bkgTree->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
    bkgTree->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
    bkgTree->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
    bkgTree->SetBranchAddress("mip_stem_length",&mip_stem_length);
    bkgTree->SetBranchAddress("mip_length_main",&mip_length_main);
    bkgTree->SetBranchAddress("mip_length_total",&mip_length_total);
    bkgTree->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
    bkgTree->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
    bkgTree->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
    bkgTree->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
    bkgTree->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
    bkgTree->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
    bkgTree->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
    bkgTree->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
    bkgTree->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
    bkgTree->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
    bkgTree->SetBranchAddress("mip_min_dis",&mip_min_dis);
    bkgTree->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
    bkgTree->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
    bkgTree->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
    bkgTree->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
    bkgTree->SetBranchAddress("gap_n_points",&gap_n_points);
    bkgTree->SetBranchAddress("gap_n_bad",&gap_n_bad);
    bkgTree->SetBranchAddress("gap_energy",&gap_energy);
    bkgTree->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
    bkgTree->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);
    bkgTree->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
    bkgTree->SetBranchAddress("cme_energy",&cme_energy);
    bkgTree->SetBranchAddress("cme_mu_length",&cme_mu_length);
    bkgTree->SetBranchAddress("cme_length",&cme_length);
    bkgTree->SetBranchAddress("cme_angle_beam",&cme_angle_beam);
    bkgTree->SetBranchAddress("anc_angle",&anc_angle);
    bkgTree->SetBranchAddress("anc_max_angle",&anc_max_angle);
    bkgTree->SetBranchAddress("anc_max_length",&anc_max_length);
    bkgTree->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
    bkgTree->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
    bkgTree->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
    bkgTree->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
    bkgTree->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
    bkgTree->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);
    bkgTree->SetBranchAddress("mgo_energy",&mgo_energy);
    bkgTree->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
    bkgTree->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
    bkgTree->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
    bkgTree->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
    bkgTree->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
    bkgTree->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
    bkgTree->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
    bkgTree->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);
    bkgTree->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
    bkgTree->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
    bkgTree->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
    bkgTree->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
    bkgTree->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
    bkgTree->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
    bkgTree->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
    bkgTree->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
    bkgTree->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
    bkgTree->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);
    bkgTree->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
    bkgTree->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
    bkgTree->SetBranchAddress("br1_1_energy",&br1_1_energy);
    bkgTree->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
    bkgTree->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
    bkgTree->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
    bkgTree->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);
    bkgTree->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
    bkgTree->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
    bkgTree->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
    bkgTree->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
    bkgTree->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
    bkgTree->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);
    bkgTree->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
    bkgTree->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
    bkgTree->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);
    bkgTree->SetBranchAddress("br3_1_energy",&br3_1_energy);
    bkgTree->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
    bkgTree->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
    bkgTree->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
    bkgTree->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
    bkgTree->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
    bkgTree->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
    bkgTree->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
    bkgTree->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
    bkgTree->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
    bkgTree->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
    bkgTree->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);
    bkgTree->SetBranchAddress("br3_4_acc_length",&br3_4_acc_length);
    bkgTree->SetBranchAddress("br3_4_total_length",&br3_4_total_length);
    bkgTree->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);
    bkgTree->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
    bkgTree->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);
    bkgTree->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
    bkgTree->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
    bkgTree->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
    bkgTree->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
    bkgTree->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
    bkgTree->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
    bkgTree->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);
    bkgTree->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
    bkgTree->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
    bkgTree->SetBranchAddress("br2_max_angle",&br2_max_angle);
    bkgTree->SetBranchAddress("br2_sg_length",&br2_sg_length);
    bkgTree->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);
    bkgTree->SetBranchAddress("stem_len_energy",&stem_len_energy);
    bkgTree->SetBranchAddress("stem_len_length",&stem_len_length);
    bkgTree->SetBranchAddress("stem_len_flag_avoid_muon_check",&stem_len_flag_avoid_muon_check);
    bkgTree->SetBranchAddress("stem_len_num_daughters",&stem_len_num_daughters);
    bkgTree->SetBranchAddress("stem_len_daughter_length",&stem_len_daughter_length);
    bkgTree->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
    bkgTree->SetBranchAddress("brm_Ep",&brm_Ep);
    bkgTree->SetBranchAddress("brm_acc_length",&brm_acc_length);
    bkgTree->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
    bkgTree->SetBranchAddress("brm_connected_length",&brm_connected_length);
    bkgTree->SetBranchAddress("brm_n_size",&brm_n_size);
    bkgTree->SetBranchAddress("brm_nacc_direct_length",&brm_acc_direct_length); // naming issue
    bkgTree->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
    bkgTree->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);
    bkgTree->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
    bkgTree->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
    bkgTree->SetBranchAddress("lem_e_charge",&lem_e_charge);
    bkgTree->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
    bkgTree->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);
    bkgTree->SetBranchAddress("br4_1_shower_main_length",&br4_1_shower_main_length);
    bkgTree->SetBranchAddress("br4_1_shower_total_length",&br4_1_shower_total_length);
    bkgTree->SetBranchAddress("br4_1_min_dis",&br4_1_min_dis);
    bkgTree->SetBranchAddress("br4_1_energy",&br4_1_energy);
    bkgTree->SetBranchAddress("br4_1_flag_avoid_muon_check",&br4_1_flag_avoid_muon_check);
    bkgTree->SetBranchAddress("br4_1_n_vtx_segs",&br4_1_n_vtx_segs);
    bkgTree->SetBranchAddress("br4_1_br4_1_n_main_segs",&br4_1_n_main_segs); // naming issue
    bkgTree->SetBranchAddress("br4_2_ratio_45",&br4_2_ratio_45);
    bkgTree->SetBranchAddress("br4_2_ratio_35",&br4_2_ratio_35);
    bkgTree->SetBranchAddress("br4_2_ratio_25",&br4_2_ratio_25);
    bkgTree->SetBranchAddress("br4_2_ratio_15",&br4_2_ratio_15);
    bkgTree->SetBranchAddress("br4_2_ratio1_45",&br4_2_ratio1_45);
    bkgTree->SetBranchAddress("br4_2_ratio1_35",&br4_2_ratio1_35);
    bkgTree->SetBranchAddress("br4_2_ratio1_25",&br4_2_ratio1_25);
    bkgTree->SetBranchAddress("br4_2_ratio1_15",&br4_2_ratio1_15);
    bkgTree->SetBranchAddress("br4_2_iso_angle",&br4_2_iso_angle);
    bkgTree->SetBranchAddress("br4_2_iso_angle1",&br4_2_iso_angle1);
    bkgTree->SetBranchAddress("br4_2_angle",&br4_2_angle);
    bkgTree->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
    bkgTree->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
    bkgTree->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
    bkgTree->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
    bkgTree->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
    bkgTree->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
    bkgTree->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
    bkgTree->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
    bkgTree->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
    bkgTree->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
    bkgTree->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
    bkgTree->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
    bkgTree->SetBranchAddress("pio_mip_id",&pio_mip_id);
    bkgTree->SetBranchAddress("pio_1_mass",&pio_1_mass);
    bkgTree->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
    bkgTree->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
    bkgTree->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
    bkgTree->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
    bkgTree->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);
    bkgTree->SetBranchAddress("stw_1_energy",&stw_1_energy);
    bkgTree->SetBranchAddress("stw_1_dis",&stw_1_dis);
    bkgTree->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
    bkgTree->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
    bkgTree->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
    bkgTree->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
    bkgTree->SetBranchAddress("spt_shower_main_length",&spt_shower_main_length);
    bkgTree->SetBranchAddress("spt_shower_total_length",&spt_shower_total_length);
    bkgTree->SetBranchAddress("spt_angle_beam",&spt_angle_beam);
    bkgTree->SetBranchAddress("spt_angle_vertical",&spt_angle_vertical);
    bkgTree->SetBranchAddress("spt_max_dQ_dx",&spt_max_dQ_dx);
    bkgTree->SetBranchAddress("spt_angle_beam_1",&spt_angle_beam_1);
    bkgTree->SetBranchAddress("spt_angle_drift",&spt_angle_drift);
    bkgTree->SetBranchAddress("spt_angle_drift_1",&spt_angle_drift_1);
    bkgTree->SetBranchAddress("spt_num_valid_tracks",&spt_num_valid_tracks);
    bkgTree->SetBranchAddress("spt_n_vtx_segs",&spt_n_vtx_segs);
    bkgTree->SetBranchAddress("spt_max_length",&spt_max_length);
    bkgTree->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
    bkgTree->SetBranchAddress("vis_1_energy",&vis_1_energy);
    bkgTree->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
    bkgTree->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
    bkgTree->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
    bkgTree->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
    bkgTree->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
    bkgTree->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
    bkgTree->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
    bkgTree->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
    bkgTree->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
    bkgTree->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
    bkgTree->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
    bkgTree->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
    bkgTree->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
    bkgTree->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
    bkgTree->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
    bkgTree->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);
    bkgTree->SetBranchAddress("hol_1_n_valid_tracks",&hol_1_n_valid_tracks);
    bkgTree->SetBranchAddress("hol_1_min_angle",&hol_1_min_angle);
    bkgTree->SetBranchAddress("hol_1_energy",&hol_1_energy);
    bkgTree->SetBranchAddress("hol_1_all_shower",&hol_1_flag_all_shower); // naming issue
    bkgTree->SetBranchAddress("hol_1_min_length",&hol_1_min_length);
    bkgTree->SetBranchAddress("hol_2_min_angle",&hol_2_min_angle);
    bkgTree->SetBranchAddress("hol_2_medium_dQ_dx",&hol_2_medium_dQ_dx);
    bkgTree->SetBranchAddress("hol_2_ncount",&hol_2_ncount);
    bkgTree->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
    bkgTree->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
    bkgTree->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
    bkgTree->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
    bkgTree->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
    bkgTree->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
    bkgTree->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_5",&mip_vec_dQ_dx_5);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_6",&mip_vec_dQ_dx_6);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_7",&mip_vec_dQ_dx_7);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_8",&mip_vec_dQ_dx_8);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_9",&mip_vec_dQ_dx_9);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_10",&mip_vec_dQ_dx_10);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_11",&mip_vec_dQ_dx_11);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_12",&mip_vec_dQ_dx_12);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_13",&mip_vec_dQ_dx_13);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_14",&mip_vec_dQ_dx_14);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_15",&mip_vec_dQ_dx_15);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_16",&mip_vec_dQ_dx_16);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_17",&mip_vec_dQ_dx_17);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_18",&mip_vec_dQ_dx_18);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_19",&mip_vec_dQ_dx_19);
    bkgTree->SetBranchAddress("br3_3_score",&br3_3_score);
    bkgTree->SetBranchAddress("br3_5_score",&br3_5_score);
    bkgTree->SetBranchAddress("br3_6_score",&br3_6_score);
    bkgTree->SetBranchAddress("pio_2_score",&pio_2_score);
    bkgTree->SetBranchAddress("stw_2_score",&stw_2_score);
    bkgTree->SetBranchAddress("stw_3_score",&stw_3_score);
    bkgTree->SetBranchAddress("stw_4_score",&stw_4_score);
    bkgTree->SetBranchAddress("sig_1_score",&sig_1_score);
    bkgTree->SetBranchAddress("sig_2_score",&sig_2_score);
    bkgTree->SetBranchAddress("lol_1_score",&lol_1_score);
    bkgTree->SetBranchAddress("lol_2_score",&lol_2_score);
    bkgTree->SetBranchAddress("tro_1_score",&tro_1_score);
    bkgTree->SetBranchAddress("tro_2_score",&tro_2_score);
    bkgTree->SetBranchAddress("tro_4_score",&tro_4_score);
    bkgTree->SetBranchAddress("tro_5_score",&tro_5_score);

    /// new output file
    TFile* round2 = new TFile("BDT_all.root", "RECREATE");
    TTree* sigR2 = sigTree->CloneTree(0);
    Float_t sigbdtscore;
    sigR2->Branch("bdt_all_score", &sigbdtscore, "bdt_all_score/F");
    TTree* bkgR2 = bkgTree->CloneTree(0);
    Float_t bkgbdtscore;
    bkgR2->Branch("bdt_all_score", &bkgbdtscore, "bdt_all_score/F");
    /////

    //// type conversion to float 
    float mip_energy_f;
    float mip_n_end_reduction_f;
    float mip_n_first_mip_f;
    float mip_n_first_non_mip_f;
    float mip_n_first_non_mip_1_f;
    float mip_n_first_non_mip_2_f;
    float mip_vec_dQ_dx_0_f;
    float mip_vec_dQ_dx_1_f;
    float mip_max_dQ_dx_sample_f;
    float mip_n_below_threshold_f;
    float mip_n_below_zero_f;
    float mip_n_lowest_f;
    float mip_n_highest_f;
    float mip_lowest_dQ_dx_f;
    float mip_highest_dQ_dx_f;
    float mip_medium_dQ_dx_f;
    float mip_stem_length_f;
    float mip_length_main_f;
    float mip_length_total_f;
    float mip_angle_beam_f;
    float mip_iso_angle_f;
    float mip_n_vertex_f;
    float mip_n_good_tracks_f;
    float mip_E_indirect_max_energy_f;
    float mip_flag_all_above_f;
    float mip_min_dQ_dx_5_f;
    float mip_n_other_vertex_f;
    float mip_n_stem_size_f;
    float mip_flag_stem_trajectory_f;
    float mip_min_dis_f;
    float gap_flag_prolong_u_f;
    float gap_flag_prolong_v_f;
    float gap_flag_prolong_w_f;
    float gap_flag_parallel_f;
    float gap_n_points_f;
    float gap_n_bad_f;
    float gap_energy_f;
    float gap_num_valid_tracks_f;
    float gap_flag_single_shower_f;
    float cme_mu_energy_f;
    float cme_energy_f;
    float cme_mu_length_f;
    float cme_length_f;
    float cme_angle_beam_f;
    float anc_angle_f;
    float anc_max_angle_f;
    float anc_max_length_f;
    float anc_acc_forward_length_f;
    float anc_acc_backward_length_f;
    float anc_acc_forward_length1_f;
    float anc_shower_main_length_f;
    float anc_shower_total_length_f;
    float anc_flag_main_outside_f;
    float mgo_energy_f;
    float mgo_max_energy_f;
    float mgo_total_energy_f;
    float mgo_n_showers_f;
    float mgo_max_energy_1_f;
    float mgo_max_energy_2_f;
    float mgo_total_other_energy_f;
    float mgo_n_total_showers_f;
    float mgo_total_other_energy_1_f;
    float mgt_flag_single_shower_f;
    float mgt_max_energy_f;
    float mgt_total_other_energy_f;
    float mgt_max_energy_1_f;
    float mgt_e_indirect_max_energy_f;
    float mgt_e_direct_max_energy_f;
    float mgt_n_direct_showers_f;
    float mgt_e_direct_total_energy_f;
    float mgt_flag_indirect_max_pio_f;
    float mgt_e_indirect_total_energy_f;
    float br1_1_shower_type_f;
    float br1_1_vtx_n_segs_f;
    float br1_1_energy_f;
    float br1_1_n_segs_f;
    float br1_1_flag_sg_topology_f;
    float br1_1_flag_sg_trajectory_f;
    float br1_1_sg_length_f;
    float br1_2_n_connected_f;
    float br1_2_max_length_f;
    float br1_2_n_connected_1_f;
    float br1_2_n_shower_segs_f;
    float br1_2_max_length_ratio_f;
    float br1_2_shower_length_f;
    float br1_3_n_connected_p_f;
    float br1_3_max_length_p_f;
    float br1_3_n_shower_main_segs_f;
    float br3_1_energy_f;
    float br3_1_n_shower_segments_f;
    float br3_1_sg_flag_trajectory_f;
    float br3_1_sg_direct_length_f;
    float br3_1_sg_length_f;
    float br3_1_total_main_length_f;
    float br3_1_total_length_f;
    float br3_1_iso_angle_f;
    float br3_1_sg_flag_topology_f;
    float br3_2_n_ele_f;
    float br3_2_n_other_f;
    float br3_2_other_fid_f;
    float br3_4_acc_length_f;
    float br3_4_total_length_f;
    float br3_7_min_angle_f;
    float br3_8_max_dQ_dx_f;
    float br3_8_n_main_segs_f;
    float stem_dir_flag_single_shower_f;
    float stem_dir_angle_f;
    float stem_dir_energy_f;
    float stem_dir_angle1_f;
    float stem_dir_angle2_f;
    float stem_dir_angle3_f;
    float stem_dir_ratio_f;
    float br2_num_valid_tracks_f;
    float br2_n_shower_main_segs_f;
    float br2_max_angle_f;
    float br2_sg_length_f;
    float br2_flag_sg_trajectory_f;
    float stem_len_energy_f;
    float stem_len_length_f;
    float stem_len_flag_avoid_muon_check_f;
    float stem_len_num_daughters_f;
    float stem_len_daughter_length_f;
    float brm_n_mu_segs_f;
    float brm_Ep_f;
    float brm_acc_length_f;
    float brm_shower_total_length_f;
    float brm_connected_length_f;
    float brm_n_size_f;
    float brm_acc_direct_length_f;
    float brm_n_shower_main_segs_f;
    float brm_n_mu_main_f;
    float lem_shower_main_length_f;
    float lem_n_3seg_f;
    float lem_e_charge_f;
    float lem_e_dQdx_f;
    float lem_shower_num_main_segs_f;
    float br4_1_shower_main_length_f;
    float br4_1_shower_total_length_f;
    float br4_1_min_dis_f;
    float br4_1_energy_f;
    float br4_1_flag_avoid_muon_check_f;
    float br4_1_n_vtx_segs_f;
    float br4_1_n_main_segs_f;
    float br4_2_ratio_45_f;
    float br4_2_ratio_35_f;
    float br4_2_ratio_25_f;
    float br4_2_ratio_15_f;
    float br4_2_ratio1_45_f;
    float br4_2_ratio1_35_f;
    float br4_2_ratio1_25_f;
    float br4_2_ratio1_15_f;
    float br4_2_iso_angle_f;
    float br4_2_iso_angle1_f;
    float br4_2_angle_f;
    float tro_3_stem_length_f;
    float tro_3_n_muon_segs_f;
    float mip_quality_energy_f;
    float mip_quality_overlap_f;
    float mip_quality_n_showers_f;
    float mip_quality_n_tracks_f;
    float mip_quality_flag_inside_pi0_f;
    float mip_quality_n_pi0_showers_f;
    float mip_quality_shortest_length_f;
    float mip_quality_acc_length_f;
    float mip_quality_shortest_angle_f;
    float mip_quality_flag_proton_f;
    float pio_mip_id_f;
    float pio_1_mass_f;
    float pio_1_pio_type_f;
    float pio_1_energy_1_f;
    float pio_1_energy_2_f;
    float pio_1_dis_1_f;
    float pio_1_dis_2_f;
    float stw_1_energy_f;
    float stw_1_dis_f;
    float stw_1_dQ_dx_f;
    float stw_1_flag_single_shower_f;
    float stw_1_n_pi0_f;
    float stw_1_num_valid_tracks_f;
    float spt_shower_main_length_f;
    float spt_shower_total_length_f;
    float spt_angle_beam_f;
    float spt_angle_vertical_f;
    float spt_max_dQ_dx_f;
    float spt_angle_beam_1_f;
    float spt_angle_drift_f;
    float spt_angle_drift_1_f;
    float spt_num_valid_tracks_f;
    float spt_n_vtx_segs_f;
    float spt_max_length_f;
    float vis_1_n_vtx_segs_f;
    float vis_1_energy_f;
    float vis_1_num_good_tracks_f;
    float vis_1_max_angle_f;
    float vis_1_max_shower_angle_f;
    float vis_1_tmp_length1_f;
    float vis_1_tmp_length2_f;
    float vis_2_n_vtx_segs_f;
    float vis_2_min_angle_f;
    float vis_2_min_weak_track_f;
    float vis_2_angle_beam_f;
    float vis_2_min_angle1_f;
    float vis_2_iso_angle1_f;
    float vis_2_min_medium_dQ_dx_f;
    float vis_2_min_length_f;
    float vis_2_sg_length_f;
    float vis_2_max_angle_f;
    float vis_2_max_weak_track_f;
    float hol_1_n_valid_tracks_f;
    float hol_1_min_angle_f;
    float hol_1_energy_f;
    float hol_1_flag_all_shower_f;
    float hol_1_min_length_f;
    float hol_2_min_angle_f;
    float hol_2_medium_dQ_dx_f;
    float hol_2_ncount_f;
    float lol_3_angle_beam_f;
    float lol_3_n_valid_tracks_f;
    float lol_3_min_angle_f;
    float lol_3_vtx_n_segs_f;
    float lol_3_shower_main_length_f;
    float lol_3_n_out_f;
    float lol_3_n_sum_f;
    float mip_vec_dQ_dx_2_f;
    float mip_vec_dQ_dx_3_f;
    float mip_vec_dQ_dx_4_f;
    float mip_vec_dQ_dx_5_f;
    float mip_vec_dQ_dx_6_f;
    float mip_vec_dQ_dx_7_f;
    float mip_vec_dQ_dx_8_f;
    float mip_vec_dQ_dx_9_f;
    float mip_vec_dQ_dx_10_f;
    float mip_vec_dQ_dx_11_f;
    float mip_vec_dQ_dx_12_f;
    float mip_vec_dQ_dx_13_f;
    float mip_vec_dQ_dx_14_f;
    float mip_vec_dQ_dx_15_f;
    float mip_vec_dQ_dx_16_f;
    float mip_vec_dQ_dx_17_f;
    float mip_vec_dQ_dx_18_f;
    float mip_vec_dQ_dx_19_f;
    float br3_3_score_f;
    float br3_5_score_f;
    float br3_6_score_f;
    float pio_2_score_f;
    float stw_2_score_f;
    float stw_3_score_f;
    float stw_4_score_f;
    float sig_1_score_f;
    float sig_2_score_f;
    float lol_1_score_f;
    float lol_2_score_f;
    float tro_1_score_f;
    float tro_2_score_f;
    float tro_4_score_f;
    float tro_5_score_f;
    //// type conversion end

    TMVA::Reader *reader = new TMVA::Reader();
    /// ordering of variables matters in TMVA

    reader->AddVariable("cme_mu_energy",&cme_mu_energy_f);
    reader->AddVariable("cme_energy",&cme_energy_f);
    reader->AddVariable("cme_mu_length",&cme_mu_length_f);
    reader->AddVariable("cme_length",&cme_length_f);
    reader->AddVariable("cme_angle_beam",&cme_angle_beam_f);
    reader->AddVariable("anc_angle",&anc_angle_f);
    reader->AddVariable("anc_max_angle",&anc_max_angle_f);
    reader->AddVariable("anc_max_length",&anc_max_length_f);
    reader->AddVariable("anc_acc_forward_length",&anc_acc_forward_length_f);
    reader->AddVariable("anc_acc_backward_length",&anc_acc_backward_length_f);
    reader->AddVariable("anc_acc_forward_length1",&anc_acc_forward_length1_f);
    reader->AddVariable("anc_shower_main_length",&anc_shower_main_length_f);
    reader->AddVariable("anc_shower_total_length",&anc_shower_total_length_f);
    reader->AddVariable("anc_flag_main_outside",&anc_flag_main_outside_f);
    reader->AddVariable("gap_flag_prolong_u",&gap_flag_prolong_u_f);
    reader->AddVariable("gap_flag_prolong_v",&gap_flag_prolong_v_f);
    reader->AddVariable("gap_flag_prolong_w",&gap_flag_prolong_w_f);
    reader->AddVariable("gap_flag_parallel",&gap_flag_parallel_f);
    reader->AddVariable("gap_n_points",&gap_n_points_f);
    reader->AddVariable("gap_n_bad",&gap_n_bad_f);
    reader->AddVariable("gap_energy",&gap_energy_f);
    reader->AddVariable("gap_num_valid_tracks",&gap_num_valid_tracks_f);
    reader->AddVariable("gap_flag_single_shower",&gap_flag_single_shower_f);
    reader->AddVariable("hol_1_n_valid_tracks",&hol_1_n_valid_tracks_f);
    reader->AddVariable("hol_1_min_angle",&hol_1_min_angle_f);
    reader->AddVariable("hol_1_energy",&hol_1_energy_f);
    reader->AddVariable("hol_1_min_length",&hol_1_min_length_f);
    reader->AddVariable("hol_2_min_angle",&hol_2_min_angle_f);
    reader->AddVariable("hol_2_medium_dQ_dx",&hol_2_medium_dQ_dx_f);
    reader->AddVariable("hol_2_ncount",&hol_2_ncount_f);
    reader->AddVariable("lol_3_angle_beam",&lol_3_angle_beam_f);
    reader->AddVariable("lol_3_n_valid_tracks",&lol_3_n_valid_tracks_f);
    reader->AddVariable("lol_3_min_angle",&lol_3_min_angle_f);
    reader->AddVariable("lol_3_vtx_n_segs",&lol_3_vtx_n_segs_f);
    reader->AddVariable("lol_3_shower_main_length",&lol_3_shower_main_length_f);
    reader->AddVariable("lol_3_n_out",&lol_3_n_out_f);
    reader->AddVariable("lol_3_n_sum",&lol_3_n_sum_f);
    reader->AddVariable("hol_1_all_shower",&hol_1_flag_all_shower_f); // naming issue
    reader->AddVariable("mgo_energy",&mgo_energy_f);
    reader->AddVariable("mgo_max_energy",&mgo_max_energy_f);
    reader->AddVariable("mgo_total_energy",&mgo_total_energy_f);
    reader->AddVariable("mgo_n_showers",&mgo_n_showers_f);
    reader->AddVariable("mgo_max_energy_1",&mgo_max_energy_1_f);
    reader->AddVariable("mgo_max_energy_2",&mgo_max_energy_2_f);
    reader->AddVariable("mgo_total_other_energy",&mgo_total_other_energy_f);
    reader->AddVariable("mgo_n_total_showers",&mgo_n_total_showers_f);
    reader->AddVariable("mgo_total_other_energy_1",&mgo_total_other_energy_1_f);
    reader->AddVariable("mgt_flag_single_shower",&mgt_flag_single_shower_f);
    reader->AddVariable("mgt_max_energy",&mgt_max_energy_f);
    reader->AddVariable("mgt_total_other_energy",&mgt_total_other_energy_f);
    reader->AddVariable("mgt_max_energy_1",&mgt_max_energy_1_f);
    reader->AddVariable("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy_f);
    reader->AddVariable("mgt_e_direct_max_energy",&mgt_e_direct_max_energy_f);
    reader->AddVariable("mgt_n_direct_showers",&mgt_n_direct_showers_f);
    reader->AddVariable("mgt_e_direct_total_energy",&mgt_e_direct_total_energy_f);
    reader->AddVariable("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio_f);
    reader->AddVariable("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy_f);
    reader->AddVariable("mip_quality_energy",&mip_quality_energy_f);
    reader->AddVariable("mip_quality_overlap",&mip_quality_overlap_f);
    reader->AddVariable("mip_quality_n_showers",&mip_quality_n_showers_f);
    reader->AddVariable("mip_quality_n_tracks",&mip_quality_n_tracks_f);
    reader->AddVariable("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0_f);
    reader->AddVariable("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers_f);
    reader->AddVariable("mip_quality_shortest_length",&mip_quality_shortest_length_f);
    reader->AddVariable("mip_quality_acc_length",&mip_quality_acc_length_f);
    reader->AddVariable("mip_quality_shortest_angle",&mip_quality_shortest_angle_f);
    reader->AddVariable("mip_quality_flag_proton",&mip_quality_flag_proton_f);
    reader->AddVariable("br1_1_shower_type",&br1_1_shower_type_f);
    reader->AddVariable("br1_1_vtx_n_segs",&br1_1_vtx_n_segs_f);
    reader->AddVariable("br1_1_energy",&br1_1_energy_f);
    reader->AddVariable("br1_1_n_segs",&br1_1_n_segs_f);
    reader->AddVariable("br1_1_flag_sg_topology",&br1_1_flag_sg_topology_f);
    reader->AddVariable("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory_f);
    reader->AddVariable("br1_1_sg_length",&br1_1_sg_length_f);
    reader->AddVariable("br1_2_n_connected",&br1_2_n_connected_f);
    reader->AddVariable("br1_2_max_length",&br1_2_max_length_f);
    reader->AddVariable("br1_2_n_connected_1",&br1_2_n_connected_1_f);
    reader->AddVariable("br1_2_n_shower_segs",&br1_2_n_shower_segs_f);
    reader->AddVariable("br1_2_max_length_ratio",&br1_2_max_length_ratio_f);
    reader->AddVariable("br1_2_shower_length",&br1_2_shower_length_f);
    reader->AddVariable("br1_3_n_connected_p",&br1_3_n_connected_p_f);
    reader->AddVariable("br1_3_max_length_p",&br1_3_max_length_p_f);
    reader->AddVariable("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs_f);
    reader->AddVariable("br3_1_energy",&br3_1_energy_f);
    reader->AddVariable("br3_1_n_shower_segments",&br3_1_n_shower_segments_f);
    reader->AddVariable("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory_f);
    reader->AddVariable("br3_1_sg_direct_length",&br3_1_sg_direct_length_f);
    reader->AddVariable("br3_1_sg_length",&br3_1_sg_length_f);
    reader->AddVariable("br3_1_total_main_length",&br3_1_total_main_length_f);
    reader->AddVariable("br3_1_total_length",&br3_1_total_length_f);
    reader->AddVariable("br3_1_iso_angle",&br3_1_iso_angle_f);
    reader->AddVariable("br3_1_sg_flag_topology",&br3_1_sg_flag_topology_f);
    reader->AddVariable("br3_2_n_ele",&br3_2_n_ele_f);
    reader->AddVariable("br3_2_n_other",&br3_2_n_other_f);
    reader->AddVariable("br3_2_other_fid",&br3_2_other_fid_f);
    reader->AddVariable("br3_4_acc_length",&br3_4_acc_length_f);
    reader->AddVariable("br3_4_total_length",&br3_4_total_length_f);
    reader->AddVariable("br3_7_min_angle",&br3_7_min_angle_f);
    reader->AddVariable("br3_8_max_dQ_dx",&br3_8_max_dQ_dx_f);
    reader->AddVariable("br3_8_n_main_segs",&br3_8_n_main_segs_f);
    reader->AddVariable("vis_1_n_vtx_segs",&vis_1_n_vtx_segs_f);
    reader->AddVariable("vis_1_energy",&vis_1_energy_f);
    reader->AddVariable("vis_1_num_good_tracks",&vis_1_num_good_tracks_f);
    reader->AddVariable("vis_1_max_angle",&vis_1_max_angle_f);
    reader->AddVariable("vis_1_max_shower_angle",&vis_1_max_shower_angle_f);
    reader->AddVariable("vis_1_tmp_length1",&vis_1_tmp_length1_f);
    reader->AddVariable("vis_1_tmp_length2",&vis_1_tmp_length2_f);
    reader->AddVariable("vis_2_n_vtx_segs",&vis_2_n_vtx_segs_f);
    reader->AddVariable("vis_2_min_angle",&vis_2_min_angle_f);
    reader->AddVariable("vis_2_min_weak_track",&vis_2_min_weak_track_f);
    reader->AddVariable("vis_2_angle_beam",&vis_2_angle_beam_f);
    reader->AddVariable("vis_2_min_angle1",&vis_2_min_angle1_f);
    reader->AddVariable("vis_2_iso_angle1",&vis_2_iso_angle1_f);
    reader->AddVariable("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx_f);
    reader->AddVariable("vis_2_min_length",&vis_2_min_length_f);
    reader->AddVariable("vis_2_sg_length",&vis_2_sg_length_f);
    reader->AddVariable("vis_2_max_angle",&vis_2_max_angle_f);
    reader->AddVariable("vis_2_max_weak_track",&vis_2_max_weak_track_f);
    reader->AddVariable("pio_1_mass",&pio_1_mass_f);
    reader->AddVariable("pio_1_pio_type",&pio_1_pio_type_f);
    reader->AddVariable("pio_1_energy_1",&pio_1_energy_1_f);
    reader->AddVariable("pio_1_energy_2",&pio_1_energy_2_f);
    reader->AddVariable("pio_1_dis_1",&pio_1_dis_1_f);
    reader->AddVariable("pio_1_dis_2",&pio_1_dis_2_f);
    reader->AddVariable("pio_mip_id",&pio_mip_id_f);
    reader->AddVariable("stem_dir_flag_single_shower",&stem_dir_flag_single_shower_f);
    reader->AddVariable("stem_dir_angle",&stem_dir_angle_f);
    reader->AddVariable("stem_dir_energy",&stem_dir_energy_f);
    reader->AddVariable("stem_dir_angle1",&stem_dir_angle1_f);
    reader->AddVariable("stem_dir_angle2",&stem_dir_angle2_f);
    reader->AddVariable("stem_dir_angle3",&stem_dir_angle3_f);
    reader->AddVariable("stem_dir_ratio",&stem_dir_ratio_f);
    reader->AddVariable("br2_num_valid_tracks",&br2_num_valid_tracks_f);
    reader->AddVariable("br2_n_shower_main_segs",&br2_n_shower_main_segs_f);
    reader->AddVariable("br2_max_angle",&br2_max_angle_f);
    reader->AddVariable("br2_sg_length",&br2_sg_length_f);
    reader->AddVariable("br2_flag_sg_trajectory",&br2_flag_sg_trajectory_f);
    reader->AddVariable("stem_len_energy",&stem_len_energy_f);
    reader->AddVariable("stem_len_length",&stem_len_length_f);
    reader->AddVariable("stem_len_flag_avoid_muon_check",&stem_len_flag_avoid_muon_check_f);
    reader->AddVariable("stem_len_num_daughters",&stem_len_num_daughters_f);
    reader->AddVariable("stem_len_daughter_length",&stem_len_daughter_length_f);
    reader->AddVariable("brm_n_mu_segs",&brm_n_mu_segs_f);
    reader->AddVariable("brm_Ep",&brm_Ep_f);
    reader->AddVariable("brm_acc_length",&brm_acc_length_f);
    reader->AddVariable("brm_shower_total_length",&brm_shower_total_length_f);
    reader->AddVariable("brm_connected_length",&brm_connected_length_f);
    reader->AddVariable("brm_n_size",&brm_n_size_f);
    reader->AddVariable("brm_n_shower_main_segs",&brm_n_shower_main_segs_f);
    reader->AddVariable("brm_n_mu_main",&brm_n_mu_main_f);
    reader->AddVariable("lem_shower_main_length",&lem_shower_main_length_f);
    reader->AddVariable("lem_n_3seg",&lem_n_3seg_f);
    reader->AddVariable("lem_e_charge",&lem_e_charge_f);
    reader->AddVariable("lem_e_dQdx",&lem_e_dQdx_f);
    reader->AddVariable("lem_shower_num_main_segs",&lem_shower_num_main_segs_f);
    reader->AddVariable("brm_nacc_direct_length",&brm_acc_direct_length_f); // naming issue
    reader->AddVariable("stw_1_energy",&stw_1_energy_f);
    reader->AddVariable("stw_1_dis",&stw_1_dis_f);
    reader->AddVariable("stw_1_dQ_dx",&stw_1_dQ_dx_f);
    reader->AddVariable("stw_1_flag_single_shower",&stw_1_flag_single_shower_f);
    reader->AddVariable("stw_1_n_pi0",&stw_1_n_pi0_f);
    reader->AddVariable("stw_1_num_valid_tracks",&stw_1_num_valid_tracks_f);
    reader->AddVariable("spt_shower_main_length",&spt_shower_main_length_f);
    reader->AddVariable("spt_shower_total_length",&spt_shower_total_length_f);
    reader->AddVariable("spt_angle_beam",&spt_angle_beam_f);
    reader->AddVariable("spt_angle_vertical",&spt_angle_vertical_f);
    reader->AddVariable("spt_max_dQ_dx",&spt_max_dQ_dx_f);
    reader->AddVariable("spt_angle_beam_1",&spt_angle_beam_1_f);
    reader->AddVariable("spt_angle_drift",&spt_angle_drift_f);
    reader->AddVariable("spt_angle_drift_1",&spt_angle_drift_1_f);
    reader->AddVariable("spt_num_valid_tracks",&spt_num_valid_tracks_f);
    reader->AddVariable("spt_n_vtx_segs",&spt_n_vtx_segs_f);
    reader->AddVariable("spt_max_length",&spt_max_length_f);
    reader->AddVariable("mip_energy",&mip_energy_f);
    reader->AddVariable("mip_n_end_reduction",&mip_n_end_reduction_f);
    reader->AddVariable("mip_n_first_mip",&mip_n_first_mip_f);
    reader->AddVariable("mip_n_first_non_mip",&mip_n_first_non_mip_f);
    reader->AddVariable("mip_n_first_non_mip_1",&mip_n_first_non_mip_1_f);
    reader->AddVariable("mip_n_first_non_mip_2",&mip_n_first_non_mip_2_f);
    reader->AddVariable("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0_f);
    reader->AddVariable("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1_f);
    reader->AddVariable("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample_f);
    reader->AddVariable("mip_n_below_threshold",&mip_n_below_threshold_f);
    reader->AddVariable("mip_n_below_zero",&mip_n_below_zero_f);
    reader->AddVariable("mip_n_lowest",&mip_n_lowest_f);
    reader->AddVariable("mip_n_highest",&mip_n_highest_f);
    reader->AddVariable("mip_lowest_dQ_dx",&mip_lowest_dQ_dx_f);
    reader->AddVariable("mip_highest_dQ_dx",&mip_highest_dQ_dx_f);
    reader->AddVariable("mip_medium_dQ_dx",&mip_medium_dQ_dx_f);
    reader->AddVariable("mip_stem_length",&mip_stem_length_f);
    reader->AddVariable("mip_length_main",&mip_length_main_f);
    reader->AddVariable("mip_length_total",&mip_length_total_f);
    reader->AddVariable("mip_angle_beam",&mip_angle_beam_f);
    reader->AddVariable("mip_iso_angle",&mip_iso_angle_f);
    reader->AddVariable("mip_n_vertex",&mip_n_vertex_f);
    reader->AddVariable("mip_n_good_tracks",&mip_n_good_tracks_f);
    reader->AddVariable("mip_E_indirect_max_energy",&mip_E_indirect_max_energy_f);
    reader->AddVariable("mip_flag_all_above",&mip_flag_all_above_f);
    reader->AddVariable("mip_min_dQ_dx_5",&mip_min_dQ_dx_5_f);
    reader->AddVariable("mip_n_other_vertex",&mip_n_other_vertex_f);
    reader->AddVariable("mip_n_stem_size",&mip_n_stem_size_f);
    reader->AddVariable("mip_flag_stem_trajectory",&mip_flag_stem_trajectory_f);
    reader->AddVariable("mip_min_dis",&mip_min_dis_f);
    reader->AddVariable("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2_f);
    reader->AddVariable("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3_f);
    reader->AddVariable("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4_f);
    reader->AddVariable("mip_vec_dQ_dx_5",&mip_vec_dQ_dx_5_f);
    reader->AddVariable("mip_vec_dQ_dx_6",&mip_vec_dQ_dx_6_f);
    reader->AddVariable("mip_vec_dQ_dx_7",&mip_vec_dQ_dx_7_f);
    reader->AddVariable("mip_vec_dQ_dx_8",&mip_vec_dQ_dx_8_f);
    reader->AddVariable("mip_vec_dQ_dx_9",&mip_vec_dQ_dx_9_f);
    reader->AddVariable("mip_vec_dQ_dx_10",&mip_vec_dQ_dx_10_f);
    reader->AddVariable("mip_vec_dQ_dx_11",&mip_vec_dQ_dx_11_f);
    reader->AddVariable("mip_vec_dQ_dx_12",&mip_vec_dQ_dx_12_f);
    reader->AddVariable("mip_vec_dQ_dx_13",&mip_vec_dQ_dx_13_f);
    reader->AddVariable("mip_vec_dQ_dx_14",&mip_vec_dQ_dx_14_f);
    reader->AddVariable("mip_vec_dQ_dx_15",&mip_vec_dQ_dx_15_f);
    reader->AddVariable("mip_vec_dQ_dx_16",&mip_vec_dQ_dx_16_f);
    reader->AddVariable("mip_vec_dQ_dx_17",&mip_vec_dQ_dx_17_f);
    reader->AddVariable("mip_vec_dQ_dx_18",&mip_vec_dQ_dx_18_f);
    reader->AddVariable("mip_vec_dQ_dx_19",&mip_vec_dQ_dx_19_f);
    reader->AddVariable("br3_3_score",&br3_3_score_f);
    reader->AddVariable("br3_5_score",&br3_5_score_f);
    reader->AddVariable("br3_6_score",&br3_6_score_f);
    reader->AddVariable("pio_2_score",&pio_2_score_f);
    reader->AddVariable("stw_2_score",&stw_2_score_f);
    reader->AddVariable("stw_3_score",&stw_3_score_f);
    reader->AddVariable("stw_4_score",&stw_4_score_f);
    reader->AddVariable("sig_1_score",&sig_1_score_f);
    reader->AddVariable("sig_2_score",&sig_2_score_f);
    reader->AddVariable("lol_1_score",&lol_1_score_f);
    reader->AddVariable("lol_2_score",&lol_2_score_f);
    reader->AddVariable("tro_1_score",&tro_1_score_f);
    reader->AddVariable("tro_2_score",&tro_2_score_f);
    reader->AddVariable("tro_4_score",&tro_4_score_f);
    reader->AddVariable("tro_5_score",&tro_5_score_f);
    reader->AddVariable("br4_1_shower_main_length",&br4_1_shower_main_length_f);
    reader->AddVariable("br4_1_shower_total_length",&br4_1_shower_total_length_f);
    reader->AddVariable("br4_1_min_dis",&br4_1_min_dis_f);
    reader->AddVariable("br4_1_energy",&br4_1_energy_f);
    reader->AddVariable("br4_1_flag_avoid_muon_check",&br4_1_flag_avoid_muon_check_f);
    reader->AddVariable("br4_1_n_vtx_segs",&br4_1_n_vtx_segs_f);
    reader->AddVariable("br4_2_ratio_45",&br4_2_ratio_45_f);
    reader->AddVariable("br4_2_ratio_35",&br4_2_ratio_35_f);
    reader->AddVariable("br4_2_ratio_25",&br4_2_ratio_25_f);
    reader->AddVariable("br4_2_ratio_15",&br4_2_ratio_15_f);
    reader->AddVariable("br4_2_ratio1_45",&br4_2_ratio1_45_f);
    reader->AddVariable("br4_2_ratio1_35",&br4_2_ratio1_35_f);
    reader->AddVariable("br4_2_ratio1_25",&br4_2_ratio1_25_f);
    reader->AddVariable("br4_2_ratio1_15",&br4_2_ratio1_15_f);
    reader->AddVariable("br4_2_iso_angle",&br4_2_iso_angle_f);
    reader->AddVariable("br4_2_iso_angle1",&br4_2_iso_angle1_f);
    reader->AddVariable("br4_2_angle",&br4_2_angle_f);
    reader->AddVariable("tro_3_stem_length",&tro_3_stem_length_f);
    reader->AddVariable("tro_3_n_muon_segs",&tro_3_n_muon_segs_f);
    reader->AddVariable("br4_1_br4_1_n_main_segs",&br4_1_n_main_segs_f); // naming issue

    //reader->BookMVA( "MyBDT", "dataset_combine/weights/BDTcombine_BDT800_3.weights.xml");
    reader->BookMVA( "MyBDT", "xgboost_set8_kaicheng.xml");

    // temp_flag = nueTag in this code
    Int_t temp_flag = 0; 

    for(int i=0; i<sigTree->GetEntries(); i++){
    //for(int i=0; i<10; i++){
        sigTree->GetEntry(i);
        //std::cout<<"Event entry: "<<i<<endl;
        if(mip_min_dis>1000) mip_min_dis = 1000.0;
        if(mip_quality_shortest_length>1000) mip_quality_shortest_length = 1000;
        if(std::isnan(mip_quality_shortest_angle)) mip_quality_shortest_angle = 0;
        if(std::isnan(stem_dir_ratio)) stem_dir_ratio = 1.0; 
        mip_energy_f = mip_energy;
        mip_n_end_reduction_f = mip_n_end_reduction;
        mip_n_first_mip_f = mip_n_first_mip;
        mip_n_first_non_mip_f = mip_n_first_non_mip;
        mip_n_first_non_mip_1_f = mip_n_first_non_mip_1;
        mip_n_first_non_mip_2_f = mip_n_first_non_mip_2;
        mip_vec_dQ_dx_0_f = mip_vec_dQ_dx_0;
        mip_vec_dQ_dx_1_f = mip_vec_dQ_dx_1;
        mip_max_dQ_dx_sample_f = mip_max_dQ_dx_sample;
        mip_n_below_threshold_f = mip_n_below_threshold;
        mip_n_below_zero_f = mip_n_below_zero;
        mip_n_lowest_f = mip_n_lowest;
        mip_n_highest_f = mip_n_highest;
        mip_lowest_dQ_dx_f = mip_lowest_dQ_dx;
        mip_highest_dQ_dx_f = mip_highest_dQ_dx;
        mip_medium_dQ_dx_f = mip_medium_dQ_dx;
        mip_stem_length_f = mip_stem_length;
        mip_length_main_f = mip_length_main;
        mip_length_total_f = mip_length_total;
        mip_angle_beam_f = mip_angle_beam;
        mip_iso_angle_f = mip_iso_angle;
        mip_n_vertex_f = mip_n_vertex;
        mip_n_good_tracks_f = mip_n_good_tracks;
        mip_E_indirect_max_energy_f = mip_E_indirect_max_energy;
        mip_flag_all_above_f = mip_flag_all_above;
        mip_min_dQ_dx_5_f = mip_min_dQ_dx_5;
        mip_n_other_vertex_f = mip_n_other_vertex;
        mip_n_stem_size_f = mip_n_stem_size;
        mip_flag_stem_trajectory_f = mip_flag_stem_trajectory;
        mip_min_dis_f = mip_min_dis;
        gap_flag_prolong_u_f = gap_flag_prolong_u;
        gap_flag_prolong_v_f = gap_flag_prolong_v;
        gap_flag_prolong_w_f = gap_flag_prolong_w;
        gap_flag_parallel_f = gap_flag_parallel;
        gap_n_points_f = gap_n_points;
        gap_n_bad_f = gap_n_bad;
        gap_energy_f = gap_energy;
        gap_num_valid_tracks_f = gap_num_valid_tracks;
        gap_flag_single_shower_f = gap_flag_single_shower;
        cme_mu_energy_f = cme_mu_energy;
        cme_energy_f = cme_energy;
        cme_mu_length_f = cme_mu_length;
        cme_length_f = cme_length;
        cme_angle_beam_f = cme_angle_beam;
        anc_angle_f = anc_angle;
        anc_max_angle_f = anc_max_angle;
        anc_max_length_f = anc_max_length;
        anc_acc_forward_length_f = anc_acc_forward_length;
        anc_acc_backward_length_f = anc_acc_backward_length;
        anc_acc_forward_length1_f = anc_acc_forward_length1;
        anc_shower_main_length_f = anc_shower_main_length;
        anc_shower_total_length_f = anc_shower_total_length;
        anc_flag_main_outside_f = anc_flag_main_outside;
        mgo_energy_f = mgo_energy;
        mgo_max_energy_f = mgo_max_energy;
        mgo_total_energy_f = mgo_total_energy;
        mgo_n_showers_f = mgo_n_showers;
        mgo_max_energy_1_f = mgo_max_energy_1;
        mgo_max_energy_2_f = mgo_max_energy_2;
        mgo_total_other_energy_f = mgo_total_other_energy;
        mgo_n_total_showers_f = mgo_n_total_showers;
        mgo_total_other_energy_1_f = mgo_total_other_energy_1;
        mgt_flag_single_shower_f = mgt_flag_single_shower;
        mgt_max_energy_f = mgt_max_energy;
        mgt_total_other_energy_f = mgt_total_other_energy;
        mgt_max_energy_1_f = mgt_max_energy_1;
        mgt_e_indirect_max_energy_f = mgt_e_indirect_max_energy;
        mgt_e_direct_max_energy_f = mgt_e_direct_max_energy;
        mgt_n_direct_showers_f = mgt_n_direct_showers;
        mgt_e_direct_total_energy_f = mgt_e_direct_total_energy;
        mgt_flag_indirect_max_pio_f = mgt_flag_indirect_max_pio;
        mgt_e_indirect_total_energy_f = mgt_e_indirect_total_energy;
        br1_1_shower_type_f = br1_1_shower_type;
        br1_1_vtx_n_segs_f = br1_1_vtx_n_segs;
        br1_1_energy_f = br1_1_energy;
        br1_1_n_segs_f = br1_1_n_segs;
        br1_1_flag_sg_topology_f = br1_1_flag_sg_topology;
        br1_1_flag_sg_trajectory_f = br1_1_flag_sg_trajectory;
        br1_1_sg_length_f = br1_1_sg_length;
        br1_2_n_connected_f = br1_2_n_connected;
        br1_2_max_length_f = br1_2_max_length;
        br1_2_n_connected_1_f = br1_2_n_connected_1;
        br1_2_n_shower_segs_f = br1_2_n_shower_segs;
        br1_2_max_length_ratio_f = br1_2_max_length_ratio;
        br1_2_shower_length_f = br1_2_shower_length;
        br1_3_n_connected_p_f = br1_3_n_connected_p;
        br1_3_max_length_p_f = br1_3_max_length_p;
        br1_3_n_shower_main_segs_f = br1_3_n_shower_main_segs;
        br3_1_energy_f = br3_1_energy;
        br3_1_n_shower_segments_f = br3_1_n_shower_segments;
        br3_1_sg_flag_trajectory_f = br3_1_sg_flag_trajectory;
        br3_1_sg_direct_length_f = br3_1_sg_direct_length;
        br3_1_sg_length_f = br3_1_sg_length;
        br3_1_total_main_length_f = br3_1_total_main_length;
        br3_1_total_length_f = br3_1_total_length;
        br3_1_iso_angle_f = br3_1_iso_angle;
        br3_1_sg_flag_topology_f = br3_1_sg_flag_topology;
        br3_2_n_ele_f = br3_2_n_ele;
        br3_2_n_other_f = br3_2_n_other;
        br3_2_other_fid_f = br3_2_other_fid;
        br3_4_acc_length_f = br3_4_acc_length;
        br3_4_total_length_f = br3_4_total_length;
        br3_7_min_angle_f = br3_7_min_angle;
        br3_8_max_dQ_dx_f = br3_8_max_dQ_dx;
        br3_8_n_main_segs_f = br3_8_n_main_segs;
        stem_dir_flag_single_shower_f = stem_dir_flag_single_shower;
        stem_dir_angle_f = stem_dir_angle;
        stem_dir_energy_f = stem_dir_energy;
        stem_dir_angle1_f = stem_dir_angle1;
        stem_dir_angle2_f = stem_dir_angle2;
        stem_dir_angle3_f = stem_dir_angle3;
        stem_dir_ratio_f = stem_dir_ratio;
        br2_num_valid_tracks_f = br2_num_valid_tracks;
        br2_n_shower_main_segs_f = br2_n_shower_main_segs;
        br2_max_angle_f = br2_max_angle;
        br2_sg_length_f = br2_sg_length;
        br2_flag_sg_trajectory_f = br2_flag_sg_trajectory;
        stem_len_energy_f = stem_len_energy;
        stem_len_length_f = stem_len_length;
        stem_len_flag_avoid_muon_check_f = stem_len_flag_avoid_muon_check;
        stem_len_num_daughters_f = stem_len_num_daughters;
        stem_len_daughter_length_f = stem_len_daughter_length;
        brm_n_mu_segs_f = brm_n_mu_segs;
        brm_Ep_f = brm_Ep;
        brm_acc_length_f = brm_acc_length;
        brm_shower_total_length_f = brm_shower_total_length;
        brm_connected_length_f = brm_connected_length;
        brm_n_size_f = brm_n_size;
        brm_acc_direct_length_f = brm_acc_direct_length;
        brm_n_shower_main_segs_f = brm_n_shower_main_segs;
        brm_n_mu_main_f = brm_n_mu_main;
        lem_shower_main_length_f = lem_shower_main_length;
        lem_n_3seg_f = lem_n_3seg;
        lem_e_charge_f = lem_e_charge;
        lem_e_dQdx_f = lem_e_dQdx;
        lem_shower_num_main_segs_f = lem_shower_num_main_segs;
        br4_1_shower_main_length_f = br4_1_shower_main_length;
        br4_1_shower_total_length_f = br4_1_shower_total_length;
        br4_1_min_dis_f = br4_1_min_dis;
        br4_1_energy_f = br4_1_energy;
        br4_1_flag_avoid_muon_check_f = br4_1_flag_avoid_muon_check;
        br4_1_n_vtx_segs_f = br4_1_n_vtx_segs;
        br4_1_n_main_segs_f = br4_1_n_main_segs;
        br4_2_ratio_45_f = br4_2_ratio_45;
        br4_2_ratio_35_f = br4_2_ratio_35;
        br4_2_ratio_25_f = br4_2_ratio_25;
        br4_2_ratio_15_f = br4_2_ratio_15;
        br4_2_ratio1_45_f = br4_2_ratio1_45;
        br4_2_ratio1_35_f = br4_2_ratio1_35;
        br4_2_ratio1_25_f = br4_2_ratio1_25;
        br4_2_ratio1_15_f = br4_2_ratio1_15;
        br4_2_iso_angle_f = br4_2_iso_angle;
        br4_2_iso_angle1_f = br4_2_iso_angle1;
        br4_2_angle_f = br4_2_angle;
        tro_3_stem_length_f = tro_3_stem_length;
        tro_3_n_muon_segs_f = tro_3_n_muon_segs;
        mip_quality_energy_f = mip_quality_energy;
        mip_quality_overlap_f = mip_quality_overlap;
        mip_quality_n_showers_f = mip_quality_n_showers;
        mip_quality_n_tracks_f = mip_quality_n_tracks;
        mip_quality_flag_inside_pi0_f = mip_quality_flag_inside_pi0;
        mip_quality_n_pi0_showers_f = mip_quality_n_pi0_showers;
        mip_quality_shortest_length_f = mip_quality_shortest_length;
        mip_quality_acc_length_f = mip_quality_acc_length;
        mip_quality_shortest_angle_f = mip_quality_shortest_angle;
        mip_quality_flag_proton_f = mip_quality_flag_proton;
        pio_mip_id_f = pio_mip_id;
        pio_1_mass_f = pio_1_mass;
        pio_1_pio_type_f = pio_1_pio_type;
        pio_1_energy_1_f = pio_1_energy_1;
        pio_1_energy_2_f = pio_1_energy_2;
        pio_1_dis_1_f = pio_1_dis_1;
        pio_1_dis_2_f = pio_1_dis_2;
        stw_1_energy_f = stw_1_energy;
        stw_1_dis_f = stw_1_dis;
        stw_1_dQ_dx_f = stw_1_dQ_dx;
        stw_1_flag_single_shower_f = stw_1_flag_single_shower;
        stw_1_n_pi0_f = stw_1_n_pi0;
        stw_1_num_valid_tracks_f = stw_1_num_valid_tracks;
        spt_shower_main_length_f = spt_shower_main_length;
        spt_shower_total_length_f = spt_shower_total_length;
        spt_angle_beam_f = spt_angle_beam;
        spt_angle_vertical_f = spt_angle_vertical;
        spt_max_dQ_dx_f = spt_max_dQ_dx;
        spt_angle_beam_1_f = spt_angle_beam_1;
        spt_angle_drift_f = spt_angle_drift;
        spt_angle_drift_1_f = spt_angle_drift_1;
        spt_num_valid_tracks_f = spt_num_valid_tracks;
        spt_n_vtx_segs_f = spt_n_vtx_segs;
        spt_max_length_f = spt_max_length;
        vis_1_n_vtx_segs_f = vis_1_n_vtx_segs;
        vis_1_energy_f = vis_1_energy;
        vis_1_num_good_tracks_f = vis_1_num_good_tracks;
        vis_1_max_angle_f = vis_1_max_angle;
        vis_1_max_shower_angle_f = vis_1_max_shower_angle;
        vis_1_tmp_length1_f = vis_1_tmp_length1;
        vis_1_tmp_length2_f = vis_1_tmp_length2;
        vis_2_n_vtx_segs_f = vis_2_n_vtx_segs;
        vis_2_min_angle_f = vis_2_min_angle;
        vis_2_min_weak_track_f = vis_2_min_weak_track;
        vis_2_angle_beam_f = vis_2_angle_beam;
        vis_2_min_angle1_f = vis_2_min_angle1;
        vis_2_iso_angle1_f = vis_2_iso_angle1;
        vis_2_min_medium_dQ_dx_f = vis_2_min_medium_dQ_dx;
        vis_2_min_length_f = vis_2_min_length;
        vis_2_sg_length_f = vis_2_sg_length;
        vis_2_max_angle_f = vis_2_max_angle;
        vis_2_max_weak_track_f = vis_2_max_weak_track;
        hol_1_n_valid_tracks_f = hol_1_n_valid_tracks;
        hol_1_min_angle_f = hol_1_min_angle;
        hol_1_energy_f = hol_1_energy;
        hol_1_flag_all_shower_f = hol_1_flag_all_shower;
        hol_1_min_length_f = hol_1_min_length;
        hol_2_min_angle_f = hol_2_min_angle;
        hol_2_medium_dQ_dx_f = hol_2_medium_dQ_dx;
        hol_2_ncount_f = hol_2_ncount;
        lol_3_angle_beam_f = lol_3_angle_beam;
        lol_3_n_valid_tracks_f = lol_3_n_valid_tracks;
        lol_3_min_angle_f = lol_3_min_angle;
        lol_3_vtx_n_segs_f = lol_3_vtx_n_segs;
        lol_3_shower_main_length_f = lol_3_shower_main_length;
        lol_3_n_out_f = lol_3_n_out;
        lol_3_n_sum_f = lol_3_n_sum;
        mip_vec_dQ_dx_2_f = mip_vec_dQ_dx_2;
        mip_vec_dQ_dx_3_f = mip_vec_dQ_dx_3;
        mip_vec_dQ_dx_4_f = mip_vec_dQ_dx_4;
        mip_vec_dQ_dx_5_f = mip_vec_dQ_dx_5;
        mip_vec_dQ_dx_6_f = mip_vec_dQ_dx_6;
        mip_vec_dQ_dx_7_f = mip_vec_dQ_dx_7;
        mip_vec_dQ_dx_8_f = mip_vec_dQ_dx_8;
        mip_vec_dQ_dx_9_f = mip_vec_dQ_dx_9;
        mip_vec_dQ_dx_10_f = mip_vec_dQ_dx_10;
        mip_vec_dQ_dx_11_f = mip_vec_dQ_dx_11;
        mip_vec_dQ_dx_12_f = mip_vec_dQ_dx_12;
        mip_vec_dQ_dx_13_f = mip_vec_dQ_dx_13;
        mip_vec_dQ_dx_14_f = mip_vec_dQ_dx_14;
        mip_vec_dQ_dx_15_f = mip_vec_dQ_dx_15;
        mip_vec_dQ_dx_16_f = mip_vec_dQ_dx_16;
        mip_vec_dQ_dx_17_f = mip_vec_dQ_dx_17;
        mip_vec_dQ_dx_18_f = mip_vec_dQ_dx_18;
        mip_vec_dQ_dx_19_f = mip_vec_dQ_dx_19;
        br3_3_score_f = br3_3_score;
        br3_5_score_f = br3_5_score;
        br3_6_score_f = br3_6_score;
        pio_2_score_f = pio_2_score;
        stw_2_score_f = stw_2_score;
        stw_3_score_f = stw_3_score;
        stw_4_score_f = stw_4_score;
        sig_1_score_f = sig_1_score;
        sig_2_score_f = sig_2_score;
        lol_1_score_f = lol_1_score;
        lol_2_score_f = lol_2_score;
        tro_1_score_f = tro_1_score;
        tro_2_score_f = tro_2_score;
        tro_4_score_f = tro_4_score;
        tro_5_score_f = tro_5_score;
        temp_flag = nueTag;
        double bdt0 = reader->EvaluateMVA("MyBDT");
        double bdt = TMath::Log10( (1+bdt0)/(1-bdt0) );
        if(bdt0 == 1) {cout<<"bdt0=1: "<<bdt<<endl;
        std::cout<<"Event entry: "<<i<<endl;}
        /* std::cout<<"BDT score: "<<bdt<<std::endl; */
        /* std::cout */
        /* <<mip_energy_f<<" " */
        /* <<mip_n_end_reduction_f<<" " */
        /* <<mip_n_first_mip_f<<" " */
        /* <<mip_n_first_non_mip_f<<" " */
        /* <<mip_n_first_non_mip_1_f<<" " */
        /* <<mip_n_first_non_mip_2_f<<" " */
        /* <<mip_vec_dQ_dx_0_f<<" " */
        /* <<mip_vec_dQ_dx_1_f<<" " */
        /* <<mip_max_dQ_dx_sample_f<<" " */
        /* <<mip_n_below_threshold_f<<" " */
        /* <<mip_n_below_zero_f<<" " */
        /* <<mip_n_lowest_f<<" " */
        /* <<mip_n_highest_f<<" " */
        /* <<mip_lowest_dQ_dx_f<<" " */
        /* <<mip_highest_dQ_dx_f<<" " */
        /* <<mip_medium_dQ_dx_f<<" " */
        /* <<mip_stem_length_f<<" " */
        /* <<mip_length_main_f<<" " */
        /* <<mip_length_total_f<<" " */
        /* <<mip_angle_beam_f<<" " */
        /* <<mip_iso_angle_f<<" " */
        /* <<mip_n_vertex_f<<" " */
        /* <<mip_n_good_tracks_f<<" " */
        /* <<mip_E_indirect_max_energy_f<<" " */
        /* <<mip_flag_all_above_f<<" " */
        /* <<mip_min_dQ_dx_5_f<<" " */
        /* <<mip_n_other_vertex_f<<" " */
        /* <<mip_n_stem_size_f<<" " */
        /* <<mip_flag_stem_trajectory_f<<" " */
        /* <<mip_min_dis_f<<" " */
        /* <<gap_flag_prolong_u_f<<" " */
        /* <<gap_flag_prolong_v_f<<" " */
        /* <<gap_flag_prolong_w_f<<" " */
        /* <<gap_flag_parallel_f<<" " */
        /* <<gap_n_points_f<<" " */
        /* <<gap_n_bad_f<<" " */
        /* <<gap_energy_f<<" " */
        /* <<gap_num_valid_tracks_f<<" " */
        /* <<gap_flag_single_shower_f<<" " */
        /* <<cme_mu_energy_f<<" " */
        /* <<cme_energy_f<<" " */
        /* <<cme_mu_length_f<<" " */
        /* <<cme_length_f<<" " */
        /* <<cme_angle_beam_f<<" " */
        /* <<anc_angle_f<<" " */
        /* <<anc_max_angle_f<<" " */
        /* <<anc_max_length_f<<" " */
        /* <<anc_acc_forward_length_f<<" " */
        /* <<anc_acc_backward_length_f<<" " */
        /* <<anc_acc_forward_length1_f<<" " */
        /* <<anc_shower_main_length_f<<" " */
        /* <<anc_shower_total_length_f<<" " */
        /* <<anc_flag_main_outside_f<<" " */
        /* <<mgo_energy_f<<" " */
        /* <<mgo_max_energy_f<<" " */
        /* <<mgo_total_energy_f<<" " */
        /* <<mgo_n_showers_f<<" " */
        /* <<mgo_max_energy_1_f<<" " */
        /* <<mgo_max_energy_2_f<<" " */
        /* <<mgo_total_other_energy_f<<" " */
        /* <<mgo_n_total_showers_f<<" " */
        /* <<mgo_total_other_energy_1_f<<" " */
        /* <<mgt_flag_single_shower_f<<" " */
        /* <<mgt_max_energy_f<<" " */
        /* <<mgt_total_other_energy_f<<" " */
        /* <<mgt_max_energy_1_f<<" " */
        /* <<mgt_e_indirect_max_energy_f<<" " */
        /* <<mgt_e_direct_max_energy_f<<" " */
        /* <<mgt_n_direct_showers_f<<" " */
        /* <<mgt_e_direct_total_energy_f<<" " */
        /* <<mgt_flag_indirect_max_pio_f<<" " */
        /* <<mgt_e_indirect_total_energy_f<<" " */
        /* <<br1_1_shower_type_f<<" " */
        /* <<br1_1_vtx_n_segs_f<<" " */
        /* <<br1_1_energy_f<<" " */
        /* <<br1_1_n_segs_f<<" " */
        /* <<br1_1_flag_sg_topology_f<<" " */
        /* <<br1_1_flag_sg_trajectory_f<<" " */
        /* <<br1_1_sg_length_f<<" " */
        /* <<br1_2_n_connected_f<<" " */
        /* <<br1_2_max_length_f<<" " */
        /* <<br1_2_n_connected_1_f<<" " */
        /* <<br1_2_n_shower_segs_f<<" " */
        /* <<br1_2_max_length_ratio_f<<" " */
        /* <<br1_2_shower_length_f<<" " */
        /* <<br1_3_n_connected_p_f<<" " */
        /* <<br1_3_max_length_p_f<<" " */
        /* <<br1_3_n_shower_main_segs_f<<" " */
        /* <<br3_1_energy_f<<" " */
        /* <<br3_1_n_shower_segments_f<<" " */
        /* <<br3_1_sg_flag_trajectory_f<<" " */
        /* <<br3_1_sg_direct_length_f<<" " */
        /* <<br3_1_sg_length_f<<" " */
        /* <<br3_1_total_main_length_f<<" " */
        /* <<br3_1_total_length_f<<" " */
        /* <<br3_1_iso_angle_f<<" " */
        /* <<br3_1_sg_flag_topology_f<<" " */
        /* <<br3_2_n_ele_f<<" " */
        /* <<br3_2_n_other_f<<" " */
        /* <<br3_2_other_fid_f<<" " */
        /* <<br3_4_acc_length_f<<" " */
        /* <<br3_4_total_length_f<<" " */
        /* <<br3_7_min_angle_f<<" " */
        /* <<br3_8_max_dQ_dx_f<<" " */
        /* <<br3_8_n_main_segs_f<<" " */
        /* <<stem_dir_flag_single_shower_f<<" " */
        /* <<stem_dir_angle_f<<" " */
        /* <<stem_dir_energy_f<<" " */
        /* <<stem_dir_angle1_f<<" " */
        /* <<stem_dir_angle2_f<<" " */
        /* <<stem_dir_angle3_f<<" " */
        /* <<stem_dir_ratio_f<<" " */
        /* <<br2_num_valid_tracks_f<<" " */
        /* <<br2_n_shower_main_segs_f<<" " */
        /* <<br2_max_angle_f<<" " */
        /* <<br2_sg_length_f<<" " */
        /* <<br2_flag_sg_trajectory_f<<" " */
        /* <<stem_len_energy_f<<" " */
        /* <<stem_len_length_f<<" " */
        /* <<stem_len_flag_avoid_muon_check_f<<" " */
        /* <<stem_len_num_daughters_f<<" " */
        /* <<stem_len_daughter_length_f<<" " */
        /* <<brm_n_mu_segs_f<<" " */
        /* <<brm_Ep_f<<" " */
        /* <<brm_acc_length_f<<" " */
        /* <<brm_shower_total_length_f<<" " */
        /* <<brm_connected_length_f<<" " */
        /* <<brm_n_size_f<<" " */
        /* <<brm_acc_direct_length_f<<" " */
        /* <<brm_n_shower_main_segs_f<<" " */
        /* <<brm_n_mu_main_f<<" " */
        /* <<lem_shower_main_length_f<<" " */
        /* <<lem_n_3seg_f<<" " */
        /* <<lem_e_charge_f<<" " */
        /* <<lem_e_dQdx_f<<" " */
        /* <<lem_shower_num_main_segs_f<<" " */
        /* <<br4_1_shower_main_length_f<<" " */
        /* <<br4_1_shower_total_length_f<<" " */
        /* <<br4_1_min_dis_f<<" " */
        /* <<br4_1_energy_f<<" " */
        /* <<br4_1_flag_avoid_muon_check_f<<" " */
        /* <<br4_1_n_vtx_segs_f<<" " */
        /* <<br4_1_n_main_segs_f<<" " */
        /* <<br4_2_ratio_45_f<<" " */
        /* <<br4_2_ratio_35_f<<" " */
        /* <<br4_2_ratio_25_f<<" " */
        /* <<br4_2_ratio_15_f<<" " */
        /* <<br4_2_ratio1_45_f<<" " */
        /* <<br4_2_ratio1_35_f<<" " */
        /* <<br4_2_ratio1_25_f<<" " */
        /* <<br4_2_ratio1_15_f<<" " */
        /* <<br4_2_iso_angle_f<<" " */
        /* <<br4_2_iso_angle1_f<<" " */
        /* <<br4_2_angle_f<<" " */
        /* <<tro_3_stem_length_f<<" " */
        /* <<tro_3_n_muon_segs_f<<" " */
        /* <<mip_quality_energy_f<<" " */
        /* <<mip_quality_overlap_f<<" " */
        /* <<mip_quality_n_showers_f<<" " */
        /* <<mip_quality_n_tracks_f<<" " */
        /* <<mip_quality_flag_inside_pi0_f<<" " */
        /* <<mip_quality_n_pi0_showers_f<<" " */
        /* <<mip_quality_shortest_length_f<<" " */
        /* <<mip_quality_acc_length_f<<" " */
        /* <<mip_quality_shortest_angle_f<<" " */
        /* <<mip_quality_flag_proton_f<<" " */
        /* <<pio_mip_id_f<<" " */
        /* <<pio_1_mass_f<<" " */
        /* <<pio_1_pio_type_f<<" " */
        /* <<pio_1_energy_1_f<<" " */
        /* <<pio_1_energy_2_f<<" " */
        /* <<pio_1_dis_1_f<<" " */
        /* <<pio_1_dis_2_f<<" " */
        /* <<stw_1_energy_f<<" " */
        /* <<stw_1_dis_f<<" " */
        /* <<stw_1_dQ_dx_f<<" " */
        /* <<stw_1_flag_single_shower_f<<" " */
        /* <<stw_1_n_pi0_f<<" " */
        /* <<stw_1_num_valid_tracks_f<<" " */
        /* <<spt_shower_main_length_f<<" " */
        /* <<spt_shower_total_length_f<<" " */
        /* <<spt_angle_beam_f<<" " */
        /* <<spt_angle_vertical_f<<" " */
        /* <<spt_max_dQ_dx_f<<" " */
        /* <<spt_angle_beam_1_f<<" " */
        /* <<spt_angle_drift_f<<" " */
        /* <<spt_angle_drift_1_f<<" " */
        /* <<spt_num_valid_tracks_f<<" " */
        /* <<spt_n_vtx_segs_f<<" " */
        /* <<spt_max_length_f<<" " */
        /* <<vis_1_n_vtx_segs_f<<" " */
        /* <<vis_1_energy_f<<" " */
        /* <<vis_1_num_good_tracks_f<<" " */
        /* <<vis_1_max_angle_f<<" " */
        /* <<vis_1_max_shower_angle_f<<" " */
        /* <<vis_1_tmp_length1_f<<" " */
        /* <<vis_1_tmp_length2_f<<" " */
        /* <<vis_2_n_vtx_segs_f<<" " */
        /* <<vis_2_min_angle_f<<" " */
        /* <<vis_2_min_weak_track_f<<" " */
        /* <<vis_2_angle_beam_f<<" " */
        /* <<vis_2_min_angle1_f<<" " */
        /* <<vis_2_iso_angle1_f<<" " */
        /* <<vis_2_min_medium_dQ_dx_f<<" " */
        /* <<vis_2_min_length_f<<" " */
        /* <<vis_2_sg_length_f<<" " */
        /* <<vis_2_max_angle_f<<" " */
        /* <<vis_2_max_weak_track_f<<" " */
        /* <<hol_1_n_valid_tracks_f<<" " */
        /* <<hol_1_min_angle_f<<" " */
        /* <<hol_1_energy_f<<" " */
        /* <<hol_1_flag_all_shower_f<<" " */
        /* <<hol_1_min_length_f<<" " */
        /* <<hol_2_min_angle_f<<" " */
        /* <<hol_2_medium_dQ_dx_f<<" " */
        /* <<hol_2_ncount_f<<" " */
        /* <<lol_3_angle_beam_f<<" " */
        /* <<lol_3_n_valid_tracks_f<<" " */
        /* <<lol_3_min_angle_f<<" " */
        /* <<lol_3_vtx_n_segs_f<<" " */
        /* <<lol_3_shower_main_length_f<<" " */
        /* <<lol_3_n_out_f<<" " */
        /* <<lol_3_n_sum_f<<" " */
        /* <<mip_vec_dQ_dx_2_f<<" " */
        /* <<mip_vec_dQ_dx_3_f<<" " */
        /* <<mip_vec_dQ_dx_4_f<<" " */
        /* <<mip_vec_dQ_dx_5_f<<" " */
        /* <<mip_vec_dQ_dx_6_f<<" " */
        /* <<mip_vec_dQ_dx_7_f<<" " */
        /* <<mip_vec_dQ_dx_8_f<<" " */
        /* <<mip_vec_dQ_dx_9_f<<" " */
        /* <<mip_vec_dQ_dx_10_f<<" " */
        /* <<mip_vec_dQ_dx_11_f<<" " */
        /* <<mip_vec_dQ_dx_12_f<<" " */
        /* <<mip_vec_dQ_dx_13_f<<" " */
        /* <<mip_vec_dQ_dx_14_f<<" " */
        /* <<mip_vec_dQ_dx_15_f<<" " */
        /* <<mip_vec_dQ_dx_16_f<<" " */
        /* <<mip_vec_dQ_dx_17_f<<" " */
        /* <<mip_vec_dQ_dx_18_f<<" " */
        /* <<mip_vec_dQ_dx_19_f<<" " */
        /* <<br3_3_score_f<<" " */
        /* <<br3_5_score_f<<" " */
        /* <<br3_6_score_f<<" " */
        /* <<pio_2_score_f<<" " */
        /* <<stw_2_score_f<<" " */
        /* <<stw_3_score_f<<" " */
        /* <<stw_4_score_f<<" " */
        /* <<sig_1_score_f<<" " */
        /* <<sig_2_score_f<<" " */
        /* <<lol_1_score_f<<" " */
        /* <<lol_2_score_f<<" " */
        /* <<tro_1_score_f<<" " */
        /* <<tro_2_score_f<<" " */
        /* <<tro_4_score_f<<" " */
        /* <<tro_5_score_f<<"\n"; */
        //// new output file
        sigbdtscore = bdt;
        sigR2->Fill();
        ////
        if(nuvtx_diff>1 || showervtx_diff>1) continue;
        hsig->Fill(bdt);
        hsigROC->Fill(bdt);

        hstotal->Fill(trueEdep);
        if(temp_flag==1) hs_cut->Fill(trueEdep);
        if(nueTag==1) hs_allcut->Fill(trueEdep);
        if(bdt>BDTcut) hs_bdt->Fill(trueEdep);

    }
    for(int i=0; i<bkgTree->GetEntries(); i++){
    //for(int i=0; i<10; i++){
        bkgTree->GetEntry(i);       
        //std::cout<<"Event entry: "<<i<<endl;
        if(mip_min_dis>1000) mip_min_dis = 1000.0;
        if(mip_quality_shortest_length>1000) mip_quality_shortest_length = 1000;
        if(std::isnan(mip_quality_shortest_angle)) mip_quality_shortest_angle = 0;
        if(std::isnan(stem_dir_ratio)) stem_dir_ratio = 1.0; 

        mip_energy_f = mip_energy;
        mip_n_end_reduction_f = mip_n_end_reduction;
        mip_n_first_mip_f = mip_n_first_mip;
        mip_n_first_non_mip_f = mip_n_first_non_mip;
        mip_n_first_non_mip_1_f = mip_n_first_non_mip_1;
        mip_n_first_non_mip_2_f = mip_n_first_non_mip_2;
        mip_vec_dQ_dx_0_f = mip_vec_dQ_dx_0;
        mip_vec_dQ_dx_1_f = mip_vec_dQ_dx_1;
        mip_max_dQ_dx_sample_f = mip_max_dQ_dx_sample;
        mip_n_below_threshold_f = mip_n_below_threshold;
        mip_n_below_zero_f = mip_n_below_zero;
        mip_n_lowest_f = mip_n_lowest;
        mip_n_highest_f = mip_n_highest;
        mip_lowest_dQ_dx_f = mip_lowest_dQ_dx;
        mip_highest_dQ_dx_f = mip_highest_dQ_dx;
        mip_medium_dQ_dx_f = mip_medium_dQ_dx;
        mip_stem_length_f = mip_stem_length;
        mip_length_main_f = mip_length_main;
        mip_length_total_f = mip_length_total;
        mip_angle_beam_f = mip_angle_beam;
        mip_iso_angle_f = mip_iso_angle;
        mip_n_vertex_f = mip_n_vertex;
        mip_n_good_tracks_f = mip_n_good_tracks;
        mip_E_indirect_max_energy_f = mip_E_indirect_max_energy;
        mip_flag_all_above_f = mip_flag_all_above;
        mip_min_dQ_dx_5_f = mip_min_dQ_dx_5;
        mip_n_other_vertex_f = mip_n_other_vertex;
        mip_n_stem_size_f = mip_n_stem_size;
        mip_flag_stem_trajectory_f = mip_flag_stem_trajectory;
        mip_min_dis_f = mip_min_dis;
        gap_flag_prolong_u_f = gap_flag_prolong_u;
        gap_flag_prolong_v_f = gap_flag_prolong_v;
        gap_flag_prolong_w_f = gap_flag_prolong_w;
        gap_flag_parallel_f = gap_flag_parallel;
        gap_n_points_f = gap_n_points;
        gap_n_bad_f = gap_n_bad;
        gap_energy_f = gap_energy;
        gap_num_valid_tracks_f = gap_num_valid_tracks;
        gap_flag_single_shower_f = gap_flag_single_shower;
        cme_mu_energy_f = cme_mu_energy;
        cme_energy_f = cme_energy;
        cme_mu_length_f = cme_mu_length;
        cme_length_f = cme_length;
        cme_angle_beam_f = cme_angle_beam;
        anc_angle_f = anc_angle;
        anc_max_angle_f = anc_max_angle;
        anc_max_length_f = anc_max_length;
        anc_acc_forward_length_f = anc_acc_forward_length;
        anc_acc_backward_length_f = anc_acc_backward_length;
        anc_acc_forward_length1_f = anc_acc_forward_length1;
        anc_shower_main_length_f = anc_shower_main_length;
        anc_shower_total_length_f = anc_shower_total_length;
        anc_flag_main_outside_f = anc_flag_main_outside;
        mgo_energy_f = mgo_energy;
        mgo_max_energy_f = mgo_max_energy;
        mgo_total_energy_f = mgo_total_energy;
        mgo_n_showers_f = mgo_n_showers;
        mgo_max_energy_1_f = mgo_max_energy_1;
        mgo_max_energy_2_f = mgo_max_energy_2;
        mgo_total_other_energy_f = mgo_total_other_energy;
        mgo_n_total_showers_f = mgo_n_total_showers;
        mgo_total_other_energy_1_f = mgo_total_other_energy_1;
        mgt_flag_single_shower_f = mgt_flag_single_shower;
        mgt_max_energy_f = mgt_max_energy;
        mgt_total_other_energy_f = mgt_total_other_energy;
        mgt_max_energy_1_f = mgt_max_energy_1;
        mgt_e_indirect_max_energy_f = mgt_e_indirect_max_energy;
        mgt_e_direct_max_energy_f = mgt_e_direct_max_energy;
        mgt_n_direct_showers_f = mgt_n_direct_showers;
        mgt_e_direct_total_energy_f = mgt_e_direct_total_energy;
        mgt_flag_indirect_max_pio_f = mgt_flag_indirect_max_pio;
        mgt_e_indirect_total_energy_f = mgt_e_indirect_total_energy;
        br1_1_shower_type_f = br1_1_shower_type;
        br1_1_vtx_n_segs_f = br1_1_vtx_n_segs;
        br1_1_energy_f = br1_1_energy;
        br1_1_n_segs_f = br1_1_n_segs;
        br1_1_flag_sg_topology_f = br1_1_flag_sg_topology;
        br1_1_flag_sg_trajectory_f = br1_1_flag_sg_trajectory;
        br1_1_sg_length_f = br1_1_sg_length;
        br1_2_n_connected_f = br1_2_n_connected;
        br1_2_max_length_f = br1_2_max_length;
        br1_2_n_connected_1_f = br1_2_n_connected_1;
        br1_2_n_shower_segs_f = br1_2_n_shower_segs;
        br1_2_max_length_ratio_f = br1_2_max_length_ratio;
        br1_2_shower_length_f = br1_2_shower_length;
        br1_3_n_connected_p_f = br1_3_n_connected_p;
        br1_3_max_length_p_f = br1_3_max_length_p;
        br1_3_n_shower_main_segs_f = br1_3_n_shower_main_segs;
        br3_1_energy_f = br3_1_energy;
        br3_1_n_shower_segments_f = br3_1_n_shower_segments;
        br3_1_sg_flag_trajectory_f = br3_1_sg_flag_trajectory;
        br3_1_sg_direct_length_f = br3_1_sg_direct_length;
        br3_1_sg_length_f = br3_1_sg_length;
        br3_1_total_main_length_f = br3_1_total_main_length;
        br3_1_total_length_f = br3_1_total_length;
        br3_1_iso_angle_f = br3_1_iso_angle;
        br3_1_sg_flag_topology_f = br3_1_sg_flag_topology;
        br3_2_n_ele_f = br3_2_n_ele;
        br3_2_n_other_f = br3_2_n_other;
        br3_2_other_fid_f = br3_2_other_fid;
        br3_4_acc_length_f = br3_4_acc_length;
        br3_4_total_length_f = br3_4_total_length;
        br3_7_min_angle_f = br3_7_min_angle;
        br3_8_max_dQ_dx_f = br3_8_max_dQ_dx;
        br3_8_n_main_segs_f = br3_8_n_main_segs;
        stem_dir_flag_single_shower_f = stem_dir_flag_single_shower;
        stem_dir_angle_f = stem_dir_angle;
        stem_dir_energy_f = stem_dir_energy;
        stem_dir_angle1_f = stem_dir_angle1;
        stem_dir_angle2_f = stem_dir_angle2;
        stem_dir_angle3_f = stem_dir_angle3;
        stem_dir_ratio_f = stem_dir_ratio;
        br2_num_valid_tracks_f = br2_num_valid_tracks;
        br2_n_shower_main_segs_f = br2_n_shower_main_segs;
        br2_max_angle_f = br2_max_angle;
        br2_sg_length_f = br2_sg_length;
        br2_flag_sg_trajectory_f = br2_flag_sg_trajectory;
        stem_len_energy_f = stem_len_energy;
        stem_len_length_f = stem_len_length;
        stem_len_flag_avoid_muon_check_f = stem_len_flag_avoid_muon_check;
        stem_len_num_daughters_f = stem_len_num_daughters;
        stem_len_daughter_length_f = stem_len_daughter_length;
        brm_n_mu_segs_f = brm_n_mu_segs;
        brm_Ep_f = brm_Ep;
        brm_acc_length_f = brm_acc_length;
        brm_shower_total_length_f = brm_shower_total_length;
        brm_connected_length_f = brm_connected_length;
        brm_n_size_f = brm_n_size;
        brm_acc_direct_length_f = brm_acc_direct_length;
        brm_n_shower_main_segs_f = brm_n_shower_main_segs;
        brm_n_mu_main_f = brm_n_mu_main;
        lem_shower_main_length_f = lem_shower_main_length;
        lem_n_3seg_f = lem_n_3seg;
        lem_e_charge_f = lem_e_charge;
        lem_e_dQdx_f = lem_e_dQdx;
        lem_shower_num_main_segs_f = lem_shower_num_main_segs;
        br4_1_shower_main_length_f = br4_1_shower_main_length;
        br4_1_shower_total_length_f = br4_1_shower_total_length;
        br4_1_min_dis_f = br4_1_min_dis;
        br4_1_energy_f = br4_1_energy;
        br4_1_flag_avoid_muon_check_f = br4_1_flag_avoid_muon_check;
        br4_1_n_vtx_segs_f = br4_1_n_vtx_segs;
        br4_1_n_main_segs_f = br4_1_n_main_segs;
        br4_2_ratio_45_f = br4_2_ratio_45;
        br4_2_ratio_35_f = br4_2_ratio_35;
        br4_2_ratio_25_f = br4_2_ratio_25;
        br4_2_ratio_15_f = br4_2_ratio_15;
        br4_2_ratio1_45_f = br4_2_ratio1_45;
        br4_2_ratio1_35_f = br4_2_ratio1_35;
        br4_2_ratio1_25_f = br4_2_ratio1_25;
        br4_2_ratio1_15_f = br4_2_ratio1_15;
        br4_2_iso_angle_f = br4_2_iso_angle;
        br4_2_iso_angle1_f = br4_2_iso_angle1;
        br4_2_angle_f = br4_2_angle;
        tro_3_stem_length_f = tro_3_stem_length;
        tro_3_n_muon_segs_f = tro_3_n_muon_segs;
        mip_quality_energy_f = mip_quality_energy;
        mip_quality_overlap_f = mip_quality_overlap;
        mip_quality_n_showers_f = mip_quality_n_showers;
        mip_quality_n_tracks_f = mip_quality_n_tracks;
        mip_quality_flag_inside_pi0_f = mip_quality_flag_inside_pi0;
        mip_quality_n_pi0_showers_f = mip_quality_n_pi0_showers;
        mip_quality_shortest_length_f = mip_quality_shortest_length;
        mip_quality_acc_length_f = mip_quality_acc_length;
        mip_quality_shortest_angle_f = mip_quality_shortest_angle;
        mip_quality_flag_proton_f = mip_quality_flag_proton;
        pio_mip_id_f = pio_mip_id;
        pio_1_mass_f = pio_1_mass;
        pio_1_pio_type_f = pio_1_pio_type;
        pio_1_energy_1_f = pio_1_energy_1;
        pio_1_energy_2_f = pio_1_energy_2;
        pio_1_dis_1_f = pio_1_dis_1;
        pio_1_dis_2_f = pio_1_dis_2;
        stw_1_energy_f = stw_1_energy;
        stw_1_dis_f = stw_1_dis;
        stw_1_dQ_dx_f = stw_1_dQ_dx;
        stw_1_flag_single_shower_f = stw_1_flag_single_shower;
        stw_1_n_pi0_f = stw_1_n_pi0;
        stw_1_num_valid_tracks_f = stw_1_num_valid_tracks;
        spt_shower_main_length_f = spt_shower_main_length;
        spt_shower_total_length_f = spt_shower_total_length;
        spt_angle_beam_f = spt_angle_beam;
        spt_angle_vertical_f = spt_angle_vertical;
        spt_max_dQ_dx_f = spt_max_dQ_dx;
        spt_angle_beam_1_f = spt_angle_beam_1;
        spt_angle_drift_f = spt_angle_drift;
        spt_angle_drift_1_f = spt_angle_drift_1;
        spt_num_valid_tracks_f = spt_num_valid_tracks;
        spt_n_vtx_segs_f = spt_n_vtx_segs;
        spt_max_length_f = spt_max_length;
        vis_1_n_vtx_segs_f = vis_1_n_vtx_segs;
        vis_1_energy_f = vis_1_energy;
        vis_1_num_good_tracks_f = vis_1_num_good_tracks;
        vis_1_max_angle_f = vis_1_max_angle;
        vis_1_max_shower_angle_f = vis_1_max_shower_angle;
        vis_1_tmp_length1_f = vis_1_tmp_length1;
        vis_1_tmp_length2_f = vis_1_tmp_length2;
        vis_2_n_vtx_segs_f = vis_2_n_vtx_segs;
        vis_2_min_angle_f = vis_2_min_angle;
        vis_2_min_weak_track_f = vis_2_min_weak_track;
        vis_2_angle_beam_f = vis_2_angle_beam;
        vis_2_min_angle1_f = vis_2_min_angle1;
        vis_2_iso_angle1_f = vis_2_iso_angle1;
        vis_2_min_medium_dQ_dx_f = vis_2_min_medium_dQ_dx;
        vis_2_min_length_f = vis_2_min_length;
        vis_2_sg_length_f = vis_2_sg_length;
        vis_2_max_angle_f = vis_2_max_angle;
        vis_2_max_weak_track_f = vis_2_max_weak_track;
        hol_1_n_valid_tracks_f = hol_1_n_valid_tracks;
        hol_1_min_angle_f = hol_1_min_angle;
        hol_1_energy_f = hol_1_energy;
        hol_1_flag_all_shower_f = hol_1_flag_all_shower;
        hol_1_min_length_f = hol_1_min_length;
        hol_2_min_angle_f = hol_2_min_angle;
        hol_2_medium_dQ_dx_f = hol_2_medium_dQ_dx;
        hol_2_ncount_f = hol_2_ncount;
        lol_3_angle_beam_f = lol_3_angle_beam;
        lol_3_n_valid_tracks_f = lol_3_n_valid_tracks;
        lol_3_min_angle_f = lol_3_min_angle;
        lol_3_vtx_n_segs_f = lol_3_vtx_n_segs;
        lol_3_shower_main_length_f = lol_3_shower_main_length;
        lol_3_n_out_f = lol_3_n_out;
        lol_3_n_sum_f = lol_3_n_sum;
        mip_vec_dQ_dx_2_f = mip_vec_dQ_dx_2;
        mip_vec_dQ_dx_3_f = mip_vec_dQ_dx_3;
        mip_vec_dQ_dx_4_f = mip_vec_dQ_dx_4;
        mip_vec_dQ_dx_5_f = mip_vec_dQ_dx_5;
        mip_vec_dQ_dx_6_f = mip_vec_dQ_dx_6;
        mip_vec_dQ_dx_7_f = mip_vec_dQ_dx_7;
        mip_vec_dQ_dx_8_f = mip_vec_dQ_dx_8;
        mip_vec_dQ_dx_9_f = mip_vec_dQ_dx_9;
        mip_vec_dQ_dx_10_f = mip_vec_dQ_dx_10;
        mip_vec_dQ_dx_11_f = mip_vec_dQ_dx_11;
        mip_vec_dQ_dx_12_f = mip_vec_dQ_dx_12;
        mip_vec_dQ_dx_13_f = mip_vec_dQ_dx_13;
        mip_vec_dQ_dx_14_f = mip_vec_dQ_dx_14;
        mip_vec_dQ_dx_15_f = mip_vec_dQ_dx_15;
        mip_vec_dQ_dx_16_f = mip_vec_dQ_dx_16;
        mip_vec_dQ_dx_17_f = mip_vec_dQ_dx_17;
        mip_vec_dQ_dx_18_f = mip_vec_dQ_dx_18;
        mip_vec_dQ_dx_19_f = mip_vec_dQ_dx_19;
        br3_3_score_f = br3_3_score;
        br3_5_score_f = br3_5_score;
        br3_6_score_f = br3_6_score;
        pio_2_score_f = pio_2_score;
        stw_2_score_f = stw_2_score;
        stw_3_score_f = stw_3_score;
        stw_4_score_f = stw_4_score;
        sig_1_score_f = sig_1_score;
        sig_2_score_f = sig_2_score;
        lol_1_score_f = lol_1_score;
        lol_2_score_f = lol_2_score;
        tro_1_score_f = tro_1_score;
        tro_2_score_f = tro_2_score;
        tro_4_score_f = tro_4_score;
        tro_5_score_f = tro_5_score;
        temp_flag = nueTag;
        double bdt0 = reader->EvaluateMVA("MyBDT");
        double bdt = TMath::Log10( (1+bdt0)/(1-bdt0) );
        if(bdt0 == -1) {cout<<"bdt0=-1: "<<bdt<<endl;
        std::cout<<"Event entry: "<<i<<endl;}
        /* std::cout<<"BDT score: "<<bdt<<std::endl; */
        /* std::cout */
        /* <<mip_energy_f<<" " */
        /* <<mip_n_end_reduction_f<<" " */
        /* <<mip_n_first_mip_f<<" " */
        /* <<mip_n_first_non_mip_f<<" " */
        /* <<mip_n_first_non_mip_1_f<<" " */
        /* <<mip_n_first_non_mip_2_f<<" " */
        /* <<mip_vec_dQ_dx_0_f<<" " */
        /* <<mip_vec_dQ_dx_1_f<<" " */
        /* <<mip_max_dQ_dx_sample_f<<" " */
        /* <<mip_n_below_threshold_f<<" " */
        /* <<mip_n_below_zero_f<<" " */
        /* <<mip_n_lowest_f<<" " */
        /* <<mip_n_highest_f<<" " */
        /* <<mip_lowest_dQ_dx_f<<" " */
        /* <<mip_highest_dQ_dx_f<<" " */
        /* <<mip_medium_dQ_dx_f<<" " */
        /* <<mip_stem_length_f<<" " */
        /* <<mip_length_main_f<<" " */
        /* <<mip_length_total_f<<" " */
        /* <<mip_angle_beam_f<<" " */
        /* <<mip_iso_angle_f<<" " */
        /* <<mip_n_vertex_f<<" " */
        /* <<mip_n_good_tracks_f<<" " */
        /* <<mip_E_indirect_max_energy_f<<" " */
        /* <<mip_flag_all_above_f<<" " */
        /* <<mip_min_dQ_dx_5_f<<" " */
        /* <<mip_n_other_vertex_f<<" " */
        /* <<mip_n_stem_size_f<<" " */
        /* <<mip_flag_stem_trajectory_f<<" " */
        /* <<mip_min_dis_f<<" " */
        /* <<gap_flag_prolong_u_f<<" " */
        /* <<gap_flag_prolong_v_f<<" " */
        /* <<gap_flag_prolong_w_f<<" " */
        /* <<gap_flag_parallel_f<<" " */
        /* <<gap_n_points_f<<" " */
        /* <<gap_n_bad_f<<" " */
        /* <<gap_energy_f<<" " */
        /* <<gap_num_valid_tracks_f<<" " */
        /* <<gap_flag_single_shower_f<<" " */
        /* <<cme_mu_energy_f<<" " */
        /* <<cme_energy_f<<" " */
        /* <<cme_mu_length_f<<" " */
        /* <<cme_length_f<<" " */
        /* <<cme_angle_beam_f<<" " */
        /* <<anc_angle_f<<" " */
        /* <<anc_max_angle_f<<" " */
        /* <<anc_max_length_f<<" " */
        /* <<anc_acc_forward_length_f<<" " */
        /* <<anc_acc_backward_length_f<<" " */
        /* <<anc_acc_forward_length1_f<<" " */
        /* <<anc_shower_main_length_f<<" " */
        /* <<anc_shower_total_length_f<<" " */
        /* <<anc_flag_main_outside_f<<" " */
        /* <<mgo_energy_f<<" " */
        /* <<mgo_max_energy_f<<" " */
        /* <<mgo_total_energy_f<<" " */
        /* <<mgo_n_showers_f<<" " */
        /* <<mgo_max_energy_1_f<<" " */
        /* <<mgo_max_energy_2_f<<" " */
        /* <<mgo_total_other_energy_f<<" " */
        /* <<mgo_n_total_showers_f<<" " */
        /* <<mgo_total_other_energy_1_f<<" " */
        /* <<mgt_flag_single_shower_f<<" " */
        /* <<mgt_max_energy_f<<" " */
        /* <<mgt_total_other_energy_f<<" " */
        /* <<mgt_max_energy_1_f<<" " */
        /* <<mgt_e_indirect_max_energy_f<<" " */
        /* <<mgt_e_direct_max_energy_f<<" " */
        /* <<mgt_n_direct_showers_f<<" " */
        /* <<mgt_e_direct_total_energy_f<<" " */
        /* <<mgt_flag_indirect_max_pio_f<<" " */
        /* <<mgt_e_indirect_total_energy_f<<" " */
        /* <<br1_1_shower_type_f<<" " */
        /* <<br1_1_vtx_n_segs_f<<" " */
        /* <<br1_1_energy_f<<" " */
        /* <<br1_1_n_segs_f<<" " */
        /* <<br1_1_flag_sg_topology_f<<" " */
        /* <<br1_1_flag_sg_trajectory_f<<" " */
        /* <<br1_1_sg_length_f<<" " */
        /* <<br1_2_n_connected_f<<" " */
        /* <<br1_2_max_length_f<<" " */
        /* <<br1_2_n_connected_1_f<<" " */
        /* <<br1_2_n_shower_segs_f<<" " */
        /* <<br1_2_max_length_ratio_f<<" " */
        /* <<br1_2_shower_length_f<<" " */
        /* <<br1_3_n_connected_p_f<<" " */
        /* <<br1_3_max_length_p_f<<" " */
        /* <<br1_3_n_shower_main_segs_f<<" " */
        /* <<br3_1_energy_f<<" " */
        /* <<br3_1_n_shower_segments_f<<" " */
        /* <<br3_1_sg_flag_trajectory_f<<" " */
        /* <<br3_1_sg_direct_length_f<<" " */
        /* <<br3_1_sg_length_f<<" " */
        /* <<br3_1_total_main_length_f<<" " */
        /* <<br3_1_total_length_f<<" " */
        /* <<br3_1_iso_angle_f<<" " */
        /* <<br3_1_sg_flag_topology_f<<" " */
        /* <<br3_2_n_ele_f<<" " */
        /* <<br3_2_n_other_f<<" " */
        /* <<br3_2_other_fid_f<<" " */
        /* <<br3_4_acc_length_f<<" " */
        /* <<br3_4_total_length_f<<" " */
        /* <<br3_7_min_angle_f<<" " */
        /* <<br3_8_max_dQ_dx_f<<" " */
        /* <<br3_8_n_main_segs_f<<" " */
        /* <<stem_dir_flag_single_shower_f<<" " */
        /* <<stem_dir_angle_f<<" " */
        /* <<stem_dir_energy_f<<" " */
        /* <<stem_dir_angle1_f<<" " */
        /* <<stem_dir_angle2_f<<" " */
        /* <<stem_dir_angle3_f<<" " */
        /* <<stem_dir_ratio_f<<" " */
        /* <<br2_num_valid_tracks_f<<" " */
        /* <<br2_n_shower_main_segs_f<<" " */
        /* <<br2_max_angle_f<<" " */
        /* <<br2_sg_length_f<<" " */
        /* <<br2_flag_sg_trajectory_f<<" " */
        /* <<stem_len_energy_f<<" " */
        /* <<stem_len_length_f<<" " */
        /* <<stem_len_flag_avoid_muon_check_f<<" " */
        /* <<stem_len_num_daughters_f<<" " */
        /* <<stem_len_daughter_length_f<<" " */
        /* <<brm_n_mu_segs_f<<" " */
        /* <<brm_Ep_f<<" " */
        /* <<brm_acc_length_f<<" " */
        /* <<brm_shower_total_length_f<<" " */
        /* <<brm_connected_length_f<<" " */
        /* <<brm_n_size_f<<" " */
        /* <<brm_acc_direct_length_f<<" " */
        /* <<brm_n_shower_main_segs_f<<" " */
        /* <<brm_n_mu_main_f<<" " */
        /* <<lem_shower_main_length_f<<" " */
        /* <<lem_n_3seg_f<<" " */
        /* <<lem_e_charge_f<<" " */
        /* <<lem_e_dQdx_f<<" " */
        /* <<lem_shower_num_main_segs_f<<" " */
        /* <<br4_1_shower_main_length_f<<" " */
        /* <<br4_1_shower_total_length_f<<" " */
        /* <<br4_1_min_dis_f<<" " */
        /* <<br4_1_energy_f<<" " */
        /* <<br4_1_flag_avoid_muon_check_f<<" " */
        /* <<br4_1_n_vtx_segs_f<<" " */
        /* <<br4_1_n_main_segs_f<<" " */
        /* <<br4_2_ratio_45_f<<" " */
        /* <<br4_2_ratio_35_f<<" " */
        /* <<br4_2_ratio_25_f<<" " */
        /* <<br4_2_ratio_15_f<<" " */
        /* <<br4_2_ratio1_45_f<<" " */
        /* <<br4_2_ratio1_35_f<<" " */
        /* <<br4_2_ratio1_25_f<<" " */
        /* <<br4_2_ratio1_15_f<<" " */
        /* <<br4_2_iso_angle_f<<" " */
        /* <<br4_2_iso_angle1_f<<" " */
        /* <<br4_2_angle_f<<" " */
        /* <<tro_3_stem_length_f<<" " */
        /* <<tro_3_n_muon_segs_f<<" " */
        /* <<mip_quality_energy_f<<" " */
        /* <<mip_quality_overlap_f<<" " */
        /* <<mip_quality_n_showers_f<<" " */
        /* <<mip_quality_n_tracks_f<<" " */
        /* <<mip_quality_flag_inside_pi0_f<<" " */
        /* <<mip_quality_n_pi0_showers_f<<" " */
        /* <<mip_quality_shortest_length_f<<" " */
        /* <<mip_quality_acc_length_f<<" " */
        /* <<mip_quality_shortest_angle_f<<" " */
        /* <<mip_quality_flag_proton_f<<" " */
        /* <<pio_mip_id_f<<" " */
        /* <<pio_1_mass_f<<" " */
        /* <<pio_1_pio_type_f<<" " */
        /* <<pio_1_energy_1_f<<" " */
        /* <<pio_1_energy_2_f<<" " */
        /* <<pio_1_dis_1_f<<" " */
        /* <<pio_1_dis_2_f<<" " */
        /* <<stw_1_energy_f<<" " */
        /* <<stw_1_dis_f<<" " */
        /* <<stw_1_dQ_dx_f<<" " */
        /* <<stw_1_flag_single_shower_f<<" " */
        /* <<stw_1_n_pi0_f<<" " */
        /* <<stw_1_num_valid_tracks_f<<" " */
        /* <<spt_shower_main_length_f<<" " */
        /* <<spt_shower_total_length_f<<" " */
        /* <<spt_angle_beam_f<<" " */
        /* <<spt_angle_vertical_f<<" " */
        /* <<spt_max_dQ_dx_f<<" " */
        /* <<spt_angle_beam_1_f<<" " */
        /* <<spt_angle_drift_f<<" " */
        /* <<spt_angle_drift_1_f<<" " */
        /* <<spt_num_valid_tracks_f<<" " */
        /* <<spt_n_vtx_segs_f<<" " */
        /* <<spt_max_length_f<<" " */
        /* <<vis_1_n_vtx_segs_f<<" " */
        /* <<vis_1_energy_f<<" " */
        /* <<vis_1_num_good_tracks_f<<" " */
        /* <<vis_1_max_angle_f<<" " */
        /* <<vis_1_max_shower_angle_f<<" " */
        /* <<vis_1_tmp_length1_f<<" " */
        /* <<vis_1_tmp_length2_f<<" " */
        /* <<vis_2_n_vtx_segs_f<<" " */
        /* <<vis_2_min_angle_f<<" " */
        /* <<vis_2_min_weak_track_f<<" " */
        /* <<vis_2_angle_beam_f<<" " */
        /* <<vis_2_min_angle1_f<<" " */
        /* <<vis_2_iso_angle1_f<<" " */
        /* <<vis_2_min_medium_dQ_dx_f<<" " */
        /* <<vis_2_min_length_f<<" " */
        /* <<vis_2_sg_length_f<<" " */
        /* <<vis_2_max_angle_f<<" " */
        /* <<vis_2_max_weak_track_f<<" " */
        /* <<hol_1_n_valid_tracks_f<<" " */
        /* <<hol_1_min_angle_f<<" " */
        /* <<hol_1_energy_f<<" " */
        /* <<hol_1_flag_all_shower_f<<" " */
        /* <<hol_1_min_length_f<<" " */
        /* <<hol_2_min_angle_f<<" " */
        /* <<hol_2_medium_dQ_dx_f<<" " */
        /* <<hol_2_ncount_f<<" " */
        /* <<lol_3_angle_beam_f<<" " */
        /* <<lol_3_n_valid_tracks_f<<" " */
        /* <<lol_3_min_angle_f<<" " */
        /* <<lol_3_vtx_n_segs_f<<" " */
        /* <<lol_3_shower_main_length_f<<" " */
        /* <<lol_3_n_out_f<<" " */
        /* <<lol_3_n_sum_f<<" " */
        /* <<mip_vec_dQ_dx_2_f<<" " */
        /* <<mip_vec_dQ_dx_3_f<<" " */
        /* <<mip_vec_dQ_dx_4_f<<" " */
        /* <<mip_vec_dQ_dx_5_f<<" " */
        /* <<mip_vec_dQ_dx_6_f<<" " */
        /* <<mip_vec_dQ_dx_7_f<<" " */
        /* <<mip_vec_dQ_dx_8_f<<" " */
        /* <<mip_vec_dQ_dx_9_f<<" " */
        /* <<mip_vec_dQ_dx_10_f<<" " */
        /* <<mip_vec_dQ_dx_11_f<<" " */
        /* <<mip_vec_dQ_dx_12_f<<" " */
        /* <<mip_vec_dQ_dx_13_f<<" " */
        /* <<mip_vec_dQ_dx_14_f<<" " */
        /* <<mip_vec_dQ_dx_15_f<<" " */
        /* <<mip_vec_dQ_dx_16_f<<" " */
        /* <<mip_vec_dQ_dx_17_f<<" " */
        /* <<mip_vec_dQ_dx_18_f<<" " */
        /* <<mip_vec_dQ_dx_19_f<<" " */
        /* <<br3_3_score_f<<" " */
        /* <<br3_5_score_f<<" " */
        /* <<br3_6_score_f<<" " */
        /* <<pio_2_score_f<<" " */
        /* <<stw_2_score_f<<" " */
        /* <<stw_3_score_f<<" " */
        /* <<stw_4_score_f<<" " */
        /* <<sig_1_score_f<<" " */
        /* <<sig_2_score_f<<" " */
        /* <<lol_1_score_f<<" " */
        /* <<lol_2_score_f<<" " */
        /* <<tro_1_score_f<<" " */
        /* <<tro_2_score_f<<" " */
        /* <<tro_4_score_f<<" " */
        /* <<tro_5_score_f<<"\n"; */
        //// new output file
        bkgbdtscore = bdt;
        bkgR2->Fill();
        ////
        if(truth_nue==1 && truth_CC==1) continue; // remove nueCC in BNB sample

        if(nueTag==0 && temp_flag==0) hbkg3->Fill(bdt); // case 3 sample
        if(nueTag==0 && temp_flag==1) hbkg2->Fill(bdt); // case 2 sample
        if(nueTag==1 && temp_flag==1) hbkg1->Fill(bdt); // case 1 sample
        hbkgROC->Fill(bdt);

        hbtotal->Fill(trueEdep);
        if(temp_flag==1) hb_cut->Fill(trueEdep);
        if(nueTag==1) hb_allcut->Fill(trueEdep);
        if(bdt>BDTcut) hb_bdt->Fill(trueEdep);

       
        /* if(bdt>BDTcut){ */
        /*     std::cout<<run<<" "<<subrun<<" "<<event<<" "<<bdt */
        /*     <<" "<<mip_energy */ 
        /*     <<" "<<mipid_score */ 
        /*     <<" "<<tro_1_score */ 
        /*     <<" "<<tro_2_score */ 
        /*     <<" "<<tro_4_score */ 
        /*     <<" "<<tro_5_score */
        /*     <<" "<<br4_tro_score */ 
        /*     <<" "<<lol_1_score */ 
        /*     <<" "<<lol_2_score */ 
        /*     <<" "<<sig_1_score */ 
        /*     <<" "<<sig_2_score */ 
        /*     <<" "<<br3_score */ 
        /*     <<" "<<br3_3_score */ 
        /*     <<" "<<br3_5_score */ 
        /*     <<" "<<br3_6_score */ 
        /*     <<" "<<gap_score */ 
        /*     <<" "<<stw_spt_score */ 
        /*     <<" "<<stw_2_score */ 
        /*     <<" "<<stw_3_score */ 
        /*     <<" "<<stw_4_score */ 
        /*     <<" "<<mgo_mgt_score */ 
        /*     <<" "<<vis_1_score */ 
        /*     <<" "<<vis_2_score */ 
        /*     <<" "<<pio_1_score */ 
        /*     <<" "<<pio_2_score */ 
        /*     <<" "<<mipquality_score */ 
        /*     <<" "<<br1_score */ 
        /*     <<" "<<stemdir_br2_score */ 
        /*     <<" "<<trimuon_score */ 
        /*     <<" "<<hol_lol_score */ 
        /*     <<" "<<cme_anc_score */ 
        /*     <<std::endl; */
        /* } */


    }
    hbkg1->Scale(hsig->Integral()/hbkg1->Integral());
    hbkg2->Scale(hsig->Integral()/hbkg2->Integral());
    hbkg3->Scale(hsig->Integral()/hbkg3->Integral());

    TGraph* roc = new TGraph(nbins);
    GetROC(hsigROC, hbkgROC, roc);

    double x0[1], y0[1];
    x0[0]=(66.11)/100;
    y0[0]=97.1/100;
    TGraph* boxcut = new TGraph(1, x0, y0);
    TCanvas* c0 = new TCanvas("c0", "ROC", 950,600);
    c0->cd();
    c0->SetGridx();
    c0->SetGridy();
    gPad->SetTicks();
    roc->GetXaxis()->SetRangeUser(-0.01,1.01);
    roc->GetXaxis()->SetNdivisions(511);
    roc->GetYaxis()->SetRangeUser(-0.01,1.01);
    roc->GetYaxis()->SetNdivisions(511);
    roc->Draw("AC");
    roc->SetLineColor(2);
    roc->SetLineWidth(4);
    roc->GetXaxis()->SetTitle("Sig Efficency");
    roc->GetYaxis()->SetTitle("Bkg Rejection");
    boxcut->Draw("Psame");
    boxcut->SetMarkerStyle(29);
    boxcut->SetMarkerSize(3.0);
    boxcut->SetMarkerColor(kRed);
    c0->SaveAs("EvalROC_type1.pdf");

    double boxeff = (66.11)/100; 
    double boxpurity = (66.11)/(66.11+78*(2.9));
    double boxeffpurity = boxeff*boxpurity;
    double sfrac = 1.0;
    double bfrac = 78; // actual weight
    double ratio = hstotal->GetEntries()/hbtotal->GetEntries(); // energy indepedent overall weight = this * bfrac
    TGraph* geff = new TGraph(nbins);
    GetEff(hsigROC, sfrac, geff);
    TGraph* gpurity = new TGraph(nbins);
    GetPurity(hsigROC, sfrac, hbkgROC, bfrac, gpurity);
    TGraph* gep = new TGraph(nbins);
    GetEffPurity(hsigROC, sfrac, hbkgROC, bfrac, gep);

    hseff_cut->Add(hs_cut);
    hseff_cut->Divide(hstotal);
    hseff_allcut->Add(hs_allcut);
    hseff_allcut->Divide(hstotal);
    hseff_bdt->Add(hs_bdt);
    hseff_bdt->Divide(hstotal);
    
    hbeff_cut->Add(hb_cut);
    hbeff_cut->Divide(hbtotal);
    hbeff_allcut->Add(hb_allcut);
    hbeff_allcut->Divide(hbtotal);
    hbeff_bdt->Add(hb_bdt);
    hbeff_bdt->Divide(hbtotal);

    /// purity weighting factor should not be energy dependent, use overall weighting factor
    for(int i=1; i<=hpurity_cut->GetNbinsX(); i++)
    {
        float sss = hs_cut->GetBinContent(i);
        float bbb = hb_cut->GetBinContent(i);
        if(sss==0 && bbb==0) hpurity_cut->SetBinContent(i, 0);
        else hpurity_cut->SetBinContent(i, sss*sfrac/(sss*sfrac+bbb*bfrac*ratio));
    }
    cout<<hs_cut->Integral()<<" "<<hb_cut->Integral()<<endl;
    for(int i=1; i<=hpurity_allcut->GetNbinsX(); i++)
    {
        float sss = hs_allcut->GetBinContent(i);
        float bbb = hb_allcut->GetBinContent(i);
        if(sss==0 && bbb==0) hpurity_allcut->SetBinContent(i, 0);
        else hpurity_allcut->SetBinContent(i, sss*sfrac/(sss*sfrac+bbb*bfrac*ratio));
    }
    cout<<hs_allcut->Integral()<<" "<<hb_allcut->Integral()<<endl;
    for(int i=1; i<=hpurity_bdt->GetNbinsX(); i++)
    {
        float sss = hs_bdt->GetBinContent(i);
        float bbb = hb_bdt->GetBinContent(i);
        if(sss==0 && bbb==0) hpurity_bdt->SetBinContent(i, 0);
        else hpurity_bdt->SetBinContent(i, sss*sfrac/(sss*sfrac+bbb*bfrac*ratio));
        //cout<<"Purity bdt: "<<sss<<" "<<bbb<<endl;
    }
    cout<<hs_bdt->Integral()<<" "<<hb_bdt->Integral()<<endl;
    cout<<hs_bdt->Integral()/hstotal->Integral()<<" "<<hs_bdt->Integral()*sfrac/(hs_bdt->Integral()*sfrac+hb_bdt->Integral()*ratio*bfrac)<<endl;


    TCanvas* c1 = new TCanvas("c1", "Cut performance", 950,600);
    c1->cd();
    gPad->SetBottomMargin(0.2);
    geff->Draw("AC");
    geff->GetXaxis()->SetRangeUser(min, max);
    geff->GetXaxis()->SetTitle("BDT score");
    geff->GetYaxis()->SetRangeUser(-0.01,1.0);
    geff->GetYaxis()->SetNdivisions(511);
    geff->SetLineColor(kBlue);
    gpurity->Draw("Csame");
    gpurity->SetLineColor(kRed);
    gep->Draw("Csame");
    gep->SetLineColor(kBlack);
    TLine* Leff = new TLine(min, boxeff, max, boxeff);
    Leff->Draw("same");
    Leff->SetLineColor(kBlue);
    Leff->SetLineStyle(kDashed);
    TLine* Lpurity = new TLine(min, boxpurity, max, boxpurity);
    Lpurity->Draw("same");
    Lpurity->SetLineColor(kRed);
    Lpurity->SetLineStyle(kDashed);
    TLine* Lep = new TLine(min, boxeffpurity, max, boxeffpurity);
    Lep->Draw("same");
    Lep->SetLineColor(kBlack);
    Lep->SetLineStyle(kDashed);
    TLegend *lg = new TLegend(0.15,0.7,0.48,0.9);
    lg->SetFillStyle(0);
    lg->AddEntry(geff, "Efficiency", "l");
    lg->AddEntry(gpurity, "Purity", "l");
    lg->AddEntry(gep, "Eff*Purity", "l");
    lg->AddEntry(Lpurity, "Cut-based purity", "l");
    lg->AddEntry(Leff, "Cut-based efficiency", "l");
    lg->AddEntry(Lep, "Cut-based efficiency*purity", "l");
    lg->Draw();
    c1->SaveAs("EvalCut_type1.pdf");

    float hmax = hsig->GetBinContent(hsig->GetMaximumBin());
    float hmax1 = hbkg1->GetBinContent(hbkg1->GetMaximumBin());
    float hmax2 = hbkg2->GetBinContent(hbkg2->GetMaximumBin());
    float hmax3 = hbkg3->GetBinContent(hbkg3->GetMaximumBin());
    if(hmax<hmax1) hmax = hmax1;
    if(hmax<hmax2) hmax = hmax2;
    if(hmax<hmax3) hmax = hmax3;

    TCanvas* c = new TCanvas("c","BDT score",950,600);
    c->cd();
    hbkg1->Draw("E1");
    hbkg1->SetLineColor(kBlack);
    hbkg1->GetYaxis()->SetRangeUser(1, hmax*1.5);
    hsig->Draw("E1 same");
    hsig->SetMarkerColor(kBlue);
    hsig->SetLineColor(kBlue);
//    hbkg2->Draw("E1 same");
//    hbkg2->SetLineColor(kGreen);
    hbkg3->Draw("E1 same");
    hbkg3->SetLineColor(kRed);
    TLegend *lg2 = new TLegend(0.15,0.7,0.48,0.9);
    lg2->SetFillStyle(0);
    lg2->AddEntry(hsig, "Signal", "l");
    lg2->AddEntry(hbkg1, "Accepted bkg by cut-based taggers", "l");
//    lg2->AddEntry(hbkg2, "Case 2 bkg", "l");
    lg2->AddEntry(hbkg3, "Rejected bkg by cut-based taggers", "l");
    lg2->Draw();
    c->SaveAs("EvalScore_type1.pdf");

    TCanvas *ca = new TCanvas("ca","Eff/Pur vs energy", 1200, 600);
    ca->Divide(2,1);
    ca->cd(1);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.15);
    hseff_cut->Draw("");
    hseff_cut->GetYaxis()->SetRangeUser(0,1.1);
    hseff_cut->GetXaxis()->SetNdivisions(505);
    hseff_cut->GetXaxis()->SetTitle("true Edep [MeV]");
    hseff_cut->GetYaxis()->SetTitle("Sig efficiency");
    hseff_cut->GetYaxis()->SetTitleOffset(1.1);
    hseff_allcut->Draw("hist same");
    hseff_allcut->SetLineColor(kBlack);
    hseff_allcut->SetLineStyle(kDashed);
    hseff_bdt->Draw("hist same");
    hseff_bdt->SetLineColor(kRed);
    TLegend *lg3 = new TLegend(0.4,0.3,0.85,0.5);
    lg3->SetFillStyle(0);
    lg3->AddEntry(hseff_cut, "Targeting tagger cut-based", "l");
    lg3->AddEntry(hseff_allcut, "All taggers", "l");
    lg3->AddEntry(hseff_bdt, "Targeting tagger BDT", "l");
    lg3->SetTextSize(0.04);
    lg3->Draw();
    ca->cd(2);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.15);
    hpurity_cut->Draw("");
    hpurity_cut->GetYaxis()->SetRangeUser(0,1.0);
    hpurity_cut->GetXaxis()->SetNdivisions(505);
    hpurity_cut->GetXaxis()->SetTitle("true Edep [MeV]");
    hpurity_cut->GetYaxis()->SetTitle("Purity");
    hpurity_cut->GetYaxis()->SetTitleOffset(1.3);
    hpurity_allcut->Draw("hist same");
    hpurity_allcut->SetLineColor(kBlack);
    hpurity_allcut->SetLineStyle(kDashed);
    hpurity_bdt->Draw("hist same");
    hpurity_bdt->SetLineColor(kRed);
    ca->SaveAs("EffPurity_type1.pdf");

    //// new output file
    round2->cd();
    sigR2->Write("", TObject::kOverwrite);
    bkgR2->Write("", TObject::kOverwrite);
    round2->Close();
    ////
    
    TFile* check = new TFile("check.root","RECREATE");
    roc->SetName("ROC");
    roc->Write();
    geff->SetName("Eff");
    geff->Write();
    gpurity->SetName("Purity");
    gpurity->Write();
    hseff_bdt->Write();
    hbeff_bdt->Write();
    hpurity_bdt->Write();
    check->Close();
}

void GetROC(TH1F* hs, TH1F* hb, TGraph* roc)
{    
    int nbins =  hs->GetNbinsX();
    double x[nbins], y[nbins];
    for(int i=0; i<nbins; i++)
    {
        //cout<<"Bin: "<<i<<endl;
        x[i]=1.0*hs->Integral(i+1, -1)/hs->Integral(0,-1);
        y[i]=hb->Integral(0, i+1)/hb->Integral(0,-1);
        //cout<<"x: "<<x[i]<<" y: "<<y[i]<<endl;
        roc->SetPoint(i, x[i], y[i]);
    }
}

void GetEffPurity(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TH1F* hb, double bfrac, TGraph* ep)
{
    int nbins =  hs->GetNbinsX();
    double x[nbins], y[nbins];
    for(int i=0; i<nbins; i++)
    {
        //cout<<"Bin: "<<i<<endl;
        x[i]=hs->GetXaxis()->GetBinCenter(i+1); // BDT score
        double eff = 1.0*hs->Integral(i+1, nbins)/hs->Integral()*sfrac;
        double bbb = 1.0*hb->Integral(i+1, nbins)/hb->Integral()*bfrac;
        y[i]=eff*eff/(eff+bbb);
        //cout<<"x: "<<x[i]<<" y: "<<y[i]<<endl;
        ep->SetPoint(i, x[i], y[i]);
    }
}

void GetEff(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TGraph* geff)
{
    int nbins =  hs->GetNbinsX();
    double x[nbins], y[nbins];
    for(int i=0; i<nbins; i++)
    {
        //cout<<"Bin: "<<i<<endl;
        x[i]=hs->GetXaxis()->GetBinCenter(i+1); // BDT score
        double eff = 1.0*hs->Integral(i+1, nbins)/hs->Integral()*sfrac;
        y[i]=eff;
        //cout<<"x: "<<x[i]<<" y: "<<y[i]<<endl;
        geff->SetPoint(i, x[i], y[i]);
    }
}

void GetPurity(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TH1F* hb, double bfrac, TGraph* gpurity)
{
    int nbins =  hs->GetNbinsX();
    double x[nbins], y[nbins];
    for(int i=0; i<nbins; i++)
    {
        //cout<<"Bin: "<<i<<endl;
        x[i]=hs->GetXaxis()->GetBinCenter(i+1); // BDT score
        double eff = 1.0*hs->Integral(i+1, nbins)/hs->Integral()*sfrac;
        double bbb = 1.0*hb->Integral(i+1, nbins)/hb->Integral()*bfrac;
        y[i]=eff/(eff+bbb);
        //cout<<"x: "<<x[i]<<" y: "<<y[i]<<endl;
        gpurity->SetPoint(i, x[i], y[i]);
    }
}


int main( int argc, char** argv )
{

    cout << "Apply xgboost BDT in ROOT TMVA" << endl;
    if(argc==6) TestEvaluate(argv[1], atof(argv[2]), atof(argv[3]), atoi(argv[4]), atof(argv[5])); // min, max, nbins
    else cout << "Incorrect input arguments!" << endl;
    return 1;

}

