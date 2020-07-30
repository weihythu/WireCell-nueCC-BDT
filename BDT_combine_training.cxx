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


void TestBDT();
void InitInput();
void InitOutput();
void InitBDT();
void TestEvaluate(double min, double max, int nbins, float BDTcut);
void GetROC(TH1F* hs, TH1F* hb, TGraph* roc);
void GetEff(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TGraph* eff);
void GetPurity(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TH1F* hb, double bfrac, TGraph* purity);
void GetEffPurity(TH1F* hs, double sfrac /*fraction of signal sample in denominator*/, TH1F* hb, double bfrac, TGraph* ep);

void TestBDT()
{
    TMVA::Tools::Instance();
    InitInput();
    InitOutput();
    InitBDT();


    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();
    
    output->Close();
    cout << "Results saved at: " << output->GetName() << std::endl;

    delete factory;
    delete dataloader;

    //TestEvaluate();

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( output->GetName() );

}

void InitBDT()
{
    factory = new TMVA::Factory( "BDTcombine", output,
        "!V:!Silent:Color:DrawProgressBar:"
        "AnalysisType=Classification" );


    dataloader = new TMVA::DataLoader("dataset_combine");
    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    //dataloader->AddVariable( "sum := var1+var2", "v1+v2", "", 'F' );
    //dataloader->AddVariable( "diff := var1-var2", "v1-v2", "", 'F' );
    //dataloader->AddVariable( "var3", "v3", "m", 'F' );
    //dataloader->AddVariable( "var4", "v4", "MeV", 'F' );
    
    dataloader->AddVariable("mip_energy","mip_energy", "", 'F');
    dataloader->AddVariable("mip_vec_dQ_dx_0","mip_vec_dQ_dx_0", "", 'F');
    dataloader->AddVariable("mip_vec_dQ_dx_1","mip_vec_dQ_dx_1", "", 'F');
    dataloader->AddVariable("mip_vec_dQ_dx_2","mip_vec_dQ_dx_2", "", 'F');
    dataloader->AddVariable("mip_vec_dQ_dx_3","mip_vec_dQ_dx_3", "", 'F');
    dataloader->AddVariable("mip_vec_dQ_dx_4","mip_vec_dQ_dx_4", "", 'F');
    /* dataloader->AddVariable("spt_angle_beam","spt_angle_beam", "", 'F'); */ 
    /* dataloader->AddVariable("spt_angle_drift","spt_angle_drift", "", 'F'); */ 
    /* dataloader->AddVariable("cme_mu_energy","cme_mu_energy", "", 'F'); */ 
    /* dataloader->AddVariable("cme_mu_length","cme_mu_length", "", 'F'); */ 
    /* dataloader->AddVariable("cme_length","cme_length", "", 'F'); */ 


    dataloader->AddVariable("mipid_score", "mipid_score", "", 'F');
    dataloader->AddVariable("gap_score", "gap_score", "", 'F');
    dataloader->AddVariable("hol_lol_score", "hol_lol_score", "", 'F');
    dataloader->AddVariable("cme_anc_score", "cme_anc_score", "", 'F');
    dataloader->AddVariable("mgo_mgt_score", "mgo_mgt_score", "", 'F');
    dataloader->AddVariable("br1_score", "br1_score", "", 'F');
    dataloader->AddVariable("br3_score", "br3_score", "", 'F');
    dataloader->AddVariable("br3_3_score", "br3_3_score", "", 'F');
    dataloader->AddVariable("br3_5_score", "br3_5_score", "", 'F');
    dataloader->AddVariable("br3_6_score", "br3_6_score", "", 'F');
    dataloader->AddVariable("stemdir_br2_score", "stemdir_br2_score", "", 'F');
    dataloader->AddVariable("trimuon_score", "trimuon_score", "", 'F');
    dataloader->AddVariable("br4_tro_score", "br4_tro_score", "", 'F');
    dataloader->AddVariable("mipquality_score", "mipquality_score", "", 'F');
    dataloader->AddVariable("pio_1_score", "pio_1_score", "", 'F');
    dataloader->AddVariable("pio_2_score", "pio_2_score", "", 'F');
    dataloader->AddVariable("stw_spt_score", "stw_spt_score", "", 'F');
    dataloader->AddVariable("vis_1_score", "vis_1_score", "", 'F');
    dataloader->AddVariable("vis_2_score", "vis_2_score", "", 'F');
    dataloader->AddVariable("stw_2_score", "stw_2_score", "", 'F');
    dataloader->AddVariable("stw_3_score", "stw_3_score", "", 'F');
    dataloader->AddVariable("stw_4_score", "stw_4_score", "", 'F');
    dataloader->AddVariable("sig_1_score", "sig_1_score", "", 'F');
    dataloader->AddVariable("sig_2_score", "sig_2_score", "", 'F');
    dataloader->AddVariable("lol_1_score", "lol_1_score", "", 'F');
    dataloader->AddVariable("lol_2_score", "lol_2_score", "", 'F');
    dataloader->AddVariable("tro_1_score", "tro_1_score", "", 'F');
    dataloader->AddVariable("tro_2_score", "tro_2_score", "", 'F');
    dataloader->AddVariable("tro_4_score", "tro_4_score", "", 'F');
    dataloader->AddVariable("tro_5_score", "tro_5_score", "", 'F');

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 78*40000./20000.);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "lowEweight" );
    //dataloader->SetBackgroundWeightExpression( "lowEweight" );
    //dataloader->SetBackgroundWeightExpression( "1.0+nueTag*10.0" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycut_b = "truth_nue==0 || truth_CC==0"; // for example: TCut mycutb = "abs(var1)<0.5";
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=40000:"
        "nTrain_Background=20000:"
        // "nTest_Signal=2000:"
        // "nTest_Background=2000:"
        "SplitMode=Random:"
        "NormMode=None:"
        //"NormMode=NumEvents:" // norm to nTrain_Signal/Background numbers
        //"NormMode=EqualNumEvents:" // norm to nTrain_Signal/Background numbers
        "!V" );

    // variations of BDTs
    // BDT: uses Adaptive Boost
    // BDTG: uses Gradient Boost
    // BDTB: uses Bagging
    // BDTD: decorrelation + Adaptive Boost
    // BDTF: allow usage of fisher discriminant for node splitting
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT800_3",
        "!H:!V:"
        "VarTransform=D,P,U,G,D,N:"
        "NTrees=800:"
        "MinNodeSize=2.5%:"
        "MaxDepth=3:"
        "BoostType=AdaBoost:"
        "AdaBoostBeta=0.5:"
        "UseBaggedBoost:"
        "BaggedSampleFraction=0.5:"
        "SeparationType=GiniIndex:"
        "nCuts=20");

    /* factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT100_20", */
    /*     "!H:!V:" */
    /*     "VarTransform=D,P,U,G,D,N:" */
    /*     "NTrees=100:" */
    /*     "MinNodeSize=2.5%:" */
    /*     "MaxDepth=20:" */
    /*     "BoostType=AdaBoost:" */
    /*     "AdaBoostBeta=0.5:" */
    /*     "UseBaggedBoost:" */
    /*     "BaggedSampleFraction=0.5:" */
    /*     "SeparationType=GiniIndex:" */
    /*     "nCuts=20"); */

    /* factory->BookMethod(dataloader, TMVA::Types::kFisher, "Fisher", */
    /*     "!H:!V:" */
    /*     "VarTransform=D,P,U,G,D,N:"); */
}


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
    
  Float_t mip_energy = 0;
  Float_t mip_vec_dQ_dx_0 = 0;
  Float_t mip_vec_dQ_dx_1 = 0;
  Float_t mip_vec_dQ_dx_2 = 0;
  Float_t mip_vec_dQ_dx_3 = 0;
  Float_t mip_vec_dQ_dx_4 = 0;
  Float_t spt_angle_beam = 0; 
  Float_t spt_angle_drift = 0; 
  Float_t cme_mu_energy = 0; 
  Float_t cme_mu_length = 0; 
  Float_t cme_length = 0; 

  Float_t mipid_score = 0;
  Float_t gap_score = 0;
  Float_t hol_lol_score = 0;
  Float_t cme_anc_score = 0;
  Float_t mgo_mgt_score = 0;
  Float_t br1_score = 0;
  Float_t br3_score = 0;
  Float_t br3_3_score = 0;
  Float_t br3_5_score = 0;
  Float_t br3_6_score = 0;
  Float_t stemdir_br2_score = 0;
  Float_t trimuon_score = 0;
  Float_t br4_tro_score = 0;
  Float_t mipquality_score = 0;
  Float_t pio_1_score = 0;
  Float_t pio_2_score = 0;
  Float_t stw_spt_score = 0;
  Float_t vis_1_score = 0;
  Float_t vis_2_score = 0;
  Float_t stw_2_score = 0;
  Float_t stw_3_score = 0;
  Float_t stw_4_score = 0;
  Float_t sig_1_score = 0;
  Float_t sig_2_score = 0;
  Float_t lol_1_score = 0;
  Float_t lol_2_score = 0;
  Float_t tro_1_score = 0;
  Float_t tro_2_score = 0;
  Float_t tro_4_score = 0;
  Float_t tro_5_score = 0;

  Int_t temp_flag = 1; // current integrated tagger cut-based result
  Int_t rest_flag = 1; // the rest tagger cut-based result

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
    sigTree->SetBranchAddress("temp_flag",&temp_flag);
    
    sigTree->SetBranchAddress("mip_energy", &mip_energy);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3);
    sigTree->SetBranchAddress("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4);
    sigTree->SetBranchAddress("spt_angle_beam",&spt_angle_beam); 
    sigTree->SetBranchAddress("spt_angle_drift",&spt_angle_drift); 
    sigTree->SetBranchAddress("cme_mu_energy",&cme_mu_energy); 
    sigTree->SetBranchAddress("cme_mu_length",&cme_mu_length); 
    sigTree->SetBranchAddress("cme_length",&cme_length); 
    sigTree->SetBranchAddress("mipid_score", &mipid_score);
    sigTree->SetBranchAddress("gap_score", &gap_score);
    sigTree->SetBranchAddress("hol_lol_score", &hol_lol_score);
    sigTree->SetBranchAddress("cme_anc_score", &cme_anc_score);
    sigTree->SetBranchAddress("mgo_mgt_score", &mgo_mgt_score);
    sigTree->SetBranchAddress("br1_score", &br1_score);
    sigTree->SetBranchAddress("br3_score", &br3_score);
    sigTree->SetBranchAddress("br3_3_score", &br3_3_score);
    sigTree->SetBranchAddress("br3_5_score", &br3_5_score);
    sigTree->SetBranchAddress("br3_6_score", &br3_6_score);
    sigTree->SetBranchAddress("stemdir_br2_score", &stemdir_br2_score);
    sigTree->SetBranchAddress("trimuon_score", &trimuon_score);
    sigTree->SetBranchAddress("br4_tro_score", &br4_tro_score);
    sigTree->SetBranchAddress("mipquality_score", &mipquality_score);
    sigTree->SetBranchAddress("pio_1_score", &pio_1_score);
    sigTree->SetBranchAddress("pio_2_score", &pio_2_score);
    sigTree->SetBranchAddress("stw_spt_score", &stw_spt_score);
    sigTree->SetBranchAddress("vis_1_score", &vis_1_score);
    sigTree->SetBranchAddress("vis_2_score", &vis_2_score);
    sigTree->SetBranchAddress("stw_2_score", &stw_2_score);
    sigTree->SetBranchAddress("stw_3_score", &stw_3_score);
    sigTree->SetBranchAddress("stw_4_score", &stw_4_score);
    sigTree->SetBranchAddress("sig_1_score", &sig_1_score);
    sigTree->SetBranchAddress("sig_2_score", &sig_2_score);
    sigTree->SetBranchAddress("lol_1_score", &lol_1_score);
    sigTree->SetBranchAddress("lol_2_score", &lol_2_score);
    sigTree->SetBranchAddress("tro_1_score", &tro_1_score);
    sigTree->SetBranchAddress("tro_2_score", &tro_2_score);
    sigTree->SetBranchAddress("tro_4_score", &tro_4_score);
    sigTree->SetBranchAddress("tro_5_score", &tro_5_score);

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
    bkgTree->SetBranchAddress("temp_flag",&temp_flag);
    
    bkgTree->SetBranchAddress("mip_energy", &mip_energy);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3);
    bkgTree->SetBranchAddress("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4);
    bkgTree->SetBranchAddress("spt_angle_beam",&spt_angle_beam); 
    bkgTree->SetBranchAddress("spt_angle_drift",&spt_angle_drift); 
    bkgTree->SetBranchAddress("cme_mu_energy",&cme_mu_energy); 
    bkgTree->SetBranchAddress("cme_mu_length",&cme_mu_length); 
    bkgTree->SetBranchAddress("cme_length",&cme_length); 
    bkgTree->SetBranchAddress("mipid_score", &mipid_score);
    bkgTree->SetBranchAddress("gap_score", &gap_score);
    bkgTree->SetBranchAddress("hol_lol_score", &hol_lol_score);
    bkgTree->SetBranchAddress("cme_anc_score", &cme_anc_score);
    bkgTree->SetBranchAddress("mgo_mgt_score", &mgo_mgt_score);
    bkgTree->SetBranchAddress("br1_score", &br1_score);
    bkgTree->SetBranchAddress("br3_score", &br3_score);
    bkgTree->SetBranchAddress("br3_3_score", &br3_3_score);
    bkgTree->SetBranchAddress("br3_5_score", &br3_5_score);
    bkgTree->SetBranchAddress("br3_6_score", &br3_6_score);
    bkgTree->SetBranchAddress("stemdir_br2_score", &stemdir_br2_score);
    bkgTree->SetBranchAddress("trimuon_score", &trimuon_score);
    bkgTree->SetBranchAddress("br4_tro_score", &br4_tro_score);
    bkgTree->SetBranchAddress("mipquality_score", &mipquality_score);
    bkgTree->SetBranchAddress("pio_1_score", &pio_1_score);
    bkgTree->SetBranchAddress("pio_2_score", &pio_2_score);
    bkgTree->SetBranchAddress("stw_spt_score", &stw_spt_score);
    bkgTree->SetBranchAddress("vis_1_score", &vis_1_score);
    bkgTree->SetBranchAddress("vis_2_score", &vis_2_score);
    bkgTree->SetBranchAddress("stw_2_score", &stw_2_score);
    bkgTree->SetBranchAddress("stw_3_score", &stw_3_score);
    bkgTree->SetBranchAddress("stw_4_score", &stw_4_score);
    bkgTree->SetBranchAddress("sig_1_score", &sig_1_score);
    bkgTree->SetBranchAddress("sig_2_score", &sig_2_score);
    bkgTree->SetBranchAddress("lol_1_score", &lol_1_score);
    bkgTree->SetBranchAddress("lol_2_score", &lol_2_score);
    bkgTree->SetBranchAddress("tro_1_score", &tro_1_score);
    bkgTree->SetBranchAddress("tro_2_score", &tro_2_score);
    bkgTree->SetBranchAddress("tro_4_score", &tro_4_score);
    bkgTree->SetBranchAddress("tro_5_score", &tro_5_score);

    /// new output file
    TFile* round2 = new TFile("BDT_all.root", "RECREATE");
    TTree* sigR2 = sigTree->CloneTree(0);
    Float_t sigbdtscore;
    sigR2->Branch("bdt_all_score", &sigbdtscore, "bdt_all_score/F");
    TTree* bkgR2 = bkgTree->CloneTree(0);
    Float_t bkgbdtscore;
    bkgR2->Branch("bdt_all_score", &bkgbdtscore, "bdt_all_score/F");
    /////

    TMVA::Reader *reader = new TMVA::Reader();

    reader->AddVariable("mip_energy",&mip_energy);
    reader->AddVariable("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
    reader->AddVariable("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
    reader->AddVariable("mip_vec_dQ_dx_2",&mip_vec_dQ_dx_2);
    reader->AddVariable("mip_vec_dQ_dx_3",&mip_vec_dQ_dx_3);
    reader->AddVariable("mip_vec_dQ_dx_4",&mip_vec_dQ_dx_4);
    /* reader->AddVariable("spt_angle_beam",&spt_angle_beam); */ 
    /* reader->AddVariable("spt_angle_drift",&spt_angle_drift); */ 
    /* reader->AddVariable("cme_mu_energy",&cme_mu_energy); */ 
    /* reader->AddVariable("cme_mu_length",&cme_mu_length); */ 
    /* reader->AddVariable("cme_length",&cme_length); */ 
    reader->AddVariable("mipid_score",&mipid_score);
    reader->AddVariable("gap_score",&gap_score);
    reader->AddVariable("hol_lol_score",&hol_lol_score);
    reader->AddVariable("cme_anc_score",&cme_anc_score);
    reader->AddVariable("mgo_mgt_score",&mgo_mgt_score);
    reader->AddVariable("br1_score",&br1_score);
    reader->AddVariable("br3_score",&br3_score);
    reader->AddVariable("br3_3_score",&br3_3_score);
    reader->AddVariable("br3_5_score",&br3_5_score);
    reader->AddVariable("br3_6_score",&br3_6_score);
    reader->AddVariable("stemdir_br2_score",&stemdir_br2_score);
    reader->AddVariable("trimuon_score",&trimuon_score);
    reader->AddVariable("br4_tro_score",&br4_tro_score);
    reader->AddVariable("mipquality_score",&mipquality_score);
    reader->AddVariable("pio_1_score",&pio_1_score);
    reader->AddVariable("pio_2_score",&pio_2_score);
    reader->AddVariable("stw_spt_score",&stw_spt_score);
    reader->AddVariable("vis_1_score",&vis_1_score);
    reader->AddVariable("vis_2_score",&vis_2_score);
    reader->AddVariable("stw_2_score",&stw_2_score);
    reader->AddVariable("stw_3_score",&stw_3_score);
    reader->AddVariable("stw_4_score",&stw_4_score);
    reader->AddVariable("sig_1_score",&sig_1_score);
    reader->AddVariable("sig_2_score",&sig_2_score);
    reader->AddVariable("lol_1_score",&lol_1_score);
    reader->AddVariable("lol_2_score",&lol_2_score);
    reader->AddVariable("tro_1_score",&tro_1_score);
    reader->AddVariable("tro_2_score",&tro_2_score);
    reader->AddVariable("tro_4_score",&tro_4_score);
    reader->AddVariable("tro_5_score",&tro_5_score);

    reader->BookMVA( "MyBDT", "dataset_combine/weights/BDTcombine_BDT800_3.weights.xml");

    for(int i=0; i<sigTree->GetEntries(); i++){
        sigTree->GetEntry(i); 
        double bdt = reader->EvaluateMVA("MyBDT");
        //// new output file
        sigbdtscore = bdt;
        sigR2->Fill();
        ////
        if(nuvtx_diff>=1 || showervtx_diff>=1) continue;
        hsig->Fill(bdt);
        hsigROC->Fill(bdt);

        hstotal->Fill(trueEdep);
        if(temp_flag==1) hs_cut->Fill(trueEdep);
        if(nueTag==1) hs_allcut->Fill(trueEdep);
        if(bdt>BDTcut) hs_bdt->Fill(trueEdep);

    }
    for(int i=0; i<bkgTree->GetEntries(); i++){
        bkgTree->GetEntry(i);        
        double bdt = reader->EvaluateMVA("MyBDT");
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
    hbkg1->GetYaxis()->SetRangeUser(0, hmax*1.5);
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
}

void GetROC(TH1F* hs, TH1F* hb, TGraph* roc)
{    
    int nbins =  hs->GetNbinsX();
    double x[nbins], y[nbins];
    for(int i=0; i<nbins; i++)
    {
        //cout<<"Bin: "<<i<<endl;
        x[i]=1.0*hs->Integral(i+1, nbins)/hs->Integral();
        y[i]=1.0-1.0*hb->Integral(i+1, nbins)/hb->Integral();
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

void InitInput()
{
    TString fname = "bdtscore_combine_training.root";
    if (!gSystem->AccessPathName( fname )) {
        input = TFile::Open( fname ); // check if file in local directory exists
    }
    else {
        TFile::SetCacheFileDir(".");
        input = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD");
    }
    if (!input) {
        std::cout << "ERROR: could not open data file" << std::endl;
        exit(1);
    }
    cout << "Using input file: " << input->GetName() << std::endl;
}

void InitOutput()
{
    TString outfileName( "TMVA.root" );
    output = TFile::Open( outfileName, "RECREATE" );
}


int main( int argc, char** argv )
{

    cout << "testing BDT in ROOT TMVA" << endl;
    if(argc==1) TestBDT();
    else if(argc==6) TestEvaluate(argv[1], atof(argv[2]), atof(argv[3]), atoi(argv[4]), atof(argv[5])); // min, max, nbins
    else cout << "Incorrect input arguments!" << endl;
    return 1;

}

