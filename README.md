 Wire-Cell-nueCC-BDT-TMVA
 * These are initial codes for integrating 30 pre-trained BDTs into a single BDT to do a final decision-making.
 * Each pre-trained BDT provides a BDT score as input to the final single BDT.
 * Some events need to be assigned a default BDT value, e.g. empty vector, or un-filled scalar variables (default value = 0).
 * "preprocess" is to prepare BDT scores from each pretrained BDT. 
 * "BDT_combine_training" is a regular TMVA code to train an evaluate a single BDT based on output from "preprocess". 
