# MultipleRegression

Files: 
Summer 
  SummerBeta.m: generates the beta values
  SummerNeural.m: generates the neural dsms, provide visualization of 2nd order RDM (MDS plots)
  SignificanceSummer.m: uses randomize_r.m, performs permutation test 
  saving_nii_to_mat.m: converts the .nii fmri file to .mat files to be able to crop
  
Sherlock: (files similar structure as Summer) 
  SherlockBeta.m
  SherlockNeural.m
  SignificanceSherlock.m
  Cropnii.m: uses sherlock_avg.m to convert to 3TR
  
CorrelationPlot.R: generates the correlation plot for the feature and neural RDMs
  Obtain correlation information from Sherlock_with action.xlsx and Summer_action.xlsx
  
featureRDM.m: generate feature RDMs and provide visualization

PublicationR.R: creates the beta value graph
Pub_graphs.xlsx: contains beta values, p values, standard deviation for all ROIs and features for both Sherlock and Summer
