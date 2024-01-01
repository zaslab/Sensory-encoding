if (!require("reticulate")) install.packages("reticulate")
if (!require("rmatio")) install.packages("rmatio")
library(reticulate)
library(rmatio)
source_python('G:/Google Drive/osm6 Manuscript/Code/ClassifierFunctions.py')
source_python('G:/Google Drive/osm6 Manuscript/Code/LinClassifierPython.py')
py_run_file("G:/Google Drive/osm6 Manuscript/Code/ClassifierFunctions.py")
py_run_file("G:/Google Drive/osm6 Manuscript/Code/LinClassifierPython.py")

RunLinClassifier = function(activations,trace_feats,vars,include_acts, include_traces,all_stim,best){
  all_confusions = list()
  
  feat_num = dim(trace_feats)[2]/11 
  if(best | !include_traces){
    feat_run = 1
  }else{
    feat_run = dim(trace_feats)[2]/11 
  }
  
  best_feats = c(1 + c(0:10)*feat_num, 2 + c(0:10)*feat_num, 3 + c(0:10)*feat_num)
  
  for (feat in c(1:feat_run)){
    
    if(best){
      neuron_feature = best_feats
    }else{
      neuron_feature = feat + c(0:10)*feat_num
    }
    all_repeats = list()
    
    data_list = GetDataIdx(activations, trace_feats, vars, all_stim, samples, include_traces,include_acts,neuron_feature)
    
    for (i in c(1:10)){  
      t = Test(RandomForestClassifier, dict(max_depth=as.integer(5), n_estimators=as.integer(200), max_features='log2'))
      t$train_and_test(data_list$wt_data, data_list$wt_y)
      all_repeats[[i]] = t$confusion_matrix
      if( i == 1){
        confusion_sum = t$confusion_matrix
      }else{
        confusion_sum = confusion_sum + t$confusion_matrix
      }
    }
    all_confusions[[feat]] = all_repeats
    message(feat)
    message(sum(diag(confusion_sum))/ sum(confusion_sum))
    
  }
  return(all_confusions)
}

GetDataIdx = function(activations, trace_feats, vars, all_stim, samples, include_traces,include_acts,neuron_feature){
  
  if (all_stim == 1){
    conc_idx = 1:samples
    wt_y = vars[conc_idx, 2]
  }else{
    start_idx = (conc - 1) * (samples/3) + 1
    end_idx = start_idx + (samples/3) - 1
    conc_idx = c(start_idx : end_idx)
    wt_y = vars[conc_idx, 1]
  }
  
  if (include_traces){
    if(include_acts){
      wt_data = cbind(activations[conc_idx, 1:11], trace_feats[conc_idx, neuron_feature]) 
    } else{  
      wt_data =   trace_feats[conc_idx, neuron_feature]
    }}else{
      wt_data =   activations[conc_idx, 1:11]
    }
  return(list(wt_y = wt_y, wt_data = wt_data))
  }

traces = read.mat("G:/Google Drive/MultiK/All matrices/DataForStatistics/WithFeatures/aravi_best_features_concatenated.mat")[[1]]
all_PCs = read.mat("G:/Google Drive/MultiK/All matrices/DataForStatistics/WithFeatures/aravi_all_PC_features_concatenated.mat")[[1]]
vars = read.mat("G:/Google Drive/MultiK/All matrices/DataForStatistics/WithFeatures/aravi_activations_by_worm_identity.mat")[[1]]
activations = read.mat("G:/Google Drive/MultiK/All matrices/DataForStatistics/WithFeatures/aravi_activations_by_worm.mat")[[1]]
samples = dim(vars)[1]

all_stim = 0
include_traces = 1
include_acts = 1
best = 3

trace_feats = all_PCs

all_confusions = RunLinClassifier(activations,trace_feats,vars,include_acts, include_traces,all_stim,best)

names(all_confusions) = letters[1:length(all_confusions)]
write.mat(all_confusions,"G:/Google Drive/MultiK/All matrices/DataForStatistics/WithFeatures/SingleFeatures/synth_aravi_106_on_3bestraw.mat")


