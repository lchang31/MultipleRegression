clear all; clc; close all;
basepath = pwd;
fmripath = basepath;
mask = fullfile(basepath,'avg152T1_gray_3mm_summer_binary.nii'); % where is your mask (mask that I sent you) and what is the name?;

config=cosmo_config();
for subj=18:18
data_path = [fmripath '/summer_movie_s' num2str(subj) '.nii']; % name of the file
data = cosmo_fmri_dataset(data_path,'mask',mask); %read nii file
data_clean=cosmo_remove_useless_data(data); % clear
save([fmripath '/summer_movie_s' num2str(subj) '.mat'], 'data_clean', '-v7.3'); % save
end

