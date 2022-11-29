%% Performs the Multiple Regression
% Finds the set of beta coefficients for each social feature for each mask
% 1/21/20, Lucy Chang
%%
load('featureRDM.mat')

%calculate the distance and normalize the values for the features
speaking_self= normalize(squareform(speaking_self));
speaking_others= normalize(squareform(speaking_others));
speaking_things=normalize(squareform(speaking_things));
social_nonsocial= normalize(squareform(social_nonsocial));
mentalization= normalize(squareform(mentalization));
amplitude= normalize(squareform(amplitude));
visual= normalize(squareform(visual));
face= normalize(squareform(face));
action = normalize(squareform(action));

%Subject list (Subject 5 removed)
subjects = {'s1','s2','s3','s4','s6','s7','s8',...
    's9','s10','s11','s12','s13','s14','s15','s16', 's17'};

%cropped and trimmed masks
 masks = {'bin_STG_without_TPJ.nii.gz', 'bin_finalMTG.nii.gz',...
    'bin_finalTPJ.nii.gz','TP.nii.gz',...
    'Precu.nii.gz','bin_aMPFCwithoutpMPFC.nii.gz','bin_finalpMPFC.nii.gz',...
    'IFG_oper.nii.gz','IPS.nii.gz','Auditory.nii',...
    'bin_VisualwithoutMTG.nii.gz',};

config=cosmo_config();
data_path=fullfile(config.tutorial_data_path);


n_subjects=numel(subjects);
n_masks=numel(masks);

targets= (1:988);
for m = 1:n_masks
    counter=0;
    msk = masks{m};
    %Print the current mask
    disp(['Mask: ' msk]);
    for s = 1:length(subjects)
        
        % Print the current subject
        sub = subjects{s};
        disp(['Subject: ' sub]);

        % Set path for subjects
        sub_path = '/Volumes/My Passport for Mac/SherlockRSA/sherlock_fmri/3TR';
        
        % load fmri data 
        ds_fn=fullfile(sub_path,[sub '_trimmed.nii']);
        
        % Set ROI path location
        % load masks
        mask_fn=fullfile(data_path, 'ROI_Masks','combinedROI',...
            'ROI_cropall', msk);
        ds_full = cosmo_fmri_dataset(ds_fn,...
                                    'mask',mask_fn,...
                                    'targets',targets);

        % compute average for each unique target
        ds=cosmo_fx(ds_full, @(x)mean(x,1), 'targets', 1);


        % remove constant features
        ds=cosmo_remove_useless_data(ds);


        % demean
        ds.samples = bsxfun(@minus, ds.samples, mean(ds.samples, 1));  
        
        %pidst and normalize the brain data
        brain=normalize(pdist(ds.samples,'correlation')); 
        
        %Run Regression Analysis
        % Using matlab function fitglm, we investigate the relative
        % contributions of the feature similarity matricies.
        
        %make a table consisting of all predictors and brain data 
        %dependent variable (here brain data) should be at the end.
        %check the dimensionlity of the variables if they match. 
        tbl = table(speaking_self', speaking_others', speaking_things', ...
            social_nonsocial', mentalization', amplitude', visual', face', action', brain',...
            'VariableNames',{'self', 'other', 'thing', 'social', 'TOM', ...
            'audio','DNN','face', 'action', 'ROI'});
        %fit glm with the provided table
        lm=fitglm(tbl);
        %lm contains lots of info but we want the beta values.
        beta=lm.Coefficients(2:height(lm.Coefficients),:)
        if counter==0
            % first dsm, allocate space
            betavalues=zeros(n_subjects,9);
        end
        
        % increase counter and store the beta values
        counter=counter+1;
        betavalues(counter,:)=beta.Estimate';
        
        %total variance (r2) explained by all predictors.
        rsquared_adj(counter,:)=lm.Rsquared.Adjusted;
    end
    
    %Output path
    save(['./SherlockBeta/Beta' msk '_cosmo.mat'], 'betavalues');
    save(['./SherlockBeta/rsquared_adj' msk '_cosmo.mat'], 'rsquared_adj');
end