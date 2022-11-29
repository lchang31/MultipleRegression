%% Finds the set of beta coefficients for each social feature for each mask
% 1/21/20, Lucy Chang

% load all the features

predictor_path = './SummerPredictors';
load([predictor_path 'featureRDM.mat'])
%calculate the distance and normalize the values before doing the for-loop.
speaking_self= normalize(squareform(speaking_self));
speaking_others= normalize(squareform(speaking_others));
speaking_things=normalize(squareform(speaking_things));
social_nonsocial= normalize(squareform(social_nonsocial));
mentalization= normalize(squareform(mentalization));
social_touch= normalize(squareform(social_touch));
visual= normalize(squareform(Visual));
amplitude = normalize(squareform(amplitude));
face= normalize(squareform(face));
action = normalize(squareform(action));

subjects = {'s1','s2','s3','s4', 's5', 's6','s7','s8',...
    's9','s10','s11','s12','s13','s14','s15','s16', 's17', 's18'};

 masks = {'bin_STG_without_TPJ.nii', 'bin_finalMTG.nii',...
    'bin_finalTPJ.nii','TP.nii',...
    'Precu.nii','bin_aMPFCwithoutpMPFC.nii','bin_finalpMPFC.nii',...
    'IFG_oper.nii','IPS.nii','Auditory.nii',...
    'bin_VisualwithoutMTG.nii',};

config=cosmo_config();


n_subjects=numel(subjects);
n_masks=numel(masks);

targets= (1:1722);
for m = 1:n_masks
    counter=0;
    msk = masks{m};
    %Print the current mask
    disp(['Mask: ' msk]);
    for s = 1:length(subjects)
        
        % Print the current subject
        sub = subjects{s};
        disp(['Subject: ' sub]);
        %data_path = '/Volumes/My Passport for Mac/SummerRSA/';

        data_path = './';
        
        % load fmri data 
        ds_fn=fullfile(data_path, 'summer_fMRI', ...
            ['summer_movie_' sub '.mat']);
        load(ds_fn);
        data_clean.samples = data_clean.samples(5:5209,:);
        movingAverageA = movmean(data_clean.samples,3,1);
        Y = movingAverageA(2:3:end,:);
        data_clean.samples = (Y(14:end,:));
        
        % load masks
        mask_fn=fullfile(data_path, 'ROI', msk);
        ds_full = cosmo_fmri_dataset(data_clean,...
                                    'mask',mask_fn,...
                                    'targets',targets);
        
        % compute average for each unique target
        ds=cosmo_fx(ds_full, @(x)mean(x,1), 'targets', 1);

        % remove constant features
        ds=cosmo_remove_useless_data(ds);

        % demean
        % Comment this out to see the effects of demeaning vs. not
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
            social_nonsocial', mentalization', social_touch', visual', ...
             amplitude', face', action', brain',...
            'VariableNames',{'self', 'other', 'thing', 'social', 'TOM', ...
            'touch','DNN', 'audio','face', 'action', 'ROI'});
        %fit glm with the provided table
        lm=fitglm(tbl);
        %lm contains lots of info but we want the beta values.
        beta=lm.Coefficients(2:height(lm.Coefficients),:)
        if counter==0
            % first dsm, allocate space
            betavalues=zeros(n_subjects,10);
        end
        
        % increase counter and store the beta values
        counter=counter+1;
        betavalues(counter,:)=beta.Estimate';
        %also we get the total variance (r2) explained by all predictors.
        rsquared_adj(counter,:)=lm.Rsquared.Adjusted;
    end
    opath = './SummerBeta/BetaROI';
    save([opath 'Beta' msk '_cosmo.mat'], 'betavalues');
    save([opath 'rsquared_adj' msk '_cosmo.mat'], 'rsquared_adj');
end