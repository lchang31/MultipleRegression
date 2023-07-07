%%
% Performs a permutation test and grooups the beta values according to 
% feature and ROI
% 1/5/21, Lucy Chang
masks = {'STG','MTG',...
    'TPJ','TP',...
    'Precu','aMPFC','pMPFC',...
    'IFG','IPS','Auditory',...
    'Visual', ...
    };
features = {'self', 'other', 'thing', 'social_nonsocial', 'mentalization', ...
            'audio','visual','face', 'action'};

%% Group data by features
i = 1;
n_masks=numel(masks);
n_features = numel(features);
for k = 1:n_features
    counter = 0;
    feature = features{k};
    for m = 1:n_masks
        msk = masks{m};
        load(['Beta' msk '.mat'])
        if counter == 0
            beta = zeros(16, n_masks);
        end
        counter = counter + 1;
        beta(:, counter) = betavalues(:,k);
    end 
    save([feature 'beta.mat'], 'beta');
end 
%%
% Perform Permutation test
counter = 0;
for k = 1:n_features
    feature = features{k};
    disp(['Analyzing ', feature]);
    load([feature 'beta.mat'])
    [pvalue, obs_stat, rand_stat, pvalue_corr] = randomize_r(beta);
    if counter==0
        pvalue_o = zeros(n_features, n_masks);
        obs_tstat = zeros(n_features, n_masks);
        pvalue_n = zeros(n_features, n_masks);
    end
    counter = counter + 1;
    pvalue_o(counter,:)= pvalue;
    obs_tstat(counter,:) = obs_stat;
    pvalue_n(counter,:)= pvalue_corr;
end


