%% Generate Neural RDMs for each subject
% 1/21/20

% Subject list (subject 5 removed)
subjects = {'s1','s2','s3','s4','s6','s7','s8',...
    's9','s10','s11','s12','s13','s14','s15','s16', 's17'};
masks = {'bin_STG_without_TPJ.nii.gz', 'bin_finalMTG.nii.gz',...
    'bin_finalTPJ.nii.gz','TP.nii.gz',...
    'Precu.nii.gz','bin_aMPFCwithoutpMPFC.nii.gz','bin_finalpMPFC.nii.gz',...
    'IFG_oper.nii.gz','IPS.nii.gz','Auditory.nii',...
    'bin_VisualwithoutMTG.nii.gz',};
%Set data path
config=cosmo_config();
data_path=fullfile(config.tutorial_data_path);

n_subjects=numel(subjects);
n_masks=numel(masks);
counter=0;

%Set target--> number of TR
targets= (1:988);
for m = 1:n_masks
    msk = masks{m};
    for s = 1:length(subjects)
        sub = subjects{s};
        % Set subject path
        sub_path = fullfile(data_path, 'sherlock_fmri', '3TR');
        ds_fn=fullfile(sub_path,['/' sub '_trimmed.nii']);
        
        % set data path    
        mask_fn=fullfile(data_path, 'ROI_Masks', 'combinedROI','ROI_cropall',...
            msk);
        ds_full = cosmo_fmri_dataset(ds_fn,...
                                    'mask',mask_fn,...
                                    'targets',targets);

        % compute average for each unique target
        ds=cosmo_fx(ds_full, @(x)mean(x,1), 'targets', 1);


        % remove constant features
        ds=cosmo_remove_useless_data(ds);


        % demean
        ds.samples = bsxfun(@minus, ds.samples, mean(ds.samples, 1));

        % compute the one-minus-correlation value for each pair of
        % targets
        dsm=cosmo_pdist(ds.samples, 'correlation');
   
        if counter==0
            % first dsm, allocate space
            n_pairs=numel(dsm);
            neural_dsms=zeros(n_subjects*n_masks,n_pairs);
        end

        % increase counter and store the dsm as the counter-th row in
        % 'neural_dsms'
        % >@@>
        counter=counter+1;
        neural_dsms(counter,:)=dsm;
        % <@@<
    end
end
save('Neural_dsm.mat', 'neural_dsms'); % save the neural dsms

%% Visualize Neural RDM
N = 16;
movingAverageA = movmean(neural_dsms,[(N-1) 0]);
average = movingAverageA(N:N:end,:);
masks = {'STG','MTG',...
    'TPJ','TP',...
    'Precu','aMPFC','pMPFC',...
    'IFG oper','IPS','Auditory',...
    'Visual', ...
    };
figure();
for m = 1:numel(masks)
    msk = masks{m};
    subplot(4,3,m);
    imagesc(squareform(average(m,:)));
    colorbar();
    colormap(jet);
    title(msk);
end

%% 2nd order RDM
cc = cosmo_corr(average');
figure();
imagesc(cc);

% MDS Visualization
[x,stress_2D] = mdscale(cc,2,'Criterion','stress','Start','random','Replicates',100);
[y,stress_3D] = mdscale(cc,3,'Criterion','stress','Start','random','Replicates',100);

% ROI categories
 social = {'STG', 'MTG'};
 mentalization = {'TPJ','TP', 'Precu','aMPFC','pMPFC'};
 action_observation = {'IFG oper','IPS'};
 sensory = {'Auditory','Visual'};

figure();

% 2D visualization plot 
subplot(1,2,1)
plot(x(1:2,1), x(1:2,2), '.',  'Color',  [1, 0, 0], 'MarkerSize', 20);
hold on
plot(x(3:7,1), x(3:7,2), '.', 'Color',  [47/256, 151/256, 102/256],  'MarkerSize', 20);
plot(x(8:9,1), x(8:9,2), '.','Color',  [52/256, 106/256, 233/256], 'MarkerSize', 20);
plot(x(10:11,1), x(10:11,2),  '.','Color',  [255/256, 188/256, 0], 'MarkerSize',20 );
hold off
text(x(1:2,1) + 0.025, x(1:2,2), social);
text(x(3:7,1) + 0.025, x(3:7,2), mentalization);
text(x(8:9,1) + 0.025, x(8:9,2), action_observation);
text(x(10:11,1) + 0.025, x(10:11,2), sensory);
legend('social', 'mentalization', 'action observation', 'sensory')

% 3D visualization
subplot(1,2,2)
plot3(y(1:2,1), y(1:2,2), y(1:2,3), 'r.', y(3:7,1), y(3:7,2), y(3:7,3), 'b.',...
y(8:9,1), y(8:9,2), y(8:9,3), 'g.', 'MarkerSize',20);
text(y(1:2,1) + 0.025, y(1:2,2), y(1:2,3), social);
text(y(3:7,1) + 0.025, y(3:7,2),y(3:7,3), mentalization);
text(y(8:9,1) + 0.025, y(8:9,2), y(8:9,3), action_observation);
legend('social', 'mentalization', 'action observation')

%% Intersubject reliability test
masks = {'STG','MTG',...
    'TPJ','TP',...
    'Precu','aMPFC','pMPFC',...
    'IFG oper','IPS','Auditory',...
    'Visual', ...
    };
n_masks = numel(masks);
k = 1;
counter = 0;
for i= 1:n_masks
    m = k + 15;
    cc = cosmo_corr((neural_dsms((k:m), :))');
    z = triu(cc, 1);  %take the upper triangular matrix
    z(z == 0) = NaN;  
    r = z(~isnan(z));   %remove all the NaN
    avgr = mean(r);
    [pvalue_corr] = randomize_r(r);
    if counter==0
        p = zeros(n_masks, 1);
        avg_r = zeros(n_masks, 1);
    end
    counter = counter + 1;
    p(counter,:)= pvalue_corr;
    avg_r(counter, :) = avgr;
    k = m + 1;
end


