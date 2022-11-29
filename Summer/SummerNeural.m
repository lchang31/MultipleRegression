%% 
% Generates the individual Neural RDMs for each ROI for the Summer fMRI data
% Generates the second order RDM and provides visualization of RDMs
% 4/13/21 Lucy Chang
%%
subjects = {'s1','s2','s3','s4', 's5', 's6','s7','s8',...
    's9','s10','s11','s12','s13','s14','s15','s16', 's17', 's18'};

masks = {'bin_STG_without_TPJ.nii', 'bin_finalMTG.nii',...
    'bin_finalTPJ.nii','TP.nii',...
    'Precu.nii','bin_aMPFCwithoutpMPFC.nii','bin_finalpMPFC.nii',...
    'IFG_oper.nii','IPS.nii','Auditory.nii',...
    'bin_VisualwithoutMTG.nii'};

%Set data path
config=cosmo_config();
data_path=fullfile(config.tutorial_data_path);

n_subjects=numel(subjects);
n_masks=numel(masks);

%Set target--> number of TR
targets= (1:1722);
for m = 1:n_masks
    counter=0;
    msk = masks{m};
    disp(['Mask: ' msk]);
    for s = 1:length(subjects)
        sub = subjects{s};
        disp(['Subject: ' sub]);
        % Set subject path
        data_path = './';
        ds_fn=fullfile(data_path, 'summer_fMRI', ...
            ['summer_movie_' sub '.mat']);
        
        load(ds_fn);
        %crop the beginning and average
        data_clean.samples = data_clean.samples(5:5209,:);
        movingAverageA = movmean(data_clean.samples,3,1);
        Y = movingAverageA(2:3:end,:);
        data_clean.samples = (Y(14:end,:));
        
        % set data path
        mask_fn=fullfile(data_path, 'ROI', msk);
        ds_full = cosmo_fmri_dataset(data_clean,...
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
            neural_dsms=zeros(n_subjects,n_pairs);
        end

        % increase counter and store the dsm as the counter-th row in
        % 'neural_dsms'
        % >@@>
        counter=counter+1;
        neural_dsms(counter,:)=dsm;
        % <@@<
    end
    %cd('/Volumes/My Passport for Mac/SummerRSA/SummerPredictors')
    save(['./SummerNeural/' msk 'neural.mat', 'neural_dsms'], 'neural_dsms', '-v7.3');
end
% save('Neural_dsm.mat', 'neural_dsms'); % save the neural dsms

%% Combine Neural RDM
counter = 1;
z = 0;
n_masks=numel(masks);
for m = 1:n_masks
    msk = masks{m};
    if counter==1
        neural_dsms=zeros(18*11,1481781);
    end
    dsm = importdata([msk 'neural.mat']);
    z = z + 18;
    neural_dsms(counter:z,:)=dsm;
    counter = counter + 18;
end

%% Visualize Neural RDM
N = 18;
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
% General 2D plot (including auditory and Visual)
subplot(2,2,1)
plot(x(:,1), x(:,2), 'r.', 'MarkerSize', 20);
text(x(:,1) + 0.025, x(:,2), masks);

% General 3D plot (including auditory and Visual)
subplot(2,2,2)
plot3(y(:,1), y(:,2), y(:,3), 'r.', 'MarkerSize', 20);
text(y(:,1)+ 0.025, y(:,2), y(:,3), masks);

% 2D visualization plot 
subplot(2,2,3)
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

% % 3D visualization
subplot(2,2,4)
plot3(y(1:2,1), y(1:2,2), y(1:2,3), 'r.', y(3:7,1), y(3:7,2), y(3:7,3), 'b.',...
y(8:9,1), y(8:9,2), y(8:9,3), 'g.', y(10:11,1), y(10:11,2),...
 y(10:11,3), 'y.', 'MarkerSize',20);
text(y(1:2,1) + 0.025, y(1:2,2), y(1:2,3), social);
text(y(3:7,1) + 0.025, y(3:7,2),y(3:7,3), mentalization);
text(y(8:9,1) + 0.025, y(8:9,2), y(8:9,3), action_observation);
text(y(10:11,1) + 0.025, y(10:11,2),y(10:11,3), sensory);
legend('social', 'mentalization', 'action observation', 'sensory')

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
    m = k + 17;
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




