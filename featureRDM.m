%% Create feature RDM
% generate feature RDM and visualize
% 10/31/20, Lucy Chang

%% Predictors for Summer
predictors = {'speaking_self', 'speaking_others', 'speaking_things',...
     'social_nonsocial', 'mentalization', 'amplitude', 'visual', 'social_touch', 'face',...
     'action'};
%% Predictors for Sherlock 
predictors = {'speaking_self', 'speaking_others', 'speaking_things',...
     'social_nonsocial', 'mentalization', 'amplitude', 'layer2', 'face',...
     'action'};
n_predictors=numel(predictors);

% indicate location of RDM files 
cd('/Volumes/My Passport for Mac/SummerRSA/SummerPredictors')

% Loop through each predictor
for i = 1:n_predictors
    x = predictors{i};
    array = importdata([x '.mat']);
    subplot(3,3,i); % create a subplot
    
    % if the feature is visual
    if i == 7
        dsm = squareform(pdist(array, 'correlation'));
        imagesc(dsm);
        title('Visual');
        colorbar();
        save('Visual_layer2_3TR.mat', 'dsm');
    % if input is action
    elseif i == 9
        imagesc(array);
        title('Action');
        colorbar();
    else
        dsm = squareform(pdist(array));
        save([predictors{i} '_3TR.mat'], 'dsm');
        imagesc(dsm);
        title(predictors{i});
        colorbar();
    end
end 

%% 2nd order RDM
% predictors = {'speaking self', 'speaking others', 'speaking things',...
%     'social nonsocial', 'mentalization', 'social_touch', 'amplitude', 'visual', 'face', 'action'};

speaking_self = squareform(speaking_self);
speaking_others= squareform(speaking_others);
speaking_things= squareform(speaking_things);
social_nonsocial= squareform(social_nonsocial);
mentalization= squareform(mentalization);
% social_touch = squareform(social_touch);
amplitude= squareform(amplitude);
visual= squareform(visual);
face= squareform(face);
action = squareform(action);

%Add predictors all in one matrix
predictor = [speaking_self; speaking_others; speaking_things; social_nonsocial;...
    mentalization; amplitude; visual; face; action];

cc = cosmo_corr(predictor');
figure();
imagesc(cc);

D = pdist(cc);
% MDS Visualization
[x,stress_2D] = mdscale(D,2,'Criterion','stress','Start','random','Replicates',100);
[y,stress_3D] = mdscale(D,3,'Criterion','stress','Start','random','Replicates',100);

figure();
% General 2D plot (including auditory and Visual)
subplot(2,2,1)
plot(x(:,1), x(:,2), 'r.', 'MarkerSize', 20);
text(x(:,1) + 0.025, x(:,2), predictors);

% General 3D plot (including auditory and Visual)
subplot(2,2,2)
plot3(y(:,1), y(:,2), y(:,3), 'r.', 'MarkerSize', 20);
text(y(:,1)+ 0.025, y(:,2), y(:,3), predictors);

%% create a function that takes the average of RDM and normalize them
function new_array = crop_TR(array, N)
    movingAverageA = movmean(array,[(N-1) 0]);
    new_array = movingAverageA(N:N:end,:);
end
    
