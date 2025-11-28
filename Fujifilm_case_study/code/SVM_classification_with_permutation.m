%% Start
%%
clear
rng default

feature_cases = [repmat({'MEP'}, 1, 1)];
phase_cases = repmat({'post'}, 1, 1);
num_samples = 22;
num_cv_repetitions = 100;
permutation = 'label'; % 'label', 'feature'
save_results = true;



for f_case = 1:length(feature_cases)
    features = feature_cases{f_case};
    phase = phase_cases{f_case};
    
    %% Load and pre-process data
    %%
    if strcmp(features, 'DoE')
        load DoE_and_offline_data.mat
        exp_data = dummyvar({DoE_and_offline_data.CLD, DoE_and_offline_data.Gene, DoE_and_offline_data.Host, ...
            DoE_and_offline_data.Plasmid, DoE_and_offline_data.T}); % plasmid = IPTG
        exp_variables = [repmat({'CLD'}, 1, 11) repmat({'Gene'}, 1, 3) repmat({'Host'}, 1, 2) ...
                repmat({'Plasmid/IPTG'}, 1, 2) repmat({'T'}, 1, 2)];
        exp_data = exp_data(:, sum(isnan(exp_data)) == 0); % remove features with NaNs
        exp_variables = exp_variables(sum(isnan(exp_data)) == 0);
        exp_data = exp_data(:, not(sum(diff(exp_data) == 0) == num_samples-1)); % remove constant features
        exp_variables = exp_variables(not(sum(diff(exp_data) == 0) == num_samples-1));
    
    elseif strcmp(features, 'concentrations')
        load DoE_and_offline_data.mat
        if strcmp(phase, 'pre')
            exp_data = [DoE_and_offline_data.AC2B_t1 DoE_and_offline_data.AC2B_t2 ...
                DoE_and_offline_data.GLYB_t1 DoE_and_offline_data.GLYB_t2 ...
                DoE_and_offline_data.NH3B_t1 DoE_and_offline_data.NH3B_t2 ...
                DoE_and_offline_data.DCW_t1 DoE_and_offline_data.DCW_t2]; % take experimental data only for the first two time points
            exp_variables = ['AC2B_t1' 'AC2B_t2' 'GLYB_t1' 'GLYB_t2' 'NH3B_t1' 'NH3B_t2' 'DCW_t1' 'DCW_t2'];
        elseif strcmp(phase, 'post')
            exp_data = [DoE_and_offline_data.AC2B_t2 DoE_and_offline_data.AC2B_t3 ...
                DoE_and_offline_data.GLYB_t2 DoE_and_offline_data.GLYB_t3 ...
                DoE_and_offline_data.NH3B_t2 DoE_and_offline_data.NH3B_t3 ...
                DoE_and_offline_data.DCW_t2 DoE_and_offline_data.DCW_t3]; % take experimental data only for the last two time points
            exp_variables = ['AC2B_t2' 'AC2B_t3' 'GLYB_t2' 'GLYB_t3' 'NH3B_t2' 'NH3B_t3' 'DCW_t2' 'DCW_t3'];
        elseif strcmp(phase, 'prepost')
            exp_data = [DoE_and_offline_data{:, 7:end}]; % take experimental data for all time points
            exp_variables = [DoE_and_offline_data.Properties.VariableNames{:, 7:end}];
        end
        exp_data = exp_data(:, sum(isnan(exp_data)) == 0); % remove features with NaNs
        exp_variables = exp_variables(sum(isnan(exp_data)) == 0);
        exp_data = exp_data(:, not(sum(diff(exp_data) == 0) == num_samples-1)); % remove constant features
        exp_variables = exp_variables(not(sum(diff(exp_data) == 0) == num_samples-1));
    
    else
        if strcmp(features, 'pFBA')
            if strcmp(phase, 'pre')
                load preinduction_pFBA_fluxes.mat
            elseif strcmp(phase, 'post')
                load postinduction_pFBA_fluxes.mat
            elseif strcmp(phase, 'end')
                load end_pFBA_fluxes.mat
            end
        elseif contains(features, 'MEP')
            if strcmp(phase, 'pre')
                load preinduction_MEP_fluxes.mat
                fluxes = [av_fluxes va_flux_stds];
                rxns = [rxns; strcat(rxns, repmat('_var', length(rxns), 1))];
            elseif strcmp(phase, 'post')
                load postinduction_MEP_fluxes.mat
                fluxes = mu_fluxes;
            end
        end
        
        % remove fluxes that are constant in all conditions
        constant_fluxes_idx = sum(diff(fluxes) == 0) == num_samples-1;
        fluxes = fluxes(:, not(constant_fluxes_idx));
        rxns = rxns(not(constant_fluxes_idx))';
        subSystems = subSystems(not(constant_fluxes_idx))';
    end
    
    % titres
    load titres.mat
    Y = titres.labels([2:5, 7:24]); % targets
    
    
    
    
    
    %% Cross-validation
    %%
    lambda_grid = logspace(-3,3,10)/21;
    num_folds_ext = 22; % for model evaluation
    num_folds_int = 7; % for hyper-parameter selection
    SVM_results.predicted = nan(num_folds_ext, fix(num_samples/num_folds_ext), num_cv_repetitions);
    SVM_results.true = nan(num_folds_ext, fix(num_samples/num_folds_ext), num_cv_repetitions);
    SVM_results.mdl = cell(num_folds_ext, num_cv_repetitions);
    SVM_results.FitInfo = cell(num_folds_ext, num_cv_repetitions);
    SVM_results.internal_accuracy = nan(num_folds_ext, num_cv_repetitions);
    if strcmp(features, 'DoE') || strcmp(features, 'concentrations')
        data = exp_data;
        SVM_results.predictors = exp_variables;
    elseif strcmp(features, 'pFBA') || strcmp(features, 'MEP')
        data = fluxes;
        SVM_results.predictors = rxns;
        SVM_results.pathways = subSystems;
    end

    disp(['features = ', features])
    disp(['phase = ', phase])
    disp(['permutation = ', permutation])
    
    for k = 1:num_cv_repetitions
        disp(num2str(k))
        if strcmp(permutation, 'label')
            permutation_idx = randperm(22);
            permuted_data = data(permutation_idx, :);
        elseif strcmp(permutation, 'feature')
            % shuffle the feature data
            permuted_data = nan(size(data));
            for i = 1:size(data,2) % for each dimension
                permutation_idx = randperm(size(data, 1)); % shuffle the data values of each axis
                permuted_data(:, i) = data(permutation_idx, i);
            end
        end
        
        for i = 1:num_folds_ext % external cv
            test_idx_ext = logical(1:22 == i)';
            Xtrain = permuted_data(~test_idx_ext, :);
            Ytrain = Y(~test_idx_ext);

            cv_idx_int = crossvalind('Kfold', length(Ytrain), num_folds_int);
            internal_accuracy = zeros(length(lambda_grid), 1);
            for l = 1:length(lambda_grid)
                for m = 1:num_folds_int
                    test_idx_int = (cv_idx_int == m);
                    xtrain = Xtrain(~test_idx_int, :);
                    ytrain = Ytrain(~test_idx_int);
                    [x, ps] = mapstd(xtrain');
                    xtrain = x';
                    mdl = fitclinear(xtrain', ytrain', 'Learner', 'svm', 'Regularization', 'ridge', ...
                        'Lambda', lambda_grid(l), 'ObservationsIn', 'columns'); % train model

                    xtest = Xtrain(test_idx_int, :);
                    ytest = Ytrain(test_idx_int);
                    xtest = mapstd('apply', xtest', ps)';
                    pred_Y = predict(mdl, xtest); % predict
                    internal_accuracy(l) = internal_accuracy(l) + acc(pred_Y, ytest);
                end
            end
            internal_accuracy = internal_accuracy / num_folds_int;

            [max_1] = ind2sub(length(lambda_grid), ...
                find(internal_accuracy == max(internal_accuracy(:))));
            if ismatrix(max_1)
                max_1 = max_1(1);
            end
            best_lambda = lambda_grid(max_1);

            [x, ps] = mapstd(Xtrain');
            Xtrain = x';
            [mdl,FitInfo] = fitclinear(Xtrain', Ytrain', 'Learner', 'svm', 'Regularization', 'ridge', ...
                        'Lambda', best_lambda, 'ObservationsIn', 'columns'); % train model

            Xtest = permuted_data(test_idx_ext, :);
            Ytest = Y(test_idx_ext);
            Xtest = mapstd('apply', Xtest', ps)';
            pred_Y = predict(mdl, Xtest); % predict

            % evaluate
            SVM_results.predicted(i, 1:length(pred_Y), k) = pred_Y;
            SVM_results.true(i, 1:length(pred_Y), k) = Ytest;
            SVM_results.mdl{i, k} = mdl;
            SVM_results.FitInfo{i, k} = FitInfo;
            SVM_results.internal_accuracy(i, k) = max(internal_accuracy(:));
        end
    end
    
    SVM_results = evaluate_classification(SVM_results);
    disp(['total LOOCV accuracy = ', num2str(SVM_results.total_accuracy)])
    disp(['total LOOCV MCC = ', num2str(SVM_results.total_MCC)])
    




    %% Save
    %%
    if save_results
        save(['..\data\results\SVM_results_' features '_' phase '_' permutation '_permutation'], 'SVM_results');
    end
end





%% Functions
%%
function s = evaluate_classification(s)
    
    num_folds = size(s.predicted, 1);
    s.bioreactor_accuracy = nan(num_folds, 1);
    s.total_accuracy = acc(s.predicted(:), s.true(:));
    s.total_MCC = mcc(s.predicted(:), s.true(:));
    
    for i = 1:num_folds
        s.bioreactor_accuracy(i) = acc(s.predicted(i, :, :), s.true(i, :, :));
    end
end

function res = acc(pred_Y, Ytest)
    res = sum(pred_Y == Ytest) / length(Ytest);
end

function res = mcc(pred_Y, Ytest)
    if length(unique([pred_Y; Ytest])) > 1
        C = confusionmat(Ytest, pred_Y);
        TN = C(1, 1);
        TP = C(2, 2);
        FN = C(2, 1);
        FP = C(1, 2);
        res = (TP*TN-FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    else
        res = 0;
    end
end
