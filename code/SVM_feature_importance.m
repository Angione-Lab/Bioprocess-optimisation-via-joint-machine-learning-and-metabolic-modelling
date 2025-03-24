
clear
load SVM_results_pFBA_end_feature_permutation.mat
SVM_results_permutation = SVM_results;
load SVM_results_pFBA_end.mat



%%
%%
% get feature importance
[num_samples, r] = size(SVM_results.mdl);
[~, k] = size(SVM_results_permutation.mdl);
predictors = SVM_results.predictors';
pathways = SVM_results.pathways';
beta_matrix = zeros(length(predictors), r);
permutation_beta_matrix = zeros(length(predictors), k);
% get LOOCV-averaged feature importance for each LOOCV repetition
for i = 1:r
    for j = 1:num_samples
        % sum normalised absolute feature scores
        beta_matrix(:, i) = beta_matrix(:, i) + abs(SVM_results.mdl{j, i}.Beta) ./ max(abs(SVM_results.mdl{j, i}.Beta));
    end
end
% average scores over CV iterations
beta_matrix = beta_matrix ./ num_samples;
% get LOOCV-averaged feature importance for each permuted LOOCV repetition
for i = 1:k
    for j = 1:num_samples
        % sum normalised absolute feature scores
        permutation_beta_matrix(:, i) = permutation_beta_matrix(:, i) + abs(SVM_results_permutation.mdl{j, i}.Beta) ./ max(abs(SVM_results_permutation.mdl{j, i}.Beta));
    end
end
% average scores over CV iterations
permutation_beta_matrix = permutation_beta_matrix ./ num_samples;

% get p-values
p_importance = zeros(length(predictors), 1);
for i = 1:length(predictors)
    p_importance(i) = ranksum(beta_matrix(i, :), permutation_beta_matrix(i, :), 'tail', 'right');
end
% adjust p-values
pFDR_importance = mafdr(p_importance, 'BHFDR', true);

figure
hold on
histogram(p_importance, 20)
histogram(pFDR_importance, 20)
hold off
legend('Unadjusted p-values', 'Adjusted p-values')

% pathway enrichment with hypergeometric test and FDR correction
enrichment_table = my_FEA(pFDR_importance <= 0.05, pathways);
for i = 2:size(enrichment_table, 1)
    for j = 1:2
        enrichment_table{i,j} = round(enrichment_table{i,j}, 3, 'significant');
    end
end

save('..\data\results\SVM_results_pFBA_end_weight_enrichment')




%%
%%
function resultCell = my_FEA(rxnSet, group)
% Significance analysis - Flux enrichment analysis using hypergeometric
% 1-sided test and FDR correction for multiple testing
%
% USAGE:
%
%    resultCell = FEA(model, rxnSet, 'subSystems')
%
% INPUTS:
%    model:           COBRA structure model
%    rxnSet:          reaction set to be enriched (vector of reaction indices e.g. 1:10)
%    group:           model.group structure e.g.
%                    'subSystems' : FEA looks for significantly enriched subsystems in rxnSet
%
% OUTPUT:
%    resultCell:    cell structure of enriched groups
%
% EXAMPLE:
%
%    load ecoli_core_model;
%    resultCell = FEA(modelEcore, 1:10, 'subSystems');
%

if nargin < 2
    error('The function FEA must be called with reaction set and group as arguments')
end
if ~isvector(rxnSet)
    error('Please provide the indices of the reactions e.g. 1:10')
end
if ~iscell(group)
    error('Please provide the group name as cell array of characters e.g. the subSystem field of any metabolic model ')
end

%Temporary Warning until FEA statistics are checked for multiple classes.
if iscell(group{1}) %Potentially multiple subSystems
    if any(cellfun(@numel, group) > 2)
        warning('Multiple subSystems detected for some reactions. FEA statistics might not be correct.\n Please consider using only one subSystem per reaction.')
    end
end

% compute frequency of enriched terms
%groups = eval(['model.' group]);
groups = group;
if iscell([groups{:}])
   [uniquehSubsystemsA] = unique([groups{:}]);
   presenceindicator = false(numel(uniquehSubsystemsA),numel(group));
   for i = 1:numel(groups)
       presenceindicator(:,i) = ismember(uniquehSubsystemsA,groups{i});
   end   
   [K,~] = find(presenceindicator);
else
    %This works only for fields which have a single entry.
    [uniquehSubsystemsA, ~, K] = unique(groups);
end
% fetch group
%enRxns = eval(['model.' group '(rxnSet)']);
enRxns = group(rxnSet);
m = length(uniquehSubsystemsA);
allSubsystems = zeros(1, m);

% look for unique occurences
if iscell([enRxns{:}])
   [uniquehSubsystems] = unique([enRxns{:}]);
   presenceindicator = false(numel(uniquehSubsystems),numel(group));
   for i = 1:numel(enRxns)       
        presenceindicator(:,i) = ismember(uniquehSubsystems,enRxns{i});
   end   
   [J,~] = find(presenceindicator);
else
    %This works only for fields which have a single entry.
    [uniquehSubsystems, ~, J] = unique(enRxns);
end

occ = histc(J, 1:numel(uniquehSubsystems));
[l, p] = intersect(uniquehSubsystemsA, uniquehSubsystems);
allSubsystems(p) = occ;

% compute total number of reactions per group
nRxns = histc(K, 1:numel(uniquehSubsystemsA));  % the number of reactions per susbsystem

% Compute p-values
% gopvalues = hygepdf(allSubsystems', max(nRxns), max(allSubsystems), nRxns);
gopvalues = hygecdf(allSubsystems'-1, repmat(sum(nRxns), length(nRxns), 1), nRxns, repmat(sum(allSubsystems), length(nRxns), 1), 'upper');

% take out the zeros for one-sided test
nonZerInd = find(allSubsystems);

% sort p-values
[m, rxnInd] = sort(gopvalues);

% intersect non zero sets with ordered pvalues
[~, nonZeroInd] = intersect(rxnInd, nonZerInd);
orderedPval = rxnInd(sort(nonZeroInd));

% Build result cell
% initilize variable
resultCell = cell(length(orderedPval) + 1, 5);
resultCell(1, :) = {'p-value', 'Adjusted p-value', 'Pathway', 'Enriched set size', 'Total set size'};

% P values
resultCell(2:end, 1) = num2cell(gopvalues(orderedPval));

% correct for multiple testing with FDR
resultCell(2:end, 2) = num2cell(mafdr(cell2mat(resultCell(2:end, 1)), 'BHFDR', true));

% Group name
resultCell(2:end, 3) = uniquehSubsystemsA(orderedPval);

% Test size
resultCell(2:end, 4) = num2cell(allSubsystems(orderedPval))';

% Total group size
resultCell(2:end, 5) = num2cell(nRxns(orderedPval));

end
