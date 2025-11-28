%% Start
%%
clear
rng default
addpath('..\..\packages\PCHA');

task = 'archetype_enrichment'; % archetypes, varexpl_curve, archetype_error, archetype_enrichment
save_results = true;



%% Task 1: explained variation curve
%%
load end_pFBA_fluxes.mat
fluxes = fluxes(1:18, :);
num_samples = size(fluxes, 1);
U = 1:size(fluxes, 1); % Entries in X used that is modelled by the AA model
I = 1:size(fluxes, 1); % Entries in X used to define archetypes
delta = 0;
opts.maxiter = 1000;
opts.conv_crit = 1e-6;

if strcmp(task, 'archetypes')
    noc = 3; % Number of archetypes
    [XC, S, C, SSE, varexpl] = PCHA(zscore(fluxes)', noc, I, U, delta, opts);

elseif strcmp(task, 'varexpl_curve')
    SSE = NaN(num_samples-1, 1);
    varexpl = NaN(num_samples-1, 1);
    for noc = 2:num_samples-1
        [XC, S, C, SSE(noc), varexpl(noc)] = PCHA(zscore(fluxes)', noc, I, U, delta, opts);
    end

    figure
    plot(1:num_samples-1, varexpl, '-o')
    xlabel('Number of archetypes')
    ylabel('Fraction of explained variance')

elseif strcmp(task, 'archetype_error')
    noc = 3; % Number of archetypes
    num_bootstraps = 100;
    xc = NaN(size(fluxes, 2), noc, num_bootstraps);
    for i = 1:num_bootstraps
        bootstrapped_fluxes = fluxes(randi(num_samples, num_samples, 1), :);
        [xc(:, :, i), S, C, SSE, varexpl] = PCHA(zscore(bootstrapped_fluxes)', noc, I, U, delta, opts);
    end
    XC = mean(xc, 3);
    XC_std = std(xc, 0, 3);
    
elseif strcmp(task, 'archetype_enrichment')
    load('AA_results_exploration_archetypes.mat', 'XC', 'noc');
    % remove constant rows
    idx = sum(diff(XC, 1, 2) == 0, 2) < 2;
    XC = XC(idx, :);
    subSystems = subSystems(idx);
    % calculate archetype centroid
    centroid = mean(XC, 2);
    enrichment_table = {};
    for i = 1:noc
        idx = ones(noc, 1);
        idx(i) = 0;
        idx = find(idx);
        idx = (XC(:, i) - XC(:, idx(1)) > 0) & (XC(:, i) - XC(:, idx(2)) > 0);
        t = my_FEA(idx, subSystems);
        enrichment_table{i} = t;
    end
end

if save_results
    save(['..\data\results\AA_results_exploration_' task]);
end
