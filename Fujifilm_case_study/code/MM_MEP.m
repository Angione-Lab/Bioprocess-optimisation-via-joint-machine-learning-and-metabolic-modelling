%% Preliminary setup
%%
clear
addpath(genpath('..\..\packages\Metabolic-EP'));
changeCobraSolver('ibm_cplex','LP');
rng default
if not(isfolder('..\data\results'))
    mkdir('..\data\results')
end
addpath('..\data\results')

% Specify parameters:
Beta = 1e6;
damping = 0.9;
precision = 1e-6;
maxit = 1000;
minvar = 1e-50;
maxvar = 1e50;
process_phase = 'post'; % 'pre', 'post' or 'end'
save_results = true;





%% Initialize data
%%
load constraint_data.mat
load DoE_and_offline_data.mat
load critical_temperatures.mat
load iEC1356_Bl21DE3.mat
load iEC1372_W3110.mat
load BL21_gene_ID_map.mat
load W3110_gene_ID_map.mat
load titres.mat

% get RNA base frequencies for the product protein
[A_freq, U_freq, C_freq, G_freq] = get_RNA_base_frequences();
A_freq = num2str(A_freq / 30);
U_freq = num2str(U_freq / 30);
C_freq = num2str(C_freq / 30);
G_freq = num2str(G_freq / 30);

% add heterologous production to the reactions
iEC1356_Bl21DE3 = addMetabolite(iEC1356_Bl21DE3, 'product[c]', 'Protein product in cytoplasm');
iEC1372_W3110 = addMetabolite(iEC1372_W3110, 'product[c]', 'Protein product in cytoplasm');
iEC1356_Bl21DE3 = addReaction(iEC1356_Bl21DE3, 'PRODUCTtranslation1', ...
    'reactionFormula', ['10 ala__L[c] + 10 arg__L[c] + 10 asn__L[c] + 10 asp__L[c] + 10 cys__L[c] + 10 glu__L[c] + 10 gln__L[c] + 10 gly[c] + 10 his__L[c] + ' ...
    '10 ile__L[c] + 10 leu__L[c] + 10 lys__L[c] + 10 met__L[c] + 10 phe__L[c] + 10 pro__L[c] + 10 ser__L[c] + 10 thr__L[c] + 10 trp__L[c] + ' ...
    '10 tyr__L[c] + 10 val__L[c] + 10 alatrna[c] + 10 argtrna[c] + 10 asntrna[c] + 10 asptrna[c] + 10 cystrna[c] + 10 glutrna[c] + 10 glntrna[c] + 10 glytrna[c] + 10 histrna[c] + ' ...
    '10 iletrna[c] + 10 leutrna[c] + 10 lystrna[c] + 10 mettrna[c] + 10 phetrna[c] + 10 protrna[c] + 10 sertrna[c] + 10 thrtrna[c] + 10 trptrna[c] + ' ...
    '10 tyrtrna[c] + 10 valtrna[c] + ' C_freq ' ctp[c] + ' G_freq ' gtp[c] + ' A_freq ' atp[c] + ' U_freq ' utp[c] + 840 atp[c] -> ' ...
    'translation_complex[c]'], ...
    'printLevel', 0);
iEC1356_Bl21DE3 = addReaction(iEC1356_Bl21DE3, 'PRODUCTtranslation2', ...
    'reactionFormula', ['translation_complex[c] -> ' ...
    '10 trnaala[c] + 10 trnaarg[c] + 10 trnaasn[c] + 10 trnaasp[c] + 10 trnacys[c] + 10 trnaglu[c] + 10 trnagln[c] + 10 trnagly[c] + 10 trnahis[c] + ' ...
    '10 trnaile[c] + 10 trnaleu[c] + 10 trnalys[c] + 10 trnamet[c] + 10 trnaphe[c] + 10 trnapro[c] + 10 trnaser[c] + 10 trnathr[c] + 10 trnatrp[c] + ' ...
    '10 trnatyr[c] + 10 trnaval[c] + ' C_freq ' ctp[c] + ' G_freq ' gtp[c] + ' A_freq ' atp[c] + ' U_freq ' utp[c] + product[c] + 840 adp[c] + 840 pi[c]'], ...
    'printLevel', 0);
iEC1372_W3110 = addReaction(iEC1372_W3110, 'PRODUCTtranslation1', ...
    'reactionFormula', ['10 ala__L[c] + 10 arg__L[c] + 10 asn__L[c] + 10 asp__L[c] + 10 cys__L[c] + 10 glu__L[c] + 10 gln__L[c] + 10 gly[c] + 10 his__L[c] + ' ...
    '10 ile__L[c] + 10 leu__L[c] + 10 lys__L[c] + 10 met__L[c] + 10 phe__L[c] + 10 pro__L[c] + 10 ser__L[c] + 10 thr__L[c] + 10 trp__L[c] + ' ...
    '10 tyr__L[c] + 10 val__L[c] + 10 alatrna[c] + 10 argtrna[c] + 10 asntrna[c] + 10 asptrna[c] + 10 cystrna[c] + 10 glutrna[c] + 10 glntrna[c] + 10 glytrna[c] + 10 histrna[c] + ' ...
    '10 iletrna[c] + 10 leutrna[c] + 10 lystrna[c] + 10 mettrna[c] + 10 phetrna[c] + 10 protrna[c] + 10 sertrna[c] + 10 thrtrna[c] + 10 trptrna[c] + ' ...
    '10 tyrtrna[c] + 10 valtrna[c] + ' C_freq ' ctp[c] + ' G_freq ' gtp[c] + ' A_freq ' atp[c] + ' U_freq ' utp[c] + 840 atp[c] -> ' ...
    'translation_complex[c]'], ...
    'printLevel', 0);
iEC1372_W3110 = addReaction(iEC1372_W3110, 'PRODUCTtranslation2', ...
    'reactionFormula', ['translation_complex[c] -> ' ...
    '10 trnaala[c] + 10 trnaarg[c] + 10 trnaasn[c] + 10 trnaasp[c] + 10 trnacys[c] + 10 trnaglu[c] + 10 trnagln[c] + 10 trnagly[c] + 10 trnahis[c] + ' ...
    '10 trnaile[c] + 10 trnaleu[c] + 10 trnalys[c] + 10 trnamet[c] + 10 trnaphe[c] + 10 trnapro[c] + 10 trnaser[c] + 10 trnathr[c] + 10 trnatrp[c] + ' ...
    '10 trnatyr[c] + 10 trnaval[c] + ' C_freq ' ctp[c] + ' G_freq ' gtp[c] + ' A_freq ' atp[c] + ' U_freq ' utp[c] + product[c] + 840 adp[c] + 840 pi[c]'], ...
    'printLevel', 0);

% rescale "infinity" bounds
iEC1356_Bl21DE3.lb(iEC1356_Bl21DE3.lb == -1000) = -100;
iEC1356_Bl21DE3.ub(iEC1356_Bl21DE3.ub == 1000) = 100;
iEC1372_W3110.lb(iEC1372_W3110.lb == -1000) = -100;
iEC1372_W3110.ub(iEC1372_W3110.ub == 1000) = 100;

% homogenise IDs of reactions shared by iEC1356_Bl21DE3 and iEC1372_W3110
iEC1372_W3110.rxns = strrep(iEC1372_W3110.rxns, 'ADNt2pp_copy1', 'ADNt2pp');
iEC1372_W3110.rxns = strrep(iEC1372_W3110.rxns, 'CYTDt2pp_copy1', 'CYTDt2pp');
iEC1372_W3110.rxns = strrep(iEC1372_W3110.rxns, 'INSt2pp_copy1', 'INSt2pp');
iEC1372_W3110.rxns = strrep(iEC1372_W3110.rxns, 'THMDt2pp_copy1', 'THMDt2pp');
iEC1372_W3110.rxns = strrep(iEC1372_W3110.rxns, 'URIt2pp_copy1', 'URIt2pp');

% convert gene names in the "genes" field into "b****" identifiers
for i = 1:length(iEC1356_Bl21DE3.genes)
    if isKey(BL21_gene_ID_map, iEC1356_Bl21DE3.genes(i))
        iEC1356_Bl21DE3.genes{i} = BL21_gene_ID_map(iEC1356_Bl21DE3.genes{i});
    end
end
for i = 1:length(iEC1372_W3110.genes)
    if isKey(W3110_gene_ID_map, iEC1372_W3110.genes(i))
        iEC1372_W3110.genes{i} = W3110_gene_ID_map(iEC1372_W3110.genes{i});
    end
end

% merge duplicate genes (1 per model)
[~, ia] = unique(iEC1356_Bl21DE3.genes);
duplicate_idxs = setdiff(1:length(iEC1356_Bl21DE3.genes), ia);
[iEC1356_Bl21DE3] = mergeModelFieldPositions(iEC1356_Bl21DE3, 'genes', strcmp(iEC1356_Bl21DE3.genes, iEC1356_Bl21DE3.genes{duplicate_idxs}));
[~, ia] = unique(iEC1372_W3110.genes);
duplicate_idxs = setdiff(1:length(iEC1372_W3110.genes), ia);
[iEC1372_W3110] = mergeModelFieldPositions(iEC1372_W3110, 'genes', strcmp(iEC1372_W3110.genes, iEC1372_W3110.genes{duplicate_idxs}));

% generate the grRules and the rxnGeneMat fields
iEC1356_Bl21DE3 = creategrRulesField(iEC1356_Bl21DE3);
iEC1356_Bl21DE3 = buildRxnGeneMat(iEC1356_Bl21DE3);
iEC1372_W3110 = creategrRulesField(iEC1372_W3110);
iEC1372_W3110 = buildRxnGeneMat(iEC1372_W3110);

% fix subsystems as they are missing from the models
iEC1356_Bl21DE3 = set_subSystems(iEC1356_Bl21DE3);
iEC1372_W3110 = set_subSystems(iEC1372_W3110);

% get all the reactions of the irreversible models
BL21_excRxns = iEC1356_Bl21DE3.rxns(findExcRxns(iEC1356_Bl21DE3));
W3110_excRxns = iEC1372_W3110.rxns(findExcRxns(iEC1372_W3110));
[rxns, idxs, ~] = unique([iEC1356_Bl21DE3.rxns; iEC1372_W3110.rxns]);
subSystems = [iEC1356_Bl21DE3.subSystems; iEC1372_W3110.subSystems];
subSystems = subSystems(idxs(idxs ~= 0));

% set nutritional constraints
[nutrient_rxns, ~] = define_nutrient_constraints();
iEC1356_Bl21DE3 = changeRxnBounds(iEC1356_Bl21DE3, nutrient_rxns, -10, 'l');
iEC1356_Bl21DE3 = changeRxnBounds(iEC1356_Bl21DE3, BL21_excRxns(not(ismember(BL21_excRxns, nutrient_rxns))), 0, 'l'); % block import of metabolites not in medium
iEC1356_Bl21DE3 = changeRxnBounds(iEC1356_Bl21DE3, BL21_excRxns, 10, 'u'); % limited export of all metabolites
iEC1372_W3110 = changeRxnBounds(iEC1372_W3110, nutrient_rxns, -10, 'l');
iEC1372_W3110 = changeRxnBounds(iEC1372_W3110, W3110_excRxns(not(ismember(W3110_excRxns, nutrient_rxns))), 0, 'l'); % block import of metabolites not in medium
iEC1372_W3110 = changeRxnBounds(iEC1372_W3110, W3110_excRxns, 10, 'u'); % limited export of all metabolites

% check that FBA works:
iEC1356_Bl21DE3 = changeObjective(iEC1356_Bl21DE3, 'BIOMASS_Ec_iJO1366_WT_53p95M');
FBAsolution = optimizeCbModel(iEC1356_Bl21DE3);
disp(['iEC1356_Bl21DE3 base biomass :  ', num2str(FBAsolution.f)])
iEC1372_W3110 = changeObjective(iEC1372_W3110, 'BIOMASS_Ec_iJO1366_WT_53p95M');
FBAsolution = optimizeCbModel(iEC1372_W3110);
disp(['iEC1372_W3110 base biomass :  ', num2str(FBAsolution.f)])

% compute data for gene set rule propagation
if isfile('..\data\processed\reaction_expression_data_BL21') && isfile('..\data\processed\reaction_expression_data_W3110')
    load reaction_expression_data_W3110.mat
    load reaction_expression_data_BL21.mat
else
    [reaction_expression_BL21, pos_genes_in_react_expr_BL21, ixs_geni_sorted_by_length_BL21] = compute_reaction_expression(iEC1356_Bl21DE3_irrev);
    save('..\data\processed\reaction_expression_data_BL21', 'reaction_expression_BL21', 'pos_genes_in_react_expr_BL21', 'ixs_geni_sorted_by_length_BL21')
    [reaction_expression_W3110, pos_genes_in_react_expr_W3110, ixs_geni_sorted_by_length_W3110] = compute_reaction_expression(iEC1372_W3110_irrev);
    save('..\data\processed\reaction_expression_data_W3110', 'reaction_expression_W3110', 'pos_genes_in_react_expr_W3110', 'ixs_geni_sorted_by_length_W3110')
end





%% Condition-specific simulations for pre-induction phase
%%
num_conditions = 22;
mu_fluxes = zeros(num_conditions, length(rxns));
s_flux_stds = zeros(num_conditions, length(rxns));
av_fluxes = zeros(num_conditions, length(rxns));
va_flux_stds = zeros(num_conditions, length(rxns));
predicted_atp = NaN(num_conditions, 1);
predicted_production = NaN(num_conditions, 1);
execution_times = NaN(num_conditions, 1);

if strcmp(process_phase, 'pre')
    for i = 1:num_conditions
        % generate condition-specific model
        if DoE_and_offline_data.Host(i) == 1
            bioreactor_model = iEC1356_Bl21DE3;
            bioreactor_model = change_enzyme_activity(bioreactor_model, critical_temperatures, 37, ...
                reaction_expression_BL21, pos_genes_in_react_expr_BL21, ixs_geni_sorted_by_length_BL21, ...
                0.00, 0.5, 0.5, 0);
        else
            bioreactor_model = iEC1372_W3110;
            bioreactor_model = change_enzyme_activity(bioreactor_model, critical_temperatures, 37, ...
                reaction_expression_W3110, pos_genes_in_react_expr_W3110, ixs_geni_sorted_by_length_W3110, ...
                0.00, 0.5, 0.5, 0);
        end

        % glycerol exchange constraints
        bioreactor_model = changeRxnBounds(bioreactor_model, 'EX_glyc_e', constraint_data.preinduction_glycerol_lb(i), 'l');
        bioreactor_model = changeRxnBounds(bioreactor_model, 'EX_glyc_e', constraint_data.preinduction_glycerol_ub(i), 'u');

        % other metabolites' exchange constraints (values +- std)
        bioreactor_model1 = bioreactor_model;
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, constraint_data.rxns, ...
            constraint_data.preinduction_exchange_rates(i, :) - constraint_data.preinduction_exchange_rates_errors(i, :), 'l');
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, constraint_data.rxns, ...
            constraint_data.preinduction_exchange_rates(i, :) + constraint_data.preinduction_exchange_rates_errors(i, :), 'u');

        % growth constraints (values +- std)
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, 'BIOMASS_Ec_iJO1366_WT_53p95M', ...
            constraint_data.preinduction_growth_rates(i) - constraint_data.preinduction_growth_rates_errors(i), 'l');
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, 'BIOMASS_Ec_iJO1366_WT_53p95M', ...
            constraint_data.preinduction_growth_rates(i) + constraint_data.preinduction_growth_rates_errors(i), 'u');
        
        % enforce high ATP production
        bioreactor_model1 = changeObjective(bioreactor_model1, 'ATPM');
        FBAsolution = optimizeCbModel(bioreactor_model1);
        bioreactor_model = changeRxnBounds(bioreactor_model, 'ATPM', 0.9*FBAsolution.f, 'l');

        % calculate flux distribution parameters
        rxns_to_constrain = [constraint_data.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M'];
        [~, b] = ismember(rxns_to_constrain, bioreactor_model.rxns);
        [~, sort_idx] = sort(b);
        exp_rates = [constraint_data.preinduction_exchange_rates(i, :), constraint_data.preinduction_growth_rates(i)];
        exp_rate_errors = [constraint_data.preinduction_exchange_rates_errors(i, :), constraint_data.preinduction_growth_rates_errors(i)];
        av_exp = exp_rates(sort_idx);
        va_exp = exp_rate_errors(sort_idx).^2;
        exp_i = find(ismember(bioreactor_model.rxns, rxns_to_constrain));
        [mu, s, ~, ~, av, va, ~, t_EP]  = MetabolicEP(full(bioreactor_model.S), bioreactor_model.b, bioreactor_model.lb, bioreactor_model.ub, ...
            Beta, damping, maxit, minvar, maxvar, precision, av_exp, va_exp, exp_i);

        [a, b] = ismember(rxns, bioreactor_model.rxns);
        mu_fluxes(i, a) = mu(b(a));
        s_flux_stds(i, a) = s(b(a));
        av_fluxes(i, a) = av(b(a));
        va_flux_stds(i, a) = va(b(a));
        execution_times(i) = t_EP;
        predicted_atp(i) = av(strcmp(bioreactor_model.rxns, 'ATPM'));
        predicted_production(i) = av(strcmp(bioreactor_model.rxns, 'PRODUCTtranslation2'));
    end
    
    % Save results:
    if save_results == true
        save('..\data\results\preinduction_MEP_fluxes.mat', 'mu_fluxes', 's_flux_stds', 'av_fluxes', 'va_flux_stds', 'rxns', 'subSystems')
        save('..\data\results\preinduction_MEP_data.mat')
    end





%% Condition-specific simulations for post-induction phase
%%
elseif strcmp(process_phase, 'post')
    % calculate [0,1]-normalised growth
    zgrowth = (constraint_data.postinduction_growth_rates - min(constraint_data.postinduction_growth_rates)) ./ ...
        (max(constraint_data.postinduction_growth_rates) - min(constraint_data.postinduction_growth_rates));

    % get product contribution to total protein
    total_protein_percent = titres.total_protein_percent([2:5 7:24]);

    % get total protein contribution to biomass
    protein_biomass_sum = abs(sum(iEC1356_Bl21DE3.S(ismember(iEC1356_Bl21DE3.mets, {'ala__L[c]','arg__L[c]','asn__L[c]','asp__L[c]', ...
        'cys__L[c]','glu__L[c]','gln__L[c]','gly[c]','his__L[c]','ile__L[c]','leu__L[c]','lys__L[c]','met__L[c]','phe__L[c]','pro__L[c]', ...
        'ser__L[c]','thr__L[c]','trp__L[c]','tyr__L[c]','val__L[c]'}), ...
        strcmp(iEC1356_Bl21DE3.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M'))));

    fluxes = zeros(num_conditions, length(rxns));
    predicted_atp_post = NaN(num_conditions, 1);
    predicted_production_post = NaN(num_conditions, 1);

    for i = 1:num_conditions
        % generate condition-specific model
        if DoE_and_offline_data.Host(i) == 1
            bioreactor_model = iEC1356_Bl21DE3;
            bioreactor_model = change_enzyme_activity(bioreactor_model, critical_temperatures, double(string(DoE_and_offline_data.T(i))), ...
                reaction_expression_BL21, pos_genes_in_react_expr_BL21, ixs_geni_sorted_by_length_BL21, ...
                0.00, 0.5, 0.5, 0);
        else
            bioreactor_model = iEC1372_W3110;
            bioreactor_model = change_enzyme_activity(bioreactor_model, critical_temperatures, double(string(DoE_and_offline_data.T(i))), ...
                reaction_expression_W3110, pos_genes_in_react_expr_W3110, ixs_geni_sorted_by_length_W3110, ...
                0.00, 0.5, 0.5, 0);
        end

        % glycerol exchange constraints
        bioreactor_model = changeRxnBounds(bioreactor_model, 'EX_glyc_e', constraint_data.postinduction_glycerol_lb(i), 'l');
        bioreactor_model = changeRxnBounds(bioreactor_model, 'EX_glyc_e', constraint_data.postinduction_glycerol_ub(i), 'u');

        % include the product in biomass
        bioreactor_model.S(strcmp(bioreactor_model.mets, 'product[c]'), strcmp(bioreactor_model.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M')) = ...
            (-0.01 - 0.39*(1-zgrowth(i)))*protein_biomass_sum;
        bioreactor_model.S(ismember(bioreactor_model.mets, {'ala__L[c]','arg__L[c]','asn__L[c]','asp__L[c]', ...
        'cys__L[c]','glu__L[c]','gln__L[c]','gly[c]','his__L[c]','ile__L[c]','leu__L[c]','lys__L[c]','met__L[c]','phe__L[c]','pro__L[c]', ...
        'ser__L[c]','thr__L[c]','trp__L[c]','tyr__L[c]','val__L[c]'}), ...
            strcmp(bioreactor_model.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M')) = ...
            bioreactor_model.S(ismember(bioreactor_model.mets, {'ala__L[c]','arg__L[c]','asn__L[c]','asp__L[c]', ...
        'cys__L[c]','glu__L[c]','gln__L[c]','gly[c]','his__L[c]','ile__L[c]','leu__L[c]','lys__L[c]','met__L[c]','phe__L[c]','pro__L[c]', ...
        'ser__L[c]','thr__L[c]','trp__L[c]','tyr__L[c]','val__L[c]'}), ...
            strcmp(bioreactor_model.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M')) * (0.6 + 0.39*zgrowth(i));

        % other metabolites' exchange constraints (values +- std)
        bioreactor_model1 = bioreactor_model;
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, constraint_data.rxns, ...
            constraint_data.postinduction_exchange_rates(i, :) - constraint_data.postinduction_exchange_rates_errors(i, :), 'l');
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, constraint_data.rxns, ...
            constraint_data.postinduction_exchange_rates(i, :) + constraint_data.postinduction_exchange_rates_errors(i, :), 'u');

        % growth constraints (values +- std)
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, 'BIOMASS_Ec_iJO1366_WT_53p95M', ...
            constraint_data.postinduction_growth_rates(i) - constraint_data.postinduction_growth_rates_errors(i), 'l');
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, 'BIOMASS_Ec_iJO1366_WT_53p95M', ...
            constraint_data.postinduction_growth_rates(i) + constraint_data.postinduction_growth_rates_errors(i), 'u');

        % enforce high ATP production
        bioreactor_model1 = changeObjective(bioreactor_model1, 'ATPM');
        FBAsolution = optimizeCbModel(bioreactor_model1);
        bioreactor_model = changeRxnBounds(bioreactor_model, 'ATPM', 0.8*FBAsolution.f, 'l');

        % calculate flux distribution parameters
        rxns_to_constrain = [constraint_data.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M'];
        [~, b] = ismember(rxns_to_constrain, bioreactor_model.rxns);
        [~, sort_idx] = sort(b);
        exp_rates = [constraint_data.postinduction_exchange_rates(i, :), constraint_data.postinduction_growth_rates(i)];
        exp_rate_errors = [constraint_data.postinduction_exchange_rates_errors(i, :), constraint_data.postinduction_growth_rates_errors(i)];
        av_exp = exp_rates(sort_idx);
        va_exp = exp_rate_errors(sort_idx).^2;
        exp_i = find(ismember(bioreactor_model.rxns, rxns_to_constrain));
        [mu, s, ~, ~, av, va, ~, t_EP]  = MetabolicEP(full(bioreactor_model.S), bioreactor_model.b, bioreactor_model.lb, bioreactor_model.ub, ...
            Beta, damping, maxit, minvar, maxvar, precision, av_exp, va_exp, exp_i);

        [a, b] = ismember(rxns, bioreactor_model.rxns);
        mu_fluxes(i, a) = mu(b(a));
        s_flux_stds(i, a) = s(b(a));
        av_fluxes(i, a) = av(b(a));
        va_flux_stds(i, a) = va(b(a));
        execution_times(i) = t_EP;
        predicted_atp(i) = av(strcmp(bioreactor_model.rxns, 'ATPM'));
        predicted_production(i) = av(strcmp(bioreactor_model.rxns, 'PRODUCTtranslation2'));
    end
    
    % Save results:
    if save_results == true
        save('..\data\results\postinduction_MEP_fluxes.mat', 'mu_fluxes', 's_flux_stds', 'av_fluxes', 'va_flux_stds', 'rxns', 'subSystems')
        save('..\data\results\postinduction_MEP_data.mat')
    end





%% Condition-specific simulations for end of process phase
%%
elseif strcmp(process_phase, 'end')
    % get product contribution to total protein
    total_protein_percent = titres.total_protein_percent([2:5 7:24]);

    for i = 1:num_conditions
        % generate condition-specific model
        if DoE_and_offline_data.Host(i) == 1
            bioreactor_model = iEC1356_Bl21DE3;
            bioreactor_model = change_enzyme_activity(bioreactor_model, critical_temperatures, double(string(DoE_and_offline_data.T(i))), ...
                reaction_expression_BL21, pos_genes_in_react_expr_BL21, ixs_geni_sorted_by_length_BL21, ...
                0.00, 0.5, 0.5, 0);
        else
            bioreactor_model = iEC1372_W3110;
            bioreactor_model = change_enzyme_activity(bioreactor_model, critical_temperatures, double(string(DoE_and_offline_data.T(i))), ...
                reaction_expression_W3110, pos_genes_in_react_expr_W3110, ixs_geni_sorted_by_length_W3110, ...
                0.00, 0.5, 0.5, 0);
        end

        % glycerol exchange constraints
        bioreactor_model = changeRxnBounds(bioreactor_model, 'EX_glyc_e', constraint_data.postinduction_glycerol_lb(i), 'l');
        bioreactor_model = changeRxnBounds(bioreactor_model, 'EX_glyc_e', constraint_data.postinduction_glycerol_ub(i), 'u');

        % include the product in biomass
        bioreactor_model.S(strcmp(bioreactor_model.mets, 'product[c]'), strcmp(bioreactor_model.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M')) = ...
            -total_protein_percent(i)*protein_biomass_sum/100;
        bioreactor_model.S(ismember(bioreactor_model.mets, {'ala__L[c]','arg__L[c]','asn__L[c]','asp__L[c]', ...
        'cys__L[c]','glu__L[c]','gln__L[c]','gly[c]','his__L[c]','ile__L[c]','leu__L[c]','lys__L[c]','met__L[c]','phe__L[c]','pro__L[c]', ...
        'ser__L[c]','thr__L[c]','trp__L[c]','tyr__L[c]','val__L[c]'}), ...
            strcmp(bioreactor_model.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M')) = ...
            bioreactor_model.S(ismember(bioreactor_model.mets, {'ala__L[c]','arg__L[c]','asn__L[c]','asp__L[c]', ...
        'cys__L[c]','glu__L[c]','gln__L[c]','gly[c]','his__L[c]','ile__L[c]','leu__L[c]','lys__L[c]','met__L[c]','phe__L[c]','pro__L[c]', ...
        'ser__L[c]','thr__L[c]','trp__L[c]','tyr__L[c]','val__L[c]'}), ...
            strcmp(bioreactor_model.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M')) * (100-total_protein_percent(i))/100;

        % other metabolites' exchange constraints (values +- std)
        bioreactor_model1 = bioreactor_model;
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, constraint_data.rxns, ...
            constraint_data.postinduction_exchange_rates(i, :) - constraint_data.postinduction_exchange_rates_errors(i, :), 'l');
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, constraint_data.rxns, ...
            constraint_data.postinduction_exchange_rates(i, :) + constraint_data.postinduction_exchange_rates_errors(i, :), 'u');

        % growth constraints (values +- std)
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, 'BIOMASS_Ec_iJO1366_WT_53p95M', ...
            constraint_data.postinduction_growth_rates(i) - constraint_data.postinduction_growth_rates_errors(i), 'l');
        bioreactor_model1 = changeRxnBounds(bioreactor_model1, 'BIOMASS_Ec_iJO1366_WT_53p95M', ...
            constraint_data.postinduction_growth_rates(i) + constraint_data.postinduction_growth_rates_errors(i), 'u');

        % enforce high ATP production
        bioreactor_model1 = changeObjective(bioreactor_model1, 'ATPM');
        FBAsolution = optimizeCbModel(bioreactor_model1);
        bioreactor_model = changeRxnBounds(bioreactor_model, 'ATPM', 0.8*FBAsolution.f, 'l');

        % calculate flux distribution parameters
        rxns_to_constrain = [constraint_data.rxns, 'BIOMASS_Ec_iJO1366_WT_53p95M'];
        [~, b] = ismember(rxns_to_constrain, bioreactor_model.rxns);
        [~, sort_idx] = sort(b);
        exp_rates = [constraint_data.postinduction_exchange_rates(i, :), constraint_data.postinduction_growth_rates(i)];
        exp_rate_errors = [constraint_data.postinduction_exchange_rates_errors(i, :), constraint_data.postinduction_growth_rates_errors(i)];
        av_exp = exp_rates(sort_idx);
        va_exp = exp_rate_errors(sort_idx).^2;
        exp_i = find(ismember(bioreactor_model.rxns, rxns_to_constrain));
        [mu, s, ~, ~, av, va, ~, t_EP]  = MetabolicEP(full(bioreactor_model.S), bioreactor_model.b, bioreactor_model.lb, bioreactor_model.ub, ...
            Beta, damping, maxit, minvar, maxvar, precision, av_exp, va_exp, exp_i);

        [a, b] = ismember(rxns, bioreactor_model.rxns);
        mu_fluxes(i, a) = mu(b(a));
        s_flux_stds(i, a) = s(b(a));
        av_fluxes(i, a) = av(b(a));
        va_flux_stds(i, a) = va(b(a));
        execution_times(i) = t_EP;
        predicted_atp(i) = av(strcmp(bioreactor_model.rxns, 'ATPM'));
        predicted_production(i) = av(strcmp(bioreactor_model.rxns, 'PRODUCTtranslation2'));
    end
    
    % Save results:
    if save_results == true
        save('..\data\results\end_MEP_fluxes.mat', 'mu_fluxes', 's_flux_stds', 'av_fluxes', 'va_flux_stds', 'rxns', 'subSystems')
        save('..\data\results\end_MEP_data.mat')
    end
end




%% Functions
%%

function [A_frequence, U_frequence, C_frequence, G_frequence] = get_RNA_base_frequences()
    % random proteic sequence of length 200 with uniform aminoacid frequency
    % http://molbiotools.com/randomsequencegenerator.html
    % corresponding most likely DNA sequence based on E. coli codon frequencies
    % http://www.bioinformatics.org/sms2/rev_trans.html

    sequence = ...
        ['auggaaugggugcaucauagccugccguuucuguuuguggaaaacaaacaggugccgggc' ...
        'cgcauguauaaacgcuggagcgugcugaugauggaauggcuggauauuuauuggguguuu' ...
        'uauguggauaacaaaguguuuauggcgcgcagcugccugauucauccguuuccggaaaug' ...
        'ugccgcauuaacugggaacagaacaugggcuuuacccuggcgcagggccuggaacaguuu' ...
        'cauaaaccgugcgauaaagugcgcgugccgaacagcgaaaaauuuuaucugcagcuguau' ...
        'uuucuggcggauacccuguuucaugaacugcuggaucgcauucuggcgcugguguuuacc' ...
        'aaaccgagcagcaaagugugguaugaagaacugcgcaacaugugcgaacauguggcgaaa' ...
        'auucaucauuaucuguaugaaaugcauaccuggaaccgcgaauuuaugugcuuuaaaccg' ...
        'auuuauuuuuuuggcaugauugugcgccugcgcgauuaucagccgggcugcggccagcug' ...
        'gaaaccugccagggcaacuuuugcuggcaucugaccauggauuggaaaaacaacgcgccg'];
    A_frequence = count(sequence, 'a');
    U_frequence = count(sequence, 'u');
    C_frequence = count(sequence, 'c');
    G_frequence = count(sequence, 'g');

end

function model = set_subSystems(model)
    % load table containing subsystem annotation from Monk et al 2016
    pathway_table = readtable('1-s2.0-S2405471216302903-mmc3.xlsx', 'Sheet', 'to_import_in_matlab');
    % load complete reaction table from BIGG to map old rxn IDs
    bigg_table = readtable('bigg_models_reactions.txt');
    % set missing pathway information
    pathway_table.Subsystem(contains(pathway_table.Reaction, 'EX_')) = {'Pseudo-reactions'};
    pathway_table.Subsystem(contains(pathway_table.Reaction, 'DM_')) = {'Pseudo-reactions'};
    pathway_table.Subsystem(strcmp(pathway_table.Reaction, 'Ec_biomass_iJO1366_WT_53p95M')) = {'Pseudo-reactions'};
    pathway_table.Subsystem(strcmp(pathway_table.Reaction, 'Ec_biomass_iJO1366_core_53p95M')) = {'Pseudo-reactions'};
    pathway_table.Subsystem(strcmp(pathway_table.Subsystem, '')) = {'Unassigned'};
    % homogenise reactions names with those in models
    for i = 1:length(pathway_table.Reaction)
        if contains(pathway_table.Reaction{i}, '_DASH_LPAREN_e_RPAREN_')
            s = strsplit(pathway_table.Reaction{i}, '_DASH_LPAREN_e_RPAREN_');
            pathway_table.Reaction{i} = [s{1} '_e'];
        end
        if contains(pathway_table.Reaction{i}, '_DASH_D_RPAREN_e_LPAREN_')
            s = strsplit(pathway_table.Reaction{i}, '_DASH_D_RPAREN_e_LPAREN_');
            pathway_table.Reaction{i} = [s{1} '__D_e'];
        end
        if contains(pathway_table.Reaction{i}, 'DASH_R_DASH_L_LPAREN_e_RPAREN_')
            s = strsplit(pathway_table.Reaction{i}, 'DASH_R_DASH_L_LPAREN_e_RPAREN_');
            pathway_table.Reaction{i} = [s{1} '_R__L_e'];
        end
        if contains(pathway_table.Reaction{i}, 'DASH_S_DASH_L_LPAREN_e_RPAREN_')
            s = strsplit(pathway_table.Reaction{i}, 'DASH_S_DASH_L_LPAREN_e_RPAREN_');
            pathway_table.Reaction{i} = [s{1} '_S__L_e'];
        end
        if contains(pathway_table.Reaction{i}, '_LPAREN_e_RPAREN_')
            s = strsplit(pathway_table.Reaction{i}, '_LPAREN_e_RPAREN_');
            pathway_table.Reaction{i} = [s{1} '_e'];
        end
        if contains(pathway_table.Reaction{i}, '_DASH_')
            s = strsplit(pathway_table.Reaction{i}, '_DASH_');
            pathway_table.Reaction{i} = [s{1} '__' s{2}];
        end
    end
    pathway_table.Reaction(strcmp(pathway_table.Reaction, 'Ec_biomass_iJO1366_WT_53p95M')) = {'BIOMASS_Ec_iJO1366_WT_53p95M'};
    pathway_table.Reaction(strcmp(pathway_table.Reaction, 'Ec_biomass_iJO1366_core_53p95M')) = {'BIOMASS_Ec_iJO1366_core_53p95M'};

    for i = 1:length(model.rxns)
        if contains(model.rxns{i}, '_copy1')
            s = strsplit(model.rxns{i}, '_copy1');
            model.rxns{i} = s{1};
        end
        if contains(model.rxns{i}, '_copy2')
            s = strsplit(model.rxns{i}, '_copy2');
            model.rxns{i} = s{1};
        end
    end
    
    % map old rxn IDs through the BIGG table
    idxs = find(not(ismember(model.rxns, pathway_table.Reaction)) & not(contains(model.rxns, 'PRODUCT')));
    for i = 1:length(idxs)
        idx = strcmp(bigg_table.bigg_id, model.rxns(idxs(i)));
        old_ids = strsplit(bigg_table.old_bigg_ids{idx}, '; ');
        [a, b] = ismember(old_ids, pathway_table.Reaction);
        if sum(a) > 1
            a = and(a, not(contains(old_ids, 'R_DM_')));
            pathway_table.Reaction(b(a)) = bigg_table.bigg_id(idx);
        else
            pathway_table.Reaction(b(a)) = bigg_table.bigg_id(idx);
        end
    end
    
    % assign subsystems to reactions
    [a, b] = ismember(model.rxns, pathway_table.Reaction);
    model.subSystems(a) = pathway_table.Subsystem(b(a));
    model.subSystems(contains(model.rxns, 'PRODUCT')) = {'Heterologous production'};
end

function [reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(model)

    genesets = model.grRules;
    genesets = regexprep(genesets,' AND ',' and '); 
    genesets = regexprep(genesets,' OR ',' or '); 

    reaction_expression = cell(length(genesets),1);
    reaction_expression(:) = {''};
    for i = 1:length(genesets)
        str_geneset = genesets{i};
        aux = associate_genes_reactions(str_geneset);
        reaction_expression{i} = aux; 
    end
    reaction_expression=strrep(reaction_expression,' ','');

    genes = model.genes;
    len = NaN(1, length(genes));
    for i = 1:length(genes)
        len(i)=length(genes{i});
    end
    [~, ixs_geni_sorted_by_length] = sort(len,'descend');
    reaction_expression_aux = reaction_expression;
    for i = 1:numel(ixs_geni_sorted_by_length)
        j = ixs_geni_sorted_by_length(i);
        matches = strfind(reaction_expression_aux,genes{j});
        pos_genes_in_react_expr{j} = find(~cellfun('isempty', matches));
        reaction_expression_aux(pos_genes_in_react_expr{j}) = strrep(reaction_expression_aux(pos_genes_in_react_expr{j}),genes{j},'');
    end

end

function str_output = associate_genes_reactions(str)
    % Sometimes grRules have no parentheses, which means that we need
    % to find a smart way to make sure that AND is solved before OR when we
    % substitute MIN and MAX respectively. This means that in the final expression, the MINs have to be
    % calculated before the MAXs.
    % To do so, we substitute the ORs first (which become MAXs), and then the ANDs inside
    % the MAXs. This is to ensure that we have an expression that first solves
    % the  ANDs (which are internal) and then solves the ORs (which are
    % external), thus respecting the common rule that AND is solved before OR
   
    while ( ~isempty(findstr(' or ', str)) ) % loops until all the AND and OR are not found because they have been substituted by MIN and MAX
    
        i = 1;
       while ( (strcmp(str(i:i+3),' or ')==0) )
            i = i+1; % while it does not find any 'and' and any 'or', keeps scrolling the array
        end
        
        str = substitute(str,i);
    end
 
    while ( ~isempty(findstr(' and ', str) ) ) % loops until all the AND and OR are not found because they have been substituted by MIN and MAX
    
        i = 1;
        
        while ( (strcmp(str(i:i+4),' and ')==0)    )
            i = i+1; % while it does not find any 'and' and any 'or', keeps scrolling the array
        end
        
        str = substitute(str,i);
    end
    
    if (isempty(findstr('max(', str)) && isempty(findstr('min(', str)) && ~strcmp(str,''))
        str_output = str;  
    else
        str_output = str; %if str is empty or is a combination of max/min of genes, we leave it as it is
    end
end

function str = substitute(str,i)
    i = i+1;
    % i is now positioned on the initial character of either 'and' or 'or'
    if (str(i)=='a')
        found_and = 1;
    else
        found_and = 0 ;
    end

    bracket_found = 0;
    j = i;
    while (  (strcmp(str(j),'(')==0) || (bracket_found~=-1) ) && (strcmp(str(j),',')==0 || (bracket_found~=0) ) && (j>1)
        j = j-1;
        if (str(j)==')')
            bracket_found = bracket_found+1; % signals further parentheses found along the path
        end
        if (str(j)=='(')
            bracket_found = bracket_found-1; % signals the closure of parentheses found along the path
        end
    end
    if (bracket_found == -1 || strcmp(str(j),',')~=0)
        j = j+1;
    end
    if (found_and == 1)
        k = i+3;
    else
        k = i+2;
    end

    bracket_found = 0;
    while ( (strcmp(str(k),')')==0) || ( bracket_found~=-1 )) && (strcmp(str(k),',')==0 || (bracket_found~=0) ) && (k<length(str))
        k = k+1;
        if (str(k)=='(')
            bracket_found = bracket_found+1; % signals further parentheses found along the path
        end
        if (str(k)==')')
            bracket_found = bracket_found-1; % signals the closure of parentheses found along the path
        end
    end
    if bracket_found == -1 ||  strcmp(str(k),',')~=0
        k = k-1;
    end

    if (found_and == 1)
        str_new = strrep( str, str(j:k), [' min(',str(j:i-1),', ',str(i+4:k), ') '] );
    else
        str_new = strrep( str, str(j:k), [' max(',str(j:i-1),', ',str(i+3:k), ') '] );
    end
    str = str_new;
end

function [uptakes, uptake_names] = define_nutrient_constraints()
    % source https://www.sigmaaldrich.com/catalog/product/sial/07533
    
    yeast_extract_uptakes = {'EX_thm_e', ... % vitamins
        'EX_nac_e', ...
        'EX_pnto__R_e', ...
        'EX_pydam_e', ...
        'EX_pydx_e', ...
        'EX_pydxn_e', ...
        'EX_btn_e', ...
        'EX_adocbl_e', ...
        'EX_cbi_e', ...
        'EX_cbl1_e', ...
        'EX_inost_e', ...
        'EX_ascb__L_e', ...
        'EX_chol_e', ...
        'EX_ala_B_e', ... % aminoacids
        'EX_ala__D_e', ...
        'EX_ala__L_e', ...
        'EX_arg__L_e', ...
        'EX_asn__L_e', ...
        'EX_asp__L_e', ...
        'EX_cys__D_e', ...
        'EX_cys__L_e', ...
        'EX_glu__L_e', ...
        'EX_gln__L_e', ...
        'EX_gly_e', ...
        'EX_his__L_e', ...
        'EX_ile__L_e', ...
        'EX_leu__L_e', ...
        'EX_lys__L_e', ...
        'EX_met__D_e', ...
        'EX_met__L_e', ...
        'EX_phe__L_e', ...
        'EX_pro__L_e', ...
        'EX_ser__D_e', ...
        'EX_ser__L_e', ...
        'EX_thr__L_e', ...
        'EX_trp__L_e', ...
        'EX_tyr__L_e', ...
        'EX_val__L_e', ...
        'EX_glc__D_e', ... % carbohydrates
        'EX_gal_e', ...
        'EX_gal_bD_e', ...
        'EX_fru_e', ...%         'EX_xyl__D_e', ...
        'EX_sucr_e', ...
        'EX_lcts_e', ...%         'EX_malt_e', ...
        'EX_tre_e', ...
        'EX_14glucan_e', ...
        'EX_pyr_e', ...
        'EX_gtp_e', ... % nucleic acids
        'EX_gsn_e', ...
        'EX_csn_e', ...
        'EX_ade_e', ...
        'EX_thym_e', ...
        'EX_gmp_e', ...
        'EX_dgmp_e', ...
        'EX_3gmp_e', ...
        'EX_cmp_e', ...
        'EX_dcmp_e', ...
        'EX_3cmp_e', ...
        'EX_amp_e', ...
        'EX_damp_e', ...
        'EX_3amp_e', ...
        'EX_dtmp_e', ...
        'EX_ump_e', ...
        'EX_dump_e', ...
        'EX_3ump_e', ...
        'EX_uri_e', ...
        'EX_ins_e', ...
        'EX_enlipa_e', ... % lipids
        'EX_hdca_e', ...
        'EX_hdcea_e', ...
        'EX_ocdca_e', ...
        'EX_ocdcea_e', ...
        'EX_colipa_e', ...
        'EX_colipap_e', ...
        'EX_etha_e', ... % other
        'EX_pheme_e', ...
        };
    
    yeast_extract_uptake_names = {'Thiamin exchange', ... % vitamins
        'Nicotinate exchange', ...
        '(R)-Pantothenate exchange', ...
        'Pyridoxamine exchange', ...
        'Pyridoxal exchange', ...
        'Pyridoxine exchange', ...
        'Biotin exchange', ...
        'Adenosylcobalamin exchange', ...
        'Cobinamide exchange', ...
        'Cob(I)alamin exchange', ...
        'Myo-Inositol exchange', ...
        'L-Ascorbate exchange', ...
        'Choline exchange', ...
        'Beta-Alanine exchange', ... % aminoacids
        'D-Alanine exchange', ...
        'L-Alanine exchange', ...
        'L-Arginine exchange', ...
        'L-Asparagine exchange', ...
        'L-Aspartate exchange', ...
        'D-Cysteine exchange', ...
        'L-Cysteine exchange', ...
        'L-Glutamate exchange', ...
        'L-Glutamine exchange', ...
        'Glycine exchange', ...
        'L-Histidine exchange', ...
        'L-Isoleucine exchange', ...
        'L-Leucine exchange', ...
        'L-Lysine exchange', ...
        'D-Methionine exchange', ...
        'L-Methionine exchange', ...
        'L-Phenylalanine exchange', ...
        'L-Proline exchange', ...
        'D-Serine exchange', ...
        'L-Serine exchange', ...
        'L-Threonine exchange', ...
        'L-Tryptophan exchange', ...
        'L-Tyrosine exchange', ...
        'L-Valine exchange', ...
        'D-Glucose exchange', ... % carbohydrates
        'D-Galactose exchange', ...
        'Beta D-Galactose exchange', ...
        'D-Fructose exchange', ...
        'D-Xylose exchange', ...
        'Sucrose exchange', ...
        'Lactose exchange', ...
        'Maltose exchange', ...
        'Trehalose exchange', ...
        '1,4-alpha-D-glucan exchange', ...
        'GTP exchange', ... % nucleic acids
        'Guanosine exchange', ...
        'Cytosine exchange', ...
        'Adenine exchange', ...
        'Thymine exchange', ...
        'GMP exchange', ...
        'DGMP exchange', ...
        '3-GMP exchange', ...
        'CMP exchange', ...
        'DCMP exchange', ...
        '3-CMP exchange', ...
        'AMP exchange', ...
        'DAMP exchange', ...
        '3-AMP exchange', ...
        'DTMP exchange', ...
        'UMP exchange', ...
        'DUMP exchange', ...
        '3-UMP exchange', ...
        'Uridine exchange', ...
        'Inosine exchange', ...
        'Phosphoethanolamine KDO(2)-lipid (A) exchange', ... % lipids
        'Hexadecanoate (n-C16:0) exchange', ...
        'Hexadecenoate (n-C16:1) exchange', ...
        'Octadecanoate (n-C18:0) exchange', ...
        'Octadecenoate (n-C18:1) exchange', ...
        'Core oligosaccharide lipid A exchange', ...
        'Core oligosaccharide lipid A diphosphate exchange', ...
        'Ethanolamine exchange', ... % other
        'Protoheme exchange', ...
        };
    
    medium_uptakes = {'EX_so4_e', ... % ions
        'EX_k_e', ...
        'EX_pi_e', ...
        'EX_na1_e', ...
        'EX_cit_e', ... % other
        'EX_glyc_e', ...
        };
    
    medium_uptake_names = {'Sulfate exchange', ... % ions
        'K+ exchange', ...
        'Phosphate exchange', ...
        'Sodium exchange', ...
        'Citrate exchange', ... % other
        'Glycerol exchange', ...
        };
    
    feed_uptakes = {'EX_nh4_e', ... % ions
        'EX_h2o_e',... % other
        };
    
    feed_uptake_names = {'Ammonia exchange', ... % ions
        'H2O exchange', ... % other
        };
    
    other_uptakes = {'EX_o2_e', ...
        'EX_ca2_e', ...
        'EX_h_e', ...
        'EX_mn2_e', ...
        'EX_cu_e', ...
        'EX_cu2_e', ...
        'EX_fe2_e', ...
        'EX_fe3_e', ...
        'EX_mg2_e', ...
        'EX_zn2_e', ...
        'EX_cl_e', ...
        'EX_co2_e', ...
        'EX_cobalt2_e', ...
        'EX_mobd_e', ...
        'EX_ni2_e', ...
        'EX_sel_e', ...
        'EX_slnt_e', ...
        'EX_tungs_e', ...
        };
    
    other_uptake_names = {'O2 exchange', ...
        'Calcium exchange', ...
        'H+ exchange', ...
        'Mn2+ exchange', ...
        'Cu+ exchange', ...
        'Cu2+ exchange', ...
        'Fe2+ exchange', ...
        'Fe3+ exchange', ...
        'Mg exchange', ...
        'Zinc exchange', ...
        'Chloride exchange', ...
        'CO2 exchange', ...
        'Co2+ exchange', ...
        'Molybdate exchange', ...
        'Ni2+ exchange', ...
        'Selenate exchange', ...
        'Selenite exchange', ...
        'Tungstate exchange', ...
        };
    
    uptakes = [yeast_extract_uptakes medium_uptakes feed_uptakes other_uptakes]';
    uptake_names = [yeast_extract_uptake_names medium_uptake_names feed_uptake_names other_uptake_names]';
end

function model = change_enzyme_activity(model, T_table, T, input_reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length, ...
    tc, tf, tm, th)
    % tc, tf, tm and th are parameters to regulate fractional activity at
    % Tc, Tf, Tm and Th respectively
    
    % sort T_table so as to have genes in model.genes order
    enzyme = model.genes;
    Tc = NaN(length(model.genes), 1);
    Tf = NaN(length(model.genes), 1);
    To = NaN(length(model.genes), 1);
    Tm = NaN(length(model.genes), 1);
    Th = NaN(length(model.genes), 1);
    crit_temp = table(enzyme, Tc, Tf, To, Tm, Th, 'RowNames', model.genes);
    [a, b] = ismember(model.genes, T_table.enzyme);
    crit_temp(a, :) = T_table(b(a), :);
    T_table = crit_temp;
    clear crit_temp
    
    T = T*ones(size(T_table, 1), 1);
    
    case1 = T <= T_table.Tc;
    case2 = T > T_table.Tc & T <= T_table.Tf;
    case3 = T > T_table.Tf & T <= T_table.To;
    case4 = T > T_table.To & T <= T_table.Tm;
    case5 = T > T_table.Tm & T <= T_table.Th;
    case6 = T > T_table.Th;
    case_missing_info = not(case1 | case2 | case3 | case4 | case5 | case6); % some T is unavailable
    case7 = case_missing_info & T > T_table.Tc & T <= T_table.To; % Tc < T <= To with no Tf
    case8 = case_missing_info & T > T_table.To & T <= T_table.Th; % To < T <= Th with no Tm
    case9 = case_missing_info & T <= T_table.Tf; % T <= Tf with no Tc
    case10 = case_missing_info & T > T_table.Tm; % T > Tm with no Th
    case11 = case_missing_info & not(case7 | case8 | case9 | case10);
        
    x = ones(length(model.genes), 1);
    x(case1) = tc;
    x(case2) = (tf*T_table.Tc(case2)+(tc-tf)*T(case2)-tc*T_table.Tf(case2)) ./ (T_table.Tc(case2)-T_table.Tf(case2));
    x(case3) = (T_table.Tf(case3)+(tf-1)*T(case3)-tf*T_table.To(case3)) ./ (T_table.Tf(case3)-T_table.To(case3));
    x(case4) = (T_table.Tm(case4)+(tm-1)*T(case4)-tm*T_table.To(case4)) ./ (T_table.Tm(case4)-T_table.To(case4));
    x(case5) = (tm*T_table.Th(case5)+(tm-th)*T(case5)-th*T_table.To(case5)) ./ (T_table.Th(case5)-T_table.Tm(case5));
    x(case6) = th;
    x(case7) = (T_table.Tc(case7)+(tc-1)*T(case7)-tc*T_table.To(case7)) ./ (T_table.Tc(case7)-T_table.To(case7));
    x(case8) = (T_table.Th(case8)+(th-1)*T(case8)-th*T_table.To(case8)) ./ (T_table.Th(case8)-T_table.To(case8));
    x(case9) = tf;
    x(case10) = tm;
    x(case11) = 1;
        
    num_reaction_expression = zeros(length(input_reaction_expression), 1);
    % searches for all the strings of the kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
    to_replace = ['min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.|min..\d*+\.+\d*.,\d*+\.+\d*.|max..\d*+\.+\d*.,\d*+\.+\d*.|' ...
        'min..\d*+\.+\d*.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.\d*+\.+\d*..|min.\d*+\.+\d*,.\d*+\.+\d*..|max.\d*+\.+\d*,.\d*+\.+\d*..'];
    to_replace_nan = ['|min.NaN,\d*+\.+\d*.|min.\d*+\.+\d*,NaN.|max.NaN,\d*+\.+\d*.|max.\d*+\.+\d*,NaN.|min..NaN.,\d*+\.+\d*.|min..\d*+\.+\d*.,NaN.|' ...
        'max..NaN.,\d*+\.+\d*.|max..\d*+\.+\d*.,NaN.|min..NaN.,.\d*+\.+\d*..|min..\d*+\.+\d*.,.NaN..|max..NaN.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.NaN..|' ...
        'min.NaN,.\d*+\.+\d*..|min.\d*+\.+\d*,.NaN..|max.NaN,.\d*+\.+\d*..|max.\d*+\.+\d*,.NaN..' ...
        '|min.NaN,NaN.|max.NaN,NaN.|min..NaN.,NaN.|max..NaN.,NaN.|min..NaN.,.NaN..|max..NaN.,.NaN..|min.NaN,.NaN..|max.NaN,.NaN..'];
    to_replace = [to_replace to_replace_nan];

    reaction_expression = input_reaction_expression;
    for i = ixs_geni_sorted_by_length
        gene_positions = pos_genes_in_react_expr{i};
        for j = 1:length(gene_positions) % for each string of reaction_expression, we replace the substring 'b****' with the number representing its enzyme activity
            reaction_expression{gene_positions(j)} = strrep(reaction_expression{gene_positions(j)}, model.genes{i}, num2str(x(i),'%.15f'));
        end
    end
    reaction_expression(cellfun(@isempty, reaction_expression)) = {'1.0'};
    
    for i = 1:length(num_reaction_expression)
        str = reaction_expression{i};
        num_parenthesis = numel(strfind(str, ')'));
        while (num_parenthesis > 32) % if there are more than 32 parentheses, matlab is unable to run EVAL, so we need to reduce these parentheses manually by starting to eval smaller pieces of the string
            substrings_to_replace = regexp(str, to_replace, 'match');
            if isempty(substrings_to_replace)
                num_parenthesis = 0;
            else
                for j = 1:numel(substrings_to_replace)
                    ss_rep = substrings_to_replace{j};
                    str = strrep(str, ss_rep, num2str(eval(ss_rep), '%.15f'));
                end
                num_parenthesis = numel(strfind(str, ')'));
            end
        end
        str = regexprep(str, '/', '');
        num_reaction_expression(i) = eval(str); % evaluates the cells like they are numerical expressions (so as to compute min and max of enzyme activities)
    end
    
    model.lb = model.lb.*num_reaction_expression;
    model.ub = model.ub.*num_reaction_expression;
end
