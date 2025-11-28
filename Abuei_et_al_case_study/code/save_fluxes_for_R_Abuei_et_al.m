load('end_pFBA_fluxes.mat')
fluxes = fluxes(1:18, :);
idx = sum(diff(fluxes) == 0) == 17;
fluxes = fluxes(:, not(idx));
rxns = rxns(not(idx));
subSystems = subSystems(not(idx));
t = array2table(fluxes);
t.Properties.VariableNames = cellfun(@strcat, repmat({'v_'}, length(rxns), 1), rxns, 'UniformOutput', false);
writetable(t, '..\data\results\end_pFBA_fluxes_exploration.xlsx')