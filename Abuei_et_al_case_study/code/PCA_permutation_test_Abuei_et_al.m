
clear
load end_pFBA_fluxes.mat
fluxes = fluxes(1:18, :);

num_permutations = 999;

lambdas = nan(num_permutations, 50);
for i = 1:num_permutations
    permuted_fluxes = nan(size(fluxes));
    for j = 1:size(fluxes, 2) % for each data dimension
        permutation_idx = randperm(size(fluxes, 1)); % shuffle the data values of each axis
        permuted_fluxes(:, j) = fluxes(permutation_idx, j);
    end
    [~,~,latent,~,~,~] = pca(zscore(permuted_fluxes));
    lambdas(i, :) = latent';
end

[~,~,latent,~,~,~] = pca(zscore(fluxes));
% Rnd-Lambda test in https://doi.org/10.1016/j.csda.2004.06.015
p_lambdas = (sum(lambdas >= repmat(latent', num_permutations, 1)) + 1) ./ (num_permutations + 1);

save('..\data\results\PCA_permutation_test')