
addpath(genpath('/data/'));

fclose all; close all; clear all; clc;

%% Prep

% Load glasser parcellation and fsaverage template space
load('fsaverage_7.1.0.mat');
load('glasser_fsaverage_7.1.0.mat');

%% Compute indices for spin-based spatial permutations

% % Read the centroid coordinates of Glasser parcellation (Vasa et al.)
% C = readmatrix('/data/Vasa_Glasser_Permutations/sphere_HCP.txt');
% CL = C(1:180,:);
% CR = C(181:end,:);
% 
% % Perform the permutations
% CP = rotate_parcellation(CL, CR, 10000);
% 
% % Save the computed permutations
% save('CP.mat', 'CP', '-v7.3');

% Load the pre-computed permutations
load('CP.mat');

%% Read the cortical maps for oxyphos and evoexp

% oxphos
table = readtable('oxphos_mean.csv');
oxphos = table.mean;

% evoexp
table = readtable('brainmaps_glasser.csv');
evoexp = table.Evolution_Expansion;
evoexp = [evoexp; evoexp]; % Mirroring on the other hemisphere

% Checking the maps
h = HyoSurfacePlotFunction(oxphos, IS, G360, flip(cbrewer2('PuOr', 256)), 'oxphos', [0.3384 0.7384], 'square'); saveas(h, 'oxphos.tif');
h = HyoSurfacePlotFunction(evoexp, IS, G360, flip(cbrewer2('PuOr', 256)), 'evoexp', [-2 2], 'square'); saveas(h, 'evoexp.tif');

%% Perform statistics

% Compute the permuted maps for oxphos
for i = 1:10000
    oxphos_perm(:,i) = oxphos(CP(:,i));
end

% Compute spatial correlation between oxphos and evoexp
r_oxphos_evoexp = corr(oxphos, evoexp);

% Perform spin test
[p_oxphos_evoexp_perm r_oxphos_evoexp_perm] = doPerms(oxphos_perm, evoexp, r_oxphos_evoexp, 'top');

% Display
h = figure('WindowState', 'maximized'); hold on;
b = histogram(r_oxphos_evoexp_perm);
b.NumBins = 200;
b.FaceColor = [0.5 0.5 0.5];
xlim([-0.5 0.5]);
ylim([0 100]);
xlabel('Cortex-wide spatial correlation (r)');
s = stem(r_oxphos_evoexp, 80, 'o', 'LineWidth', 1.25, 'MarkerSize', 8)
s.Color = [0 0 0];
s.MarkerFaceColor = [0 0 0];
legend({'Null', ['r = ' num2str(r_oxphos_evoexp) '(p\_spin = ' num2str(p_oxphos_evoexp_perm) ')']});
saveas(h, 'spin_test_oxphos_evoexp.tif');
