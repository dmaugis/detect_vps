% Benchmark VP detection on Eurasian Cities Dataset
%
% This program is written by Jose Lezama <jlezama@gmail.com> and
% distributed under the terms of the GPLv3 license.
%
% Eurasian Cities Dataset (ECD) by Olga Barinova
% http://graphics.cs.msu.ru/en/research/projects/msr/geometry
%
% Reported performance should be close to  89.8% for the training set 
% (first 25 images) and 89.6% for the test set

clear
close all

tic

% point the following path to your local copy of ECD
DATASET_PATH = '/Users/jose/Documents/cmla/datasets/20130829_EurasianCitiesBase';

dirlist = dir(sprintf('%s/*.jpg',DATASET_PATH));

params.REFINE_THRESHOLD =  .62;% 0.375;
params.VARIATION_THRESHOLD = .3;%.305;
params.DUPLICATES_THRESHOLD = .0001;%.0039;
params.DIST_THRESHOLD = .165;
params.INFINITY_THRESHOLD = 3.58; % (lambda)

focal_ratio = 1.08;% 1.05;


all_errors = zeros(size(dirlist));

PLOT = 0;

for f = 26:103%1:length(dirlist)
    % load image and ground truth
    img_filename = sprintf('%s/%s',DATASET_PATH, dirlist(f).name);
    mat_filename = sprintf('%s/%shor.mat',DATASET_PATH, dirlist(f).name(1:end-4));
    
    img = imread(img_filename);
    [H, W, ~] = size(img);
    
    if PLOT
        img = imread(img_filename);
        [H, W, ~] = size(img);
        figure, imagesc(img)
    end
    
    manhattan = 0;
    acceleration = 1;    
    
    % compute  horizon line detection
    horizon_line_ours = detect_vps(img_filename, '.', manhattan, acceleration, focal_ratio, params);
    
    
    % load ground truth
    load(mat_filename)
    
    % compute horizon line estimation error
    X = [1 W];
    P_gt = [-horizon(1)/horizon(2) -horizon(3)/horizon(2)];
    Y_gt = polyval(P_gt,X);
    Y_ours = [horizon_line_ours(1) horizon_line_ours(W)];
    err = max(abs(Y_gt-Y_ours))/H;
    
    if PLOT
        hold on
        plot(X, Y_gt);
    end
    
    fprintf('----------------------------------------\n')
    fprintf('File %i, error: %f\n', f, err);
    fprintf('----------------------------------------\n')
    
    all_errors(f) = err;
end

% compute AUC score for training  and testing set
auc_train = calc_auc(all_errors(1:25));
auc_test = calc_auc(all_errors(26:end));

fprintf('AUCtrain: %2.4f. AUCtest:%2.4f\n', auc_train, auc_test);

% should output something close to AUCtrain: 0.8940. AUCtest:0.8970

toc
