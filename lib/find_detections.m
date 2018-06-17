function [detections, m, b] = find_detections(points, params);
% runs alignment detector
% add mixtures code to path
if ~isdeployed
    addpath mixtures/
    addpath mex_files/
end


%% check that points exist
if isempty(points)
    detections = [];
    m = [];
    b = [];
    return
end

%% check if acceleration is needed, and activated
if size(points,1)> params.MAX_POINTS_ACCELERATION && params.ACCELERATION
    disp('- USE ACCELERATED VERSION -');
    ALGORITHM_VERSION =2;
else
    ALGORITHM_VERSION = 1;
end

%% normalize points

M = max(points(:));
m = min(points(:));

points = (points-m)/(M-m)*512; % this version of the alignment detector expects a 512x512 domain


if ALGORITHM_VERSION==1
    detections = alignments_slow(points);
    detections = detections';
else
    candidates = run_mixtures(points, params.GMM_Ks, '');
    detections = alignments_fast(points,candidates);
    detections = detections';
end

% convert detections to line coordinates
if ~isempty(detections)
    dets = detections(:, [1:4]);
    dets = dets/512*(M-m)+m;
    
    detections(:,1:4) = dets;
    
    detections(:,5) = detections(:,5)/512*(M-m);
    
    x1 = dets(:,1);
    y1 = dets(:,2);
    x2 = dets(:,3);
    y2 = dets(:,4);
    
    dy = y2 - y1;
    dx = x2 - x1;
    
    m = dy./dx;
    b = (y1.*x2 -y2.*x1)./dx;
else
    dets=[];
    m=[];
    b=[];
end

