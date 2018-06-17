function mvp_all = refine_detections(mvp_all, lines_lsd, params);
% Refines VP detections using lines from LSD
%
% Version 0.6, July 2015
%
% This program is written by Jose Lezama <jlezama@gmail.com>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

D = size(mvp_all,2);

mvp_refined = zeros(D,2);
for i =1:D
    vp = mvp_all(:,i)';
    vp = refine_vp(lines_lsd, vp, params);
    mvp_refined(i,:) = vp;
end

mvp_all = mvp_refined';
end

%% iterator function
function  vp = refine_vp( lines,  vp, params);

% Given a cluster of line segments, and two segments indicated by p1 and
% p2, obtain the main vanishing point determined by the segments

THRESHOLD = params.REFINE_THRESHOLD;

H = params.H;
W = params.W;

vp_orig = vp;
vp_tmp = vp;

[vp] = refine_vp_iteration(lines,  vp, THRESHOLD,H, W);


variation = norm(vp-vp_orig)/norm(vp_orig);
fprintf('variation: %f\n', variation);
if (variation > params.VARIATION_THRESHOLD)
    fprintf('vp changed too much (%f)... probably unstable conditions, keeping initial vp\n', variation)
    vp = vp_orig;
end

end

%% iteration function
function [vp] = refine_vp_iteration(lines,  vp, THRESHOLD, H,  W)
% finds intersection of each line in cluster with vanishing point segments

z=1:length(lines);

vp_orig = vp;


z2= [];

X = [0 W];

mp =[lines(:,1)+(lines(:,3)-lines(:,1))/2, lines(:,2)+(lines(:,4)-lines(:,2))/2];

L = size(lines,1);
O = ones(L,1);
Z = zeros(L,1);
vpmat = my_repmat2(vp,[L 1]);


VP = my_cross([mp O], [vpmat O]);
VP3 = my_repmat(VP(:,3),[1 3]);
VP = VP./VP3;

mp_vp = [VP(:,1) VP(:,2)];


a = VP(:,1);
b = VP(:,2);

% get angle betwen lines
lt2 = [Z -1./b W*O -W*a./b-1./b];

A = lines(:,3:4)-lines(:,1:2);
B = lt2(:,3:4)-lt2(:,1:2);

normA = sqrt(A(:,1).^2+A(:,2).^2);
normB = sqrt(B(:,1).^2+B(:,2).^2);

A = A./my_repmat(normA,[1 2]);
B = B./my_repmat(normB,[1 2]);

angle = acos(dot(A',B')');
angle = real(angle); % numerical errors
angle = min(angle, pi-angle);

angle = angle*180/pi;

z2 = find(angle<THRESHOLD);

Z = length(z);



%% obtain a refined VP estimate from sub-cluster z2
lengths = sqrt(sum(([lines(:,1) lines(:,2)] - [lines(:,3) lines(:,4)]).^2,2));
weights = lengths/max(lengths);
lis=line_to_homogeneous(lines);

Z2 = length(z2);

Q=zeros(3);
Is = [1 0 0; 0 1 0; 0 0 0];

%
l2 = lis(z2,:)';
w2 = weights(z2).^2;
w2 = my_repmat(w2,[1 3])';

b = dot((l2'*Is)',l2)'; %diag(l2'*Is*l2);
b = my_repmat(b,[1 3])';
Q = (l2./b.*w2)*l2';


p = [0 0 1]';

A = [2*Q -p];

vp = null(A);


vp = vp(1:2,1)/vp(3,1);
vp = vp';

end



