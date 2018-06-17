
function [Y, vpoimg] = compute_horizon_line_manhattan(mvp_all, NFAs, lines_lsd, params);
% computes horizontal line from vps and using the NFA values to apply
% orthogonality constraints. saves data to output image and output text
% file
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

H = params.H;
W = params.W;

%% York Urban parameters (given)
% focal = 6.05317058975369;
% pixelSize = 0.00896875;
% pp = [307.551305282635, 251.454244960136];

pp = [W/2.08095361329015 H/1.90889599050553]; % parameters from YUD
FOCAL_RATIO = params.FOCAL_RATIO; % in YUD = 1.05455933619402


vp_all_unit = [mvp_all(1,:)-pp(1);(H- mvp_all(2,:))-(H-pp(2)); ones(size(mvp_all(2,:)))*W*FOCAL_RATIO];


% normalize vectors;
for i=1:size(vp_all_unit,2)
    vp_all_unit(:,i)=vp_all_unit(:,i)/norm(vp_all_unit(:,i));% * sign(vp_all_unit(3,i));
end

vp_all_unit(isnan(vp_all_unit))=1;

my_vps = vp_all_unit;

%% impose orthogonality
my_orthogonal_vps = orthogonal_triplet(my_vps, NFAs, params.ORTHOGONALITY_THRESHOLD);

if size(my_orthogonal_vps,2) == 1
    my_orthogonal_vps = repmat(my_orthogonal_vps,[1 3]);
end
%

if size(my_orthogonal_vps,2) == 2
    
    fprintf('obtained only 2 VP... estimating the third one...\n')
    estimated_vp = cross(my_orthogonal_vps(:,2),my_orthogonal_vps(:,1));
    
    %% refine estimation
    % *** evpimg = 1/pixelSize*[focal*estimated_vp(1)/estimated_vp(3)+pp(1)*pixelSize; focal*estimated_vp(2)/estimated_vp(3) + (H-pp(2))*pixelSize];
    evpimg = [W*FOCAL_RATIO*estimated_vp(1)/estimated_vp(3)+pp(1); W*FOCAL_RATIO*estimated_vp(2)/estimated_vp(3) + (H-pp(2))];
    
    evpimg(2) = H - evpimg(2);
    
    % evpimg_refined = refine_vp(lines_lsd ,  evpimg', params);
    evpimg_refined = refine_detections(evpimg, lines_lsd , params);
    
    
    evpimg_refined = evpimg_refined';
    
    evpimg_refined = evpimg;
    
    % *** estimated_vp1 = pixelSize*[evpimg_refined(1)-pp(1);(H- evpimg_refined(2))-(H-pp(2)); focal/pixelSize];
    estimated_vp1 = [evpimg_refined(1)-pp(1);(H- evpimg_refined(2))-(H-pp(2)); W*FOCAL_RATIO];
    
    estimated_vp1 = estimated_vp1/norm(estimated_vp1);
    
    if  abs(dot(estimated_vp1, my_orthogonal_vps(:,1)))< params.ORTHOGONALITY_THRESHOLD && abs(dot(estimated_vp1, my_orthogonal_vps(:,2))) < params.ORTHOGONALITY_THRESHOLD
        estimated_vp = estimated_vp1;
    end
    
    
    %% add to current list of vps
    my_orthogonal_vps = [my_orthogonal_vps estimated_vp];
    
    
    
end

mvp=my_orthogonal_vps;

mvp1=mvp(:,1);
mvp2=mvp(:,2);
mvp3=mvp(:,3);
%
% guarantee norm 1
mvp1 = mvp1/norm(mvp1);
mvp2 = mvp2/norm(mvp2);
mvp3 = mvp3/norm(mvp3);

vpo1img = [W*FOCAL_RATIO*mvp1(1)/mvp1(3)+pp(1); W*FOCAL_RATIO*mvp1(2)/mvp1(3) + (H-pp(2))];
vpo2img = [W*FOCAL_RATIO*mvp2(1)/mvp2(3)+pp(1); W*FOCAL_RATIO*mvp2(2)/mvp2(3) + (H-pp(2))];
vpo3img = [W*FOCAL_RATIO*mvp3(1)/mvp3(3)+pp(1); W*FOCAL_RATIO*mvp3(2)/mvp3(3) + (H-pp(2))];


vpoimg = [vpo1img vpo2img vpo3img];


vpoimg(2,:) = H-vpoimg(2,:);


%% which one is the vertical vanishing point?

% calculate angle with image
vpoimg_centered = vpoimg - [W/2 W/2 W/2; H/2 H/2 H/2];
angles = atan2(vpoimg_centered(2,:), vpoimg_centered(1,:));
angles = mod(angles,2*pi);

[V I_vert] = min(abs(cos(angles)))

I_hor = setdiff([1:3], I_vert);


vpoimg = vpoimg(:, [I_hor(1) I_vert I_hor(2)]);



%% get horizon line
P_ours = polyfit([vpoimg(1,1) vpoimg(1,3)], [vpoimg(2,1) vpoimg(2,3)],1);

X = [1:W];
Y = polyval(P_ours,X);

%% plot
img = imread(params.img_in);

% draw segments with colors
if params.PRINT

    
    imorig=img;
    
    img = draw_segments(img,vpoimg, lines_lsd, params);


    
    % draw horizon line
    linewidth = 4;
    step = 10;
    for i=1:length(X);
        if mod(i,20)>=10
            col = [255 128 0];
        else
            col = [128 0 255];
        end
        
        for j=-round(linewidth/2):round(linewidth/2)
            try % sometimes falls outside image
                img(round(Y(i)+j),round(X(i)),:)=col;
            end
        end
    end
    
    if ~isdeployed && params.PLOT
        
        figure, imagesc(img);
    end
    
    % write output
    imwrite(img,sprintf('%s/horizon_line.png',params.folder_out));
end

fileID = fopen(sprintf('%s/out.txt',params.folder_out),'w');
fprintf(fileID,'horizontal vps:\n');
fprintf(fileID,'(%f, %f)\n',vpoimg(1,1), vpoimg(2,1));
fprintf(fileID,'(%f, %f)\n',vpoimg(1,3), vpoimg(2,3));
fprintf(fileID,'vertical vp:\n');
fprintf(fileID,'(%f, %f)\n',vpoimg(1,2), vpoimg(2,2));
fprintf(fileID,'horizon line:\n');
fprintf(fileID,'(%f, %f), (%f, %f)\n',X(1), Y(1), X(end), Y(end));
fclose(fileID);

