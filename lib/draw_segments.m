function img2 = draw_segments(img,vpimg, lines, params)
% Overlays line segments and vanishing directions on image
% Author: Jose Lezama <jlezama@gmail.com>

DEBUG = 0;

W = params.W;
H = params.H;


THRESHOLD = params.REFINE_THRESHOLD;
Z = size(lines,1);
V = size(vpimg,2);

colors = hsv(V);


zi = cell(1,V);
for i=1:V
    zi{i}=zeros(H,W);
end

zi2 = cell(1,V);
for i=1:V
    zi2{i}=zeros(H,W);
end


assign = ones(Z,V)*inf;

figure, imagesc(img)
hold on

all_js = [];
for i=1:V
    vp = vpimg(:,i)';
    %%%%%
    %         for j=1:Z
    %             lt = lines(j,:);
    %
    %             % get line from segment middle point to vp
    %             mp =[lt(1)+(lt(3)-lt(1))/2, lt(2)+(lt(4)-lt(2))/2];
    %
    %             mp_vp = get_intersection([mp 1], [vp 1]);
    %
    %             a = mp_vp(1);
    %             b = mp_vp(2);
    %
    %             % get angle betwen lines
    %             lt2 = [0 -1/b W -W*a/b-1/b];
    %
    %
    %             A = lt(3:4)-lt(1:2);
    %             B = lt2(3:4)-lt2(1:2);
    %
    %
    %             angle = acos(dot(A,B)/norm(A)/norm(B));
    %             angle = min(angle, pi-angle);
    %
    %             angle = angle*180/pi;
    %
    %              if angle< THRESHOLD
    %                                   assign(j,i) = angle;
    %                                   all_js=[all_js j];
    %
    %                                   plot([lt(1) lt(3)],[lt(2) lt(4)],'Color',colors(i,:),'LineWidth',2)
    %
    %                                   a = (lt(4)-lt(2))/(lt(3)-lt(1));
    %                                   b = lt(2)-a*lt(1);
    %                                   x=[1 W];
    %                                   plot(x,x*a+b,':','Color',colors(i,:));
    % %
    % %                 assign(j,i) = angle;
    % %
    % %                 lt(1) = max(1,lt(1));
    % %                 lt(2) = max(1,lt(2));
    % %                 lt(3) = max(1,lt(3));
    % %                 lt(4) = max(1,lt(4));
    % %
    % %                 z = drawline([lt(2), lt(1)],[lt(4) lt(3)],[H W 3]);
    % %
    % %                % make lines fatter
    % %                zi{i}(z)=1;
    %              end
    %
    %         end
    %%%%%%%
    
    
    %     %%%%%%%%%%%%
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
    
    assign(:,i) = angle;
    
    %%%%%%%%%%
    
end

% find best fitting vp for each line, but only those under the threshold
[angles I] = min(assign',[],1);
z = find(angles<THRESHOLD);

if DEBUG
    figure, imagesc(img)
    title('debug')
    hold on
end

for i=1:length(z)
    lt = lines(z(i),:);
    
    if DEBUG
        plot([lt(1) lt(3)],[lt(2) lt(4)],'Color',colors(I(z(i)),:),'LineWidth',2)
    end
    
    try
        l = drawline([lt(2), lt(1)],[lt(4) lt(3)],[H W 3]);
        zi{I(z(i))}(l)=1;
    end
    
    % draw dashed lines
    a = ((lt(4)-lt(2)))/(lt(3)-lt(1));
    b = lt(2)-a*lt(1);
    ang=atan(a);
    
    if abs(ang*180/pi)>45 | isnan(ang)
        % inverse stuff
        
        Y=[1:H];
        S = 3;
        z0 = find(mod(Y,S)<1);
        Y = round(Y(z0)+mod(lt(2),S));
        X = round((Y-b)/a);
    else
        
        X=[1:W];
        S = 5;
        z0 = find(mod(X,S)<1);
        X = round(X(z0)+mod(lt(1),S));
        Y = round(X*a+b);
    end
    z1 = find(Y>=1 & Y<=H & X>=1 & X<=W);
    
    z2 = sub2ind([H W], Y(z1),X(z1));
    zi2{I(z(i))}(z2)=1;
    %     z2 = sub2ind(size(img), Y(z1),X(z1),ones(size(Y(z1)))*2);
    %     img(z2)=0;
    %     z2 = sub2ind(size(img), Y(z1),X(z1),ones(size(Y(z1)))*3);
    %     img(z2)=0;
    
    
end


if size(img,3)==1;
    img = repmat(img,[1 1 3]);
end

red = img(:,:,1);
green = img(:,:,2);
blue = img(:,:,3);


% make lines fatter
se = strel('disk',1);
for i=1:V
    zi{i}=logical(imdilate(zi{i},se));
    red(zi{i})= colors(i,1)*255;
    green(zi{i})= colors(i,2)*255;
    blue(zi{i})= colors(i,3)*255;
    zi2{i}=logical(zi2{i});
    red(zi2{i})= colors(i,1)*255;
    green(zi2{i})= colors(i,2)*255;
    blue(zi2{i})= colors(i,3)*255;
end


img2 = zeros(size(img));
img2(:,:,1) = red;
img2(:,:,2) = green;
img2(:,:,3) = blue;
img2 = uint8(img2);