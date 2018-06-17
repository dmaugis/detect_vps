function [L] = line_to_homogeneous(l);
% converts line in x1 y1 x2 y2 format to homogeneous coordinates

x1 = l(:,1);
y1 = l(:,2);
x2 = l(:,3);
y2 = l(:,4);


dx = x2-x1;
dy = y2-y1;
a = -dy;
b = dx;
c = -a.*x1 -b.*y1;

L = [a, b, c];