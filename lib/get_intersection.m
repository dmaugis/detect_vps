function vp = get_intersection(L1,L2);
% gets point of intersection between lines in homogeneous coordinates

VP = cross(L1,L2);
VP = VP/VP(3);

vp = [VP(1) VP(2)];