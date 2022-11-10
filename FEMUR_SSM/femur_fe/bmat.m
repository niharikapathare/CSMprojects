function [B J] = bmat(c)
% c is a vector of x,y coordinates [x1 y1 x2 y2 x3 y3]
% params is a vector of r,s natural coordinates
% B is the strain-displacement matrix
% J is the determinant of the Jacobian

J = norm(cross([c(3)-c(1) c(4)-c(2) 0],[c(5)-c(3) c(6)-c(4) 0]));
B = [c(4)-c(6) 0 c(6)-c(2) 0 c(2)-c(4) 0;...
     0 c(5)-c(3) 0 c(1)-c(5) 0 c(3)-c(1);...
     c(5)-c(3) c(4)-c(6) c(1)-c(5) c(6)-c(2) c(3)-c(1) c(2)-c(4)]/J;
