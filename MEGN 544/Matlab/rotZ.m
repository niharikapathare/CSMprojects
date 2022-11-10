%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function name: rotZ
%Returns rotation matrix describing rotation about Z axis

%[R] = rotZ(theta)

%R = the rotation matrix describing rotation about Z axis for a given theta

%theta = input angle rotated about the Z axis in radians


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function R = rotZ(theta)
R =[cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0;0 0 1];
end
