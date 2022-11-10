function [Dmag_peak,vM_peak] = femur_fe_v4(infile,mod,Fx,Fy,fig)
% Inputs:
% infile - string holding name of input file containint triangle mesh
% mod - modulus of elasticity for the bone
% Fx - x-comp of the force applied on the femoral head
% Fy - y-comp of the force applied on the femoral head
% fig - Boolean to turn output figure on (1) or off (0)
% Outputs:
% Dmag_peak - peak value of the displacement magnitude
% vM_peak - peak value of the von Mises stress magnitude

% Poisson ratio
nu = 0.35;
% constant thickness in (mm)
t = 5;
% plane strain constitutive matrix
E = (mod/(1+nu)/(1-2*nu))*[1-nu nu 0;nu 1-nu 0;0 0 (1-2*nu)/2];

% node coordinates
c = dlmread(infile,',',[1 1 771 2]);
e = dlmread(infile,',',[773 1 2111 3]);
[dof,col] = size(c); dof = dof*2;
[ele,col] = size(e);

% assembly
K = zeros(dof,dof);
for i=1:ele,
    cel = [c(e(i,1),:) c(e(i,2),:) c(e(i,3),:)];
    k = kmat(cel,E,t);
    K = expand(K,k,[2*e(i,1)-1 2*e(i,1) 2*e(i,2)-1 2*e(i,2) 2*e(i,3)-1 2*e(i,3)]);
end

% displacement boundary conditions
Kbig = 10e6*max(max(K));
Kmod = K;
% what nodes to fix in x-dir
xfix = [68:112];
% what nodes to fix in y-dir
yfix = [68:112];
for i=xfix,
    Kmod(2*i-1,2*i-1) = Kmod(2*i-1,2*i-1) + Kbig;
end
for i=yfix,
    Kmod(2*i,2*i) = Kmod(2*i,2*i) + Kbig;
end

% load boundary conditions
xnode = [184];
xforc = Fx;
ynode = [184];
yforc = Fy;
Pmod = zeros(dof,1); 
for i=1:length(xnode),
    Pmod(2*xnode(i)-1) = xforc(i);
end
for i=1:length(ynode),
    Pmod(2*ynode(i)) = yforc(i);
end

% solve
D = inv(Kmod)*Pmod;
R = K*D;
% deformed coordinates
Dnodal = reshape(D,2,dof/2)';
cd = c + Dnodal;
% normalized displacement magnitude at each node
for i=1:length(Dnodal),
    Dmag(i,1) = norm(Dnodal(i,:));
end
Dnorm = normc(Dmag);
% find the max location
[val,dindex] = max(Dmag); Dmag_peak = val;

% constant stresses
for i=1:ele,
    % strains: exx, eyy, exy
    ee(:,i) = bmat([c(e(i,1),:) c(e(i,2),:) c(e(i,3),:)])*...
        D([2*e(i,1)-1 2*e(i,1) 2*e(i,2)-1 2*e(i,2) 2*e(i,3)-1 2*e(i,3)]);
    % stresses: sxx, syy, sxy
    S(i,:) = E*ee(:,i);
    % stress in the z-dir
    Sz(i,1) = (mod/(1+nu)/(1-2*nu))*nu*(sum(ee(:,i)));
    % von Mises stress
    vM(i,1) = sqrt(S(i,1)^2 + S(i,2)^2 + Sz(i)^2 - S(i,1)*S(i,2) -Sz(i)*(S(i,1)+S(i,2)) + 3*S(i,3)^2);
end

% average von Mises stress at each node
vMcat = cell(dof/2,1);
for i=1:ele,
    vMcat{e(i,1)} = [vMcat{e(i,1)} vM(i)];
    vMcat{e(i,2)} = [vMcat{e(i,2)} vM(i)];
    vMcat{e(i,3)} = [vMcat{e(i,3)} vM(i)];
end
for i=1:dof/2,
    vMave(i,1) = mean(vMcat{i});
end
% normalized vM at each node
vMnorm = normc(vMave);
% find the max location
[val,vmindex] = max(vMave); vM_peak = val;

if fig ~= 0,
    % plot nominal and deformed shape with contours of displacement magnitude
    subplot(1,2,1); triplot(e,c(:,1),c(:,2),'k--');
    hold on; axis equal; axis off; axis([-120 90 -220 240]);
    title(['Displacement Magnitude, max = ',num2str(max(Dmag)),' mm']);
    p.Vertices = cd; p.Faces = e;
    p.FaceColor = 'interp'; p.FaceVertexCData = Dnorm;
    patch(p);
    % mark the max location
    plot(cd(dindex,1),cd(dindex,2),'g.','MarkerSize',25);
    
    % plot nominal and deformed shape with contours of vM stress
    subplot(1,2,2); triplot(e,c(:,1),c(:,2),'k--');
    hold on; axis equal; axis off; axis([-120 90 -220 240]);
    title(['von Mises Stress, max = ',num2str(max(vMave)),' MPa']);
    p.Vertices = cd; p.Faces = e;
    p.FaceColor = 'interp'; p.FaceVertexCData = vMnorm;
    patch(p);
    % mark the max location
    plot(cd(vmindex,1),cd(vmindex,2),'g.','MarkerSize',25);
end

