clear vars; close all; clc;

% load all nodal data for 20 files
dat_raw = [];
for i=1:20,
    infile = ['mesh_files/femur_',num2str(i),'.inp'];
    tmp = dlmread(infile,',',[1 1 771 2]);
    dat_raw = [dat_raw [tmp(:,1);tmp(:,2)]];
end
[datr,datc] = size(dat_raw);
% only need to load element data once because topology is identical
e = dlmread(infile,',',[773 1 2111 3]);
[eler,elec] = size(e);

% register all shapes to first shape: remove translation and rotation
dat = register2D(dat_raw);

% plot the first unregistered shape
color = 'grcmywgrcmywgrcmywgrcmyw';
subplot(1,2,1); triplot(e,dat_raw(1:datr/2,1),dat_raw(datr/2+1:datr,1),'k--');
hold on; axis equal; axis off; title('Unregistered Raw Data');
% plot the first registered shape
subplot(1,2,2); triplot(e,dat(1:datr/2,1),dat(datr/2+1:datr,1),'k--');
hold on; axis equal; axis off; title('Registered Data');
% plot all remaining shapes for inspection
for j=2:datc,
    % add specimen to the unregistered plot
    subplot(1,2,1); triplot(e,dat_raw(1:datr/2,j),dat_raw(datr/2+1:datr,j),color(j));
    % add specimen to the registered plot
    subplot(1,2,2); triplot(e,dat(1:datr/2,j),dat(datr/2+1:datr,j),color(j));
end

% compute mean shape
dat_mean = mean(dat,2);
% subtract mean shape from each specimen
for k=1:datc,
    dat_mod(:,k) = dat(:,k) - dat_mean;
end

% covariance and eigen analysis

C = dat_mod'*dat_mod;
[V,D] = eig(C);
V = normc(dat_mod*V); % normc() normalizes each column (eigenvector)
D = D/(datc-1);

% extract eigenvalues above 60
[d,v,cvar] = extract_eig(D,V,60);
cvar

% find principal components
b = v'*dat_mod;
b_std = std(b');

% Mode Shape Plots
figure;
% make overlay plot of mode 1
mode1_hi = dat_mean + b_std(1)*v(:,1);
mode1_lo = dat_mean - b_std(1)*v(:,1);
subplot(1,4,1); triplot(e,mode1_hi(1:datr/2,1),mode1_hi(datr/2+1:datr,1),'r--');
hold on; axis equal; axis off; title('Mode 1: Blue(-), Red(+)');
triplot(e,mode1_lo(1:datr/2,1),mode1_lo(datr/2+1:datr,1),'b--');
% show mean shape
triplot(e,dat_mean(1:datr/2,1),dat_mean(datr/2+1:datr,1),'k--');

% make overlay plot of mode 2
mode2_hi = dat_mean + b_std(2)*v(:,2);
mode2_lo = dat_mean - b_std(2)*v(:,2);
subplot(1,4,2); triplot(e,mode2_hi(1:datr/2,1),mode2_hi(datr/2+1:datr,1),'r--');
hold on; axis equal; axis off; title('Mode 2: Blue(-), Red(+)');
triplot(e,mode2_lo(1:datr/2,1),mode2_lo(datr/2+1:datr,1),'b--');
% show mean shape
triplot(e,dat_mean(1:datr/2,1),dat_mean(datr/2+1:datr,1),'k--');

% make overlay plot of mode 3
mode3_hi = dat_mean + b_std(3)*v(:,3);
mode3_lo = dat_mean - b_std(3)*v(:,3);
subplot(1,4,3); triplot(e,mode3_hi(1:datr/2,1),mode3_hi(datr/2+1:datr,1),'r--');
hold on; axis equal; axis off; title('Mode 3: Blue(-), Red(+)');
triplot(e,mode3_lo(1:datr/2,1),mode3_lo(datr/2+1:datr,1),'b--');
% show mean shape
triplot(e,dat_mean(1:datr/2,1),dat_mean(datr/2+1:datr,1),'k--');

% make overlay plot of mode 4
mode4_hi = dat_mean + b_std(4)*v(:,4);
mode4_lo = dat_mean - b_std(4)*v(:,4);
subplot(1,4,4); triplot(e,mode4_hi(1:datr/2,1),mode4_hi(datr/2+1:datr,1),'r--');
hold on; axis equal; axis off; title('Mode 4: Blue(-), Red(+)');
triplot(e,mode4_lo(1:datr/2,1),mode4_lo(datr/2+1:datr,1),'b--');
% show mean shape
triplot(e,dat_mean(1:datr/2,1),dat_mean(datr/2+1:datr,1),'k--');

% save SSM for use in Prob simulation
save SSM_femur_data e dat_mean b b_std v







