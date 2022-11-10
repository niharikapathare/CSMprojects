function outfile = SSM_femur_instance(fname)
load SSM_femur_data
[datr,~] = size(dat_mean);
[eler,~] = size(e);
rng('shuffle');
std_vals = randn(4,1);
randb = b_std'.*std_vals;
dat = dat_mean + v*randb;
outfile = fname;
dlmwrite(outfile,'*Node','delimiter','');
dlmwrite(outfile,[(1:datr/2)' reshape(dat,datr/2,2)],'-append');
dlmwrite(outfile,'*Element','-append','delimiter','');
dlmwrite(outfile,[(1:eler)' e],'-append');