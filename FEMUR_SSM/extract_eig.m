function [d,v,cvar] = extract_eig(D,V,min)
[row,col] = size(D);
ev = D(row,col);
i = 1;
while ev > min,
    d(i,1) = ev;
    v(:,i) = V(:,col);
    row = row - 1;
    col = col - 1;
    i = i + 1;
    if row == 0,break;end
    ev = D(row,col);
end
% collect the cumulative variance info for each mode
for i=1:length(d),
    cvar(i,1) = i;
    cvar(i,2) = sum(d(1:i))/trace(D);
end

