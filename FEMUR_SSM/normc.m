function y = normc(x)
[r,c] = size(x);
for ii=1:c,
    y(:,ii) = x(:,ii)/norm(x(:,ii));
end