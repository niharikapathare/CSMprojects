function [R,v] = align(p,a)

pm = (sum(p,1)/length(p))';
am = (sum(a,1)/length(a))';
pa = zeros(3,3);
for i=1:length(p),
    pa = pa + p(i,:)'*a(i,:);
end
M = (pa/length(p)) - pm*am';
[v,d] = eig(M'*M);
D11 = sqrt(d(3,3));
D22 = sqrt(d(2,2));
V = -[v(:,3) v(:,2) v(:,1)];
m = M*V;
R = [m(:,1)/D11 m(:,2)/D22 cross(m(:,1),m(:,2))/D11/D22]*V';
v = pm - R*am;
