function k = kmat(c,E,t)

% integrate each element in [k] with 1-pt quadrature
[b J] = bmat(c);
k = t*(b.'*E*b*J*0.5);