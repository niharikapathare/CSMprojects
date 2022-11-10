function expanded = expand(K,k,dof)
K(dof,dof) = K(dof,dof) + k;
expanded = K;