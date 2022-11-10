function dat_reg = register2D(dat_raw)
[datr,datc] = size(dat_raw);
dat_reg = dat_raw;
a = reshape(dat_raw(:,1),datr/2,2);
a = [a zeros(datr/2,1)];
for j=2:datc,
    % find fit for specimen
    p = reshape(dat_raw(:,j),datr/2,2);
    p = [p zeros(datr/2,1)];
    [R,v] = align(p,a);
    % transform specimen to align with first shape
    q = zeros(length(p),3);
    for i=1:length(p),
        q(i,:) = (R'*(p(i,:)'-v))';
    end
    dat_reg(:,j) = [q(:,1);q(:,2)];
end
