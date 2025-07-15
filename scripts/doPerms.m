function [p_perm r_perm] = doPerm(permMap, map, r_true, mode)

% This function computes the permutation statistics based on Pearson
% correlation.

% Check input
if (size(map,1) ~= 1) & (size(map,2) ~= 1)
    error('Variable map needs to be an 1d array');
end
nPerm = size(permMap,2);

% Reshape
if size(map,2) > size(map,1)
    map = map';
end

% Calculate correlation
for i = 1:nPerm
    r_perm(i) = corr(permMap(:,i), map);
end

% Calculate p
switch mode
    case 'top'
        p_perm = sum(r_true < r_perm) ./ nPerm;
    case 'bottom'
        p_perm = sum(r_true > r_perm) ./ nPerm;
end

end
