function [penalty] = dvlp_penalty_char(An,B)
modAn = compute_Amod(An);
Xminus = compute_Xminus(An);

if norm(full(B)) == 0
    penalty = Xminus * 0;
else
    penalty = 0.5 * (An-modAn) * Xminus * inv(B*Xminus);
end
end

function modA = compute_Amod(A)

% eig does not support sparse matrices
[V,D] = eig(full(A));
modA = V * abs(D)/V;
end

% compute the eigenvector corresponding to the negative eigenvalues
function Xminus = compute_Xminus(A)
[V,D] = eig(full(A));
D= D(sub2ind(size(D),1:size(D,1),1:size(D,2)));
loc_neg = find(D<-(1e-10));
Xminus = V(:,loc_neg);
end
