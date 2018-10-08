function penalty = dvlp_penalty_kinetic(A,M)

Amod_kinetic = compute_Amod_kinetic(M);
penalty = (A - Amod_kinetic)/2;

end

function Amod_kinetic = compute_Amod_kinetic(M)
Amod_1D = dvlp_Amod_kinetic_1D(M+1);   

all_idx = cell(M+1,1);

%% collect the multi indices
for i = 0 : M
    [all_idx{i+1},~,~] = IDX_Full(i);
end

%% exploit orthogonality in the other two directions
% extract the Hermite polynomial degrees along the x direction
all_idx_x_dir = cell2mat(cellfun(@(a) a(:,1),all_idx,'Un',0));

Amod_kinetic = Amod_1D(all_idx_x_dir + 1,all_idx_x_dir + 1);

end

function modA_kinetic = dvlp_Amod_kinetic_1D(M)
    
    poly_degree = 0:1:M;
    mat_half_int = zeros(length(poly_degree),length(poly_degree));
    
    for i = 1 : length(poly_degree) + 1
        for j = 1 : length(poly_degree) + 1
            poly_test = i - 1;
            poly_sol = j - 1;
            mat_half_int(i,j) = HermiteHalfSpace(poly_test,poly_sol);
        end
    end
    
    modA_kinetic = zeros(length(poly_degree),length(poly_degree));

    for i = 1 : length(poly_degree)
        for j = 1 : length(poly_degree)
            poly_degree_test = i-1;
            poly_degree_sol = j-1;
            
            % else the entry is zero
            if(mod(poly_degree_test-poly_degree_sol,2) == 0 && j > 1)
                modA_kinetic(i,j) = 2 * (sqrt(poly_degree_sol) * mat_half_int(i,j-1)+ ...
                                         sqrt(poly_degree_sol + 1) * mat_half_int(i,j+1));
            end
            
            if(mod(poly_degree_test-poly_degree_sol,2) == 0 && j == 1)
                modA_kinetic(i,j) = 2 * sqrt(poly_degree_sol + 1) * mat_half_int(i,j+1);
            end
        end
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
loc_neg = D<-(1e-10);
Xminus = V(:,loc_neg);
end