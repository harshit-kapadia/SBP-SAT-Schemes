% write the triple product matrices to files
% M_max is the maximum value of Hermite polynomial we store
function [] = write_triple_prod(M_max)

result_Dxi = cell(M_max + 1,1);
result_D2xi = cell(M_max + 1,1);

for i = 0 : M_max
    
    if i >= 1 % no non-trival entries else
        result_Dxi{i+1} = compute_triple_prod_Dxi(i,M_max);
    end
    
    if i >= 2 % no non-trivial entries else
        result_D2xi{i+1} = compute_triple_prod_D2xi(i,M_max);
    end
    
    filename_Dxi = strcat('triple_prod/triple_prod_Dxi_He',num2str(i),'.txt');
    filename_D2xi = strcat('triple_prod/triple_prod_D2xi_He',num2str(i),'.txt');
    
    [row,col,v] = find(result_Dxi{i+1});
    dlmwrite(filename_Dxi,[col row v], 'delimiter', '\t','precision',16);
    [row,col,v] = find(result_D2xi{i+1});
    dlmwrite(filename_D2xi,[col row v], 'delimiter', '\t','precision',16);
end

end

function output = compute_triple_prod_Dxi(i,M_max)
output = zeros(M_max+1,M_max+1);

for j = 0 : M_max
    for k = 0 : M_max
        output(j+1,k+1) = He_triple_product_Dxi(i,j,k);
    end
end

output = sparse(output);

end

function output = compute_triple_prod_D2xi(i,M_max)
output = zeros(M_max+1,M_max+1);

for j = 0 : M_max
    for k = 0 : M_max
        output(j+1,k+1) = He_triple_product_D2xi(i,j,k);
    end
end

output = sparse(output);

end
