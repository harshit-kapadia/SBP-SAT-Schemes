function error = compute_error(result_dvm,result_mom,PX)

num_mom = length(result_mom);
num_macros = size(result_dvm,1)-2;

% num of mom times the num of macroscopic quantitites
error = zeros(num_mom,num_macros);

for i = 1 : num_mom
    for j = 1 : num_macros
        temp = result_dvm(j + 2,:) - result_mom{i}(j  + 2, :);
        temp = reshape(temp,size(PX));
        int_x = dot(transpose(temp),transpose(PX * temp),2);  % integral along x.
        int_xy = sum(PX*int_x);  % integral along xy.
        error(i,j) = sqrt(int_xy);
    end    
end

end