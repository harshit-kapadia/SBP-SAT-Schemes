% Gives set of multi-index of the basis with magnitude M
function [index, total] = multi_index(M)
count = 1;
for i = 0 : M
    for j = 0 : i
        index(count,:) = [M - i, (i - j), j];
        count = count + 1;
    end
end
total = count - 1 ; 
end