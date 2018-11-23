function [convg_odd,convg_even] = convg_rate_macro(error,M_values)

num_quantities = size(error,2);
convg_odd = zeros(num_quantities,1);
convg_even = zeros(num_quantities,1);

for i = 1 : num_quantities
    
    tempOdd = polyfit(log(M_values(1:2:end)),log(error(1:2:end,i))',1);
    convg_odd(i) = abs(tempOdd(1));
    
    tempEven = polyfit(log(M_values(2:2:end)),log(error(2:2:end,i))',1);
    convg_even(i) = abs(tempEven(1));
end

end

