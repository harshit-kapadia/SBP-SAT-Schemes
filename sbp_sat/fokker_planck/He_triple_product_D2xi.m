% computes \int He_m(\xi) He''_{Dxi1}(\xi) He_n(\xi) f0 d\xi

function [output] = He_triple_product_D2xi(D2xi_index,m,n)
l = D2xi_index-2;

if (l < 0)
    output = 0;
else
    output = sqrt(D2xi_index * (D2xi_index-1)) * He_triple_prod(l,m,n);
end

end