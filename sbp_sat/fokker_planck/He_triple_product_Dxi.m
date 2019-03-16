% computes \int He_m(\xi) He'_{Dxi1}(\xi) He_n(\xi) f0 dxi
function [output] = He_triple_product_Dxi(Dxi_index,m,n)
l = Dxi_index-1;

if (l < 0)
    output = 0;
else
    output = sqrt(Dxi_index) * He_triple_prod(l,m,n);
end

end

