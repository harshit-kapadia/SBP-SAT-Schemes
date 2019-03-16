function [output] = He_triple_product_D2(D2_index,m,n)
output = 0 ;
l = D2_index-2 ;
if(mod(l+m+n,2)==0)
    if(l+m>=n && m+n>=l && n+l>=m)
        output = D2_index * (D2_index-1) * sqrt(factorial(l)*factorial(m)*factorial(n)) / ...
            ( factorial((l+m+n)/2) * factorial((m+n-l)/2) * factorial((n+l-m)/2) ) ;
    end
end
end