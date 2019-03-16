function [output] = He_triple_product_D1(D1_index,m,n)
output = 0 ;
l = D1_index-1 ;
if(mod(l+m+n,2)==0)
    if(l+m>=n && m+n>=l && n+l>=m)
        output = D1_index * sqrt(factorial(l)*factorial(m)*factorial(n)) / ...
            ( factorial((l+m+n)/2) * factorial((m+n-l)/2) * factorial((n+l-m)/2) ) ;
    end
end
end