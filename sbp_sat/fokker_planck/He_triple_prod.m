% computes \int He_m(\xi) He_l(\xi) He_n(\xi) f0 d\xi. He is normalised
% such that \int He_m He_n f0 d\xi = delta_{mn}

function output = He_triple_prod(l,m,n)
output = 0;

if(mod(l+m+n,2)==0)
    if(l+m>=n && m+n>=l && n+l>=m)
        output = sqrt(factorial(l)*factorial(m)*factorial(n)) / ...
            ( factorial_diff(l,m,n) * factorial_diff(l,n,m) * factorial_diff(n,m,l));
    end
end

end

function f = factorial_diff(l,m,n)
f = factorial((l+m-n)/2);
end