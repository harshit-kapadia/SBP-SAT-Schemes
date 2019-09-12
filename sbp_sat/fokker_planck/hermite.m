% computes normalized probabilists' Hermite polynomial of order n at x
function output = hermite(n,x)
output = hermiteH(n,x/sqrt(2)) / (2^(n/2) * sqrt(factorial(n))) ;
end