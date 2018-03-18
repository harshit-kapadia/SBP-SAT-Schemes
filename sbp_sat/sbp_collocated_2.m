function[D, P] = sbp_collocated_2(N, h)
% N - no. of grid cells & h - grid size

d = ones(N+1, 1); % for diagonal entries of D matrix; length = len of diag
D = spdiags([-d d], [-1 1], N+1, N+1);
D(1,1) = -2    ; D(N+1,N) = -2 ;
D(N+1,N+1) = 2 ; D(1,2) = 2 ; % efficient -> modify values at -1, 0 and 1
% diagonals and then create sparse matrix
D = (1/(2*h))*D;

P = spdiags([d], [0], N+1, N+1);
P(1,1) = 1/2 ; P(N+1,N+1) = 1/2 ;
P = h*P;

end