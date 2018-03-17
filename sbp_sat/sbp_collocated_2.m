function[D, P, g?, S?] = get_system_data(n_equ)

d = ones(N+1, 1); % for diagonal entries of D matrix; length = len of diag
D = spdiags([-d d], [-1 1], N+1, N+1);
D(1,1) = -2    ; D(N+1,N) = -2 ;
D(N+1,N+1) = 2 ; D(1,2) = 2 ; % efficient -> modify values at -1, 0 and 1
% diagonals and then create sparse matrix
D = (1/(2*h))*D;

P = spdiags([d], [0], N+1, N+1);
P(1,1) = 1/2 ; P(N+1,N+1) = 1/2 ;
P = h*P;

% boundary at x=0 during t=0
g = zeros(1, m);
g(1,1) = solution(grid(1), 0);
g(1,2) = solution2(grid(1), 0);

S = zeros(N+1, m);
S(1,1) = -1 * (u0(1,1)-g(1,1));
S(1,2) = -1 * (u0(1,1)-g(1,1));