function[system_data] = get_system_data(n_equ)

value_x = ones(n_equ)
system_data.Ax = spdiags([value], [0], n_equ, n_equ);

value_y = ones(n_equ)
system_data.Ay = spdiags([value], [0], n_equ, n_equ);

% system_data.Ax = full(system_data.Ax);
% system_data.Ay = full(system_data.Ay);

end