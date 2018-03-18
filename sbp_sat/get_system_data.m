function[system_data] = get_system_data(n_equ)

value_x = ones(n_equ);
system_data.Ax = spdiags([value_x], [0], n_equ, n_equ);

value_y = ones(n_equ);
system_data.Ay = spdiags([value_y], [0], n_equ, n_equ);

end