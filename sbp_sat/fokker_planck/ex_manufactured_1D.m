function [] = ex_manufactured_1D(M)

par = struct(...
    'name','manufactured_1D',... % name of example
    'compute_density',@compute_density,...
    'compute_ux',@compute_ux,...
    'compute_uy',@compute_uy,...
    'compute_theta',@compute_theta,...
    'compute_sigma_xx',@compute_sigma_xx,...
    'compute_qx',@compute_qx,...
    'compute_qy',@compute_qy,...
    );

par.M = M ;

% yet to complete

end

%========================================================================
% Problem specific functions
%========================================================================

function output = compute_density(coefficients)
output = 0 ; % yet to complete
end