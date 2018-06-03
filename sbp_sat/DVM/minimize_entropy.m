% we minize the entropy given the moments to find the maxwellian
function f = minimize_entropy(Ax,Ay,mass_matrix,inv_mass_matrix,...
                              rho,ux,uy,theta,id_sys,value_f0,id)

vx = Ax(id,id);
vy = Ay(id,id);

points_x = size(ux,1);
points_y = size(ux,2);
lagrange = cell(4,1);

lagrange{1} = (inv_mass_matrix(1,1) * ux + inv_mass_matrix(1,2) * uy ...
               +inv_mass_matrix(1,3) * rho + inv_mass_matrix(1,4) * theta);
                     
lagrange{2} = (inv_mass_matrix(2,1) * ux + inv_mass_matrix(2,2) * uy ...
                         +inv_mass_matrix(2,3) * rho + inv_mass_matrix(2,4) * theta);
                     
                     
lagrange{3} = (inv_mass_matrix(3,1) * ux + inv_mass_matrix(3,2) * uy ...
               +inv_mass_matrix(3,3) * rho + inv_mass_matrix(3,4) * theta);
           
lagrange{4} =  (inv_mass_matrix(4,1) * ux + inv_mass_matrix(4,2) * uy ...
                +inv_mass_matrix(4,3) * rho + inv_mass_matrix(4,4) * theta);


switch id_sys
    case 1
        f = (vx * lagrange{1} + vy * lagrange{2} + lagrange{3}...
                + (He2(vx)+He2(vy))/sqrt(2) * lagrange{4})*value_f0(id);
    case 2
        lagrange_2 = theta /(sqrt(2) * mass_matrix(3,3));
        f = lagrange_2 * value_f0(id);
end


end

function f = He2(x)
f = (x.^2 - 1)/sqrt(2);
end


