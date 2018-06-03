function f = store_minimized_entropy(Ax,Ay, ...
                            mass_matrix,inv_mass_matrix,rho,ux,uy,theta,value_f0)

    neqn = length(diag(Ax));
    f = cell(2,neqn);
   
    
    for id_sys = 1 : 2
        for m = 1 : neqn
            f{id_sys,m} = minimize_entropy(Ax,Ay,...
                mass_matrix,inv_mass_matrix,rho, ...
                ux,uy,theta,id_sys,value_f0,m);
        end
        
    end
end