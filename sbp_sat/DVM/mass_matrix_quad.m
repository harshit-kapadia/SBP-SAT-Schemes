function mass_matrix = mass_matrix_quad(Ax,Ay,all_weights)

mass_matrix = zeros(4,4);

vx = diag(Ax);
vy = diag(Ay);

value_f0 = arrayfun(@(x,y) f0(x,y),vx,vy);

scaled_weights = all_weights.*value_f0;
mass_matrix(1,1) = sum((He1(vx).*He1(vx)).*scaled_weights);
mass_matrix(1,2) = sum((He1(vx).*He1(vy)).*scaled_weights);
mass_matrix(1,3) = sum((He1(vx).* He0(vx)).*scaled_weights);
mass_matrix(1,4) = sum((He1(vx).*(He2(vx)+He2(vy))).*scaled_weights)/sqrt(2);

mass_matrix(2,2) = sum((He1(vy).*He1(vy)).*scaled_weights);
mass_matrix(2,3) = sum((He1(vy).*He0(vx)).*scaled_weights);
mass_matrix(2,4) = sum((He1(vy).*(He2(vx)+He2(vy))).*scaled_weights)/sqrt(2);

mass_matrix(3,3) = sum((He0(vx).*He0(vx)).*scaled_weights);
mass_matrix(3,4) = sum((He0(vx).*(He2(vx)+He2(vy))).*scaled_weights)/sqrt(2);

mass_matrix(4,4) = sum(((He2(vx)+He2(vy)).*(He2(vx)+He2(vy))).*scaled_weights)/2;

% make the whole thing symmetric
mass_matrix = mass_matrix+mass_matrix'-diag(diag(mass_matrix));
end

function f = He1(x)
f = x;
end

function f = He0(x)
f = ones(length(x),1);
end

function f = He2(x)
f = (x.^2 - 1)/sqrt(2);
end

function f = f0(vx,vy)
f = exp(-(vx^2+vy^2)/2)/(2 * pi);
end