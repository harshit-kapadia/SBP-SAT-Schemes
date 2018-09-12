function [vect,value,value_vt,pos_U] = compute_rhoW_prep_simple(Ax,Ay,all_weights)

vect = cell(1,4);
value = cell(1,4);
value_vt = cell(1,4);
pos_U = cell(1,4);

n = [1,0;0,1;-1,0;0,-1]; t = [0,1;-1,0;0,-1;1,0];% normals and tangential vectors for the walls
num_bc = 4; % total number of boundaries

for i = 1 : num_bc
    An = Ax * n(i,1) + Ay * n(i,2);
    At = Ax * t(i,1) + Ay * t(i,2);
    [vect{i},value{i},value_vt{i},pos_U{i}] = compute_rhoW_at_wall(An,At,all_weights);
end

end

% function for a particular wall. 
% An is Ax * nx + Ay * ny and At = Ax * tx + Ay * ty. n=normal to the wall
% and t = tangential direction to the wall. 
function [vect,value,value_vt,pos_U] = compute_rhoW_at_wall(An,At,all_weights)

vn = diag(An);
vt = diag(At);

neg_vn = vn(vn<0); 
pos_vn = vn(vn>0);
vt_neg_vn = vt(vn<0);

%weights multiplied by value of f0 at the negative vns
value_f0 = arrayfun(@(x,y) f0(x,y),neg_vn,vt_neg_vn);
w_times_f0 = value_f0.*all_weights(vn<0);
denominator = dot(neg_vn,w_times_f0);

He2_times_w_times_f0 = (He_2(neg_vn) + He_2(vt_neg_vn)).*w_times_f0/sqrt(2);% He_2(xi_1) + He_2(xi_2) * f0 * weights
value = dot(neg_vn,He2_times_w_times_f0)/denominator;% integral of xi_n  * (He_2(xi_1) + He_2(xi_2) * f0)

He1_times_w_times_f0 = vt_neg_vn.*w_times_f0;
value_vt = dot(neg_vn,He1_times_w_times_f0)/denominator; % integral of xi_n xi_t f0

vect = (pos_vn.*all_weights(vn>0))/denominator;% the vector which hits the kinetic solution
pos_U = find(vn > 0);% position of the positive normal velocities

end

% second degree hermite polynomial
function value = He_2(xi)
value = (xi.^2-1)/sqrt(2);
end

function f = f0(vx,vy)
f = exp(-(vx^2+vy^2)/2)/(2 * pi);
end