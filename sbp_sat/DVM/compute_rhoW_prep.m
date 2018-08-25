function [vect,value,pos_U] = compute_rhoW_prep(Ax,Ay,all_weights)

all_vx = diag(Ax);
all_vy = diag(Ay);

pos_vx_p = find(all_vx > 0); % position of positive x velocities
pos_vx_m = find(all_vx < 0); % position of negative x velocites

% same as above but for the y direction. 
pos_vy_p = find(all_vy > 0);
pos_vy_m = find(all_vy < 0);

vx_p_x = all_vx(pos_vx_p);
vx_m_x = all_vx(pos_vx_m);
vy_p_x = all_vy(pos_vx_p); % y-vel corresponding to positive x-vel
vy_m_x = all_vy(pos_vx_m); % y-vel corresponding to negative x-vel

% same as above but for the y direction
vx_p_y = all_vx(pos_vy_p);
vx_m_y = all_vx(pos_vy_m);
vy_p_y = all_vy(pos_vy_p);
vy_m_y = all_vy(pos_vy_m);

num_pos_x = length(vx_p_x);
num_neg_x = length(vx_m_x);

num_pos_y = length(vy_p_y);
num_neg_y = length(vy_m_y);

weight_vx_m = all_weights(pos_vx_m);
weight_vx_p = all_weights(pos_vx_p);

weight_vy_m = all_weights(pos_vy_m);
weight_vy_p = all_weights(pos_vy_p);

vect = cell(1,4);
value = cell(1,4);

int_f0 = 0; % actually will store int_f0_vn (vn is normal velocity at boundary)
int_f0_vx_sq = 0; % integral(f0 * vx_sq * vn)
int_f0_vy_sq = 0; % integral(f0 * vy_sq * vn)
% boundary at x = 1
% integrate the maxwellian with vx_m_x (incoming normal velocity at boundary)
for i = 1 : num_neg_x
    int_f0 = int_f0 + vx_m_x(i) * f0(vx_m_x(i),vy_m_x(i)) * weight_vx_m(i);
    int_f0_vx_sq = int_f0_vx_sq + vx_m_x(i) * vx_m_x(i)^2 * ...
        f0(vx_m_x(i),vy_m_x(i)) * weight_vx_m(i);
    int_f0_vy_sq = int_f0_vy_sq + vx_m_x(i) * vy_m_x(i)^2 * ...
        f0(vx_m_x(i),vy_m_x(i)) * weight_vx_m(i);
end
value{1} = ((int_f0_vx_sq + int_f0_vy_sq)/2-int_f0)/int_f0;
% prepration to
% integrate the kinetic solution with vx_p_x (outgoing normal velocity at boundary)
% int_f_vn = sum((vx_p_x.*U(pos_vx_p)).*weight_vx_p); (actual integration)
vect{1} = (vx_p_x.*weight_vx_p)/int_f0;


int_f0 = 0;
int_f0_vx_sq = 0;
int_f0_vy_sq = 0;
% integrate the maxwellian with vy_m_y (incoming normal velocity at boundary)
for i = 1 : num_neg_y
    int_f0 = int_f0 + vy_m_y(i) * f0(vx_m_y(i),vy_m_y(i)) * weight_vy_m(i);
    int_f0_vx_sq = int_f0_vx_sq + vy_m_y(i) * vx_m_y(i)^2 * ...
        f0(vx_m_y(i),vy_m_y(i)) * weight_vy_m(i);
    int_f0_vy_sq = int_f0_vy_sq + vy_m_y(i) * vy_m_y(i)^2 * ...
        f0(vx_m_y(i),vy_m_y(i)) * weight_vy_m(i);
end

value{2} = ((int_f0_vx_sq + int_f0_vy_sq)/2-int_f0)/int_f0;
% prepration to
% integrate the kinetic solution with vy_p_y (outgoing normal velocity at boundary)
vect{2} = (vy_p_y.*weight_vy_p)/int_f0;


int_f0 = 0;
int_f0_vx_sq = 0;
int_f0_vy_sq = 0;
% integrate the maxwellian with vx_p_x (incoming normal velocity at boundary)
for i = 1 : num_pos_x
    int_f0 = int_f0 + vx_p_x(i) * f0(vx_p_x(i),vy_p_x(i)) * weight_vx_p(i);
    int_f0_vx_sq = int_f0_vx_sq + vx_p_x(i) * vx_p_x(i)^2 * ...
        f0(vx_p_x(i),vy_p_x(i)) * weight_vx_p(i);
    int_f0_vy_sq = int_f0_vy_sq + vx_p_x(i) * vy_p_x(i)^2 * ...
        f0(vx_p_x(i),vy_p_x(i)) * weight_vx_p(i);
end
value{3} = ((int_f0_vx_sq + int_f0_vy_sq)/2-int_f0)/int_f0;
% prepration to
% integrate the kinetic solution with vx_m_x (outgoing normal velocity at boundary)
vect{3} = vx_m_x.*weight_vx_m/int_f0;


int_f0 = 0;
int_f0_vx_sq = 0;
int_f0_vy_sq = 0;
% integrate the maxwellian with vy_p_y (incoming normal velocity at boundary)
for i = 1 : num_pos_y
    int_f0 = int_f0 + vy_p_y(i) * f0(vx_p_y(i),vy_p_y(i)) * weight_vy_p(i);
    int_f0_vx_sq = int_f0_vx_sq + vy_p_y(i) * vx_p_y(i)^2 * ...
        f0(vx_p_y(i),vy_p_y(i)) * weight_vy_p(i);
    int_f0_vy_sq = int_f0_vy_sq + vy_p_y(i) * vy_p_y(i)^2 * ...
        f0(vx_p_y(i),vy_p_y(i)) * weight_vy_p(i);
end
value{4} = ((int_f0_vx_sq + int_f0_vy_sq)/2-int_f0)/int_f0;
% prepration to
% integrate the kinetic solution with vy_m_y (outgoing normal velocity at boundary)
vect{4} = vy_m_y.*weight_vy_m/int_f0;

% % final target
% rhoW = -int_f_vn-thetaW * ((int_f0_vx_sq + int_f0_vy_sq)/2-int_f0);
% rhoW = rhoW/int_f0;

pos_U{1} = pos_vx_p ;
pos_U{2} = pos_vy_p ;
pos_U{3} = pos_vx_m ;
pos_U{4} = pos_vy_m ;
end

function f = f0(vx,vy)
f = exp(-(vx^2+vy^2)/2)/(2 * pi);
end