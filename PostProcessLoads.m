function [clockface_table, extreme_loads, unit_vec_i] = PostProcessLoads(input, PLOT, SF)

span_loc = input.span_loc;

M = input.dat([1,2],:,:);

theta = linspace(0,330,12);
unit_vec = [0;1;0]; % start is flapwise toward tower
clockface_table = zeros(length(span_loc), length(theta));
extreme_loads = zeros(length(span_loc), length(input.variables), length(theta));
for i_th = 1:length(theta)
    theta_i = theta(i_th);
    Rz =  [cosd(theta_i)    -sind(theta_i)  0;...
           sind(theta_i)    cosd(theta_i)   0;...
           0                0               1];
    unit_vec_i(:,i_th) = Rz*unit_vec;
    
    
    % scalar projection of the load vector (F/M) = dot product with the
    % load vector and the unit vector to be projected on to
    M_proj = sum(unit_vec_i(1:2,i_th).*M, 1);
    
    [max_M, t_ID] = max(M_proj, [], 3);
    for iCS = 1:length(span_loc)
        extreme_loads(iCS,:,i_th) = input.dat(:, iCS, t_ID(iCS))*SF;
    end
    clockface_table(:,i_th) = max_M*SF;     
end

if PLOT
    Plot_clocktable(clockface_table, 1)
end
end