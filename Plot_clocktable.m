function Plot_clocktable(clockface_table, iCS)

    theta = linspace(0, 330, size(clockface_table,2));
    unit_vec = [0;1;0]; % start is flapwise toward tower
    for i = 1:length(theta)
        Rz =  [cosd(theta(i))    -sind(theta(i))  0;...
               sind(theta(i))    cosd(theta(i))   0;...
               0                0               1];
        unit_vec_i(:,i) = Rz*unit_vec;
    end
    figure, hold on
    for i = 1:size(clockface_table,3)
        clock_loads = clockface_table(iCS,:,i).*unit_vec_i(1:2,:);
        clock_loads = [clock_loads,clock_loads(:,1)];
        plot(clock_loads(2,:), clock_loads(1,:),'-o')
    end
    axis equal, grid on
    xlabel('Flapwise load'), ylabel('Edgewise ;load'), title(['Cross section ID: ', num2str(iCS)])
end