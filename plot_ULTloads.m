function plot_ULTloads(span_loc, ULT_loads)

figure

components = [{'Mx'}, {'Nm'};
              {'My'}, {'Nm'};
              {'Mxy'}, {'Nm'};
              {'Mz'}, {'Nm'};
              {'Fx'}, {'N'};
              {'Fy'}, {'N'};
              {'Fxy'}, {'N'};
              {'Fz'}, {'N'};];
angle = linspace(0,330,12);

for i_comp = 1:size(ULT_loads,2)
    subplot(2,4,i_comp)
    hold on
    for ith = 1:12
        plot(span_loc, squeeze(ULT_loads(:,i_comp,ith)),'DisplayName', ['Angle = ', num2str(angle(ith))])
    end
    xlabel('Spanwise location (m)')
    ylabel(['Load (',components{i_comp,2},')'])
    title(components{i_comp,1})
    if i_comp ==1
        legend
    end
end


end