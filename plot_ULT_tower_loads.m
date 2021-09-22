function plot_ULT_tower_loads(span_loc, ULT_loads)

figure

components = [{'Mx'}, {'Nm'};
              {'My'}, {'Nm'};
              {'Mxy'}, {'Nm'};
              {'Mz'}, {'Nm'};
              {'Fx'}, {'N'};
              {'Fy'}, {'N'};
              {'Fxy'}, {'N'};
              {'Fz'}, {'N'};];

for i_comp = 1:size(ULT_loads,2)/2
    subplot(2,4,i_comp)
    hold on
    plot(span_loc, squeeze(ULT_loads(:,[i_comp*2-1,i_comp*2])))
    xlabel('Spanwise location (m)')
    ylabel(['Load (',components{i_comp,2},')'])
    title(components{i_comp,1})
    if i_comp ==1
        legend
    end
    grid on
end


end