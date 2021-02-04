function [clock_table_out, ULT_loads_out, names_out] = shrink_loads(clock_table_in, ULT_loads_in)

[clock_table_out, ID_b] = max(clock_table_in,[],3);


for ith = 1:size(ULT_loads_in,3)
    for iCS = 1:size(ULT_loads_in,1)
        ULT_loads_out(iCS,:,ith) = ULT_loads_in(iCS,:,ith,ID_b(iCS,ith));
%         names_out(1,
    end
end

end