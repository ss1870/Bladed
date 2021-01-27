
clear all




root_path = "E:\SamScott\Runs\IEA_15MW\";

outputs_to_open = {'62'; % 'Blade 1 Loads: User axes'
                   '63'; % 'Blade 2 Loads: User axes'
                   '64'; % 'Blade 3 Loads: User axes'
                   '22'; % 'Hub loads: rotating GL coordinates'
                   '18'; % 'Blade 1 Deflections'
                   '19'; % 'Blade 2 Deflections'
                   '20'; % 'Blade 3 Deflections'
                   };


% Assemble available dlcs in Bladed runs folder
root_folders = dir(root_path);
root_folders(~contains({root_folders.name},'dlc')) = [];

out = [];

%% loop through dlc folders
for i_dlc = [1,3,4] %:length(root_folders)
    dlc_path = root_path + root_folders(i_dlc).name;
    run_folders = dir(dlc_path);
    run_folders = run_folders(3:end);
    
    rt_fname = erase(root_folders(i_dlc).name,'.');
    out.(rt_fname) = [];
    %% loop through runs in dlc folder
    for i_run = 1:length(run_folders)
        run_name = run_folders(i_run).name;
        run_path = dlc_path + '\' + run_name;
        out.(rt_fname).(run_name) = [];
        
        %% loop through desired output files
        for i_out = 1:length(outputs_to_open)
            % define filename paths
            filepath_binary = run_path + '\powprod.$' + outputs_to_open(i_out);
            filepath_header = run_path + '\powprod.%' + outputs_to_open(i_out);
            
            %% Extract useful info from header file
            if isfile(filepath_header)
                fileID_h = fopen(filepath_header,'r');
                tline = fgetl(fileID_h); tline_split = strsplit(tline);
                while ischar(tline)
                   if any(strcmp(tline_split, 'NDIMENS')),    n_dimens = str2double(tline_split(2)); end
                   if any(strcmp(tline_split, 'DIMENS')),     dimens = str2double(tline_split(2:end)); end
                   if any(strcmp(tline_split, 'AXIVAL')),     span_loc = str2double(tline_split(2:end)); end
                   if any(strcmp(tline_split, 'GENLAB'))
                       tline_split = strsplit(tline,'''');
                       out_name = tline_split(2); 
                   end
                   if any(strcmp(tline_split, 'MIN')),        min_t = str2double(tline_split(2)); end
                   if any(strcmp(tline_split, 'STEP')),       dt = str2double(tline_split(2)); end
                   if any(strcmp(tline_split, 'VARIAB'))    
                       tline_split = strsplit(tline,'''');
                       variables = tline_split(2:2:end); 
                   end

                   tline = fgetl(fileID_h);
                   if ischar(tline), tline_split = strsplit(tline); end
                end
                fclose(fileID_h);
            else
                error(['Header file does not exist:',filepath_header])
            end
            
            %% Extract data from binary file
            if isfile(filepath_binary)
                fileID_b = fopen(filepath_binary);
                dat = fread(fileID_b,inf,'float32=>double');
                fclose(fileID_b);    
                dat = reshape(dat,dimens);
                t = min_t:dt:(dimens(end)-1)*dt;
                
                out_name = erase(out_name{:},[" ",":"]);
                out.(rt_fname).(run_name).(out_name) = [];
                out.(rt_fname).(run_name).(out_name).dat = dat;
                out.(rt_fname).(run_name).(out_name).t = t;
                out.(rt_fname).(run_name).(out_name).span_loc = span_loc;
                out.(rt_fname).(run_name).(out_name).variables = variables;
            else
                error(['Binary file does not exist:',filepath_binary])
            end
            
            
            %% Post-processing of loads
            if isfield(out.(rt_fname).(run_name).(out_name),'dat') && ...
                    contains(out_name,'LoadsUseraxes')
                
                
            end
            
            %% Post-process displacements
            if isfield(out.(rt_fname).(run_name).(out_name),'dat') ...
                && contains(out_name,'Deflections')
                if isfield(out.(rt_fname).(run_name),'minmaxDefl')
                    out.(rt_fname).(run_name).minmaxDefl(:,:,end+1) = [...
                            min(out.(rt_fname).(run_name).(out_name).dat(:,end,:), [], 3),...
                            max(out.(rt_fname).(run_name).(out_name).dat(:,end,:), [], 3)];
                else
                    out.(rt_fname).(run_name).minmaxDefl = [...
                            min(out.(rt_fname).(run_name).(out_name).dat(:,end,:), [], 3),...
                            max(out.(rt_fname).(run_name).(out_name).dat(:,end,:), [], 3)];
                end
                out.(rt_fname).(run_name) = rmfield(out.(rt_fname).(run_name), out_name);
            end
            
            
        end % end loop through desired outputs
        
    end % end loop through run folders
    
end % end loop through dlc folders


%% Post-processing
% probably best to keep the post-processing separate from the data reading

dlcs = fieldnames(out);
clock_table_d = [];
for i_dlc = 1:length(dlcs)
    dlc_i = dlcs{i_dlc};
    runs = fieldnames(out.(dlc_i));
    clock_table_r = [];
    for i_run = 1:length(runs)
        run_i = runs{i_run};
        outputs_i = fieldnames(out.(dlc_i).(run_i));
        
        %% Post-process blade loads in user axes
        ID_user_axes = contains(outputs_i,'LoadsUseraxes');
        if any(ID_user_axes)
            clock_table_b = [];
            for ib = 1:sum(ID_user_axes)
                [clock_table_b(:,:,ib), ULT_loads_b(:,:,:,ib)] = PostProcessLoads(out.(dlc_i).(run_i).(outputs_i{ib}), 0);
            end
            [clock_table_r(:,:,i_run), ULT_loads_r(:,:,:,i_run)] = shrink_loads(clock_table_b, ULT_loads_b); % Plot_clocktable(clock_table_r(:,:,i_run), 1)
            
        end
    end % end loop through runs
    [clock_table_d(:,:,i_dlc), ULT_loads_d(:,:,:,i_dlc)] = shrink_loads(clock_table_r, ULT_loads_r); % Plot_clocktable(clock_table_r, 1)
    % Plot_clocktable(clock_table_d, 1)
    
end

[clock_table, ULT_loads] = shrink_loads(clock_table_d, ULT_loads_d);



Plot_clocktable(clock_table_d, 1)            
plot_ULTloads(span_loc, ULT_loads)


