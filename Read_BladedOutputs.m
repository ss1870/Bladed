
% clear all




root_path = "E:\SamScott\Runs\IEA_15MW\";

outputs_to_open = {'62'; % 'Blade 1 Loads: User axes'
                   '63'; % 'Blade 2 Loads: User axes'
                   '64'; % 'Blade 3 Loads: User axes'
                   '22'; % 'Hub loads: rotating GL coordinates'
                   '18'; % 'Blade 1 Deflections'
                   '19'; % 'Blade 2 Deflections'
                   '20'; % 'Blade 3 Deflections'
                   };
               
safety_factors = {
                    'dlc11',    1.25;
                    'dlc12',    1.0;
                    'dlc13',    1.35;
                    'dlc14',    1.35;
                    'dlc15',    1.35;
                    'dlc61',    1.35;
                    'dlc62',    1.1;
                    'dlc63',    1.35;
                    'dlc64',    1.0;
                    };

Weibull_params = [9.767, 2.12];
if 0
% Assemble available dlcs in Bladed runs folder
root_folders = dir(root_path);
root_folders(~contains({root_folders.name},'dlc')) = [];

out = [];

%% loop through dlc folders
for i_dlc = [1:8] %:length(root_folders)
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
end

%% Post-processing
% probably best to keep the post-processing separate from the data reading

dlcs = fieldnames(out);
clock_table_d = [];
names_all_dlc = {};
header_data = {};
for i_dlc = [1,8]%1:length(dlcs)
    dlc_i = dlcs{i_dlc};
    runs = fieldnames(out.(dlc_i));
    safety_factor = safety_factors{contains(safety_factors(:,1), dlc_i), 2};
%     clock_table_allr = [];
    names_allr = {};
    for i_run = 1:length(runs)
        run_i = runs{i_run};
        outputs_i = fieldnames(out.(dlc_i).(run_i));
        
        %% Post-process blade loads in user axes
        ID_user_axes = contains(outputs_i,'LoadsUseraxes');
        if any(ID_user_axes)
            clock_table_b = [];
            for ib = 1:sum(ID_user_axes)
                [clock_table_allb(:,:,i_dlc,i_run,ib), ULT_loads_b(:,:,:,ib)] = PostProcessLoads(out.(dlc_i).(run_i).(outputs_i{ib}), 0, safety_factor);
                names_allb{1, i_dlc, i_run, ib} = [dlc_i, '_', run_i, '_b', num2str(ib)];
            end
            [clock_table_allr(:,:,i_dlc,i_run), ULT_loads_r(:,:,:,i_run)] = shrink_loads(squeeze(clock_table_allb(:,:,i_dlc,i_run,:)), ULT_loads_b); % Plot_clocktable(clock_table_r(:,:,i_run), 1)
            names_allr(i_run,1) = {[dlc_i, '_', run_i]};
            
            %% Post-process fatigue loads
            if contains(dlc_i, {'dlc11','dlc64'})
%             if contains(dlc_i, {'dlc64'})
                ib = 1;
                dat = permute(out.(dlc_i).(run_i).(outputs_i{ib}).dat([1,2,4,5,6,8],:,:), [3, 1, 2]);
                dat = reshape(dat, size(dat,1), []);
                dat = [out.(dlc_i).(run_i).(outputs_i{ib}).t', dat];
                dat = [[0, repelem(out.(dlc_i).(run_i).(outputs_i{ib}).span_loc, 1, 6)]; dat];
                if contains(dlc_i, {'dlc64'})
                    writematrix(dat,['X:\Corporate\3. INNOVATION PROGRAMMES\01 Projects\02 Projects\01 Live\PN000450 SIF4Ax\10 Engineering\01 Data\Fatigue Load CAses\', [dlc_i, '_', run_i], '.csv']) 
                end
                name_split = split(run_i,'_');
                mean_wind_speed = str2num(erase(name_split{1},'U'));
                header_data = [header_data; {[dlc_i, '_', run_i], mean_wind_speed}];
            end
        end
    end % end loop through runs
    [clock_table_d(:,:,i_dlc), ULT_loads_d(:,:,:,i_dlc)] = shrink_loads(squeeze(clock_table_allr(:,:,i_dlc,1:length(runs))), ULT_loads_r); % Plot_clocktable(clock_table_r, 1)
%     Plot_clocktable(squeeze(clock_table_allr(:,:,i_dlc,:)), 1, names_allr)
    names_all_dlc{i_dlc,1} = names_allr;
end

[clock_table, ULT_loads] = shrink_loads(clock_table_d(:,:,1:5), ULT_loads_d(:,:,:,1:5));
Plot_clocktable(clock_table_d(:,:,1:5), 1, dlcs)            
plot_ULTloads(span_loc, ULT_loads)




% Array Prep
wind_speeds_all = cell2mat(header_data(:,2));
wind_speeds = unique(cell2mat(header_data(:,2)));
wind_bin_edges = 3:2:35;
% probability of each bin occurring
pdf  = WeibullPDF(wind_speeds, Weibull_params(1), Weibull_params(2), 'cdf', 0);
[Nreps, V_uniq] = hist(wind_speeds_all, wind_speeds);
Nrepts = [V_uniq, Nreps'];

% loop through bins
pdf_all = zeros(length(wind_speeds_all(:,1)),1);
for iv = 1:length(wind_bin_edges)-1
    ID_bin_i = find(wind_speeds_all(:,1)<wind_bin_edges(iv+1) & wind_speeds_all(:,1)>wind_bin_edges(iv));
    dlc = wind_speeds_all(ID_bin_i, 2);
    n11 = 0; n64 = 0;
    dlc11 = dlc==1.1;
    if any(dlc11)
        n11 = sum(dlc11);
        pdf_all(ID_bin_i(dlc11)) = pdf(iv)*0.975/n11;
    end
    dlc64 = dlc==6.4;
    if any(dlc64)
        n64 = sum(dlc64);
        pdf_all(ID_bin_i(dlc64)) = pdf(iv)*0.025/n64;
    end
    
end


for i = 1:length(pdf)
   pdf(i) = pdf(i)/Nrepts(Nrepts(:,1) == wind_speeds_all(i,1), 2);
end
lifetime = 25;
yrs     = pdf_all*lifetime;
repeats = yrs*365.25*24*60*60./600;
ChkYrs = sum(repeats)/6/8766
header_data = [header_data, num2cell(repeats)];
header_data = [{'RunName','WindSpeed','Repeats'}; header_data];
T = cell2table(header_data(2:end,:),'VariableNames',header_data(1,:));
 
% Write the table to a CSV file
writetable(T,'X:\Corporate\3. INNOVATION PROGRAMMES\01 Projects\02 Projects\01 Live\PN000450 SIF4Ax\10 Engineering\01 Data\Fatigue Load CAses\Header.csv')

