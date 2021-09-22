
% clear all




root_path = "E:\SamScott\Runs\IEA_15MW\";

outputs_to_open = {
%                    '62'; % 'Blade 1 Loads: User axes'
%                    '63'; % 'Blade 2 Loads: User axes'
%                    '64'; % 'Blade 3 Loads: User axes'
%                    '22'; % 'Hub loads: rotating GL coordinates'
%                    '18'; % 'Blade 1 Deflections'
%                    '19'; % 'Blade 2 Deflections'
%                    '20'; % 'Blade 3 Deflections'
                   '25'; % 'Tower loads GL coordinates'
                   };
               
safety_factors = {
                    'dlc11',    1.25;
                    'dlc12',    1.0;
                    'dlc13',    1.35;
                    'dlc14',    1.35;
                    'dlc15',    1.35;
                    'dlc61',    1.35;
                    'dlc61_wSeaState',    1.35;
                    'dlc62',    1.1;
                    'dlc63',    1.35;
                    'dlc64',    1.0;
                    };

Weibull_params = [9.767, 2.12];

if 1
% Assemble available dlcs in Bladed runs folder
root_folders = dir(root_path);
root_folders(~contains({root_folders.name},'dlc')) = [];

out = [];

%% loop through dlc folders
for i_dlc = [5:6] %:length(root_folders)
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
tow_load_table_all = [];
blade_load_flag = false;
tow_load_flag = false;
dlc_ids = [2];%1:length(dlcs)
max_tow_meta = {};
max_tow = [];
tow_fatigue_loads = {};
for i_dlc = dlc_ids
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
            blade_load_flag = true;
            clock_table_b = [];
            for ib = 1:sum(ID_user_axes)
                [clock_table_allb(:,:,i_dlc,i_run,ib), ULT_loads_b(:,:,:,ib)] = PostProcessLoads(out.(dlc_i).(run_i).(outputs_i{ib}), 0, safety_factor);
                names_allb{1, i_dlc, i_run, ib} = [dlc_i, '_', run_i, '_b', num2str(ib)];
            end
            [clock_table_allr(:,:,i_dlc,i_run), ULT_loads_r(:,:,:,i_run)] = shrink_loads(squeeze(clock_table_allb(:,:,i_dlc,i_run,:)), ULT_loads_b); % Plot_clocktable(clock_table_r(:,:,i_run), 1)
            names_allr(i_run,1) = {[dlc_i, '_', run_i]};
            
            %% Post-process fatigue loads
            if contains(dlc_i, {'dlc11', 'dlc64'})
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
        
        %% Post-process tower loads
        ID_tow_loads = contains(outputs_i,'TowerloadsGLcoordinates');
        if any(ID_tow_loads)
            tow_load_flag = true;
            for iCS = 1:length(out.(dlc_i).(run_i).(outputs_i{1}).span_loc)
                [maxtow_irun, tID] = max(out.(dlc_i).(run_i).(outputs_i{1}).dat(:,iCS,:),[],3);
                maxtow_all_components = squeeze(out.(dlc_i).(run_i).(outputs_i{1}).dat(:,iCS,tID));
                [mintow_irun, tID] = min(out.(dlc_i).(run_i).(outputs_i{1}).dat(:,iCS,:),[],3);
                mintow_all_components = squeeze(out.(dlc_i).(run_i).(outputs_i{1}).dat(:,iCS,tID));
                
                max_min_tow_loads(iCS,:,i_dlc,i_run) = [maxtow_irun', mintow_irun']*safety_factor;
                tow_load_table_all(:,:,iCS,i_dlc,i_run) = reshape([maxtow_all_components; mintow_all_components],8,[])'*safety_factor;
                
                % Save fatigue loads
                if contains(dlc_i, {'dlc11', 'dlc64'})
                    name_split = split(run_i,'_');
                    mean_wind_speed = str2num(erase(name_split{1},'U'));
                    dlc_temp = str2num(erase(dlc_i,'dlc'));
                    
                    out.(dlc_i).(run_i).(outputs_i{1}).dat(:,iCS,:);
                    for j = [1,2,4,5,6,8]
                        Loadsj = squeeze(out.(dlc_i).(run_i).(outputs_i{1}).dat(j,iCS,:)); % Loads(:,j);
                        % Find peaks and valleys
                        PV = [true;(Loadsj(2:end-1,1)>Loadsj(1:end-2,1) & Loadsj(2:end-1,1) > Loadsj(3:end,1)) | ...
                            (Loadsj(2:end-1,1)<Loadsj(1:end-2,1) & Loadsj(2:end-1,1) < Loadsj(3:end,1));true];
                        Extrema = Loadsj(PV,1);
                        [Cmax,Imax] = max(Extrema);
                        if mod(length(Extrema),2)
                            if xor(Extrema(1) > Extrema(2), Extrema(1) >= Extrema(end))
                                Extrema(1) = [];
                                Extrema = circshift(Extrema, 2 - Imax);
                            else
                                Extrema(end) = [];
                                Extrema = circshift(Extrema, 1 - Imax);
                            end
                        else
                            if xor(Extrema(1) > Extrema(2), Extrema(1) >= Extrema(end))
                                Extrema(1) = [];
                                Extrema(end) = [];
                                Extrema = circshift(Extrema, 2 - Imax);
                            else
                                Extrema = circshift(Extrema, 1 - Imax);
                            end
                        end
                        Extrema(end + 1) = Cmax; 
                        CycleData = rainflow_mex(Extrema)';
                        CycleData = [CycleData, ones(size(CycleData,1),1)*mean_wind_speed, ones(size(CycleData,1),1)*dlc_temp];
                        CycleData(:,3) = CycleData(:,3)*i_run;
                        
                        if size(tow_fatigue_loads,1)>=iCS && size(tow_fatigue_loads,2)>=j
                            tow_fatigue_loads{iCS,j} = [tow_fatigue_loads{iCS,j}; CycleData];
                        else
                            tow_fatigue_loads{iCS,j} = CycleData;
                        end
                    end

                    
                end
                
            end
        end
    end % end loop through runs
    
    %% Shrink loads for the present dlc
    % i.e. take the max from each run
    if blade_load_flag
        [clock_table_d(:,:,i_dlc), ULT_loads_d(:,:,:,i_dlc)] = shrink_loads(squeeze(clock_table_allr(:,:,i_dlc,1:length(runs))), ULT_loads_r); % Plot_clocktable(clock_table_r, 1)
    %     Plot_clocktable(squeeze(clock_table_allr(:,:,i_dlc,:)), 1, names_allr)
    end
    if tow_load_flag
        for iCS = 1:size(max_min_tow_loads,1)
            for i_c = 1:size(max_min_tow_loads,2)/2
                [max_tow(iCS, i_c, i_dlc), idtemp] = max(squeeze(max_min_tow_loads(iCS, i_c, i_dlc, 1:length(runs))));
                tow_load_table_dlc(2*i_c-1,:,iCS,i_dlc) = tow_load_table_all(2*i_c-1,:,iCS,i_dlc,idtemp);
                var = erase(out.(dlc_i).(run_i).(outputs_i{1}).variables{i_c},'Tower ');
                max_tow_meta(2*i_c-1,:,iCS,i_dlc)  = [{var}, {'Max'}, {[dlcs{i_dlc},'_', runs{idtemp}]}, {safety_factor}];
                
                [max_tow(iCS, i_c+8, i_dlc), idtemp] = min(squeeze(max_min_tow_loads(iCS, i_c+8 ,i_dlc, 1:length(runs))));
                tow_load_table_dlc(2*i_c,:,iCS,i_dlc) = tow_load_table_all(2*i_c,:,iCS,i_dlc,idtemp);
                max_tow_meta(2*i_c,:,iCS,i_dlc)  = [{var}, {'Min'}, {[dlcs{i_dlc},'_', runs{idtemp}]}, {safety_factor}];
            end
        end
    end
    
    names_all_dlc{i_dlc,1} = names_allr;
end

%% Shrink loads from all dlcs
if blade_load_flag
    [clock_table, ULT_loads] = shrink_loads(clock_table_d(:,:,1:5), ULT_loads_d(:,:,:,1:5));
    Plot_clocktable(clock_table_d(:,:,1:5), 1, dlcs)            
    plot_ULTloads(span_loc, ULT_loads)
end
if tow_load_flag
    %% Process tower ultimate loads
    max_tow_meta_final = {};
    for iCS = 1:size(max_tow,1)
        for i_c = 1:size(max_tow,2)/2
            [max_tow_final(iCS, 2*i_c-1), idtemp] = max(squeeze(max_tow(iCS, i_c, dlc_ids)));
            tow_load_table_final(2*i_c-1,:,iCS) = tow_load_table_dlc(2*i_c-1, :, iCS, dlc_ids(idtemp));
            max_tow_meta_final(2*i_c-1,:,iCS)  = max_tow_meta(2*i_c-1, :, iCS, dlc_ids(idtemp));
            
            [max_tow_final(iCS, 2*i_c), idtemp] = min(squeeze(max_tow(iCS, i_c+8, dlc_ids)));
            tow_load_table_final(2*i_c,:,iCS) = tow_load_table_dlc(2*i_c, :, iCS, dlc_ids(idtemp));
            max_tow_meta_final(2*i_c,:,iCS)  = max_tow_meta(2*i_c, :, iCS, dlc_ids(idtemp));
        end
    end
    tower_load_table = [max_tow_meta_final(:,1:3,:),num2cell(tow_load_table_final),max_tow_meta_final(:,4,:)];
    tower_load_table = [repmat([{'Height:', ''}, {'Load Case'}, erase(out.(dlc_i).(run_i).(outputs_i{1}).variables,'Tower '), {'Saftey Factor'}], 1, 1, size(max_tow,1)); tower_load_table];
    tower_load_table(1,2,:) = num2cell(out.(dlc_i).(run_i).(outputs_i{1}).span_loc);
    plot_ULT_tower_loads(out.(dlc_i).(run_i).(outputs_i{1}).span_loc, max_tow_final)
    if 0
        writecell(reshape(permute(tower_load_table, [2,1,3]), 12, [])', 'Tower_ultimate_loads_dlc61_wSeaState.csv')
    end
    
    %% Process tower fatigue loads    
    wind_bin_edges = 3:2:35;
    lifetime = 25;
    tower_DELs = [];
    NFat = 10^7;
    for iCS = 1:size(tow_fatigue_loads,1)
        for iL = [1,2,4,5,6,8]
            wind_speeds_all = tow_fatigue_loads{iCS,iL}(:,[4:5,3]);
            wind_speeds = unique(wind_speeds_all(:,1));
            % probability of each bin occurring
            pdf  = WeibullPDF(wind_speeds, Weibull_params(1), Weibull_params(2), 'cdf', 0);
            [Nreps, V_uniq] = hist(wind_speeds_all(:,1), wind_speeds);
            Nrepts = [V_uniq, Nreps'];
            % loop through bins
            pdf_all = zeros(length(wind_speeds_all(:,1)),1);
            for iv = 1:length(wind_bin_edges)-1
                ID_bin_i = find(wind_speeds_all(:,1)<wind_bin_edges(iv+1) & wind_speeds_all(:,1)>wind_bin_edges(iv));
                dlc = wind_speeds_all(ID_bin_i, :);
                n11 = 0; n64 = 0;
                dlc11_ID = dlc(:,2)==11;
                dlc11 = unique(dlc(dlc11_ID,:),'rows');
                if any(dlc11_ID)
                    n11 = size(dlc11, 1);
                    pdf_all(ID_bin_i(dlc11_ID)) = pdf(iv)*0.975/n11;
                end
                dlc64_ID = dlc(:,2)==64;
                dlc64 = unique(dlc(dlc64_ID,:),'rows');
                if any(dlc64_ID)
                    n64 = size(dlc64, 1);
                    pdf_all(ID_bin_i(dlc64_ID)) = pdf(iv)*0.025/n64;
                end

            end
            yrs     = pdf_all*lifetime;
            repeats = yrs*365.25*24*60*60./600;
%             ChkYrs = sum(repeats)/6/8766

            for iSN = 3:14
                N = (tow_fatigue_loads{iCS,iL}(:,1)./5e8).^-iSN;
                D = sum(repeats./N);
                S = 5e8*(NFat/D)^(-1/iSN);
                tower_DELs(iSN,iL,iCS) = S;
            end
        end
    end
    tower_DELs(:,[3,7],:) = [];
    tower_DELs(1:2,:,:) = [];
    
    tower_DELs_2print = [num2cell(repmat((3:14)',1,1,size(tower_DELs,3))), num2cell(tower_DELs)];
    tower_DELs_2print = [repmat([cell(1,7); {'SN slope'}, erase(out.(dlc_i).(run_i).(outputs_i{1}).variables([1,2,4,5,6,8]),'Tower ')],1,1,size(tower_DELs_2print,3)); ...
                    tower_DELs_2print];
    tower_DELs_2print(1,1:2,:) = [repmat({'Span = '},1,1,size(tower_DELs,3)), permute(num2cell(out.(dlc_i).(run_i).(outputs_i{1}).span_loc),[1,3,2])];
    writecell(reshape(permute(tower_DELs_2print, [2,1,3]), 7, [])', 'Tower_DELs.csv')
end



fr

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

