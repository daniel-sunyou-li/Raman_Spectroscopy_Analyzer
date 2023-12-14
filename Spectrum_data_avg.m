%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    SPECTRUM DATA AVERAGER v 1.0                    %%%
%%%              Author:  Daniel Li (daniel_li2@brown.edu)             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                 
%%% Take in spectra data and compute the average, variance, and points
%%% used at each bin.  The binning is based upon the smallest binning 
%%% resolution of the raw data and the minimum and maximum value of the
%%% ultimate averaged spectrum is based on the largest minimum value and
%%% the smallest maximum value from the entire data set.  A diagnostic
%%% plot is produced containing the raw data, averaged spectrum, variance 
%%% at each bin and points used for computing the average per bin. Lastly,
%%% a .txt file can be saved for the averaged spectrum.

%%% USER-INPUT FOR DATA FILES TO BE USED
num_files_cell = inputdlg('Data files to use:  ');
num_files = str2double(num_files_cell{1});

while (num_files < 0) || (round(num_files,0) ~= num_files)
    num_files = inputdlg(['Decimals and negative numbers forbidden.',...
        'Type ''0'' to exit. ','Data files to use:  ']);
    if num_files == 0
       return; 
    end
end

%%% USER-INPUT FOR X-UNITS
Unit_choice = menu('Choose the units','cm-1','nm','eV');
Units = [];
if Unit_choice == 1.0
   Units = '(cm-1)';
elseif Unit_choice == 2.0
    Units = '(nm)';
elseif Unit_choice == 3.0
    Units = '(eV)';
end

%%% INITIALIZE STRUCT TO HOLD RAW DATA, AND VARIABLE-HOLDERS
file_list = struct('raw_x',{},'raw_y',{},'path',{});
res_min = 0;
rngmin_min = 0;
rngmax_min = 0;

%%% FOR-LOOP FOR CONVERTING .TXT FILES TO ARRAYS
for i = 1:num_files
    [File,path] = uigetfile('*.txt');
    File_path = strcat(path,File); 
    
    Raw_data = dlmread(File_path,'\t');
    
    % APPENDING DATA FROM .TXT FILE TO STRUCT 'file_list'
    file_list(i).raw_x = Raw_data(:,1);
    file_list(i).raw_y = Raw_data(:,2);
    file_list(i).path = File;
    
    % DETERMINE THE MINIMUM RESOLUTION AND RANGE PARAMETERS
    res_cur = file_list(i).raw_x(2) - file_list(i).raw_x(1);
    rngmin_cur = min(file_list(i).raw_x);
    rngmax_cur = max(file_list(i).raw_x);
    
    if res_cur < res_min || res_min == 0
        res_min = res_cur;
    end
    
    if rngmin_cur < rngmin_min || rngmin_min == 0
       rngmin_min = rngmin_cur; 
    end
    
    if rngmax_cur < rngmax_min || rngmax_min == 0
       rngmax_min = rngmax_cur;
    end
end

%%% PLOT RAW DATA
fig1 = figure;
hold on;
for i = 1:num_files
    plot(file_list(i).raw_x,file_list(i).raw_y);
end
title('Raw Data Used');
xlabel(['Frequency ',Units]);
ylabel('Counts');
legend(file_list.path,'Location','NorthWest');
grid on;
grid minor;
hold off;

%%% USER-INPUT TO PROCEED WITH SELECTED RAW DATA
cont = menu('Proceed with selected data?','Yes','No');
if cont == 2.0
   return; 
end

%%% CREATE NEW ARRAYS BASED ON MINIMUM RESOLUTION AND RANGE PARAMETERS
comb_bins = ceil((rngmax_min - rngmin_min) / res_min);
comb_ls = linspace(rngmin_min,rngmax_min,comb_bins);

%%% HOLDER ARRAYS
comb_mean = zeros(0,comb_bins); 
comb_var = zeros(0,comb_bins);
comb_used = zeros(0,comb_bins);

%%% FOR-LOOP THAT COMPUTES THE MEAN FOR EACH BIN IN AVERAGED DATA
for i = 1:comb_bins
    data_used = 0;
    data_used_arr = []; % TRACKS POINTS USED PER BIN
    
%%% NESTED FOR-LOOP SCANNING FOR EACH ENTRY OF file_list
    for j = 1:num_files
        length_file = length(file_list(j).raw_x);
        
%%% NESTED FOR-LOOP SCANNING THROUGH EACH ENTRY OF file_list(j)
        for k = 1:length_file
            kval = file_list(j).raw_x(k);
            
%%% CONDITIONAL STATEMENTS FOR DETERMINING IF ENTRY FALLS WITHIN BIN RANGE
%%% AND TRACKS NUMBER OF POINTS FALLING WITHIN BIN RANGE      
            if i ~= comb_bins
                if kval >= comb_ls(i) && kval < comb_ls(i+1)
                    data_used = data_used + 1;
                    data_used_arr(data_used) = file_list(j).raw_y(k); 
                end
            elseif i == comb_bins
                if kval >= comb_ls(i) && kval < (rngmax_min + res_min)
                    data_used = data_used + 1;
                    data_used_arr(data_used) = file_list(j).raw_y(k);
                end
            end
        end
    end
    
    sum_count = sum(data_used_arr);
    
%%% EVALUATE MEAN, VARIANCE AND ASSIGN POINTS USED
    comb_mean(i) = mean(data_used_arr);
    comb_var(i) = var(data_used_arr);
    comb_used(i) = data_used;
end

comb_norm_std = (comb_var.^0.5)./comb_mean;
mean_comb_norm_std = mean(comb_norm_std,'omitnan');
mcnSE_round = num2str(round(mean_comb_norm_std,3));

%%% Plot 1: Raw Data Compiled
%%% Plot 2: Averaged Data
%%% Plot 3: Normalized Standard Error
%%% Plot 4: Data Points Used per Bin

figure('pos',[800 0 1200 1000]);

subplot(2,2,1);
hold on;
for i = 1:num_files
    plot(file_list(i).raw_x,file_list(i).raw_y);
end
title(['Raw Data for ',num2str(num_files),' Data Set(s)']);
xlabel(['Frequency ',Units]);
ylabel('Counts (AU)');
grid on;
grid minor;
legend(file_list.path,'Location','NorthWest');
hold off;

subplot(2,2,2);
hold on;
mean_plot = plot(comb_ls,comb_mean,'LineWidth',2,'Color',[0.2 0.6 0.9]);
title(['Average Data for ',num2str(num_files),' Data Set(s)']);
xlabel(['Frequency ',Units]);
ylabel('Avg. Counts (AU)');
grid on;
grid minor;
annotation('textbox',[0.58 0.805 .1 .1],'String',['Range: ',...
    num2str(rngmin_min),' - ',num2str(rngmax_min),' ',Units]);
hold off;

subplot(2,2,3);
hold on;
std_plot = plot(comb_ls,comb_norm_std,'LineStyle','none','Marker','o',...
    'MarkerSize',2, 'Color',[0.1 0.6 0.33]);
title(['Normalized Standard Deviation for ',num2str(num_files),...
    ' Data Set(s)']);
xlabel(['Frequency ',Units]);
ylabel('Norm. SE (AU)');
grid on;
grid minor;
annotation('textbox',[0.14 0.34 .1 .1],'String',['Avg. Norma. SE:  ',...
    mcnSE_round],'FitBoxToText','on');
hold off;

subplot(2,2,4);
use_plot = plot(comb_ls,comb_used, 'LineStyle','none','Marker','o',...
    'MarkerSize',2,'Color',[1 .4 .4]);
title('Data Points Used to Determine Average Count');
xlabel(['Frequency ',Units]);
ylabel('Points Used');
annotation('textbox',[0.58 0.34 .1 .1],'String',['Min Res: ',...
    num2str(res_min),' ',Units]);
ylim([0 num_files*2.5]);
grid on;
hold off;

%%% SAVE DATA
data = flip(rot90(cat(1,comb_ls,comb_mean,comb_var,comb_used),-3));

save = menu('Save data as a .txt?','Yes','No');
if save == 2.0
   return; 
elseif save == 1.0
    folder_path = uigetdir;
    file_name = inputdlg('Write file name (exclude extension): ');
    while isempty(file_name{1})
        file_name = inputdlg('Please write a file name: ');
    end
    file = strcat(folder_path,'/',file_name,'.txt');
    dlmwrite(file{1},data,'delimiter','\t');
end

return;






