%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                Multi-Spectra Scaling and Ratio v 1.1               %%%
%%%                      Last Updated: 7/12/2018                       %%%
%%%              Author:  Daniel Li (daniel_li2@brown.edu)             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plots multiple spectra, allows scaling of recently added plot 
%%% then computes the ratio of the secondary, and on, plots to the 
%%% initial plot. 

done = false;

%%% Two final products
f1 = figure('Name','Raw Data','NumberTitle','off','Visible','off');
f2 = figure('Name','Data Ratio','NumberTitle','off','Visible','off');

%%% Ask for user input on units
Units_x = inputdlg('What are the x-units?');
Units_y = inputdlg('What are the y-units?');

%%% First go through data and plot raw data, with the option to scale the
%%% data

uiwait(msgbox('The first data file used will be the reference.'));

data_sets = 0;

data = struct('raw_x',{},'raw_y',{},'scale',{},'path',{},'name',{});

while done == false
    option = menu('Choose new data spectrum?','Yes','No','Exit');
    
    if option == 1.0
        data_sets = data_sets + 1;
        
% Select data and convert to type double and append to list
        
        new = true; % if new is true, will change current data       
        new_scale = false; % if new_scale is true, choose new scale
        new_file = true; % if new_file is true, choose new file
                
        while new == true
            if new_file == true
                [File,path] = uigetfile('*.txt');
                File_path = strcat(path,File); 

                Raw_data = dlmread(File_path,'\t');

                data(data_sets).raw_x = Raw_data(:,1);
                data(data_sets).raw_y = Raw_data(:,2);
                data(data_sets).path = File;
                data(data_sets).scale = [1.0];
                name = inputdlg('Name this data set','s');
                name_scale = strcat(name,' x',...
                    num2str(data(data_sets).scale));
                data(data_sets).name = name_scale;
                clf(f1);
            end
            
            if new_scale == true
                data(data_sets).scale = str2double(inputdlg...
                    ('Choose the scale: '));
                name_scale = strcat(name,' x',...
                    num2str(data(data_sets).scale));
                data(data_sets).name = name_scale;
                clf(f1);
            end
            % Plot each data set to view shape and decide to scale or select
            % new data
            
            figure(f1);
            hold on;
            legend_array = strings(data_sets);

            for n = 1:data_sets
                plot(data(n).raw_x,data(n).raw_y * data(n).scale);
                legend_array(n) = char(data(n).name);
            end
            
            title('Raw Data');
            xlabel(Units_x);
            ylabel(Units_y);
            legend(legend_array,'Location','NorthWest');
            grid on;
            grid minor;
            hold off;
            f1.Visible = 'on';
            
            % Option to proceed
            option = menu('Choose an option: ','Continue',...
                'Choose new file','Scale','Exit');
            
            if option == 1.0 % Continue program without changing current
                new = false;
                
            elseif option == 2.0 % Replace current data
                new = true;
                new_file = true;
                new_scale = false;
                
            elseif option == 3.0 % Scale current data
                new = true;
                new_file = false;
                new_scale = true;
                
            elseif option == 4.0
                return;
            end
        end
        
    elseif option == 2.0 % If no new spectrum is chosen, continue
        if data_sets == 0 % ensures there is a data set available
           return; 
        end
        
        f1.Visible = 'on'; 
        done = true;
        
    elseif option == 3.0 % If exiting program is chosen
        return;
    end
end

%%% Produce the plot for the data ratio

%%% # Data sets - 1 ratios computed
%%% 1.) Determine minimum range
%%% 2.) Allocate appropriate points according to min resolution
%%% 3.) Calculate the mean of the newly binned data
%%% 4.) Divide the newly binned data by the reference data
%%% 5.) Produce Individual plots for each ratio and one that contains all

%%% First declare the n-1 struct arrays to hold relevant information
comp_data = struct('x',{},'ratio',{},'min',{},'index',{},...
    'max',{},'bins',{},'name',{});

%%% Determine the minimum and maximum values for each ratio computed
%%% Using the reference for computing the ratio
res = data(1).raw_x(2) - data(1).raw_x(1);

%%% Loops through each non-reference data set
for n = 1:(data_sets - 1)
    comp_data(n).min = max(min(data(1).raw_x),min(data(1+n).raw_x));
    comp_data(n).max = min(max(data(1).raw_x),max(data(1+n).raw_x));
    comp_data(n).bins = ceil((comp_data(n).max - comp_data(n).min) / res);  
    comp_data(n).x = zeros(0,comp_data(n).bins);
    comp_data(n).ratio = zeros(0,comp_data(n).bins);
    % For loop to determine the minimum index of the reference data  
    for j = 2:(length(data(1).raw_x) - 1)
        prev = data(1).raw_x(j-1);
        curr = data(1).raw_x(j);
        foll = data(1).raw_x(j+1);
        
        dprev = abs(comp_data(n).min - prev);
        dcurr = abs(comp_data(n).min - curr);
        dfoll = abs(comp_data(n).min - foll);
        
        if (dcurr < dprev) && (dcurr < dfoll)
           comp_data(n).index = j;
           break
        elseif dprev < dcurr
            comp_data(n).index = j - 1;
            break
        end
    end
    
    % Re-bin each raw data set (not the reference)
    % Loops through the number of bins for each ratio
    for m = 1:comp_data(n).bins
       data_used = 0;
       sum = 0;
       
       % Bins through each raw data set (except reference) to determine
       % if assigned to current bin of the ratio data
       for k = 1:length(data(n+1).raw_x)
           min_x = comp_data(n).min;
           current = data(n+1).raw_x(k);
           
           % terminates loop prematurely if observed raw data exceeds the
           % maximum value of x for the ratio data
           if data(n+1).raw_x(k) > comp_data(n).max
              break 
           end
           
           if current >= min_x + ((m - 1)*res) && current < min_x + m * res
                sum = sum + data(n+1).raw_y(k);
                data_used = data_used + 1;
           end
       end
       comp_data(n).x(m) = comp_data(n).min + ((m - 1)*res);
       bin_value = sum / data_used; % Mean value of raw data bins
       index = comp_data(n).index - 1 + m;
       
       % Calculating the ratio for a single index
       comp_data(n).ratio(m) = bin_value / data(1).raw_y(index);
    end
    comp_data(n).name = strcat(data(n+1).name);
end

%%% Plot figures collectively if more than 1 data set
if data_sets > 1
    legend_array2 = strings(data_sets-1);
    f3 = figure('Name','All Ratios','NumberTitle','off','Visible','off');
    hold on;
    for n = 1:(data_sets-1)
        plot(comp_data(n).x,comp_data(n).ratio);
        legend_array2(n) = char(comp_data(n).name);
    end
    title(['Ratio of Spectra to Reference: ',data(1).name]);
    xlabel(Units_x);
    ylabel(Units_y);
    grid on;
    grid minor;
    leg1 = legend(legend_array2,'Location','NorthWest');
    title(leg1,'Ratios');
    hold off;
    f3.Visible = 'on';
end

%%% Ask user if they would like to plot any of the plots individually
done = false; % when true, the loop will terminate

if data_sets > 1
    while done == false
        option = menu('Would you like to plot any ratios individually?',...
            'Yes','No');
        if option == 2.0
           done = true; % terminates program
           return; % redundant statement
        elseif option == 1.0
            % Ask which dataset to plot
            option = menu('Choose a ratio dataset to plot',comp_data.name);
            ratnum = round(option);
            f4 = figure('Name','Individual Plot','NumberTitle','off',...
                'Visible','on');
            plot(comp_data(ratnum).x,comp_data(ratnum).ratio);
            title([comp_data(ratnum).name,' vs. ',data(1).name]);
            xlabel(Units_x);
            ylabel(Units_y);
            grid on;
            grid minor;
            
            % Ask if they want to save data as a .txt file
            option = menu('Save as .txt file?','Yes','No');
            if option == 1.0
               data = vertcat(comp_data(ratnum).x,comp_data(ratnum).ratio);
               folder_path = uigetdir;
               file_name = inputdlg(['Write file name (exclude ',...
                   'extension): '],'s');
               while isempty(file_name)
                   file_name = inputdlg('Write file name: ','s');
               end
               file = strcat(folder_path,'/',file_name,'.txt');
               dlmwrite(file{1},data,'delimiter','\t');
               clf(f4);
               f4.Visible = 'off';
            elseif option == 2.0
               clf(f4);
               f4.Visible = 'off';
            end
        end
    end
end

