%%% Plots the Polynomial Fitted Center vs. the number of fit points for
%%% different Estimated centers.  Additionally plots the adjusted chi-
%%% square values for the right-hand side axis to qualify the Polynomial
%%% Fitted Center.

uiwait(msgbox(['Files must be .txt with three columns in the order ',...
    'of ''Fit Points'', ''Center Estimate'', ''Adj. Chi-Squared''']));

f1 = figure('Name','Center Estimate from Polynomial Fit',...
    'NumberTitle','off',...
    'Visible','off');

data = struct('Points',[],...
    'Center',[],...
    'Chi',[],...
    'Name',{});

done = false; % Done when user does not want to include new plots

file_num = 0;

while done == false
    option = menu('Choose an option.',...
        'Add New Data',...
        'Plot Existing Data',...
        'Exit Program');
    
    if option == 1.0
        file_num = file_num + 1;
        
        [File,path] = uigetfile('*.txt');
        File_path = strcat(path,File);

        Raw_data = dlmread(File_path,'\t');
        while size(Raw_data,2) ~= 3 % Ensures correct data file 
            uiwait(msgbox('File is invalid--be sure it has 3 columns.'));
            [File,path] = uigetfile('*.txt');
            File_path = strcat(path,File);

            Raw_data = dlmread(File_path,'\t');         
        end
        
        data(file_num).Points = Raw_data(:,1);
        data(file_num).Center = Raw_data(:,2);
        data(file_num).Chi = Raw_data(:,3);
        data(file_num).Name = inputdlg('Name this data set: ');
        
    elseif option == 2.0
        figure(f1);
        hold on
        
        for i = 1:file_num % Plot the center estimates
            cval = file_num*2;
            yyaxis left
            plot(data(i).Points,data(i).Center,'-o',...
                'MarkerFaceColor',[i/cval i/cval i/cval],...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerSize',8);
            ylabel('Center Estimate');
        end
        
        for i = 1:file_num % Plot the adj. chi-squared values
           yyaxis right
           plot(data(i).Points,data(i).Chi,'-s',...
               'MarkerFaceColor',[i/cval i/cval i/cval],...
               'MarkerEdgeColor',[0 0 0],...
               'MarkerSize',8);
           ylabel('\chi^{2}_{Adj}');
                  
        end
        
        xlabel('Fit Points Used');
        grid on;
        grid minor;
        
        f1.Visible = 'on'; % Show graph
        
        title(inputdlg('Name the title of the plot: '));
        
        option = menu('Show legend?',...
        'Yes',...
        'No');
        if option == 1.0
            legend([data.Name]);
        end
        
        done = true; % Ends program
             
    elseif option == 3.0
        return;
    end
end