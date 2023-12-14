%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Optical Spectrum Visualization and Frequency ID v 1.0        %%%
%%%              Author:  Daniel Li (daniel_li2@brown.edu)             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                 
%%% READS IN A SINGLE DATA FILE AND THEN REQUESTS USER TO INPUT LINES THAT
%%% INDICATE INTERESTING FEATURES AT THE GIVEN X-VALUE.  WORKS FOR ANY
%%% OPTICAL/EM SPECTRA (I.E. RAMAN, TRANSMISSION, ABSORPTION, ETC.)

%%% SELECT FILE AND LOAD IT AS DATA
[File,path] = uigetfile('*.txt');
File_path = strcat(path,File); 

Raw_data = dlmread(File_path,'\t');
Raw_data_x = Raw_data(:,1);
Raw_data_y = Raw_data(:,2);

%%% PRODUCES PLOT
figure;
plot(Raw_data_x,Raw_data_y);
title(File);
grid on;
grid minor;
xlabel('Frequency');
ylabel('Counts');
ylim([0 max(Raw_data_y)*1.2])

%%% RECURSIVELY ALLOW USER TO INPUT LINES
pass = false;
num = 0;
while pass == false
    option = menu('Plot a new line?','Yes','No');
    num = num + 1;
    if option == 1.0
        loc = inputdlg('Location: ');
        hold on;
        peak = str2double(loc);
        line([peak; peak],[0; max(Raw_data_y)*1.2],'color','r',...
            'LineStyle','--','LineWidth',0.5,'Clipping','on');
        text(peak+20,max(Raw_data_y)*0.05*num,loc,'FontSize',8);
        hold off;
    elseif option == 2.0
        pass = true;
        return;
    end
    
end
    