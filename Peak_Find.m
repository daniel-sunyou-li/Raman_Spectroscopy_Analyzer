%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   POLYNOMIAL PEAK FITTING v 1.2                    %%%
%%%                      Last Updated: 7/12/2018                       %%%          
%%%              Author:  Daniel Li (daniel_li2@brown.edu)             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                 
%%% Reads in a single data set to be used for polynomial peak fitting.  
%%% First, the raw data center is chosen; second, the number of points 
%%% centered around the raw data center is selected.  The user is shown
%%% the highlighted region of the raw data to be used in the fit, and is
%%% given an option to re-select the number of points.  Third, the user
%%% can select from a polynomial order of 1 (linear) to 4.  Lastly, plots
%%% are produced of the fit computed superimposed on the raw data with the
%%% center(s) marked by a line.  Additionally, the adjusted chi-square is
%%% included as a metric for the fit.

%%% USER-INPUT SELECT DATA FILE
[File,path] = uigetfile('*.txt');
File_path = strcat(path,File);

Raw_data = dlmread(File_path,'\t');
Raw_data_x = Raw_data(:,1);
Raw_data_y = Raw_data(:,2);

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

%%% USER-INPUT FOR THE NUMBER OF DATA SETS TO BE AVERAGED ACROSS
multiple_data = menu('Is this data set an average of multiple sets?',...
    'Yes','No');
mScale = 0;
if multiple_data == 1.0
    mScale = 1;
end

%%% PLOT DATA FOR USER REFERENCE WHEN DETERMINING PEAK LOCATIONS
pos1 = [200 500 500 300];
f1 = figure('Name','Raw Data','NumberTitle','off','Position',pos1);
plot(Raw_data_x,Raw_data_y);
title('Raw Data');
xlabel(Units);
ylabel('Counts');
grid on;
grid minor;

done = false;

%%% FIGURE TO HAVE SELECTED DATA FOR FIT BE PLOTTED ON
pos2 = [500 0 500 500];
f2 = figure('Name','Selected Data','NumberTitle','off','Visible','off',...
    'Position',pos2);

while done == false
%%% USER-INPUT ESTIMATED PEAK LOCATION
    Peak_est_cell = inputdlg('Estimated Peak Location:  ');
    Peak_est = str2double(Peak_est_cell{1});
    while (Peak_est < min(Raw_data_x)) || (Peak_est > max(Raw_data_x))
       Peak_est_cell = inputdlg(['Input out of bounds, choose again.',...
           ' Type ''-1'' to terminate:  ']);
       Peak_est = str2double(Peak_est_cell{1});
       if Peak_est == -1
          return; 
       end
    end

%%% USER-INPUT FIT POINTS USED
    Fit_points_c = inputdlg('Points for fit? ');
    Fit_points = str2double(Fit_points_c{1});
    x_res = Raw_data_x(2) - Raw_data_x(1);
    range = Fit_points * x_res;
    
    while (round(Fit_points,0) ~= Fit_points) || (Fit_points < 0) ||...
            (Peak_est - (range/2)) < min(Raw_data_x)
        if (round(Fit_points,0) ~= Fit_points) || (Fit_points < 0)
            Fit_points_c = inputdlg(['Dec./Neg. not allowed, choose again.',...
                ' Type ''-1'' to terminate: ']);
            Fit_points = str2double(Fit_points_c{1});
            range = Fit_points * x_res;
            if Fit_points == -1
                return;
            end
        elseif (Peak_est - (range/2)) < min(Raw_data_x)
            Fit_points_c = inputdlg(['Range out of bounds, choose again. ',...
               ' Type ''-1'' to terminate:  ']);
            Fit_points = str2double(Fit_points_c{1});
            range = Fit_points * x_res;
            if Fit_points == -1
              return; 
            end
        end
    end

%%% CREATE NEW ARRAY OF TRIMMED RAW DATA TO USE FOR FITTING
    x_ind = [];

    for n = 2:length(Raw_data_x)
%%% COMPUTE n-1, n, n+1 TO DETERMINE IF n IS THE CLOSEST VALUE TO APP.
        prev = Raw_data_x(n-1);
        curr = Raw_data_x(n);
        foll = Raw_data_x(n+1);
        
        dprev = abs(Peak_est - prev);
        dcurr = abs(Peak_est - curr);
        dfoll = abs(Peak_est - foll);
        
        if (dcurr < dprev) && (dcurr < dfoll)
            x_ind = n;
            break 
        end
    end

%%% CREATES NEW FIT DATA SET
    Fit_data_x = zeros(0,Fit_points);
    Fit_data_y = zeros(0,Fit_points);
    x_init = [];

    if mod(Fit_points,2) == 1
        x_init = x_ind - ((Fit_points - 1) / 2);
    elseif mod(Fit_points,2) == 0
        x_init = x_ind - (Fit_points/2);  
    end

%%% POPULATE Fit_data_x AND Fit_data_y
    for m = 1:Fit_points
        Fit_data_x(m) = Raw_data_x(x_init + m);
        Fit_data_y(m) = Raw_data_y(x_init + m);
    end

%%% VISUALIZE FIT POINTS
    plot(Fit_data_x, Fit_data_y,'r','LineWidth',3);
    hold on;
    plot(Raw_data_x,Raw_data_y,'k','LineWidth',0.5,'LineStyle',':');
    hold off;
    xlabel(['Frequency ', Units]);
    ylabel('Counts (AU)');
    title('Fit Data vs Raw Data');
    grid on;
    grid minor;
    legend('Fit Data','Raw Data','Location', 'northwest');
    movegui(f2,'east');
%%% USER INPUT TO CHOOSE NEW PARAMETERS FOR THE FIT
    cont_choice = menu('Proceed with analysis?','Yes',...
        'Choose Different Peak/Range','Exit');
    if cont_choice == 3.0
        return;
    elseif cont_choice == 2.0
        clf(f2);
        f2.Visible = 'off';
    elseif cont_choice == 1.0
        done = true;
    end
end

%%% USER INPUT ON POLYNOMIAL FIT ORDER
poly_order = fix(menu('Choose polynomial fit order','1','2','3','4'));

%%% EVALUATE POLYNOMIAL COEFFICIENTS IN DESCENDING ORDER
[poly,S,mu] = polyfit(Fit_data_x,Fit_data_y,poly_order);

%%% POPULATE ARRAY HOLDING COEFFICIENT DATA POINTS
poly_linspace = linspace(Fit_data_x(1) - Fit_points * x_res...
    , Fit_data_x(Fit_points) + Fit_points * x_res);
[poly_y,delta] = polyval(poly,poly_linspace,S,mu);

%%% PLOT POLYNOMIAL FIT DATA SUPERIMPOSED ON SELECTED DATA FOR FIT
f3 = figure('Name','Polynomial Fit','NumberTitle','off');
poly_plot = plot(poly_linspace,poly_y,'Color',[0.1,0.5,0.2],'LineWidth',1);
hold on;
fit_plot = plot(Fit_data_x,Fit_data_y,'Color',[0.1,0.1,0.1],...
    'Marker','o','LineStyle','none');
title('Polynomial Fit vs. Fit Data');
xlabel(['Frequency ', Units]);
ylabel('Counts (AU)');
grid on;
grid minor;
hold off;

%%% EVALUATING THE MAXMIMUM OF THE POLYNOMIAL
%%% DIFFERENTIATE THE POLYNOMIAL
poly_der = polyder(poly);

%%% IDENTIFY ROOTS OF THE DERIVATIVE OF THE POLYNOMIAL
peak_loc = (roots(poly_der).*mu(2))+mu(1);

%%% ADD THE EVALUATED MAXIMUM AS A LINE ON PREVIOUS PLOT
hold on;

%%% PLOTS MULTIPLE CENTERS
for o = 1:length(peak_loc)
    line([peak_loc(o); peak_loc(o)],[0; max(poly_y)*1.5],...
        'Color','k','LineStyle','--','LineWidth',0.5);
    peak_str = num2str(peak_loc(o));
    text(peak_loc(o),max(poly_y)*1.4,...
        ['  ',peak_str],'FontSize',9);
end

hold off;

%%% PLOT MODIFICATION
fit_legend = legend([poly_plot,fit_plot],'Polynomial Fit','Fit Data',...
    'Location', 'northwest');
xlim([Fit_data_x(1) - Fit_points * x_res,...
    Fit_data_x(Fit_points) + Fit_points * x_res]);
ylim([0,max(poly_y)*1.5]);

%%% GOODNESS OF FIT EVALUATION USING ADJ CHI-SQUARE
[poly_ref,delta2] = polyval(poly,Fit_data_x,S,mu);

%%% DOF = # DATA POINTS - # COEFFICIENTS IN POLYNOMIAL FIT
chisq = 0;

for p = 1:length(Fit_data_y)
   chisq = chisq + (((Fit_data_y(p) - poly_ref(p))^2)/poly_ref(p));
end

adj_chisq = round(chisq / (S.df - mScale),3);
adj_chisq = num2str(adj_chisq);
DOF = num2str(S.df - mScale);

title(fit_legend,['\chi^2_{adj}',' :',adj_chisq, ' DOF: ',DOF],...
    'interpreter','tex');

movegui(f3,'center');

msgbox('Reference ''chi sq p value.png'' to qualify fit.');

return;


