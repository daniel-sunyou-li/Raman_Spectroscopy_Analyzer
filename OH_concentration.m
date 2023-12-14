%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        OH Concentration Analysis (Raman & Absorption) v 1.0        %%%
%%%              Author:  Daniel Li (daniel_li2@brown.edu)             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                 
%%% OH Concentration Determination:  reads in two data files--one with
%%% known OH concentration, the other with unknown concentration--of either
%%% Raman or FTIR spectra and will compute the OH concentration of the
%%% unknown sample based upon a user-input "known" OH concentration.  There
%%% will also be an option to select a standard reference sample if a 
%%% known sample is not available.

%%% Select between Raman or FTIR determination
option = menu('Choose either Raman or Absorption analysis: ',...
    'Raman Spectroscopy','FT Infrared Absorption');
analysis_type = 0;
if option == 1.0
   analysis_type = 0; % Raman
elseif option == 2.0
    analysis_type = 1; % Absorption
end

%%% Ask user if a reference file will be used or not for Raman
%%% Absorption does not need reference because of Beer-Lambert relation
if analysis_type == 0
    option = menu('Will you be using reference data?','Yes','No','Exit');
    reference = 0;
    if option == 1.0
        reference = 1; % Yes reference
    elseif option == 2.0
        reference = 0; % No referenece
    elseif option == 3.0
        return;
    end
end

%%% Ask user for x-units if using absorption
Units_x = 'cm-1';

if analysis_type == 1
   Units_x = 'nm';
   option = menu('What units are the x-values?','nm','eV','Exit');
   if option == 1.0
       Units_x = 'nm'; % nm
   elseif option == 2.0
       Units_x = 'eV'; % eV
   elseif option == 3.0
       return;
   end
end
    
%%% Upload data files
%%% Selecting reference data
if reference == 1 % Checks if reference data is in appropriate range
    pass = false;
    if analysis_type == 0 % Checks raman data
        while pass == false
            uiwait(msgbox('Choose your reference data.'));
            [file1,path1] = uigetfile('*.txt');
            file_path1 = strcat(path1,file1);

            Ref_data = dlmread(file_path1,'\t');
            Ref_data_x = Ref_data(:,1);
            Ref_data_y = Ref_data(:,2);

            %%% Assumes cm-1 for x-values
            if max(Ref_data_x) < 3650 || min(Ref_data_x) > 400
                uiwait(msgbox(['Data must be available in the range',...
                    '400 to 4000']));
            else
                figure;
                plot(Ref_data_x,Ref_data_y,'k','LineWidth',0.5,...
                    'LineStyle',':');
                title(file1);
                xlabel(['Frequency ',Units_x]);
                ylabel('Units (AU)');
                grid on;
                grid minor;
                
                option = menu('Use this as reference data?','Yes','No');
                if option == 1.0
                    pass = true;
                end
            end
        end
    elseif analysis_type == 1 % Checks absorption data
        while pass == false
            msgbox('Choose your reference data.');
            [file1,path1] = uigetfile('*.txt');
            file_path1 = strcat(path1,file1);

            Ref_data = dlmread(file_path1,'\t');
            Ref_data_x = Ref_data(:,1);
            Ref_data_y = Ref_data(:,2);
            
            if strcmp('nm',Units_x)
                if max(Ref_data_x) < 1300 || min(Ref_data_x) > 3000
                    msgbox(['Data must be available between ',...
                    '1300 and 3000 nm']);
                else
                    figure;
                    plot(Ref_data_x,Ref_data_y,'k','LineWidth',0.5,...
                        'LineStyle',':');
                    title(file1);
                    xlabel(['Frequency ',Units_x]);
                    ylabel('Units (AU)');
                    grid on;
                    grid minor;

                    option = menu('Use this as reference data?','Yes','No');
                    if option == 1.0
                        pass = true;
                    end
                end
            elseif strcmp('eV',Units_x)
                if max(Ref_data_x) < 0.41 || min(Ref_data_x) > 1.0
                    msgbox(['Data must be available between ',...
                    '0.4 and 1.0 eV']);
                else
                    figure;
                    plot(Ref_data_x,Ref_data_y,'k','LineWidth',0.5,...
                        'LineStyle',':');
                    title(file1);
                    xlabel(['Frequency ',Units_x]);
                    ylabel('Units (AU)');
                    grid on;
                    grid minor;

                    option = menu('Use this as reference data?','Yes','No');
                    if option == 1.0
                        pass = true;
                    end
                end
            end
        end
    end
    if pass == false
        option = menu('Data selected was not valid, choose again?',...
            'Yes','No');
        if option == 2.0
            return;
        end
    end
end

%%% Selecting unknown data set
uiwait(msgbox(['Choose your unknown data.',...
    ' Make sure to use the same units']));

pass = false;
while pass == false
    [file2,path2] = uigetfile('*.txt');
    file_path2 = strcat(path2,file2);

    Unk_data = dlmread(file_path2,'\t');
    Unk_data_x = Unk_data(:,1);
    Unk_data_y = Unk_data(:,2);
    
    if analysis_type == 0 % Checks raman data
        if max(Unk_data_x) < 3650 || min(Unk_data_x) > 400
            msgbox('Data must be available in the range 400 to 3650');
        else
            figure;
            plot(Unk_data_x,Unk_data_y,'k','LineWidth',0.5,...
                'LineStyle',':');
            title(file2);
            xlabel(['Frequency ',Units_x]);
            ylabel('Units (AU)');
            grid on;
            grid minor;

            option = menu('Use this as unknown data?','Yes','No');
            if option == 1.0
                pass = true;
            end
        end
    elseif analysis_type == 1 % Checks absorption data
        if strcmp('nm',Units_x)
            if max(Unk_data_x) < 1300 || min(Unk_data_x) > 3000
                msgbox(['Data must be available between ',...
                '1300 and 3000 nm']);
            else
                figure;
                plot(Unk_data_x,Unk_data_y,'k','LineWidth',0.5,...
                    'LineStyle',':');
                title(file2);
                xlabel(['Frequency ',Units_x]);
                ylabel('Units (AU)');
                grid on;
                grid minor;
                
                option = menu('Use this as unknown data?','Yes','No');
                if option == 1.0
                    pass = true;
                end
            end
        elseif strcmp('eV',Units_x)
            if max(Unk_data_x) < 0.41 || min(Unk_data_x) > 1.0
                msgbox(['Data must be available between ',...
                '0.4 and 1.0 eV']);
            else
                figure;
                plot(Ref_data_x,Ref_data_y,'k','LineWidth',0.5,...
                    'LineStyle',':');
                title(file1);
                xlabel(['Frequency ',Units_x]);
                ylabel('Units (AU)');
                grid on;
                grid minor;
                
                option = menu('Use this as reference data?','Yes','No');
                if option == 1.0
                    pass = true;
                end
            end
        end
    end
end

%%% If user selects a reference file, then ask for OH concentration.
%%% Give option to choose known concentration in either 'ppm' or 'cm-3'

if reference == 1
   Units_y = 'ppm';
       option = menu('What are the units of the reference data','ppm',...
           'cm-3');
   if option == 2.0
       Units_y = 'cm-3';
   end
   input = inputdlg(['What is the concentration in: ',Units_y],...
       'Reference Sample',[1 40]);
   ref_conc = str2double(input);
   while ref_conc < 0 
       input = inputdlg('Negative values forbidden. Input new value: ',...
           'Reference Sample',[1 40]);
       ref_conc = str2double(input);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                      Raman Spectroscopy                      %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute the maximum intensities of 440/1600 peak to 3650 peak using 
%%% polynomial fit of the data and ask the user how many points should be
%%% used for the polynomial fit
if analysis_type == 0

    Peak_est = [440,1600,3650];
    peak_ind = [0,0,0]; % index numbers of estimated peak locations

    for n = 1:3
        for m = 2:length(Unk_data_x)
            prev = Unk_data_x(m-1);
            curr = Unk_data_x(m);
            foll = Unk_data_x(m+1);

            % Determine the abs. difference of prev, curr and foll to
            % relevant points
            dprev = abs(Peak_est(n) - prev);
            dcurr = abs(Peak_est(n) - curr);
            dfoll = abs(Peak_est(n) - foll);

            % Conditional to determine if curr is approx peak_est
            if (dcurr < dprev) && (dcurr < dfoll)
                peak_ind(n) = m;
                break
            end
        end
    end

    %%% Initialize arrays that will contain data used for fitting
    peak_structs = struct('peak',{},'x',{},'y',{},'coeff',{},'xmax',{},...
        'ymax',{});
    peak_structs(1).peak = peak_ind(1);
    peak_structs(2).peak = peak_ind(2);
    peak_structs(3).peak = peak_ind(3);

    %%% User-input how many points used to fit polynomial
    pass = false;

    while pass == false
        fit_points = inputdlg({'440 Peak','1600 Peak','3650 Peak'},...
            'Polynomial Fit Points (Unknown)',[1 55]);

        for n = 1:3
            num = str2double(fit_points(n));
            % Checks if polynomial fit points are valid values
            while round(num,0) ~= num 
                fit_points(n) = inputdlg(['Choose a different quantity',...
                    ' of fit points.  Type ''0'' to exit']);
                num = fit_points(n);
                if num == 0
                    return;
                end
            end

            % Define starting value for peak_structs(n)
            x_init = [];
            if mod(num,2) == 1
                x_init = peak_structs(n).peak - ((num - 1) / 2);
            elseif mod(num,2) == 0
                x_init = peak_structs(n).peak - (num/2);
            end

            % Populates the peak structs with x- and y-values
            for m = 1:num
                peak_structs(n).x(m) = Unk_data_x(x_init + m - 1);
                peak_structs(n).y(m) = Unk_data_y(x_init + m - 1);
            end
        end

        f1 = figure;
        plot(Unk_data_x,Unk_data_y,'k','LineWidth',0.5,'LineStyle',':');
        for n = 1:3
            hold on;
            plot(peak_structs(n).x,peak_structs(n).y,'r','LineWidth',3);
            hold off;
        end
        xlabel(['Frequency ', Units_x]);
        ylabel('Counts (AU)');
        title('Fit Data vs Raw Data');
        grid on;
        grid minor;
        legend('Raw Data','Fit Data','Location','northwest');

        option = menu('Proceed with this data?','Yes',...
            'Select Again','Exit');
        if option == 1.0
           pass = true; 
        elseif option == 2.0
            pass = false;
        elseif option == 3.0
            return;
        end
    end

    %%% Subtract off the baseline for 3650 peak using 3400 and 4000 values
    %%% User-input for what the baseline range should be
    pass = false;

    base_struct = struct('index',{},'x',{},'y',{});

    while pass == false
        base_est = inputdlg({'Baseline (left): ','Baseline (right): ',...
            'Points to use: '},'Baseline Parameters',[1 50]);
        num = str2double(base_est);
        while num(1) > num(2) 
            base_est = inputdlg({'Baseline (left): ','Baseline (right): ',...
                'Points to use: '},['Baseline "left" must be less than'...
                ' baseline "right."  Choose again.'],[1 85]);
            num = str2double(base_est);
        end

        %%% Populate baseline arrays
        for n = 1:2
            % Identify baseline indeces
            for m = 2:length(Unk_data_x)
                prev = Unk_data_x(m-1);
                curr = Unk_data_x(m);
                foll = Unk_data_x(m+1);

                dprev = abs(num(n) - prev);
                dcurr = abs(num(n) - curr);
                dfoll = abs(num(n) - foll);

                if (dcurr < dprev) && (dcurr < dfoll)
                    base_struct(n).index = m;
                    break
                end
            end
        end

        for n = 1:num(3)
            base_struct(1).x(n) = Unk_data_x(base_struct(1).index - ...
                num(3) + n - 1);
            base_struct(1).y(n) = Unk_data_y(base_struct(1).index - ...
                num(3) + n - 1);
            base_struct(2).x(n) = Unk_data_x(base_struct(2).index + ...
                n - 1);
            base_struct(2).y(n) = Unk_data_y(base_struct(2).index + ...
                n - 1);
        end

        %%% Compute baseline value
        baseline = (mean(base_struct(1).y) + mean(base_struct(2).y))/2;
        baseline_array = struct('x',[],'y',[]);
        count = 1;
        for n = base_struct(1).index:base_struct(2).index
            baseline_array.x(count) = Unk_data_x(n);
            baseline_array.y(count) = baseline;
            count = count + 1;
        end

        %%% Visualize baseline
        figure;
        plot(Unk_data_x,Unk_data_y,'k','LineWidth',0.5,'LineStyle',':');
        for n = 1:2
            hold on;
            plot(base_struct(n).x,base_struct(n).y,'r','LineWidth',3);
            hold off;
        end
        hold on;
        plot(baseline_array.x, baseline_array.y,'b','LineWidth',3);
        xlabel(['Frequency ', Units_x]);
        ylabel('Counts (AU)');
        title('Baseline Fit');
        grid on;
        grid minor;
        leg = legend('Raw Data','Fit Data','Location','northwest');
        title(leg,['Baseline = ',num2str(baseline)]);
        hold off;

        option = menu('Proceed with this base?','Yes',...
            'Select Again','Exit');
        if option == 1.0
           pass = true; 
        elseif option == 2.0
            pass = false;
        elseif option == 3.0
            return;
        end
    end

%%% Normalize the 3650 peak values by subtracting off the baseline value
    for n = 1:str2double(fit_points(3))
        peak_structs(3).y(n) = peak_structs(3).y(n) - baseline; 
    end

%%% After obtaining the desired parameters and arrays containing x- and
%%% y-values for the raw data, fit the 440, 1600, and 3650 peaks to
%%% polynomials and evaluate each at their respective maximum intensities

%%% Using polynoial of order n = 2
%%% Use the 'coeff' field in the peak_structs structure to hold the 
%%% polynomial coefficients
    for n = 1:3
        x_cent = mean(peak_structs(n).x);
        peak_structs(n).coeff = polyfit(peak_structs(n).x - x_cent,...
            peak_structs(n).y,2);
        poly_der = polyder(peak_structs(n).coeff);
        peak_structs(n).xmax = roots(poly_der);
        peak_structs(n).ymax = polyval(peak_structs(n).coeff,...
            peak_structs(n).xmax);
    end

    %%% Evaluate the ratios of the intensities
    unk_ratio_1 = peak_structs(3).ymax / peak_structs(1).ymax;
    unk_ratio_2 = peak_structs(3).ymax / peak_structs(2).ymax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%               Raman: Without Reference Sample                %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% If not using reference sample:
%%% From Mikkelsen, Galeener, and Mosby (1981), 1000 ppm for 490 peak
%%% From Walrafen and Stone (1975), 0.1% (or 1000 ppm) for 1600 peak
%%% To convert from cm-3 to ppm:  
%%% 1.) Determine Fused Silica density in g/cm3
%%% 2.) Convert into cm-3 by dividing by MW for Fused Silica
%%% 3.) Divide [I3650/I440 * 1.3E22 cm-3] by [FS] cm-3
    if reference == 0
        ref490_ppm = 1000.00;
        ref1600_ppm = 1000.00;
        rat_ref490 = 6 * 10^(-3);
        OH490_ppm = ( unk_ratio_1 / rat_ref490 ) * ref490_ppm;
        rho_sio2 = inputdlg('Input density of fused silica sample (g/cm3): ');
        OH490_cm3 = OH490_ppm *  str2double(rho_sio2) * 60.08;
        ref490_cm3 = ref490_ppm * str2double(rho_sio2) * 60.08;

        % Plot the concentration as a bar graph in which there are two
        % categories for standard reference and unknown
        option = menu('Plot in cm-3 or in ppm?','cm-3','ppm');

        if option == 1.0
            figure;
            bar0 = bar(categorical({'Mikkelsen et al (1981)','Unknown'}),...
                [ref490_cm3 OH490_cm3],0.4,'FaceColor',[0.6 0.1 0.1]);
            text(1:2,[ref490_cm3 OH490_cm3],...
                num2str([ref490_cm3 OH490_cm3]),'vert','bottom',...
                'horiz','center');
            title('OH Concentration (Without Reference)');
            ylabel('Concentration (cm-3)');
            set(gca,'YMinorTick','on','YGrid','on','GridLineStyle','-');
        elseif option == 2.0
            bar0 = bar(categorical({'Reference','Unknown'}),...
                [ref490_ppm OH490_ppm],0.4,'FaceColor',[0.1 0.6 0.1]);
            text(1:2,[ref490_ppm OH490_ppm],...
                num2str([ref490_ppm OH490_ppm]'),'vert','bottom',...
                'horiz','left');
            title('OH Concentration (Without Reference)');
            ylabel('Concentration (ppm wt)');
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                Raman: Using Reference Sample                 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 1.) Take ratio of I3650/I440 = S, for both Ref. and Unk.
%%% 2.) Divide S(Unk.) by S(Ref.) and then multiply by User-input conc.
%%% 3.) Convert cm-3 -> ppm or ppm -> cm-3

%%% First must repeat the peak intensity ratio evaluation of the reference
    if reference == 1
        Peak_est = [440,1600,3650];
        peak_ind = [0,0,0];

        for n = 1:3
            for m = 2:length(Ref_data_x)
                prev = Ref_data_x(m-1);
                curr = Ref_data_x(m);
                foll = Ref_data_x(m+1);

                dprev = abs(Peak_est(n) - prev);
                dcurr = abs(Peak_est(n) - curr);
                dfoll = abs(Peak_est(n) - foll);

                if (dcurr < dprev) && (dcurr < dfoll)
                    peak_ind(n) = m;
                    break
                end
            end
        end

        peak_structs = struct('peak',{},'x',{},'y',{},'coeff',{},'xmax',{},...
            'ymax',{});
        peak_structs(1).peak = peak_ind(1);
        peak_structs(2).peak = peak_ind(2);
        peak_structs(3).peak = peak_ind(3);

        pass = false;

        while pass == false
            input = inputdlg({'440 Peak','1600 Peak','3650 Peak'},...
                'Polynomial Fit Points (Reference)',[1 55]);
            fit_points = str2double(input);
            for n = 1:3

                while round(fit_points(n),0) ~= fit_points(n) 
                    input(n) = inputdlg(['Choose a different quantity',...
                        ' of fit points.  Type ''0'' to exit']);
                    fit_points(n) = str2double(input(n));
                    if fit_points(n) == 0
                        return;
                    end
                end

                x_init = [];
                if mod(fit_points(n),2) == 1
                    x_init = peak_structs(n).peak - ((fit_points(n) - 1) / 2);
                elseif mod(fit_points(n),2) == 0
                    x_init = peak_structs(n).peak - (fit_points(n)/2);
                end

                for m = 1:fit_points(n)
                    peak_structs(n).x(m) = Ref_data_x(x_init + m - 1);
                    peak_structs(n).y(m) = Ref_data_y(x_init + m - 1);
                end
            end

            figure;
            plot(Ref_data_x,Ref_data_y,'k','LineWidth',0.5,'LineStyle',':');
            for n = 1:3
                hold on;
                plot(peak_structs(n).x,peak_structs(n).y,'r','LineWidth',3);
                hold off;
            end
            xlabel(['Frequency ', Units_x]);
            ylabel('Counts (AU)');
            title('Raw Data vs Ref Data');
            grid on;
            grid minor;
            legend('Raw Data','Ref Data','Location','northwest');

            option = menu('Proceed with this data?','Yes',...
                'Select Again','Exit');
            if option == 1.0
               pass = true; 
            elseif option == 2.0
                pass = false;
            elseif option == 3.0
                return;
            end
        end

        pass = false;

        base_struct = struct('index',{},'x',{},'y',{});

        while pass == false
            input = inputdlg({'Baseline (left): ','Baseline (right): ',...
                'Points to use: '},'Baseline Parameters');
            base_par = str2double(input);
            while base_par(1) > base_par(2) 
                input = inputdlg({'Baseline (left): ','Baseline (right): ',...
                    'Points to use: '},['Baseline "left" must be less than'...
                    ' baseline "right."  Choose again']);
                base_par = str2double(input);
            end

            for n = 1:2

                for m = 2:length(Ref_data_x)
                    prev = Ref_data_x(m-1);
                    curr = Ref_data_x(m);
                    foll = Ref_data_x(m+1);

                    dprev = abs(base_par(n)- prev);
                    dcurr = abs(base_par(n) - curr);
                    dfoll = abs(base_par(n) - foll);

                    if (dcurr < dprev) && (dcurr < dfoll)
                        base_struct(n).index = m;
                        break
                    end
                end
            end

            for n = 1:base_par(3)
                base_struct(1).x(n) = Ref_data_x(base_struct(1).index - ...
                    base_par(3) + n - 1);
                base_struct(1).y(n) = Ref_data_y(base_struct(1).index - ...
                    base_par(3) + n - 1);
                base_struct(2).x(n) = Ref_data_x(base_struct(2).index + ...
                    n - 1);
                base_struct(2).y(n) = Ref_data_y(base_struct(2).index + ...
                    n - 1);
            end

            %%% Compute baseline value
            baseline = (mean(base_struct(1).y) + mean(base_struct(2).y))/2;
            baseline_array = struct('x',[],'y',[]);
            count = 1;
            for n = base_struct(1).index:base_struct(2).index
                baseline_array.x(count) = Ref_data_x(n);
                baseline_array.y(count) = baseline;
                count = count + 1;
            end

            %%% Visualize baseline
            figure;
            plot(Ref_data_x,Ref_data_y,'k','LineWidth',0.5,'LineStyle',':');
            for n = 1:2
                hold on;
                plot(base_struct(n).x,base_struct(n).y,'r','LineWidth',3);
                hold off;
            end
            hold on;
            plot(baseline_array.x, baseline_array.y,'b','LineWidth',3);
            xlabel(['Frequency ', Units_x]);
            ylabel('Counts (AU)');
            title('Baseline Fit of Reference Data');
            grid on;
            grid minor;
            leg = legend('Raw Data','Reference Data','Location','northwest');
            title(leg,['Baseline = ',num2str(baseline)]);
            hold off;

            option = menu('Proceed with this base?','Yes',...
                'Select Again','Exit');
            if option == 1.0
               pass = true; 
            elseif option == 2.0
                pass = false;
            elseif option == 3.0
                return;
            end
        end

    %%% Normalize the 3650 peak values by subtracting off the baseline value
        for n = 1:fit_points(3)
            peak_structs(3).y(n) = peak_structs(3).y(n) - baseline; 
        end

    %%% After obtaining the desired parameters and arrays containing x- and
    %%% y-values for the raw data, fit the 440, 1600, and 3650 peaks to
    %%% polynomials and evaluate each at their respective maximum intensities

    %%% Using polynoial of order n = 2
    %%% Use the 'coeff' field in the peak_structs structure to hold the 
    %%% polynomial coefficients
        for n = 1:3
            x_cent = mean(peak_structs(n).x);
            peak_structs(n).coeff = polyfit(peak_structs(n).x - x_cent,...
                peak_structs(n).y,2);
            poly_der = polyder(peak_structs(n).coeff);
            peak_structs(n).xmax = roots(poly_der);
            peak_structs(n).ymax = polyval(peak_structs(n).coeff,...
                peak_structs(n).xmax);
        end

        % Evaluate the ratios of the intensities
        ref_ratio_1 = peak_structs(3).ymax / peak_structs(1).ymax;
        ref_ratio_2 = peak_structs(3).ymax / peak_structs(2).ymax;

        input = inputdlg('Input density of fused silica sample (g/cm3): ');
        rho_sio2 = str2double(input);

        if strcmp(Units_y,'ppm')
            ref_conc_ppm = ref_conc;
            ref_conc_cm3 = ref_conc_ppm * str2double(rho_sio2) * 60.08;
        elseif strcmp(Units_y,'cm-3')
            ref_conc_cm3 = ref_conc;
            ref_conc_ppm = ref_conc_cm3 / str2double(rho_sio2) / 60.08;
        end

        OH490_ppm = ( unk_ratio_1 / ref_ratio_1 ) * ref_conc_ppm;
        OH490_cm3 = ( unk_ratio_1 / ref_ratio_1 ) * ref_conc_cm3;

        % Plot the concentration as a bar graph in which there are two
        % categories for standard reference and unknown

        % Ask to plot in ppm or cm-3
        option = menu('Plot in cm-3 or in ppm?','cm-3','ppm');

        if option == 1.0
            figure;
            bar0 = bar(categorical({'Reference','Unknown'}),...
                [ref_conc_cm3 OH490_cm3],0.4,'FaceColor',[0.6 0.1 0.1]);
            text(1:2,[ref_conc_cm3 OH490_cm3],...
                num2str([ref_conc_cm3 OH490_cm3]),'vert','bottom',...
                'horiz','center');
            title('OH Concentration (With Reference)');
            ylabel('Concentration (cm-3)');
            set(gca,'YMinorTick','on','YGrid','on','GridLineStyle','-');
        elseif option == 2.0
            bar0 = bar(categorical({'Reference','Unknown'}),...
                [ref_conc_ppm OH490_ppm],0.4,'FaceColor',[0.1 0.6 0.1]);
            text(1:2,[ref_conc_ppm OH490_ppm],...
                num2str([ref_conc_ppm OH490_ppm]'),'vert','bottom',...
                'horiz','left');
            title('OH Concentration (With Reference)');
            ylabel('Concentration (ppm wt)');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                    FT Infrared Absorption                    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analysis_type == 1

%%% First identify how many absorption peaks will be used in computing
%%% the concentrations.  Will be checking for the peaks: 1400 nm (0.86 eV),
%%% 2200 nm (0.56 eV), 2750 nm (0.45 eV)
    peaks_str = {};

    if strcmp(Units_x,'eV')
        if min(Unk_data_x) < 0.45 && max(Unk_data_x) < 0.56
            peak_est = [0.45];
            peaks_str = {'0.45 eV Peak'};
        elseif min(Unk_data_x) > 0.45 && max(Unk_data_x) < 0.86
            peak_est = [0.56];
            peaks_str = {'0.56 eV Peak'};
        elseif min(Unk_data_x) > 0.56
            peak_est = [0.86];
            peaks_str = {'0.86 eV Peak'};
        elseif min(Unk_data_x) < 0.45 && max(Unk_data_x) < 0.86
            peak_est = [0.45 0.56];
            peaks_str = {'0.45 eV Peak','0.56 Peak'};
        elseif min(Unk_data_x) > 0.45 && max(Unk_data_x) > 0.86
            peak_est = [0.56 0.86];
            peaks_str = {'0.56 eV Peak','0.86 Peak'};
        elseif min(Unk_data_x) < 0.45 && max(Unk_data_x) > 0.86
            peak_est = [0.45 0.56 0.86];
            peaks_str = {'0.45 eV Peak','0.56 eV Peak','0.86 eV Peak'};
        end
    elseif strcmp(Units_x,'nm')
        if min(Unk_data_x) < 1400 && max(Unk_data_x) < 2200
            peak_est = [1400];
            peaks_str = {'1400 nm Peak'};
        elseif min(Unk_data_x) > 1400 && max(Unk_data_x) < 2750
            peak_est = [2200];
            peaks_str = {'2200 nm Peak'};
        elseif min(Unk_data_x) > 2200
            peak_est = [2750];
            peaks_str = {'2750 nm Peak'};
        elseif min(Unk_data_x) < 1400 && max(Unk_data_x) < 2750
            peak_est = [1400 2200];
            peaks_str = {'1400 nm Peak','2200 nm Peak'};
        elseif min(Unk_data_x) > 1400 && max(Unk_data_x) > 2750
            peak_est = [2200 2750];
            peaks_str = {'2200 nm Peak','2750 nm Peak'};
        elseif min(Unk_data_x) < 1400 && max(Unk_data_x) > 2750
            peak_est = [1400 2200 2750];
            peaks_str = {'1400 nm Peak','2200 nm Peak','2750 nm Peak'};
        end       
    end
    
    peaks = length(peak_est);
    
    peak_ind = [0,0,0]; % index numbers of estimated peak locations
    peak_structs = struct('peak',{},'x',{},'y',{},'coeff',{},'xmax',{},...
        'ymax',{});
    
    for n = 1:peaks
        for m = 2:length(Unk_data_x)
            prev = Unk_data_x(m-1);
            curr = Unk_data_x(m);
            foll = Unk_data_x(m+1);

            % Determine the abs. difference of prev, curr and foll to
            % relevant points
            dprev = abs(peak_est(n) - prev);
            dcurr = abs(peak_est(n) - curr);
            dfoll = abs(peak_est(n) - foll);

            % Conditional to determine if curr is approx peak_est
            if (dcurr < dprev) && (dcurr < dfoll)
                peak_ind(n) = m;
                peak_structs(n).peak = peak_ind(n);
                break
            end
        end
    end

%%% Select range of points to fit the absorption peaks
    pass = false;

    while pass == false
        fit_points = inputdlg(peaks_str,'Polynomial Fit Points (Unknown)',[1 55]);

        for n = 1:peaks
            num = str2double(fit_points(n));
            % Checks if polynomial fit points are valid values
            while round(num,0) ~= num 
                fit_points(n) = inputdlg(['Choose a different quantity',...
                    ' of fit points.  Type ''0'' to exit']);
                num = fit_points(n);
                if num == 0
                    return;
                end
            end

            % Define starting value for peak_structs(n)
            x_init = [];
            if mod(num,2) == 1
                x_init = peak_structs(n).peak - ((num - 1) / 2);
            elseif mod(num,2) == 0
                x_init = peak_structs(n).peak - (num/2);
            end

            % Populates the peak structs with x- and y-values
            for m = 1:num
                peak_structs(n).x(m) = Unk_data_x(x_init + m - 1);
                peak_structs(n).y(m) = Unk_data_y(x_init + m - 1);
            end
        end

        f1 = figure;
        plot(Unk_data_x,Unk_data_y,'k','LineWidth',0.5,'LineStyle',':');
        for n = 1:3
            hold on;
            plot(peak_structs(n).x,peak_structs(n).y,'r','LineWidth',3);
            hold off;
        end
        xlabel(['Frequency ', Units_x]);
        ylabel('Absorption (AU)');
        title('Fit Points');
        grid on;
        grid minor;
        legend('Raw Data','Fit Points','Location','northwest');

        option = menu('Proceed with this data?','Yes',...
            'Select Again','Exit');
        if option == 1.0
           pass = true; 
        elseif option == 2.0
            pass = false;
        elseif option == 3.0
            return;
        end
    end    

%%% Polynomial fit to the absorption peaks and determine maximum

    for n = 1:peaks
        x_cent = mean(peak_structs(n).x);
        peak_structs(n).coeff = polyfit(peak_structs(n).x - x_cent,...
            peak_structs(n).y,2);
        poly_der = polyder(peak_structs(n).coeff);
        peak_structs(n).xmax = roots(poly_der);
        peak_structs(n).ymax = polyval(peak_structs(n).coeff,...
            peak_structs(n).xmax);
    end

%%% Evaluate concentration using Beer-Lambert Law
    input = inputdlg(['Over how many centimeters did the transmission',...
        ' occur? Type ''0'' if unsure.']);
    length = str2double(input);
    
%%% Absorption index for: 1400 nm is 8E-8, 2200 nm is 5E-7, 2750 nm
%%% is 2E06
end




