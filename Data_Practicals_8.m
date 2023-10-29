clear all
close all
%
Filename = 'C:\Users\alexamin\Dropbox (UiO)\practical-data\workshop-sandbox\sasha\mcmlter-clim-caam_airt-daily-20230706.csv';
%Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 18);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["metlocid", "date_time", "avg_airt2m", "std_airt2m", "min_airt2m", "max_airt2m", "n_airt2m", "avg_airt1m", "std_airt1m", "min_airt1m", "max_airt1m", "n_airt1m", "avg_airt3m", "std_airt3m", "min_airt3m", "max_airt3m", "n_airt3m", "n_comments"];
opts.VariableTypes = ["categorical", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Specify variable properties
opts = setvaropts(opts, "metlocid", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "date_time", "InputFormat", "MM/dd/yy");
opts = setvaropts(opts, ["avg_airt3m", "std_airt3m", "min_airt3m", "max_airt3m"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["avg_airt3m", "std_airt3m", "min_airt3m", "max_airt3m"], "ThousandsSeparator", ",");
% Import the data
Data = readtable(Filename, opts);
% Convert datetime to numeric (e.g., days since the first date)
t = days(Data.date_time - Data.date_time(1));
T0 = mean(Data.avg_airt3m(~isnan(Data.avg_airt3m)));
A0 = 2*std(Data.avg_airt3m(~isnan(Data.avg_airt3m)));
t0 = -55;
t_fun = A0*sin(2*pi*(t-t0)/365)+T0;
ii = find(t>6500 & t<9000);
%
% figure(1),clf
% plot(t,Data.avg_airt3m,'.k'), hold on,
% errorbar(t,Data.avg_airt3m,Data.std_airt3m,'ok')
% plot(t,t_fun,'r')
% xlim([6000 9000])

%%
x1 = t(ii);
x2 = Data.avg_airt1m(ii);
x2std = Data.std_airt1m(ii); % Standard error

% Weights (inverse of standard errors)
weights = 1 ./ x2std;
%
A = [x1,  ones(size(x1)), sin(2*pi*(x1-t0)/365)];
%Iterative procedure
max_iter = 100; % Maximum number of iterations
tolerance = 1e-6;  % Convergence tolerance
prev_m = zeros(3,1);  %Initialize previous coefficients to zeros
%
for it = 1:max_iter
    it
    A_weighted = sqrt(weights(:)) .* A; % weighted system matrix A 
    x2_weighted = sqrt(weights(:)) .* x2(:); % weighted data vector   
    % solve normal equations
    m = (A_weighted'*A_weighted)\ A_weighted'*x2_weighted ;
    %Check for convergence
    if all(abs(m - prev_m) < tolerance)
        disp('Converged.')
        break
    end
    %Update weights based on new residuals
    residuals = x2 - A*m;
    weights = 1 ./ (residuals.^2 + 1e-8); %least squares norm
    %Update previous coefficients
    prev_m = m;
end
% predict
p1 = linspace(6000,9000,100);
xtrend = p1*m(1)+m(2);
xp = xtrend+m(3)*sin(2*pi*(p1-t0)/365);
% display
figure(2),clf
plot(x1, x2, 'ok'), hold on
errorbar(x1, x2, x2std, '.k')
plot(p1,xp,'r')
plot(p1,xtrend,'b--','LineWidth',2)
title(sprintf('Estimated Model: x_2 = %6.4fx_1 + %6.4f + %6.4f sin(w x_1)\n', m(1), m(2), m(3)))
xlabel('x_1'), ylabel('x_2')

