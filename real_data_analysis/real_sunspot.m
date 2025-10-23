close all; clear; clc;

addpath('../basics');
warning('off', 'all')   % turn off all warnings

sunspot = readmatrix('real_data_analysis/sunspot_year.csv', ...
                     'Delimiter', ';', 'CommentStyle', '#');
%years = sunspot(:,1);
%months = sunspot(:,2);
%x=datetime(years,months, ones(size(years)));
x = sunspot(4:end,1);
y = sunspot(4:end,2);

%y = sunspot(10:end,2); signal=y;
%load sunspot
%x = sunspot(:,1);
%y = sunspot(:,2);

% train
n=length(y);
signal=y;%-mean(y);

period=2:15;
par=par_est(n,period,signal)

Yhat_general = (pred_element(n, par.p, par.w, par.Kappa, signal))';
rmse_general=sqrt(mean((signal((par.p+1):n)-Yhat_general((par.p+1):n)).^2))

a=matern_sig2_nu_ell_est(par.p,par.Kappa);
sig2_est=a(1); nu_est=a(2); ell_est=a(3);
Kappa_matern=sig2_est *toeplitz(matern_cov(0:(par.p-1), nu_est, ell_est, par.p));
Kappa_matern(logical(eye(size(Kappa_matern)))) = sig2_est;
pred=pred_element(n, par.p, par.w, Kappa_matern, signal);
Yhat_matern=pred';
rmse_matern=sqrt(mean((signal((par.p+1):n)-Yhat_matern((par.p+1):n)).^2))

b=mackay_th_sig2_est(par.p,par.Kappa);
sig2_est_mackay=b(2);th_est_mackay=b(1);
Kappa_mackay=sig2_est_mackay*toeplitz( exp(-th_est_mackay^2 * (sin(pi * (0:(par.p-1)) / par.p)).^2));
pred=pred_element(n, par.p, par.w, Kappa_mackay, signal);
Yhat_mackay=pred';
rmse_mackay=sqrt(mean((signal((par.p+1):n)-Yhat_mackay((par.p+1):n)).^2))

c=cosine_sig2_est(par.p,par.Kappa);
sig2_est_cosine=c;
Kappa_cosine=sig2_est_cosine*toeplitz(cos(2*pi * (0:(par.p-1)) / par.p));
pred=pred_element(n, par.p, par.w, Kappa_cosine, signal);
Yhat_cosine=pred';
rmse_cosine=sqrt(mean((signal((par.p+1):n)-Yhat_cosine((par.p+1):n)).^2))

se=sqrt(pred_var(n, par.p, par.w, par.Kappa)); %bootstrap_pi(signal, par.p);

x = x(:);            
%signal = signal(:)+mean(y);  
%Yhat_general = Yhat_general(:)+mean(y);    

% Average coverage probability
%avg_cov_prob = mean(signal(t) >= lower_pi & signal(t) <= upper_pi);

% Average interval length
%avg_length = mean(upper_pi - lower_pi);

%fprintf('Coverage probability of the interval: %.2f%%\n', avg_cov_prob*100);
%fprintf('Avg prediction interval length: %.2f\n', avg_length);



lower_pi = Yhat_general-1.96*(se);%lower3(:);
upper_pi = Yhat_general+1.96*(se);%upper3(:);
t = 1:length(signal);
upper_pi = upper_pi(:)'; 
lower_pi = lower_pi(:)';



x_plot = x(1:end);             % length matches upper_pi etc.
signal_plot = signal(1:end);
Yhat_plot = Yhat_matern(1:end);
upper_plot = upper_pi;
lower_plot = lower_pi;

xv = x_plot(:);
yu = upper_plot(:);
yl = lower_plot(:);
signal_v = signal_plot(:);
Yhat_v = Yhat_plot(:);


% --- IEEE double-column figure setup ---
figure('Units','inches','Position',[1 1 7.2 3.0]); hold on;
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 7.2 3.0]);
set(gca, 'FontSize', 8);

% --- Fill the region between lower and upper (light grey) ---
fill([xv; flipud(xv)], [yu; flipud(yl)], [0.9 0.9 0.9], ...
    'EdgeColor', 'none', 'FaceAlpha', 1);

% --- Plot the upper and lower interval boundaries (optional, grey dotted) ---
plot(xv, yu, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.7);
plot(xv, yl, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.7);

% --- Plot observed and predicted series ---
h_obs = plot(xv, signal_v, 'k', 'LineWidth', 0.8);
h_fit = plot(xv, Yhat_v, 'r--', 'LineWidth', 0.8);

xlabel('Year'); ylabel('Number of Sunspots');


legend([h_obs, h_fit], {'Observed Sunspot Numbers','Fitted QPGP with periodic Maternal kernel'}, ...
       'Location','northeast');

ylim([-200,550]);
xlim([x_plot(1) x_plot(end)]);
grid on;

% Number of ticks
num_ticks = 5;

% Define a fraction of the series to skip at start and end
skip_frac_start = 0.05;  % skip first 5% of observations
skip_frac_end   = 0.05;  % skip last 5% of observations

n = length(x_plot);

% Compute indices for ticks, leaving margin at start and end
start_idx = round(1 + n*skip_frac_start);
end_idx   = round(n - n*skip_frac_end);

% Generate 5 evenly spaced indices between start_idx and end_idx
tick_indices = round(linspace(start_idx, end_idx, num_ticks));

% Tick positions in x-axis units
tick_positions = xv(tick_indices);

% Integer labels (remove decimals)
tick_labels = floor(tick_positions);

% Set ticks and labels
xticks(tick_positions);
xticklabels(arrayfun(@(x) sprintf('%d', x), tick_labels, 'UniformOutput', false));

% Ensure integer formatting
ax = gca;
ax.XAxis.TickLabelFormat = '%.0f';
box on;
exportgraphics(gcf, 'sunspot.pdf', 'ContentType','vector','BackgroundColor','None');
