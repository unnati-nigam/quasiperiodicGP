close all; clear; clc;

addpath('basics');
warning('off', 'all')

data = readmatrix('real_data_analysis/co2.xlsx'); % Assumes numerical data (excluding first column)
data_vector = reshape(data(:, 2:end)', [], 1);
y = data_vector;
y(y < 0) = NaN;
y =replace_na_with_average(y);
n = length(y);
start_date = datetime(1958,1,1); 
x = start_date + calmonths(0:n-1);
xx = (1:n)';                 

X = [ones(n,1), xx, xx.^2];  
b = X \ y;                  

y_fitted = X * b;          
signal = y - y_fitted;   

period = (2:1:20)';
par=par_est(n,period,signal)

Yhat_general = (pred_element(n, par.p, par.w, par.Kappa, signal))';
rmse_general=sqrt(mean((signal((par.p+1):n)-Yhat_general((par.p+1):n)).^2))

a=matern_sig2_nu_ell_est(par.p,par.Kappa);
sig2_est_matern=a(1);nu_est_matern=a(2);ell_est_matern=a(3);
Kappa_matern=sig2_est_matern *toeplitz( matern_cov(0:(par.p-1), nu_est_matern, ell_est_matern, par.p));
Kappa_matern(logical(eye(size(Kappa_matern)))) = sig2_est_matern;
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

se=sqrt(pred_var(n, par.p, par.w, par.Kappa));

x = x(:);            
signal = signal(:);  
Yhat_general = Yhat_general(:);    
%lower_pi = lower3(:);
%upper_pi = upper3(:);
t = (1):length(signal);

% Average coverage probability
%avg_cov_prob = mean(signal(t) >= lower_pi & signal(t) <= upper_pi);

% Average interval length
%avg_length = mean(upper_pi - lower_pi);

%fprintf('Coverage probability of the interval: %.2f%%\n', avg_cov_prob*100);
%fprintf('Avg prediction interval length: %.2f\n', avg_length);

lower_pi = Yhat_general-1.96*(se);%lower3(:);
upper_pi = Yhat_general+1.96*(se);%upper3(:);

% --- Match plotting vectors ---
x_plot      = x(t);
signal_plot = signal(t);
Yhat_plot   = Yhat_general(t);
upper_plot  = upper_pi;
lower_plot  = lower_pi;

% --- Ensure column vectors ---
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
plot(xv, yu, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7);
plot(xv, yl, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7);

% --- Plot observed and predicted series ---
h_obs = plot(xv, signal_v, 'k', 'LineWidth', 0.8);
h_fit = plot(xv, Yhat_v, 'r--', 'LineWidth', 0.8);

% --- Axis labels and legend ---
xlabel('Year');
ylabel('Carbon Dioxide Emission (ppm)');
legend([h_obs, h_fit], ...
    {'Observed CO2 Emission', ...
     'Fitted QPGP with general kernel'}, ...
    'Location','northeast');

% --- Grid, ticks, and limits ---
grid on;
ax = gca;
ax.XAxis.TickLabelFormat = 'yyyy';

% --- Set 5 ticks evenly spaced along xv ---
ax.XTick = linspace(xv(par.p+1), xv(end), 5);

ylim([-10, 15]);
box on;


% --- Export as vector PDF (IEEE-ready) ---
exportgraphics(gcf, 'co2.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none');




function x = replace_na_with_average(x)
    for i = 1:length(x)
        if isnan(x(i))
            if i == 1
                x(i) = x(i + 1); % Replace first NaN with the next value
            elseif i == length(x)
                x(i) = x(i - 1); % Replace last NaN with the previous value
            else
                x(i) = mean([x(i - 1), x(i + 1)], 'omitnan'); % Average neighbors
            end
        end
    end
end