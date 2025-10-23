% Sunspot prediction result with QPGP

close all; clear; clc;

addpath('../basics');

tide_data = readmatrix('real_data_analysis/tide.csv'); 
tide=(tide_data(:,3))';
%data_matrix = reshape(tide, 12, []);
y = tide';
n=length(y);
signal=y-mean(y);

% train
period =132:156;
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

%[lower3,upper3]=bootstrap_pi(signal, par.p);
se=sqrt(pred_var(n, par.p, par.w, par.Kappa));

% --- Define datetime vector ---
startTime = datetime(2016,1,1,0,0,0);
endTime   = datetime(2016,4,9,23,50,0);
x = (startTime:minutes(10):endTime).';
signal = signal(:);
Yhat_general = Yhat_general(:);
lower_pi = Yhat_general-1.96*(se);%lower3(:);
upper_pi = Yhat_general+1.96*(se);%upper3(:);

% --- Select the valid range ---
t = (1):length(signal);

% --- Coverage Probability and Average Length ---
avg_cov_prob = mean(signal(t) >= lower_pi & signal(t) <= upper_pi);
avg_length = mean(upper_pi - lower_pi);
fprintf('Coverage probability of the interval: %.2f%%\n', avg_cov_prob*100);
fprintf('Avg prediction interval length: %.2f\n', avg_length);

% --- Match plotting vectors ---
x_plot      = x(t);
signal_plot = signal(t)+mean(y);
Yhat_plot   = Yhat_matern(t)+mean(y);
upper_plot  = upper_pi+mean(y);
lower_plot  = lower_pi+mean(y);

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
fill([xv; flipud(xv)], [yu; flipud(yl)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 1);

% --- Plot the upper and lower interval boundaries (optional, grey dotted) ---
plot(xv, yu, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
plot(xv, yl, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);

% --- Plot observed and predicted series ---
h_obs = plot(xv, signal_v, 'k', 'LineWidth', 0.6);
h_fit = plot(xv, Yhat_v, 'r--', 'LineWidth', 0.6);

% --- Axis labels and legend ---
xlabel('Year 2016');
ylabel('Water Level (m)');
legend([h_obs, h_fit], ...
    {'Observed Water Level', ...
     'Fitted QPGP with periodic Matérn kernel'}, ...
    'Location','northeast');

% --- Grid, ticks, and limits ---
grid on;
ax = gca;
ax.XAxis.TickValues = datetime(2016,1:4,1);
ax.XAxis.TickLabelFormat = 'dd MMM'; % Show only day and month
ax.XAxis.Exponent = 0;              % Avoid scientific notation
datetick('x', 'dd mmm', 'keepticks'); % Optional for older MATLAB versions

xlim([datetime(2016,1,1) xv(end)]);
ylim([0, 5]);
box on;

% --- Export as vector PDF (IEEE-ready) ---
exportgraphics(gcf, 'tide.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

