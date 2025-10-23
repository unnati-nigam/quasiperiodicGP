close all; clear; clc;

addpath('../basics');

tide_data = readmatrix('real_data_analysis/tide.csv'); 
tide=(tide_data(:,3))';
%data_matrix = reshape(tide, 12, []);
y = tide';
start=1;
n=length(y);
signal=y;

% train
period =132:156;
par=par_est(n,period,signal)

%[SE_w, SE_Kappa, w_upper, w_lower, Kappa_upper, Kappa_lower] = bootstrap_par_se(y, par.p)
[w_boot, Kappa_general_boot, Kappa_mackay_boot,Kappa_matern_boot, Kappa_cosine_boot]=bootstrap_par_se(signal, par.p);

% Example: A is m x n (rows = samples, columns = variables)
lower_w = quantile(w_boot, 0.025)   % 2.5% quantile (row vector of size 1 x n)
upper_w = quantile(w_boot, 0.975)   % 97.5% quantile (row vector of size 1 x n)
w_se = std(w_boot)

lower_Kappa_general = squeeze(quantile(Kappa_general_boot(1,:,:), 0.025, 3))'  
upper_Kappa_general = squeeze(quantile(Kappa_general_boot(1,:,:), 0.975, 3))'
Kappa_general_se = squeeze(std(Kappa_general_boot(1,:,:),0,3))'

b=mackay_th_sig2_est(par.p,par.Kappa)
lower_Kappa_mackay = (quantile(Kappa_mackay_boot, 0.025, 2))'  
upper_Kappa_mackay = (quantile(Kappa_mackay_boot, 0.975, 2))'  
Kappa_se_mackay = std(Kappa_mackay_boot,0,2)'

a=matern_sig2_nu_ell_est(par.p,par.Kappa)
lower_Kappa_matern = (quantile(Kappa_matern_boot, 0.025, 2))'  
upper_Kappa_matern = (quantile(Kappa_matern_boot, 0.975, 2))'  
Kappa_se_matern = std(Kappa_matern_boot,0,2)'

c=cosine_sig2_est(par.p,par.Kappa)
lower_Kappa_cosine = (quantile(Kappa_cosine_boot, 0.025, 2))'  
upper_Kappa_cosine = (quantile(Kappa_cosine_boot, 0.975, 2))'  
Kappa_se_cosine = std(Kappa_cosine_boot,0,2)'



% ----------------- Prepare x-axis -----------------
x = 0:(par.p)/2;       % lag values for x-axis
idx = 1:length(x);       % MATLAB-safe indices for array/matrix access

% ----------------- Create Figure -----------------
figure('Units','inches','Position',[1 1 3.5 2.5]); hold on; % IEEE single-column

% ----------------- Plot Bootstrap Quantiles -----------------
h_bounds_upper = plot(x, upper_Kappa_general(idx), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
h_bounds_lower = plot(x, lower_Kappa_general(idx), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);

% ----------------- Plot Estimated Kappa -----------------
h_est = plot(x, par.Kappa(1, idx), 'k-', 'LineWidth', 1.5);

% ----------------- Axis Labels -----------------
xlabel('Lag','FontSize',7,'FontWeight','bold');
ylabel('Estimate of generasl $\kappa_p$','FontSize',10,'FontWeight','bold','Interpreter','latex');

% ----------------- Legend -----------------
lgd = legend(h_est, {'Estimated general $\kappa_p$ for Sunspot dataset'}, ...
             'Location','best','FontSize',6,'Box','on','Interpreter','latex');
lgd.FontName = 'Times';

% ----------------- Axis Formatting -----------------
grid on; box on;
set(gca,'FontSize',6,'LineWidth',0.75);

% Set x-axis limits and ticks exactly at x values
xlim([min(x), max(x)]);

% Optionally remove extra white space around data
% axis tight;  % Not needed since xlim already matches data

% ----------------- Export Figure -----------------
exportgraphics(gca, 'tide_general_kappa.pdf', 'ContentType','vector','BackgroundColor','None');



% ----------------- Create Figure -----------------

figure('Units','inches','Position',[1 1 3.5 2.5]); hold on; % IEEE single-column


% Plot estimated Kappa SE in red solid xxxline
h_est = plot(x, Kappa_general_se(idx), 'k-', 'LineWidth', 1.5);

% ----------------- Axis Labels ----------------
xlabel('Lag','FontSize',7,'FontWeight','bold');
ylabel('Standard Error of general $\kappa_p$','FontSize',7,'FontWeight','bold','Interpreter','latex');

% ----------------- Legend -----------------
lgd = legend(h_est, 'SE of \kappa_p for Sunspot dataset', ...
             'Location','best','Box','on');
% Force font to match axes
lgd.FontName = 'Times';        % same as axes
lgd.FontSize = 6;              % match axes font size
lgd.FontWeight = 'normal';     % optional: bold/normal
lgd.TextColor = [0 0 0];  

% ----------------- Axis Formatting -----------------
grid on; box on;
set(gca,'FontSize',6,'LineWidth',0.75);
axis tight;  % remove extra white space around data

% ----------------- Export Figure -----------------
exportgraphics(gca, 'tide_kappa_se.pdf', 'ContentType','vector','BackgroundColor','none');



 