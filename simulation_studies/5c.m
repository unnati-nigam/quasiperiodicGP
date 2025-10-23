close all;
clear;
clc;
addpath('basics');
rng(1);

p = 100;
n = 10000;
theta = 1;
sig2 = 1;
omega = 0.5;
k=n/p;

runs = 1000;
m = 1000;
param = nan(runs, 3);
ww = nan(runs, m);
thth = nan(runs, m);
sig2sig2 = nan(runs, m);

for run = 1:runs
    disp(run);
    disp(datetime);

    y = QPGPsim(n, p, theta, omega, sig2);
    par = par_est(n, p, y);

    w_EST = par.w;
    th_sig2 = mackay_th_sig2_est(par.p, par.Kappa);
    th_EST = th_sig2(1);
    sig2_EST = th_sig2(2);

    param(run, :) = [w_EST, th_EST, sig2_EST];

    % residuals
    residuals = zeros(k-1, p);
    for blk = 1:(k-1)   % block index
        y_prev = y((blk-1)*p + 1 : blk*p);
        y_curr = y(blk*p + 1 : (blk+1)*p);
        residuals(blk,:) = y_curr(:)' - par.w * y_prev(:)';
    end
    residuals = residuals - mean(residuals, 'all');

    % bootstrap
    for j = 1:m
        rng(run + j);   % seed depends on run & j
        idx = randi(k-1, k-1, 1); % sample k-1 blocks
        z_star = residuals(idx, :); % (k-1, p)

        yb = zeros(1, n);
        yb(1:p) = y(1:p);

        for blk = 1:(k-1)  % again use blk, not i
            yb((blk*p + 1): ((blk+1)*p)) = par.w * yb(((blk-1)*p+1) : (blk*p)) + z_star(blk, :);
        end

        yb = yb';
        par_b = par_est(n, p, yb);

        ww(run, j) = par_b.w;
        th_sig2 = mackay_th_sig2_est(p, par_b.Kappa);
        thth(run, j) = th_sig2(1);
        sig2sig2(run, j) = th_sig2(2);
    end
end

% Compute statistics
disp('Averages:');
disp(std(param(:,1))); %w
disp(std(param(:,2)));% theta
disp(std(param(:,3)));%sig2

disp('Quartiles:')
www=std(ww,0,2);
ththth=std(thth,0,2);
sig2sig2sig2=std(sig2sig2,0,2);

quantile(www,[0.25,0.5,0.75])
quantile(ththth,[0.25,0.5,0.75])
quantile(sig2sig2sig2,[0.25,0.5,0.75])


%disp(median(std(ww, 0, 2)));
%disp(median(std(thth, 0, 2)));
%disp(median(std(sig2sig2, 0, 2)));

 %Define the Excel file name
excelFileName = 'n10000_p100.xlsx';

% Create or overwrite the Excel file
writematrix(ww, excelFileName, 'Sheet', 'w');          % Write ww to 'w' sheet
writematrix(thth, excelFileName, 'Sheet', 'theta');    % Write thth to 'theta' sheet
writematrix(sig2sig2, excelFileName, 'Sheet', 'sig2'); % Write sig2sig2 to 'sig2' sheet
writematrix(param, excelFileName, 'Sheet', 'param'); % Write sig2sig2 to 'sig2' sheet

% Display success message
disp(['Excel file "', excelFileName, '" created successfully with 5 sheets.']);
beep; 
