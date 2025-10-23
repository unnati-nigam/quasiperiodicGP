close all;
clear;
clc;
addpath('basics');
rng(1);

p = 10;
n = 2500;
theta = 1;
sig2 = 1;
omega = 0.5;
k=n/p;

runs = 1000;
m = 100;
param = nan(runs, 3);
ww = nan(runs, m);
thth = nan(runs, m);
sig2sig2 = nan(runs, m);

for i = 1:runs

    y = QPGPsim(n, p, theta, omega, sig2);
    par = par_est(k, p, y);
    w_EST = par.w;
    Kappa=sig2 * toeplitz(exp(-theta^2 * (sin(pi * (0:(p-1)) / p)).^2));
    Kappa_EST = par.Kappa;
    th_sig2 = mackay_th_sig2_est(par.p, par.Kappa);
    th_EST = th_sig2(1);
    sig2_EST = th_sig2(2);
    param(i, :) = [w_EST, th_EST, sig2_EST];
    x=nan(m,1);y=nan(m,1);
    for j = 1:m
        %j
        rng(i+j);
        Y_star=QPGPsim(k*p, p,th_EST,w_EST, sig2_EST);
        w_Kappa=w_Kappa_est(k, p, Y_star);
        ww(i, j) = w_Kappa.w;
        th_sig2 = mackay_th_sig2_est(p, w_Kappa.Kappa);
        thth(i, j) = th_sig2(1);
        sig2sig2(i, j) = th_sig2(2);
    end
    disp(i);
    disp(datetime);
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

k

%disp(median(std(ww, 0, 2)));
%disp(median(std(thth, 0, 2)));
%disp(median(std(sig2sig2, 0, 2)));

 %Define the Excel file name
%excelFileName = 'k500_final.xlsx';

% Create or overwrite the Excel file
%writematrix(ww, excelFileName, 'Sheet', 'w');          % Write ww to 'w' sheet
%writematrix(thth, excelFileName, 'Sheet', 'theta');    % Write thth to 'theta' sheet
%writematrix(sig2sig2, excelFileName, 'Sheet', 'sig2'); % Write sig2sig2 to 'sig2' sheet
%writematrix(param, excelFileName, 'Sheet', 'param'); % Write sig2sig2 to 'sig2' sheet

% Display success message
%disp(['Excel file "', excelFileName, '" created successfully with 5 sheets.']);
beep; 
