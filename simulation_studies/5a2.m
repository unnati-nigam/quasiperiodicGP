clear; close all; clc;
addpath('routine', 'basics')

n=10000; 
p=10; w=0.5;theta=1;sig2=1;
Kappa=sig2 * toeplitz(exp(-theta^2 * (sin(pi * (0:(p-1)) / p)).^2));

rng(1);
y=QPGPsim(n, p, theta, w, sig2);

%proposed
par=par_est(n,p,y);
a=mackay_th_sig2_est(p, par.Kappa);
th_est=a(1); sig2_est=a(2);
Kappa_est=sig2_est *toeplitz( exp(-th_est^2 * (sin(pi * (0:(par.p-1)) / par.p)).^2));
tic;
pred=pred_element(n, par.p, par.w, Kappa_est, y);
yhat_me=pred';
time_me=toc
ipse_me=mean((y(2:n)-yhat_me(2:n)).^2)


%base 
tic;
yhat_base=actual_prediction(n, p, w, Kappa, y);
time_base=toc
ipse_base=mean((y(2:n)-yhat_base(2:n)).^2)

