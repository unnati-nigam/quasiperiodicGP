% Simulation study

close all;
clear;
clc;
p=10;
sigma=1;
delta=0;
theta=5;
omega=0.5;
addpath('routine', 'routine/fastF0Nls');
datetime('now')
run=4;
iter=100;
%w=zeros(run,iter);
%pp=zeros(run,iter);
%th=zeros(run,iter);
%sig=zeros(run,iter);
time_est=zeros(run,iter);
time_pred=zeros(run,iter);

for i =1:run
    N=10*(10.^i);
    for j = 1:iter                                                               
        rng(i+j);
        period=10;
        lob = [0 0.1 -0.9999];      upb = [0 15 0.9999];
        y=randQPGP(N,p,sigma,delta,theta,omega,@period_sin_gauss_cov);
        tic;
        % estimation time 
        model=fit_QPGP(period,y,@regpoly0,@period_sin_gauss_cov,lob,upb);
        time_est(i,j)=toc;

        tic;
        y_pred=pred_QPGP(model,y);
        time_pred(i,j)=toc;
        %w(i,j)=model.thetahat(3);
        %pp(i,j)=model.period;
        %sig(i,j)=(model.sigmahat);
        %th(i,j)=model.thetahat(2);
    end
    i
   % datetime('now')
end

mean_est_time=(mean(time_est,2))'
mean_pred_time=(mean(time_pred,2))'


%me=(mean(w,2))'
%sd=(std(w,0,2))'
%me=(mean(th,2))'
%sd=(std(th,0,2))'
%me=(mean(sig,2))'
%sd=(std(sig,0,2))'
