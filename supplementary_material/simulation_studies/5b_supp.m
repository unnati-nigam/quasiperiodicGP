close all;
clear;
clc;
addpath('basics');

p=100;
n=600;
sigma=1;
theta=1;
omega=0.5;
iter=1000; 

w_base=zeros(1,iter);
th_base=zeros(1,iter);
sig2_base=zeros(1,iter);
time_est_base=zeros(1,iter);

w_proposed=zeros(1,iter);
th_proposed=zeros(1,iter);
sig2_proposed=zeros(1,iter);
time_est_proposed=zeros(1,iter);

rng(1);

for j = 1:iter 
        j
        period=p;
        y=QPGPsim(n,p,theta,omega,sigma^2);
        
        %base
        tic; 
        result=par3_grid_opti(n, p, y);
        w_base(j)=result(1);
        th_base(j)=result(2);
        sig2_base(j)=result(3);
        time_est_base(j)=toc;

        tic;
       %proposed 
        % estimation time 
        par=par_est(n,p,y);
        w_proposed(j)=par.w;
        a=mackay_th_sig2_est(period,par.Kappa);
        th_proposed(j)=a(1);
        sig2_proposed(j)=a(2);
    
       time_est_proposed(j)=toc;
        

end
 

rmse_w_base=sqrt((mean(w_base)-omega)^2+var(w_base))
rmse_w_proposed=sqrt((mean(w_proposed)-omega)^2+var(w_proposed))
rmse_th_base=sqrt((mean(th_base)-theta)^2+var(th_base))
rmse_th_proposed=sqrt((mean(th_proposed)-theta)^2+var(th_proposed))
rmse_sig2_base=sqrt((mean(sig2_base)-sigma^2)^2+var(sig2_base))
rmse_sig2_proposed=sqrt((mean(sig2_proposed)-sigma^2)^2+var(sig2_proposed))


time_base=mean(time_est_base)
time_proposed=mean(time_est_proposed)

