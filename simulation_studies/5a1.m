clear; close all; clc
addpath('routine', 'basics')

n=10000; 
p=10; omega=0.5;theta=1;sig2=1;
Kappa=sig2 * toeplitz(exp(-theta^2 * (sin(pi * (0:(p-1)) / p)).^2));

rng(100);
y=QPGPsim(n, p, theta, omega, sig2);

tic;
proposed=logL_exact(n, p, omega, Kappa, y)
time_proposed=toc;

tic;
base=log_L_TSP(n, [0 theta omega], y,@period_sin_gauss_cov,@regpoly0_zeromean,p)
time_base=toc;

time_base
time_proposed
