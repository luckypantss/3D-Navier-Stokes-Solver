clear 
clc
close all
isPlot = true;      % Turn to false to disable plots outside of simulation
%% Physical properties:
% For room temperature at 20deg celsius. Data taken from "Fysika".
rho     = 997;              % Density                   [kg/m^3]
mu      = 1057e-6;          % Dynamic Viscocity         [Pa*s]
nu      = mu/rho;           % Kinematic Viscocity       [m^2/s]
F_z     = - 9.823;          % Gravitatonal force        [m/s^2]

%% DISC. the Domain
L_x = 1;
L_y = 1;
L_z = 1;

n_x = 100;
n_y = 100;
n_z = 100; 

[X, Y, Z, dx, dy, dz] = SpatialDisc(L_x, L_y, L_z, n_x, n_y, n_z,isPlot);

%% Initiating the solver 
t       = 20;                                 % Simulation time           [s]
D_h     = 4 * L_y*L_x/(2*(L_x+L_y));          % Hydraluic diameter        [m]
beta    = 50;          
dt      = 0.01;

tic
NSSolver(X, Y, Z, dx, dy, dz, t, dt, rho, nu,beta, F_z, D_h, isPlot);
toc


%%
[X,Y,Z] = meshgrid(linspace(0,99,10));

a=NSSolverTest(X, Y, Z, 1, 1, 1, 1000, 1, 1, 0.05 ,0.0125, 0, 1, 1);
% NSSolver(X, Y, Z, dx, dy, dz, t, dt, rho, nu, beta, F_z, D_h, isPlot)


