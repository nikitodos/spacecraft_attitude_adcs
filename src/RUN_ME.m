%% Spacecraft Attitude Dynamics - Group 34 / Project 145
% 6U cubesat by
%
% Stefano Diambri
% Ludovico Drioli
% Michele Pellizzer
% Giovanni Nicola D'Aloisio
clear; clc; close all;

%% Sims data
%
% P L E A S E   N O T E  
% In the report the number of Simulation was set to 30. 
% For the convinience of prof. Bernelli, it is reduced to 6
%
simulation_total = 6; 

% Arrays for Global Statistic 
error_list = [];
saturation_data = [];
control_torque_data = [];

%% Parameter

% other parameter
R_E=astroConstants(23);
T_E=astroConstants(53)*3600;
w_E=2*pi/T_E;
epsilon=astroConstants(8);
mu=astroConstants(13);
r_sun=astroConstants(2);
c=astroConstants(5)*10^3;
AU=astroConstants(2);
R_sun=astroConstants(3);
Fe=1358;
alpha_G0=0;

% orbit
a=3000+R_E;
e=0.2;
i=deg2rad(30);
th0=deg2rad(0);
n=sqrt(mu/a^3);
T=2*pi/n;

% Principal inertia axis of the spacecraft
z_b = 0.10;  % [m]
y_b = 0.20;  % [m]
x_b = 0.30;  % [m]

m_b = 12;    % [kg]

% Cuboid formula
Ix = (m_b/12)*((y_b^2)+(x_b^2));
Iy = (m_b/12)*((z_b^2)+(x_b^2));
Iz = (m_b/12)*((y_b^2)+(z_b^2)); % this is intended direction of angular rate

I=[Ix 0 0; 0 Iy 0; 0 0 Iz];

% magnetic field data
j_B=[0.01; 0.05; 0.01];

g_IGRF = [
 -29350.0 -1410.3     0       0       0       0       0       0       0       0       0       0       0       0;
 -2556.2   2950.9   1648.7    0       0       0       0       0       0       0       0       0       0       0;
 1360.9  -2404.2   1243.8   453.4    0       0       0       0       0       0       0       0       0       0;
 894.7     799.6     55.8  -281.1    12.0    0       0       0       0       0       0       0       0       0;
 -232.9    369.0    187.2  -138.7  -141.9   20.9    0       0       0       0       0       0       0       0;
 64.3      63.8     76.7  -115.7   -40.9   14.9  -60.8    0       0       0       0       0       0       0;
 79.6     -76.9     -8.8    59.3    15.8    2.5  -11.2   14.3    0       0       0       0       0       0;
 23.1      10.9    -17.5     2.0   -21.8   16.9   14.9  -16.8    1.0    0       0       0       0       0;
 4.7        8.0      3.0    -0.2    -2.5  -13.1    2.4    8.6   -8.7  -12.8    0       0       0       0;
 -1.3      -6.4      0.2     2.0    -1.0   -0.5   -0.9    1.5    0.9   -2.6   -3.9    0       0       0;
 3.0       -1.4     -2.5     2.4    -0.6    0.0   -0.6   -0.1    1.1   -1.0   -0.1    2.6    0       0;
 -2.0      -0.1      0.4     1.2    -1.2    0.6    0.5    0.5   -0.1   -0.5   -0.2   -1.2   -0.7    0;
 0.2       -0.9      0.6     0.7    -0.2    0.5    0.1    0.7    0.0    0.3    0.2    0.4   -0.5   -0.4
];

h_IGRF = [
 0      4545.5       0       0       0       0       0       0       0       0       0       0       0       0;
 0     -3133.6    -814.2     0       0       0       0       0       0       0       0       0       0       0;
 0      -56.9      237.6   -549.6    0       0       0       0       0       0       0       0       0       0;
 0      278.6     -134.0    212.0  -375.4    0       0       0       0       0       0       0       0       0;
 0       45.3      220.0   -122.9    42.9   106.2    0       0       0       0       0       0       0       0;
 0      -18.4       16.8     48.9   -59.8    10.9    72.8    0       0       0       0       0       0       0;
 0      -48.9      -14.4     -1.0    23.5    -7.4   -25.1   -2.2    0       0       0       0       0       0;
 0        7.2      -12.6     11.5    -9.7    12.7     0.7    -5.2    3.9    0       0       0       0       0;
 0      -24.8       12.1      8.3    -3.4    -5.3     7.2    -0.6    0.8    9.8    0       0       0       0;
 0        3.3        0.1      2.5     5.4    -9.0     0.4    -4.2    -3.8    0.9    -9.0    0       0       0;
 0        0.0        2.8     -0.6     0.1     0.5    -0.3    -1.2    -1.7    -2.9    -1.8    -2.3    0       0;
 0       -1.2        0.6      1.0    -1.5     0.0     0.6    -0.2     0.8     0.1    -0.9     0.1     0.2    0;
 0       -0.9        0.7      1.2    -0.3    -1.3    -0.1     0.2    -0.2     0.5     0.6    -0.6    -0.3   -0.5
];

% Processor timestep
Ts = 0.1;

% Non linear Control gains
k1_tilde=5e-4; % De tumbling
k1=20e-4; % Re pointing (angular rate)
k2=2e-4; % Re pointing (attitude)


%% PD - LINEAR CONTROLLER

% Linear controller setting
settling_time=50;
damping=0.7;

% Gravitational parameters
Kp = (Iy-Ix)/Iz;
Ky = (Iz-Iy)/Ix;
Kr = (Iz-Ix)/Iy;

% State Space formulation (Gravity Gradient only)
A11 = [0 (1-Ky)*n 0; (Kr-1)*n 0 0; 0 0 0];
A12 = -n^2*diag([Ky; 4*Kr; 3*Kp]);
A21 = diag([1; 1; 1]);
A22 = zeros(3);

A = [A11 A12; A21 A22];

B1 = diag([1/Ix; 1/Iy; 1/Iz]);
B2 = zeros(3);

B = [B1; B2];

C1 = zeros(3);
C2 = diag([1 1 1]);
C = [C1 C2];

D = zeros(3);

a_p=eig(A);

% Poles exclusion zone on the complex plane
sigma=5/settling_time;
w_n=sigma/damping;

R_restriction = w_n;
h_restriction = sigma;
th_restriction = acos(damping);

% K computing (matrix of gains)

k = ComputeMinK(R_restriction, h_restriction, th_restriction, a_p);

k=min(k);
a_controlled=a_p+3*k*ones(length(a_p),1); % pole placing

K = place(A, -B, a_controlled);

xx = linspace(-3,0,100) ;
rr = linspace(-R_restriction,0,1000);

% plotting
figure
grid on
hold on
xline(-h_restriction, "LineWidth",2, "LineStyle", "--")
plot(xx, tan(th_restriction)*xx, "LineWidth",2, "LineStyle", ":")
plot(xx, -tan(th_restriction)*xx,  "LineWidth",2, "LineStyle", ":")
plot(rr, sqrt(R_restriction^2-rr.^2), "LineWidth",2, "LineStyle", "-.")
plot(rr, -sqrt(R_restriction^2-rr.^2), "LineWidth",2, "LineStyle", "-.")
plot(real(a_controlled), imag(a_controlled),"o", "LineWidth", 3)
% fontsize(36, "points")

%% Stability Margin and Phase Margin Computation
controlled_sys = ss(A+B*K, B, C+D*K, D);

figure
bodeplot(controlled_sys)
grid on
% fontsize(36, "points")


figure
nyquist(controlled_sys)
% fontsize(36, "points")

S = allmargin(controlled_sys);
phase_margin = S.PhaseMargin;
gain_margin = S.GainMargin;
disp("phase margin = " + phase_margin + "°")
disp("gain margin = " + gain_margin)


%% SIMULATIONS : ITERATIVE SECTION
% run the simulation N time. The cycle includes:
% - Random variables (sensors, initial conditions)
% - Oustaloup approx (for need of simulation)
% - Extraction of data for Statistical analysis
% - Complete plotting of the LAST run

for simulation_number = 1:simulation_total
    disp("Simulation n° " + simulation_number + "/" + simulation_total)

% inital angular velocity
wx_0=rand(1)*pi/6;
wy_0=rand(1)*pi/6;
wz_0=rand(1)*pi/6; % this is intended direction of angular rate

w_0=[wx_0; wy_0; wz_0];

% inital DCM
epsilon1=rand(1)*pi/6;
epsilon2=rand(1)*pi/6;

A_0=[cos(th0+epsilon1) sin(th0+epsilon1) 0; -sin(th0+epsilon1) cos(th0+epsilon1) 0; 0 0 1]*...
   [1 0 0; 0 cos(i+epsilon2) sin(i+epsilon2); 0 -sin(i+epsilon2) cos(i+epsilon2)];

% sun sensor data
sigma_sun=deg2rad(0.5);
b_sun=deg2rad(0.1*rand(1));
Ts_sun=0.1;
Mis_sun=deg2rad(0.1*rand(3,1));

% magnetometer data
noise_density_mag=16e-9;
precision_mag=8e-9;
b_mag=15e-9;
saturation_mag=60e-6;
Ts_mag=0.1; %1/25;
Mis_mag=deg2rad(0.1*rand(3,1));
Nonorthog_err_mag=deg2rad(1)*ones(3,1);

% gyro data % CRH02-100
Mis_gyro=deg2rad(0.1*rand(3,1));
Nonorthog_err_gyro=deg2rad(0.1*rand(3,1));
gyro_mounting = eye(3);
gyro_delay = 1;
b_gyro = deg2rad(0.01*rand(1));
Ts_gyro = 0.1;
SF=deg2rad(20*5.25e-3);
SFn=0.02/100;
ARW=deg2rad(0.017/sqrt(3600));
sigma_ARW=ARW/sqrt(Ts_gyro);
RRW=0;
sigma_RRW=RRW/sqrt(Ts_gyro);
bias_instability=deg2rad(0.12/3600);
sigma_bias=bias_instability/sqrt(Ts_gyro);
saturation_gyro = deg2rad(100);
c_time = 100;

% attitude data
alpha1=0.72;

% Actuators
RW.A = eye(3);
RW.A_inv = eye(3);
RW.omega_act_0 = 0;
RW.max_torque = 0.0015; % N*m
RW.max_torque_start_descending = 3500*2*pi/60; % rad/s (conversion from rpm)
RW.max_speed = 6000*2*pi/60; % rad/s
RW.max_torque_descending_rate = -RW.max_torque/(RW.max_speed - RW.max_torque_start_descending);
RW.min_torque = 0.000015; % N*m
RW.I = 0.3; % kg*m^2
RW.mis=deg2rad(0.1*rand(3,1));

%% Oustaloup approximation for 1/sqrt(s)

N=8;
v=zeros(N,1);
mu_o=zeros(N,1);
alpha=0.5;
w_l=0.001;
w_h=1000;
syms s

gamma=(w_h/w_l)^(alpha/N);
eta=(w_h/w_l)^((1-alpha)/N);
mu_o(1)=w_l*sqrt(eta);
G=1;

for index=1:N-1
    v(index)=mu_o(index)*gamma;
    mu_o(index+1)=v(index)*eta;
    G=G.*(1+s./v(index))./(1+s./mu_o(index));
end

v(N)=mu_o(N)*gamma;
G=G*(1+s/v(N))./(1+s/mu_o(N));
gain=subs(G,s,0);
[num, den] = numden(G);

num8 = expand(num);
den8 = expand(den);

num8_coeff = sym2poly(num8);
den8_coeff = sym2poly(den8);

FdT=tf([num8_coeff(1) num8_coeff(2) num8_coeff(3) num8_coeff(4)...
    num8_coeff(5) num8_coeff(6) num8_coeff(7) num8_coeff(8) num8_coeff(9)],...
    [den8_coeff(1) den8_coeff(2) den8_coeff(3) den8_coeff(4)...
    den8_coeff(5) den8_coeff(6) den8_coeff(7) den8_coeff(8) den8_coeff(9)]);

if simulation_number== simulation_total
    figure
    bode(FdT)
    grid on
    % fontsize(36, "points")
end

syms z
G_z_sun=subs(G,s,2/Ts_sun*(z-1)/(z+1));
G_z_sun=simplify(G_z_sun);

[num_z_sun, den_z_sun] = numden(G_z_sun);

num8_z_sun = expand(num_z_sun);
den8_z_sun = expand(den_z_sun);

num_cz_sun = double(sym2poly(num8_z_sun));
den_cz_sun = double(sym2poly(den8_z_sun));

G_z_mag=subs(G,s,2/Ts_mag*(z-1)/(z+1));
G_z_mag=simplify(G_z_mag);

[num_z_mag, den_z_mag] = numden(G_z_mag);

num8_z_mag = expand(num_z_mag);
den8_z_mag = expand(den_z_mag);

num_cz_mag = double(sym2poly(num8_z_mag));
den_cz_mag = double(sym2poly(den8_z_mag));

%% Call the Simulink and returns data
simout = sim("Project_sim.slx");
sim_time = simout.tout;
sim_angular_rate = simout.angular_rate;
sim_error = simout.control_error();
sim_switch2nadir = simout.switch_to_nadir;
sim_switch2repointing = simout.switch_to_repointing;
sim_eclipse = simout.eclipse_pin;
sim_act_torque = simout.actuator_torque;

%% Plot the last run completely
if(simulation_number == simulation_total)

    idx_start_repointing = find(sim_switch2repointing>=1, 1);
    idx_start_nadir_pointing = find(sim_switch2nadir>=2, 1);
    
    figure
    hold on
    yyaxis left
    plot(sim_time, rad2deg(sim_error))
    xline(idx_start_repointing/10, "--", "Start Repoiting")
    xline(idx_start_nadir_pointing/10, "--", "Start Nadir Pointing")
    xlabel("Timestamp")
    ylabel("1-axis Error [°]")
    yyaxis right
    plot(sim_time, sim_eclipse)
    ylabel("Eclipse")
    title("Error of pointing (1 simulation)")
    legend("Error", "Eclipse")
    grid on
    % fontsize(36, "points")

    flat_sim_act_torque = reshape(sim_act_torque, [size(sim_act_torque,3),3]);
    figure
    hold on
    yyaxis left
    plot(sim_time, flat_sim_act_torque(:,1) , "r" )
    plot(sim_time, flat_sim_act_torque(:,2) , "b" )
    plot(sim_time, flat_sim_act_torque(:,3) , "black" )
    xline(idx_start_repointing/10, "--", "Start Repoiting")
    xline(idx_start_nadir_pointing/10, "--", "Start Nadir Pointing")
    xlabel("Timestamp")
    ylabel("Actuator torque [Nm]")
    yyaxis right
    plot(sim_time, sim_eclipse)
    ylabel("Eclipse")
    title("Actuator torque 3-axis (1 simulation)")
    legend("RW-x", "RW-y", "RW-z", "Eclipse")
    grid on
    % fontsize(36, "points")
    
    figure
    hold on
    yyaxis left
    plot(sim_time, rad2deg( sim_angular_rate(:,1)' ), "r")
    plot(sim_time, rad2deg( sim_angular_rate(:,2)' ), "b")
    plot(sim_time, rad2deg( sim_angular_rate(:,3)' ), "black")
    xline(idx_start_repointing/10, "--", "Start Repoiting")
    xline(idx_start_nadir_pointing/10, "--", "Start Nadir Pointing")
    xlabel("Timestamp")
    ylabel("Angular Rate [°]/s]")
    yyaxis right
    plot(sim_time, sim_eclipse)
    ylabel("Eclipse")
    title("Angualar rate (1 simulation)")
    legend("wx", "wy", "wz", "Eclipse")
    grid on
    % fontsize(36, "points")

end


%% Saves data used in the statistical analysis of linear controller

idx_start_repointing = find(sim_switch2repointing>=1, 1);
idx_start_nadir_pointing = find(sim_switch2nadir>=2, 1);

sim_error_repointing = sim_error(idx_start_repointing:end);
sim_error_linear = sim_error(idx_start_nadir_pointing:end);

idx_start_stat_period = 200+find(sim_eclipse(1:floor(T/Ts))==0, 1, "last");
idx_end_stat_period = floor(T/Ts) + find(sim_eclipse(floor(T/Ts):2*floor(T/Ts))==0, 1, "first");

sim_error_stat = rad2deg(sim_error(idx_start_stat_period:idx_end_stat_period));

error_list(:,end+1) = sim_error_stat;

%% Data extraction: RW staturation estimation

RW.omega = simout.RW_omega;

first_sample = idx_start_nadir_pointing + 200;

coef1 = polyfit(sim_time(first_sample:end), RW.omega(1, first_sample:end), 1);
coef2 = polyfit(sim_time(first_sample:end), RW.omega(2, first_sample:end), 1);
coef3 = polyfit(sim_time(first_sample:end), RW.omega(3, first_sample:end), 1);

max_slope = max(abs([coef1(1), coef2(1), coef3(1)]));
saturation_time = RW.max_torque_start_descending/max_slope;
saturation_data(end+1) = saturation_time;
disp("The saturation time is " + saturation_time/3600/24/365 + " years")

end

%% Statistical analysis Plots

figure
plot(error_list)
xlabel("Timestamp")
ylabel("Error on pointing [°]")
title("Error of Linear Nadir Pointing using " + simulation_number + " simulations")
grid on
% fontsize(36, "points")


figure
histogram(error_list, "Normalization", "probability")
xlabel("Error [°]")
ylabel("Probability")
title("Histogram of error of Linear Nadir Pointing using " + simulation_number + " simulations")
grid on
% fontsize(36, "points")


figure
hold on
cdfplot(reshape(error_list, [size(error_list,1)*size(error_list,2),1]))
xlabel("Error [°]")
ylabel("Cumulative Probability")
title("Cumulative probability of pointing error using " + simulation_number + " simulations")
yline(0.95)
yline(0.997)
legend("Cumulative function", "2-sigma (Normal dist)", "3-sigma (Normal dist)")
% fontsize(36, "points")


figure
qqplot(reshape(error_list, [size(error_list,1)*size(error_list,2),1]))
% fontsize(36, "points")


% Saturation
figure
histogram(saturation_data/3600/24/365, "Normalization", "probability")
xlabel("Time to Saturation of 1 RW [years]")
ylabel("Probability")
title("Histogram of Saturation of RW using " + simulation_number + " simulations")
grid on
% fontsize(36, "points")


figure
hold on
cdfplot(saturation_data/3600/24/365)
xlabel("Time to Saturation of 1 RW [years]")
ylabel("Cumulative Probability")
title("Cumulative probability of Saturation of RW using " + simulation_number + " simulations")
yline(0.95)
yline(0.997)
legend("Cumulative function", "2-sigma (Normal dist)", "3-sigma (Normal dist)")
% fontsize(36, "points")

%% IC

% IC mean of error
meanError = mean(error_list)'; % dataset (delle medie)

IC(0.01, meanError, "mean of linear error", "°")

% IC 95% of error
error95 = prctile(error_list, 95)';

IC(0.05, error95, "95-percentile of linear error", "°")

% IC max of error
error100 = max(error_list)';

IC(0.05, error100, "max error of linear error", "°")

% IC mean saturation time

IC(0.01, saturation_data/3600/24/365, "mean of saturation time", "years")

%% Test d'ipotesi

MeanHypothesisTest(11, 0.01, saturation_data, "mean of saturation time", "years", "<")

MeanHypothesisTest(1, 0.05, error95, "95-percentile of linear error", "°", ">")

MeanHypothesisTest(1, 0.05, error100, "max error of linear error", "°", ">")

%% Auxillary for STATISTICS

function MeanHypothesisTest(mu0, alfa, data, name, unit, option)
    N_data = length(data);
    mean_data = mean(data);
    s_data = std(data);
    s_data_c = s_data*sqrt(N_data/(N_data-1));
    
    z_gamma = norminv(1-alfa);

    if(option == "<") % H0: mu < mu0
        test = (mean_data-mu0)/(sqrt(s_data_c^2/N_data));
        refusal = test > z_gamma;
    end
    if(option == ">") % H0: mu > mu0
        test = (mean_data-mu0)/(sqrt(s_data_c^2/N_data));
        refusal = test < - z_gamma;
    end
    
    if refusal == 0
        text = "CANNOT be refused";
    else
        text = "CAN be refused";
    end

    disp("The hypotesis regarding " + name + " is H0 = mean" + ...
        " " + option + " " + mu0 + unit + " " + text + "with a confidence of " + ...
        (1-alfa) +"%");
    disp("  ");
end

function IC(alfa, data, name, unit)

    N_data = length(data);
    % if N_data < 30
    %     disp("The number is too low, numerosity hypothesis is not met (N>=30)")
    % end

    mean_data = mean(data);
    s_data = std(data);
    s_data_c = s_data*sqrt(N_data/(N_data-1));

    z_gamma = norminv(1-alfa/2);
    IC_mean = mean_data + [-1 1]*(z_gamma*sqrt(s_data_c^2/N_data));

    disp("The IC for the mean of " + name + " is " + IC_mean(1) + unit + " - " + IC_mean(2)  + unit + " with a confidence of " + (1-alfa)*100 + "%");
    disp("  ");
end

%% Auxillary for Linear K gain matrix computing
% computes the minimum gain in a proportional controller
function k = ComputeMinK(R, h, th, a_p)

    k = [];
    re_a = real(a_p)';
    im_a = imag(a_p)';
    
    h_limit = -h;
    xth_limit = @(y) - 1/tan(th).*y.*(y>=0) + 1/tan(th).*y.*(y<0);
    xR_limit = @(y) - sqrt(R^2-y.^2);
    
    for index = 1:length(a_p)
        limit = [h_limit, xth_limit(im_a(index)), xR_limit(im_a(index))];
    
        effective_limit = min(limit, [], "all");
    
        if re_a(index) <=  effective_limit
            k(index) = 0;
        else
            k(index) = effective_limit - re_a(index);
        end 
    end
end