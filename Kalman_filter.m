close all
clear all
clc
addpath('functions')

%% Loading data
dataname = 'F16traindata_CMabV_2020';
load(dataname, 'Z_k', 'U_k', 'Cm')
Z_k = Z_k';
U_k = U_k';

%% Simulation parameters
doIEKF      = 0;            % IEKF will be used if set to 1
max_it      = 100;
dt          = 0.01;
[m,n]       = size(Z_k);
eps         = 1e-10;

%% Initialisation
% States
Ex_0    = [Z_k(3,1); 0; 0; 0]; % Initial guess for optimal value

% Covariance matrix estimate
stdx_0  = [0.1 0.1 0.1 0.1];
P_0     = diag(stdx_0.^2);

% Process noise matrix
sigma_w     = [1e-3, 1e-3, 1e-3, 0];
Q           = diag(sigma_w.^2);

% Sensor noise matrix
sigma_v     = [0.035, 0.013, 0.110];
R           = diag(sigma_v.^2);

% System equation matrices
B   = eye(m+1,m);
G   = zeros(m+1);
Fx  = zeros(m+1);

%% Kalman filter
% Create space for all data
XX_k1k1     = zeros(m+1, n);
PP_k1k1     = zeros(m+1, n);
STDx_cor    = zeros(m+1, n);
z_pred      = zeros(m, n);

% Initial estimate
x_k_1k_1 = Ex_0; % x(0|0)=E{x_0}
P_k_1k_1 = P_0; % P(0|0)=P(0)

% time set up
t0 = 0;
t1 = dt;

for i = 1:n
    % One step ahead prediction
    [t, x_kk_1] = rk4(@calc_f, x_k_1k_1, U_k(:,i), [t0 t1]);
    z_kk_1 = calc_h(x_kk_1);
    z_pred(:,i) = z_kk_1;
    
    % Continous to discrete time transformations
    [~, Psi] = c2d(Fx, B, dt);   
    [Phi, Gamma] = c2d(Fx, G, dt);
    
    % Covariance matrix prediction
    P_kk_1 = Phi*P_k_1k_1*Phi' + Gamma*Q*Gamma'; 
    P_pred = diag(P_kk_1);
    stdx_pred = sqrt(diag(P_kk_1));
    
    % Run iterive part if doIEKf = 1
    if (doIEKF)
        eta2    = x_kk_1;
        err     = 2*eps;
        
        it = 0;
        while (err > eps)
            if (it >= max_it)
                fprintf('max iterations exceeded \n');
                break
            end
            it      = it + 1;
            eta1    = eta2;
            
            % Calculate jacobian of H
            Hx = calc_hx(eta1);
            
            % Check observability of state
            if (i == 1 && it == 1)
                rankHF = calcObsRank(Hx, Fx);
                if (rankHF < (m))
                    warning('Current state is not observable')
                end
            end
            
            % Calculate Kalman gain
            K = P_kk_1*Hx'/(Hx*P_kk_1*Hx' + R);
            
            % New observation state
            z_p = calc_h(eta1);
            eta2    = x_kk_1 + K * (Z_k(:,i) - z_p - Hx*(x_kk_1 - eta1));
            err     = norm((eta2 - eta1), inf) / norm(eta1, inf);
        end
        x_k_1k_1 = eta2;
    else
        % perturbation of h(x,u,t)
        Hx = calc_hx(x_kk_1);
        
        % Calculate Kalman gain
        K = P_kk_1*Hx'/(Hx*P_kk_1*Hx' + R);
        
        x_k_1k_1 = x_kk_1 + K * (Z_k(:,i) - z_kk_1); 
    end
    
    % Covariance matrix of state estimation error
    P_k_1k_1 = (eye(m+1) - K*Hx) * P_kk_1 * (eye(m+1) - K*Hx)' + K*R*K';
    P_cor = diag(P_k_1k_1);
    stdx_cor = sqrt(diag(P_k_1k_1));
    
    % Next step
    t0 = t1; 
    t1 = t1 + dt;
    
    % Store results
    XX_k1k1(:,i) = x_k_1k_1;
    STDx_cor(:,i) = stdx_cor;
end

%% Calculate alpha_true
alpha_t = atan(XX_k1k1(3,:)./XX_k1k1(1,:));

%% Save outcome
alpha = alpha_t';
beta = Z_k(2,:)'; 
save('data/reconstructed_flight_data.mat', 'alpha', 'beta', 'Cm')

%% Plotting results
time = 0:dt:100;

figure(1)
subplot(3,2,1);
hold on
plot(time,rad2deg(Z_k(1,:)),time,rad2deg(z_pred(1,:)))
legend('Measured state', 'Estimated', 'Location', 'northwest')
ylabel('Angle of attack [degrees]')
xlabel('Time [sec]')

subplot(3,2,2);
plot(time,rad2deg(z_pred(1,:) - Z_k(1,:)))
ylabel('Angle of attack error [degrees]')
xlabel('Time [sec]')

subplot(3,2,3);
hold on
plot(time,rad2deg(Z_k(2,:)),time,rad2deg(z_pred(2,:)))
legend('Measured state', 'Estimated', 'Location', 'northwest')
ylabel('Side slip angle [degrees]')
xlabel('Time [sec]')

subplot(3,2,4);
plot(time,rad2deg(z_pred(2,:) - Z_k(2,:)))
ylabel('Side slip angle error [degrees]')
xlabel('Time [sec]')

subplot(3,2,5);
hold on
plot(time,Z_k(3,:),time,z_pred(3,:))
legend('Measured state', 'Estimated', 'Location', 'northwest')
ylabel('Velocity [m/s]')
xlabel('Time [sec]')

subplot(3,2,6);
plot(time,z_pred(3,:) - Z_k(3,:))
ylabel('Velocity error [m/s]')
xlabel('Time [sec]')

figure(2)
subplot(3,1,1);
plot(time,XX_k1k1(1,:))
ylabel('u [m/s]')
xlabel('Time [sec]')

subplot(3,1,2);
plot(time,XX_k1k1(2,:))
ylabel('v [m/s]')
xlabel('Time [sec]')

subplot(3,1,3);
plot(time,XX_k1k1(3,:))
ylabel('w [m/s]')
xlabel('Time [sec]')

figure(3)
subplot(2,1,1);
plot(time,XX_k1k1(4,:))
ylabel('C_a [-]')
xlabel('Time [sec]')

subplot(2,1,2);
plot(time,rad2deg(alpha_t),time,rad2deg(Z_k(1,:)))
legend('True state', 'Measured state', 'Location', 'northwest')
ylabel('Angle of attack [degrees]')
xlabel('Time [sec]')