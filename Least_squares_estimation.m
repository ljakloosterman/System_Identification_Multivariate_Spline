close all; clear all; clc
addpath('functions')
rng default

iterate = 0;        % Set to one if iterations over different orders is desired 
if iterate
    max_order = 25; % Set max order for iterations
else
    order = 9;      % set order
end

%% Loading data
load('data/reconstructed_flight_data')
xeval = [alpha beta];

%% Splitting data
cv = cvpartition(size(xeval,1),'HoldOut',0.5);
idx = cv.test;
X_id = xeval(~idx,:);       X_val = xeval(idx,:);   
Y_id = Cm(~idx,:);          Y_val = Cm(idx,:);

%% find order with highest accuracy of fit
if iterate
    mse_list = zeros(1,max_order);

    for order = 1:max_order
        % Calculate estimate
        P = poly2(X_id, X_val, Y_id, order);

        % Calculate error and append to list
        mse = immse(P, Y_val);
        mse_list(order) = mse;
        fprintf('Order: %2.0f, mean squared error: %5.4d \n', order, mse)
    end
    [mse, order] = min(mse_list);
end

%% Calculate estimate with highest accuracy of fit
[P, A] = poly2(X_id, X_val, Y_id, order);
mse = immse(P, Y_val);

%% Residual Validation
E = Y_val - P;
E_mean = mean(E);
conf = 1.96 / sqrt(length(E));
[acx,lags] = xcorr(E-mean(E), 'coeff');
perc_inside_bounds = length(acx(abs(acx) < conf)) / length(acx) * 100;

%% Statistical model quality analyses
cov = pinv(A' * A);
sig = (E' * E) / (size(A, 1) - size(A, 2));
VAR = sig * diag(cov);

%% Plotting the results
figure(1)
hold on
plot3(rad2deg(alpha), rad2deg(beta), Cm, '.k')
grid on
view(45, 45)
xlabel('Angle of attack [degrees]')
ylabel('Side slipe angle [degrees]')
zlabel('C_m[-]')
title('Input data')

if iterate
    figure(2)
    hold on
    plot(order, mse, 'd', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
    bar(mse_list)
    ylabel('Accuracy of fit error')
    xlabel('Order of polynominal')
    legend('Best accuracy of fit')
    title('Mean squared error for different orders')
end

figure(3)
hold on
trisurf(delaunayn(X_val), rad2deg(X_val(:,1)), rad2deg(X_val(:,2)), P, 'EdgeColor', 'none')
plot3(rad2deg(alpha), rad2deg(beta), Cm, '.k')
grid on
view(45, 45)
xlabel('Angle of attack [degrees]')
ylabel('Side slipe angle [degrees]')
zlabel('C_m[-]')
title('Polynominal')

figure(4)
plot(E)
xlabel("Data point")
ylabel('Residual [-]')
title('Model error')

figure(5);
hold on;
line([lags(1), lags(end)], [conf, conf], 'Color','red')
line([lags(1), lags(end)], [-conf, -conf], 'Color','red')
plot(lags, acx)
legend('95% confidence bounds')
xlabel('Number of lags')
ylabel('Auto-correlation')
title('Model error auto-correlation')

figure(6)
hold on
plot(1:length(VAR), VAR)
xlabel('Coefficient index')
ylabel('Coefficient variance')
title('Coefficient variances')

fprintf('Final mean squared error: %5.4d \n', mse)
fprintf('mean residual: %5.4d \n', E_mean)
fprintf('Percentage inside confidence bounds: %3.2f \n', perc_inside_bounds)