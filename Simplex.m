close all; clear all; clc
addpath('functions')
rng default

iterate = 0;        % Set to one if iterations over different orders is desired
if iterate
    max_order = 25; % Set max order for iterations
else
    order = 13;     % set order
end

%% Loading data
load('data/reconstructed_flight_data')
xeval = [alpha beta];

%% Splitting data
cv = cvpartition(size(xeval,1),'HoldOut',0.5);
idx = cv.test;
X_id = xeval(~idx,:);       X_val = xeval(idx,:);   
Y_id = Cm(~idx,:);          Y_val = Cm(idx,:);

%% Define vertices
[Vx,Vy] = minboundtri(alpha, beta); % Function to define smallest traingle to enclose all datapoints
V = [Vx(1:end-1) Vy(1:end-1)];

%% Construct triangulation
tri = delaunayTriangulation(V);

%% Find barycoordinates and the corresponding triangles
[Tn_id, Bcor_id] = tsearchn(V, tri, X_id);
[Tn_val, Bcor_val] = tsearchn(V, tri, X_val);

%% Find order with highest accuracy of fit
if iterate
    mse_list = zeros(1,max_order);

    for order = 1:max_order
        % Define kappa
        kappa = sortrows(partitions(order,[1 1 1]), 'descend');

        % Constuct B-form regression matrix
        B_id = x2fx(Bcor_id, kappa);
        B_val = x2fx(Bcor_val, kappa);

        % Formulate ordinary least squares estimator
        c_ols = pinv(B_id'*B_id)*B_id'*Y_id;

        % Calculate estimate
        Y_est = B_val*c_ols;

        % Calculate error and append to list
        mse = immse(Y_est, Y_val);
        mse_list(order) = mse;
        fprintf('Order: %2.0f, mean squared error: %5.4d \n', order, mse)
    end
    % Find order with highest accuray of fit
    [mse, order] = min(mse_list);
end

%% Calculate estimate with highest accuracy of fit
% Define kappa
kappa = sortrows(partitions(order,[1 1 1]), 'descend');

% Construct Bnet
Bnet = bsplinen_bary2cart(V, kappa./order);

% Constuct B-form regression matrix
B_id = x2fx(Bcor_id, kappa);
B_val = x2fx(Bcor_val, kappa);

% Formulate ordinary least squares estimator
c_ols = pinv(B_id'*B_id)*B_id'*Y_id;

% Calculate estimate
Y_est = B_val*c_ols;
mse = immse(Y_est, Y_val);

%% Residual Validation
E = Y_val - Y_est;
E_mean = mean(E);
conf = 1.96/sqrt(length(E));
[acx,lags] = xcorr(E-mean(E), 'coeff');
perc_inside_bounds = length(acx(abs(acx) < conf)) / length(acx) * 100;

%% Statistical model quality analyses
cov = pinv(B_val'*B_val);
sig = (E'*E)/(size(B_val,1) - size(B_val,2));
VAR = sig*diag(cov);

%% Plotting the results
figure(1)
plot3(rad2deg(alpha), rad2deg(beta), Cm, '.k')
view(45, 45)
xlabel('Angle of attack [degrees]')
ylabel('Side slipe angle [degrees]')
zlabel('C_m[-]')
title('Input data')

figure(2)
hold on
plot(rad2deg(xeval(:,1)), rad2deg(xeval(:,2)), '.', 'markerSize', 5)
triplot(delaunayTriangulation(rad2deg(V)), 'linewidth', 2, 'Color', 'r')
xlabel('Angle of attack [degrees]')
ylabel('Side slipe angle [degrees]')
title('Simplex and data')

if iterate
    figure(3)
    hold on
    plot(order, mse, 'd', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
    bar(mse_list)
    ylabel('Mean squared error')
    xlabel('Order of simplex')
    legend('Best accuracy of fit')
    title('Mean squared error for different orders')
end

figure(4)
hold on
triplot(delaunayTriangulation(V), 'linewidth', 2, 'Color', 'r')
scatter(Bnet(:,1), Bnet(:,2), 50, 'filled', 'b')
text(Bnet(:,1)+0.001, Bnet(:,2)-0.01, cellstr(num2str(kappa)));
xlabel('Angle of attack [degrees]')
ylabel('Side slipe angle [degrees]')
title('B-net')

figure(5)
hold on
plot3(rad2deg(X_val(:,1)), rad2deg(X_val(:,2)), Y_val, '.', 'markerSize', 5, 'Color', 'k')
triplot(delaunayTriangulation(rad2deg(V)), 'linewidth', 2, 'Color', 'b')
trisurf(delaunayn([rad2deg(X_val(:,1)) rad2deg(X_val(:,2))]), rad2deg(X_val(:,1)), rad2deg(X_val(:,2)), Y_est, 'EdgeColor', 'none')
view(45, 45)
xlabel('Angle of attack [degrees]')
ylabel('Side slipe angle [degrees]')
zlabel('C_m[-]')
title('Simplex polynominal')

figure(6)
plot(E)
xlabel("Data point")
ylabel('Residual [-]')
title('Model residual')

figure(7);
hold on;
line([lags(1), lags(end)], [conf, conf], 'Color','red')
line([lags(1), lags(end)], [-conf, -conf], 'Color','red')
plot(lags, acx)
xlabel('Number of lags')
ylabel('Auto-correlation')
title('Model error auto-correlation')

figure(8)
hold on
plot(1:length(VAR), VAR)
xlabel('Coefficient index')
ylabel('Coefficient variance')
title('Coefficient variances')

fprintf('Final mean squared error: %5.4d \n', mse)
fprintf('mean residual: %5.4d \n', E_mean)
fprintf('Percentage inside confidence bounds: %3.2f \n', perc_inside_bounds)