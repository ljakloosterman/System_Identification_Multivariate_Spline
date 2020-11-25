close all; clear all; clc
addpath('functions')
rng default

WLS = 0;                    % For weighted least squares change to 1
iterate = 0;                % change to 1 if iterations over different...
if iterate                  % number of triangles and orders is desired.
    max_simplices = 18;     % Iterating takes a long time with high max numbers
    max_order = 10;
else
    Nsimplices = 2;
    order = 10;
end

%% Loading data
load('data/reconstructed_flight_data')
xeval = [alpha beta];

%% Splitting data
cv = cvpartition(size(xeval,1),'HoldOut',0.5);
idx = cv.test;
X_id = xeval(~idx,:);       X_val = xeval(idx,:);   
Y_id = Cm(~idx,:);          Y_val = Cm(idx,:);

%% Find number of simplices and order with highest accuracy of fit
if iterate
    mse_mat = zeros(max_simplices/2,max_order);

    for Nsimplices = 2:2:max_simplices

        % define vertices
        V = vertices(xeval, Nsimplices);

        % Construct triangulation
        tri = delaunayTriangulation(V);
        tri = triangulation(sort(tri.ConnectivityList, 2), V);

        % Find Barycoordinates and the correspomding triangle
        [Tn_id, Bcor_id] = tsearchn(V, tri, X_id);
        [Tn_val, Bcor_val] = tsearchn(V, tri, X_val);

        % Sort X and Y based on corresponding triangles
        X_id = sortrows([Tn_id X_id],1);    X_val = sortrows([Tn_val X_val],1);
        Y_id = sortrows([Tn_id Y_id],1);    Y_val = sortrows([Tn_val Y_val],1);

        X_id = X_id(:,2:end);   X_val = X_val(:,2:end);
        Y_id = Y_id(:,2:end);   Y_val = Y_val(:,2:end);

        for order = 1:max_order

            % Construct B-net
            kappa = sortrows(partitions(order,[1 1 1]), 'descend');
            Bnet = [];
            for i = 1:Nsimplices
                Bnet = [Bnet; bsplinen_bary2cart(V(tri.ConnectivityList(i,:),:), kappa/order)];
            end

            % Constuct B-form regression matrix
            B_id = sparse([]);
            B_val = sparse([]);
            for n = 1:Nsimplices
                B_id = blkdiag(B_id, x2fx(Bcor_id(Tn_id == n,:), kappa));
                B_val = blkdiag(B_val, x2fx(Bcor_val(Tn_val == n,:), kappa));
            end

            % Construct smoothness matrix
            H = sparse(zeros(0,size(B_id,2)));
            n = 1;
            for i = 1:length(Bnet)
                for j = i:length(Bnet)
                    if round(Bnet(i,:),8) == round(Bnet(j,:),8) & i ~= j
                        H(n,i) = 1;
                        H(n,j) = -1;
                        n = n + 1;
                    end
                end
            end
            
            % Compute Karush-Kuhn-Tucker matrix
            KKT = [B_id'*B_id H'; H zeros(size(H,1))];

            % Construct global B-coefficient vector
            C_id = B_id'*Y_id;

            % Formulate ordinary least squares B-coefficint/Langrangian estimator
            c_ols = pinv(full(KKT))*[C_id; zeros(size(H,1),1)];
            c_ols = c_ols(1:size(B_val,2),1);
            

            % Calculate estimate
            Y_est = B_val*c_ols;
            
            if WLS
                % Formulate weighted least squares B-coefficint estimator
                p_ols = B_id*c_ols;
                W = diag((p_ols - Y_id).^2);
                c_ols = pinv(B_id'*inv(W)*B_id)*B_id'*inv(W)*Y_id;

                % Calculate estimate
                Y_est = B_val*c_ols;
                mse = immse(Y_est, Y_val);
            end

            % Calculate error and append to list
            mse = immse(Y_est, Y_val);
            mse_mat(Nsimplices/2,order) = mse;
            fprintf('Number of simplices %2.0f, order: %2.0f, mean squared error: %5.4d \n', Nsimplices, order, mse)
        end
    end
    % Find order with highest accuray of fit
    mse = min(min(mse_mat));
    [Nsimplices, order] = find(mse_mat == mse);
    Nsimplices = Nsimplices*2;
end

%% Calculate estimate with highest accuracy of fit
% define vertices
V = vertices(xeval, Nsimplices);

% Construct triangulation
tri = delaunayTriangulation(V);
tri = triangulation(sort(tri.ConnectivityList, 2), V);

% Find Barycoordinates and the correspomding triangle
[Tn_id, Bcor_id] = tsearchn(V, tri, X_id);
[Tn_val, Bcor_val] = tsearchn(V, tri, X_val);

% Sort X and Y based on corresponding triangles
X_id = sortrows([Tn_id X_id],1);    X_val = sortrows([Tn_val X_val],1);
Y_id = sortrows([Tn_id Y_id],1);    Y_val = sortrows([Tn_val Y_val],1);

X_id = X_id(:,2:end);   X_val = X_val(:,2:end);
Y_id = Y_id(:,2:end);   Y_val = Y_val(:,2:end);

% Construct B-net
kappa = sortrows(partitions(order,[1 1 1]), 'descend');
Bnet = [];
for i = 1:Nsimplices
    Bnet = [Bnet; bsplinen_bary2cart(V(tri.ConnectivityList(i,:),:), kappa/order)];
end

 % Constuct B-form regression matrix
B_id = sparse([]);
B_val = sparse([]);
for n = 1:Nsimplices
    B_id = blkdiag(B_id, x2fx(Bcor_id(Tn_id == n,:), kappa));
    B_val = blkdiag(B_val, x2fx(Bcor_val(Tn_val == n,:), kappa));
end

% Construct smoothness matrix
H = sparse(zeros(0,size(B_id,2)));
n = 1;
for i = 1:length(Bnet)
    for j = i:length(Bnet)
        if round(Bnet(i,:),8) == round(Bnet(j,:),8) & i ~= j
            H(n,i) = 1;
            H(n,j) = -1;
            n = n + 1;
        end
    end
end

% Compute Karush-Kuhn-Tucker matrix
KKT = [B_id'*B_id H'; H zeros(size(H,1))];

% Construct global B-coefficient vector
C_id = B_id'*Y_id;

% Formulate ordinary least squares B-coefficint/Langrangian estimator
c_ols = pinv(full(KKT))*[C_id; zeros(size(H,1),1)];
c_ols = c_ols(1:size(B_val,2),1);

% Calculate estimate
Y_est = B_val*c_ols;
mse = immse(Y_est, Y_val);

if WLS
    % Formulate weighted least squares B-coefficint estimator
    p_ols = B_id*c_ols;
    W = diag((p_ols - Y_id).^2);
    c_ols = pinv(B_id'*inv(W)*B_id)*B_id'*inv(W)*Y_id;
    
    % Calculate estimate
    Y_est = B_val*c_ols;
    mse = immse(Y_est, Y_val);
end

%% Residual Validation
E = Y_val - Y_est;
E_mean = mean(E);
conf = 1.96/sqrt(length(E));
[acx,lags] = xcorr(E-mean(E), 'coeff');
perc_inside_bounds = length(acx(abs(acx) < conf)) / length(acx) * 100;

%% Statistical model quality analyses
cov = inv(B_val'*B_val);
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

if iterate
    figure(2)
    bar = bar3(mse_mat);
    for i = 1:size(mse_mat,2)
        zdata = ones(6*size(mse_mat,1),4);
        k = 1;
        for j = 0:6:(6*size(mse_mat,1)-6)
          zdata(j+1:j+6,:) = mse_mat(k,i);
          k = k+1;
        end
        set(bar(i),'Cdata',zdata)
    end
    xlabel('Simplex order')
    ylabel('Number of simplices')
    zlabel('mean squared error')
    yticks((2:2:max_simplices)/2)
    yticklabels(split(num2str(2:2:max_simplices)))
    title('Mean squared error for different number of simplices and orders')
end


figure(3)
hold on
plot(xeval(:,1), xeval(:,2), '.')
triplot(delaunayTriangulation(V), 'linewidth', 2, 'Color', 'r')
xlabel('Angle of attack [degrees]')
ylabel('Side slipe angle [degrees]')
title('Simplex and data')

figure(4)
hold on
triplot(delaunayTriangulation(V), 'linewidth', 2, 'Color', 'r')
scatter(Bnet(:,1), Bnet(:,2), 50, 'filled', 'b')
for i = 1:length(Bnet)/length(kappa)
    if mod(i,2) == 1
        hor = 'left'; ver = 'bottom';
    else
        hor = 'right'; ver = 'top';
    end
    text(Bnet(length(kappa)*(i-1)+1:i*length(kappa),1),...
         Bnet(length(kappa)*(i-1)+1:i*length(kappa),2),...
         cellstr(num2str(kappa)),...
         'HorizontalAlignment', hor, 'VerticalAlignment', ver)
end
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

if WLS == 0
    figure(6)
    spy(KKT)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('Karush-Kuhn-Tucker matrix spy plot')
end

figure(7)
plot(E)
xlabel("Data point")
ylabel('Residual [-]')
title('Model residual')

figure(8);
hold on;
line([lags(1), lags(end)], [conf, conf], 'Color','red')
line([lags(1), lags(end)], [-conf, -conf], 'Color','red')
plot(lags, acx)
xlabel('Number of lags')
ylabel('Auto-correlation')
title('Model error auto-correlation')

figure(9)
hold on
plot(1:length(VAR), VAR)
xlabel('Coefficient index')
ylabel('Coefficient variance')
title('Coefficient variances')

fprintf('Final mean squared error: %5.4d \n', mse)
fprintf('mean residual: %5.4d \n', E_mean)
fprintf('Percentage inside confidence bounds: %3.2f \n', perc_inside_bounds)