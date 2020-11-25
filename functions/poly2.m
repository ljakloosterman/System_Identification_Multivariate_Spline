function [Polynominal, A_val] = poly2(X_in, X_val, Y_in, order)
    A_in = zeros(size(X_in,1),(order+1)^size(X_in,2));
    idx = 0;
    for i = 0:order
        for j = 0:order
            idx = idx + 1;
            A_in(:,idx) = (X_in(:,1).^i).*(X_in(:,2).^j);
        end
    end
    theta_ols = pinv(A_in'*A_in)*A_in'*Y_in;
    
    % evaluation
    A_val = zeros(size(X_val,1),(order+1)^size(X_val,2));
    idx = 0;
    for i = 0:order
        for j = 0:order
            idx = idx + 1;
            A_val(:,idx) = (X_val(:,1).^i).*(X_val(:,2).^j);
        end
    end
    Polynominal = A_val*theta_ols;
end