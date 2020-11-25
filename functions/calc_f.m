function xdot = calc_f(~, x, u)
    xdot = zeros(size(x,1),1);
    xdot(1:3) = u;
end