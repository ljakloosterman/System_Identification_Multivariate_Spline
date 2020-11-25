function zpred = calc_h(x)
    u = x(1);
    v = x(2);
    w = x(3);
    C_alpha = x(4);
    zpred = [atan(w/u)*(1 + C_alpha); atan(v/sqrt(u^2 + w^2)); sqrt(u^2 + v^2 + w^2)];
end