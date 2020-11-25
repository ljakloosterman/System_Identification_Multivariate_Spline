function Hx = calc_hx(x)
    u = x(1);
    v = x(2);
    w = x(3) ;
    C_alpha = x(4);
    
    Hx = [-C_alpha*w/(u^2*(1 + w^2/u^2)) - w/(u^2*(1 + w^2/u^2)) 0, C_alpha/(u*(1 + w^2/u^2)) + 1/(u*(1 + w^2/u^2)), atan(w/u);...
        -u*v/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)), 1/(sqrt(u^2 + w^2)*(v^2/(u^2 + w^2) + 1)), -v*w/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)), 0;...
        u/sqrt(u^2 + v^2 + w^2), v/sqrt(u^2 + v^2 + w^2), w/sqrt(u^2 + v^2 + w^2), 0];
end