function r = calcObsRank(H, Fx)

    nstates = size(Fx,1);

    F = eye(size(Fx));
    Rank = [];
    for i = 1:nstates-1
       Rank = [ Rank; H*F ];
       F = F*Fx;
    end
    Rank = [ Rank; H*F ];
    r    = rank(Rank);
end