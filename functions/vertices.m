function V = vertices(input, Ntriangles)

    [x,y] = minboundrect(input(:,1), input(:,2));
    V = [x y];
    
    Ncenters = Ntriangles/2 - 1;
    if Ncenters > 0
        [~,C] = kmeans(input, Ntriangles/2 - 1, 'Maxiter', 1000);
        V = [V; C];
    end 
    V = unique(V, 'rows');
end