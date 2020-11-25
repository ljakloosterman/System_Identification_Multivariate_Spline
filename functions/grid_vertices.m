function V = vertices(input, Ntriangles)
    % Append 1 traingle if desired number is uneven
    if mod(Ntriangles,2) ~= 0
        Ntriangles = Ntriangles + 1;
    end

    % Find most optimal meshgrid dimension
    div = divisors(Ntriangles/2);
    list = zeros(1,round(length(div)/2));
    for i = 1:length(div)/2
        list(i) = div(end-i+1) - div(i);
    end
    [~, idx] = min(list);
    
    % Construct meshgrid
    [x,y] = minboundrect(input(:,1), input(:,2));
    X = linspace(min(x), max(x), div(end-idx+1) + 1);
    Y = linspace(min(y), max(y), div(idx) + 1);
    
    % Find vertices
    [X_grid,Y_grid] = meshgrid(X,Y);
    V = [X_grid(:) Y_grid(:)];
end