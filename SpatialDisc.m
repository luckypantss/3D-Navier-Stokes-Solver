function [X, Y, Z, dx, dy, dz] = SpatialDisc(L_x, L_y, L_z, n_x, n_y, n_z, isPlot) 
% Defining the domain of x and y
    x       = linspace(0, L_x, n_x);
    y       = linspace(0, L_y, n_y);
    z       = linspace(0,L_z, n_z);
    % Spacing between the gridpoints
    dx      = L_x / (n_x-1);
    dy      = L_y / (n_y-1);
    dz      = L_z / (n_z-1);
    
    % Seetting the domain into gridpoints
    [X, Y, Z] = meshgrid(x,y,z);

    numOfGrids = length(X) * length(Y) * length(Z);
    
    % Plotting the gridpoints
    if isPlot == true 
        figure
        scatter3(X(:), Y(:), Z(:), 'filled')
        xlabel('X')
        ylabel('Y')
        title(['Discretized points of x, y and z with ', num2str( 3), ' grid points']);
        grid minor
    end
end