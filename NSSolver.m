function [] = NSSolver(X, Y, Z, dx, dy, dz, t, dt, rho, nu, beta, F_z, D_h, isPlot)
    format long
    % Define sizes based on input grids
    gridSize = size(X);
    U_size  = [t, gridSize, 3];

    U       = zeros(U_size);
    U_temp  = zeros(U_size);
    P       = zeros([t, gridSize]);
    P_temp  = zeros([t, gridSize]);

    % Variables for convective and diffusive terms
    u_conv      = zeros(3, 1);
    u_diff      = zeros(3, 1);
    % Variables for pressure gradient and velocity divergence
    p_force     = 0;
    p_moment    = 0;

    % Initializing the Turbulent Viscosity
    C_s     = 0.15;             % Smagorinsky constant 
    Delta   = sqrt((dx^2 + dy^2 + dz^2)); 
    nu_t    = 0;

    % Reynolds Number
    u_mag           = 0;
    u_tot           = 0;
    Re_map          = zeros(t-1,1);
    t_map           = zeros(t-1,1);

    % Initialize the figures for the 3D quiver plot and pressure contour plot if plotting is enabled
    if isPlot
        quiverFig = figure('Name', '3D Quiver Plot');
        contourFig = figure('Name', 'Pressure Contour Plot');
    end

    for t_n = 1:t-1
        u_tot = 0; % Reset the total velocity magnitude for each time step
        for i = 2:gridSize(1) - 1
            for j = 2:gridSize(2) - 1
                for k = 2:gridSize(3) - 1

                    % Vector components of U
                    u   = U(t_n, i, j, k, 1);
                    v   = U(t_n, i, j, k, 2);
                    w   = U(t_n, i, j, k, 3);

                    for A = 1:3
                        % Convective and diffusive terms
                        u_conv(A, 1) = -0.5 * ( ...
                            u * (U(t_n, i+1, j, k, A) - U(t_n, i-1, j, k, A)) / dx + ...
                            v * (U(t_n, i, j+1, k, A) - U(t_n, i, j-1, k, A)) / dy + ...
                            w * (U(t_n, i, j, k+1, A) - U(t_n, i, j, k-1, A)) / dz ...
                        );

                        u_diff(A, 1) =  ( ...
                            (U(t_n, i+1, j, k, A) - 2 * U(t_n, i, j, k, A) + U(t_n, i-1, j, k, A)) / dx^2 + ...
                            (U(t_n, i, j+1, k, A) - 2 * U(t_n, i, j, k, A) + U(t_n, i, j-1, k, A)) / dy^2 + ...
                            (U(t_n, i, j, k+1, A) - 2 * U(t_n, i, j, k, A) + U(t_n, i, j, k-1, A)) / dz^2 ...
                        );
                    end

                    % Pressure force and momentum terms
                    p_force = 0.5 * beta * ( ...
                        (U(t_n, i+1, j, k, 1) - U(t_n, i-1, j, k, 1)) / dx + ...
                        (U(t_n, i, j+1, k, 2) - U(t_n, i, j-1, k, 2)) / dy + ...
                        (U(t_n, i, j, k+1, 3) - U(t_n, i, j, k-1, 3)) / dz ...
                    );
                    p_moment = -0.5 * ( ...
                        u * (P(t_n, i+1, j, k) - P(t_n, i-1, j, k)) / dx + ...
                        v * (P(t_n, i, j+1, k) - P(t_n, i, j-1, k)) / dy + ...
                        w * (P(t_n, i, j, k+1) - P(t_n, i, j, k-1)) / dz ...
                    );

                    % Update pressure field
                    P_temp(t_n + 1, i, j, k) = P(t_n, i, j, k) + dt * (p_force + p_moment);

                    % Smear the pressure field a bit in 3D
                    P_temp(t_n + 1, i, j, k) = (1/3) * P_temp(t_n + 1, i, j, k) + ...
                           (1/12) * (P_temp(t_n + 1, i-1, j, k) + P_temp(t_n + 1, i+1, j, k) + ...
                                   P_temp(t_n + 1, i, j-1, k) + P_temp(t_n + 1, i, j+1, k) + ...
                                   P_temp(t_n + 1, i, j, k-1) + P_temp(t_n + 1, i, j, k+1));

                    % Turbulent Viscosity 
                    nu_t = addTurbulentViscosity(U, t_n, i, j, k, C_s, Delta);

                    % Update velocity components
                    U_temp(t_n + 1, i, j, k, 1) = u + dt * ( ...
                        u_conv(1, 1) + (nu + nu_t) * u_diff(1, 1) - ...
                        1/(2*rho*dx) * (P(t_n, i+1, j, k) - P(t_n, i-1, j, k)) + ...
                        1 / rho * F_z ...
                    );

                    U_temp(t_n + 1, i, j, k, 2) = v + dt * ( ...
                        u_conv(2, 1) + (nu + nu_t) * u_diff(2, 1) - ...
                        1/(2*rho*dy) * (P(t_n, i, j+1, k) - P(t_n, i, j-1, k)) + ...
                        1 / rho * F_z ...
                    );

                    U_temp(t_n + 1, i, j, k, 3) = w + dt * ( ...
                        u_conv(3, 1) + (nu + nu_t) * u_diff(3, 1) - ...
                        1/(2*rho*dz) * (P(t_n, i, j, k+1) - P(t_n, i, j, k-1)) + ...
                        1 / rho * F_z ...
                    );

                    % Taking the resulting vector magnitude and summing all for t_n
                    u_mag = sqrt(u^2 + v^2 + w^2);
                    u_tot = u_tot + u_mag;
                end
            end
        end

        % Apply no-slip boundary conditions
        U_temp = applyNoSlipBoundary(U_temp, t_n);

        % Updating the matrices:
        U = U_temp;
        P = P_temp; 

        % Calculating Re & updating the plot
        u_mean          = u_tot / ((size(X, 1) - 2) * (size(Y, 2) - 2) * (size(Z, 3) - 2));
        Re              = D_h * u_mean / nu;
        t_map(t_n)      = t_n;
        Re_map(t_n)     = Re;

        % Update the 3D quiver plot and pressure contour plot
        if isPlot
            figure(quiverFig);
            plotQuiver3D(U, X, Y, Z, t_n);
            figure(contourFig);
            plotPressureContour(P, X, Y, Z, t_n, round(gridSize(3) / 2)); % Plot at mid z-level
            pause(0.1); % Pause to control the update speed
        end
    end

    if isPlot
        % Final plot at the last time step
        figure(quiverFig);
        plotQuiver3D(U, X, Y, Z, t);
        
        % Plot Reynolds number change over time
        figure;
        plot(t_map, Re_map);
        xlabel('Time [t]');
        ylabel('Reynolds Number [Re]');
        title('Change of Reynolds number over time steps t');
        grid on;
    end
end

function [nu_t] = addTurbulentViscosity(U_field, t_n, i, j, k, C_s, Delta)
    % Computing the strain rate tensor components for S in 3D
    S_11 = (U_field(t_n, i+1, j, k, 1) - U_field(t_n, i-1, j, k, 1)) / (2 * Delta); 
    S_22 = (U_field(t_n, i, j+1, k, 2) - U_field(t_n, i, j-1, k, 2)) / (2 * Delta); 
    S_33 = (U_field(t_n, i, j, k+1, 3) - U_field(t_n, i, j, k-1, 3)) / (2 * Delta); 
    
    S_12 = 0.5 * ((U_field(t_n, i, j+1, k, 1) - U_field(t_n, i, j-1, k, 1)) / (2 * Delta) ...
                + (U_field(t_n, i+1, j, k, 2) - U_field(t_n, i-1, j, k, 2)) / (2 * Delta)); 
    S_21 = S_12;

    S_13 = 0.5 * ((U_field(t_n, i, j, k+1, 1) - U_field(t_n, i, j, k-1, 1)) / (2 * Delta) ...
                + (U_field(t_n, i+1, j, k, 3) - U_field(t_n, i-1, j, k, 3)) / (2 * Delta)); 
    S_31 = S_13;

    S_23 = 0.5 * ((U_field(t_n, i, j, k+1, 2) - U_field(t_n, i, j, k-1, 2)) / (2 * Delta) ...
                + (U_field(t_n, i, j+1, k, 3) - U_field(t_n, i, j-1, k, 3)) / (2 * Delta)); 
    S_32 = S_23;

    % Calculating the magnitude of S in 3D
    S_mag = sqrt(S_11.^2 + S_22.^2 + S_33.^2 + 2 * (S_12.^2 + S_13.^2 + S_23.^2));

    % Calculating the turbulent Viscosity
    nu_t = (C_s * Delta).^2 * S_mag;
end

% Apply no-slip boundary conditions
function [U_bound] = applyNoSlipBoundary(U_field, t_n)
    % Top and Bottom boundaries (along x direction)
    U_field(t_n + 1, 1, :, :, :)        = 0.1; % Top boundary
    U_field(t_n + 1, end, :, :, :)      = 0; % Bottom boundary

    % Left and Right boundaries (along y direction)
    U_field(t_n + 1, :, 1, :, :)        = 0; % Left boundary
    U_field(t_n + 1, :, end, :, :)      = 0; % Right boundary

    % Front and Back boundaries (along z direction)
    U_field(t_n + 1, :, :, 1, :)        = 0; % Front boundary
    U_field(t_n + 1, :, :, end, :)      = 0; % Back boundary

    U_bound = U_field;
end

% Quiver3D plot function
function [] = plotQuiver3D(U, X, Y, Z, t_n)
    % Function to create a 3D quiver plot of the velocity field at time step t_n
    % U: Velocity field (4D array: time x x-dim x y-dim x z-dim x velocity-components)
    % X, Y, Z: Coordinate grids (3D arrays)
    % t_n: Time step to plot (scalar)

    % Extract the 3D velocity components at time step t_n
    Ux = squeeze(U(t_n, :, :, :, 1));
    Uy = squeeze(U(t_n, :, :, :, 2));
    Uz = squeeze(U(t_n, :, :, :, 3));

    % Avoid plotting if all vectors are zero
    if all(Ux(:) == 0) && all(Uy(:) == 0) && all(Uz(:) == 0)
        return;
    end

    % Create the quiver3 plot
    cla; % Clear the current axes
    quiver3(X, Y, Z, Ux, Uy, Uz, 'AutoScale', 'on');
    title(['3D Quiver plot of velocity field at t = ', num2str(t_n)]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    grid on;
    drawnow; % Update the figure
end

% Contour plot function
function plotPressureContour(P, X, Y, Z, t_n, z_level)
    % Extract the pressure field at the given time step and z-level
    P_slice = squeeze(P(t_n, :, :, z_level));
    X_slice = squeeze(X(:, :, z_level));
    Y_slice = squeeze(Y(:, :, z_level));

    % Check if pressure data is constant
    if all(P_slice(:) == P_slice(1))
        warning('Pressure data is constant; contour plot may not render correctly.');
    end

    % Create the contour plot
    cla; % Clear the current axes
    contourf(X_slice, Y_slice, P_slice, 20); % 20 contour levels
    colorbar;
    title(['Pressure Contour at t = ', num2str(t_n), ', z = ', num2str(Z(1, 1, z_level))]);
    xlabel('X');
    ylabel('Y');
    axis equal;
    grid on;
    drawnow; % Update the figure
end
