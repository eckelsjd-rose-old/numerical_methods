%% hw10.m
% Class:       ME427 - Intro to CFD; HW10
% Author:      Joshua Eckels
% Description: Function to evaluate non-linear, 1D, steady Burgers
%              equation with a non-uniform mesh and finite difference
% Strategy:    1. Compute a non-uniform mesh
%              2. Discretize with non-uniform mesh equations
%              4. Define system of N residual eqns for N discrete points
%              5. Solve system with multi-dimensional Newton-Raphson method
%              6. Task 1: converge on solution for nu=0.00001
%              7. Task 2: converge on solution at lowest possible nu
% Date:        4/9/2020
clc
clear all
close all
%% Task 1
fprintf("TASK 1:\n");
fprintf("Enter 'q' at any time to quit.\n");
fprintf("Enter 'c' at any time to continue to TASK 2.\n");
fprintf("Use decimal notation for fractions.\n");
fprintf("Press enter to accept default values.\n");

% Main user input loop
while (1)
    isNumeric = @(S) ~isnan(str2double(S));
    N_DEFAULT = 400;
    A_DEFAULT = -0.05;
    b_DEFAULT = 20;
    c_DEFAULT = 0.5;
    
    % Get number of points N
    fprintf("\nRunning Task 1: find a solution with v=0.00001\n");
    N_str = input('Number of points : ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif N_str == 'c'
        break;
    elseif isNumeric(N_str)
        N1 = str2double(N_str);
    elseif isempty(N_str)
        N1 = N_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
    % Get tuning parameter A
    N_str = input('fGSR amplitutde "A" parameter: ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif N_str == 'c'
        break;
    elseif isNumeric(N_str)
        A = str2double(N_str);
    elseif isempty(N_str)
        A = A_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
     % Get tuning parameter b
    N_str = input('fGSR steepness "b" parameter: ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif N_str == 'c'
        break;
    elseif isNumeric(N_str)
        b = str2double(N_str);
    elseif isempty(N_str)
        b = b_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
    % Get tuning parameter c
    N_str = input('fGSR width "c" parameter: ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif N_str == 'c'
        break;
    elseif isNumeric(N_str)
        c = str2double(N_str);
    elseif isempty(N_str)
        c = c_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
    % generate the mesh
    [s,x_mat,fGSR] = SOLVE.generate_mesh(N1,A,b,c);
    figure() % plot fGSR(s)
    subplot(2,1,1) 
    plot(s,fGSR,'-r');
    xlabel('Computational domain [s]');
    ylabel('Grid-stretching ratio [fGSR]');
    title(sprintf("X mesh for N = %d",N1));
    subplot(2,1,2) 
    plot(s,x_mat(:,1),'-r');
    xlabel('Computational domain [s]');
    ylabel('Physical domain [x]');
    SOLVE.display_mesh(x_mat(1,1), x_mat(end,1), -1, 1, x_mat(:,1), [-1, 1]);

    % wait for user input
    N_str = input('Continue? [y/n] : ','s');
    if N_str == 'y'
        % do nothing
    elseif N_str == 'n'
        fprintf("Trying a new mesh. . .\n");
        clear all
        close all
        continue;
    else
        fprintf("Invalid input. Try again. . .\n");
        clear all
        close all
        continue;
    end
    
    % Solve and plot Burger's equation for nu=0.00001
    tol = 0.001;                % tolerance
    solve = @(N,nu)SOLVE.solve_burgers(N,nu,A,b,c);
    
    nu = 0.00001;               % desired diffusivity
    N2 = floor(1.5*N1);                % second set of points (larger than N1)
    fprintf("Searching for solution to Burgers equation with diffusivity nu = %7.5f. . .\n",nu);
    
    % compute target discretization using N1 points
    [x1, urt1] = solve(N1,nu);
    dx1_avg = mean(x1(2:end)-x1(1:end-1));

    % compute target discretization using N2 points
    [x2, urt2] = solve(N2,nu);
    dx2_avg = mean(x2(2:end)-x2(1:end-1));

    % interpolate N1 and N2 solutions to common domain
    x_common = linspace(x1(1),x1(end),(N1+N2)/2);
    u1c = interp1(x1,urt1,x_common);
    u2c = interp1(x2,urt2,x_common);

    % approximate error in N1 discretization
    n = dx1_avg/dx2_avg;
    diff = u2c - u1c;
    
    % Richardson's Extrapolation (assuming 2nd-order accuracy)
    error = (n^2/(n^2-1))*(diff);
    norminf = max(abs(error));
    fprintf("Error at N=%d is estimated to be %6.4f m/s.\n",N1,norminf);
    
    % Check if we are below tolerance
    if (norminf < tol)
        fprintf("Solution converged at N=%d. Error: %6.4e < %6.4e\n\n",N1,norminf,tol);
        fprintf("Parameters: A=%5.3f, b=%d, c=%5.3f, N=%d\n",A,b,c,N1);
        
        % Estimate the thickness of shock region
        [~,left_idx] = min(abs(urt1-0.9));
        [~,right_idx] = min(abs(urt1+0.9));
        shock = x1(right_idx) - x1(left_idx);
        fprintf("Shock region: %6.4E < x < %6.4E. Thickness = %6.4E\n",x1(left_idx),x1(right_idx),shock);
        
        % wait for user input
        N_str = input('Save results? [y/n] : ','s');
        if N_str == 'y'
            % plot results
            figure()
            plot(x1,urt1,'or');
            xlabel('X coordinate [m]');
            ylabel("Function u [m/s]");
            title(sprintf('Solution to 1D Burgers equation; nu=%7.5f, N=%d',nu,N1));
            xlim([-1 1]);
            ylim([-1 1]);
            grid on
            save task1.mat x1 urt1 A b c N1 shock
            return;
        elseif N_str == 'n'
            fprintf("Continuing to next task. . .\n");
            close all
            clear all
            break;
        else
            fprintf("Invalid input. Closing. . .\n");
            return;
        end
    end
    
    % plot results
    figure()
    plot(x1,urt1,'or');
    hold on
    plot(x2,urt2,'-k');
    xlabel('X coordinate [m]');
    ylabel("Function u [m/s]");
    title(sprintf('Solution to 1D Burgers equation; nu=%7.5f',nu));
    legend(sprintf("N=%d",N1),sprintf("N=%d",N2));
    xlim([-1 1]);
    ylim([-1 1]);
    grid on
    
    % wait for user input
    N_str = input('Try again? [y/n] : ','s');
    if N_str == 'y'
        close all
        clear all
        continue;
    elseif N_str == 'n'
        fprintf("Continuing to next task. . .\n");
        break;
    else
        fprintf("Invalid input. Closing. . .\n");
        return;
    end
end

%% Task 2
fprintf("\nTASK2:\n");
fprintf("Enter 'q' at any time to quit.\n");
fprintf("Use decimal notation for fractions.\n");
fprintf("Press enter to accept default values.\n");

% Main user input loop
while (1)
    isNumeric = @(S) ~isnan(str2double(S));
    nu_DEFAULT = 0.00001;
    N_DEFAULT = 400;
    A_DEFAULT = -0.05;
    b_DEFAULT = 20;
    c_DEFAULT = 0.5;
    
    % Get diffusivity value nu
    fprintf("\nRunning Task 2: find lowest resolvable diffusivity\n");
    N_str = input('Diffusivity : ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif isNumeric(N_str)
        nu = str2double(N_str);
    elseif isempty(N_str)
        nu = nu_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
    % Get number of points N
    N_str = input('Number of points : ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif isNumeric(N_str)
        N1 = str2double(N_str);
    elseif isempty(N_str)
        N1 = N_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
    % Get tuning parameter A
    N_str = input('fGSR amplitutde "A" parameter: ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif isNumeric(N_str)
        A = str2double(N_str);
    elseif isempty(N_str)
        A = A_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
     % Get tuning parameter b
    N_str = input('fGSR steepness "b" parameter: ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif isNumeric(N_str)
        b = str2double(N_str);
    elseif isempty(N_str)
        b = b_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
    % Get tuning parameter c
    N_str = input('fGSR width "c" parameter: ','s');
    if N_str == 'q'
        fprintf("Quitting. . .\n");
        return;
    elseif isNumeric(N_str)
        c = str2double(N_str);
    elseif isempty(N_str)
        c = c_DEFAULT;
    else
        fprintf("Invalid input: %s\n",N_str);
        continue;
    end
    
    % generate the mesh
    [s,x_mat,fGSR] = SOLVE.generate_mesh(N1,A,b,c);
    figure() % plot fGSR(s)
    subplot(2,1,1) 
    plot(s,fGSR,'-r');
    xlabel('Computational domain [s]');
    ylabel('Grid-stretching ratio [fGSR]');
    title(sprintf("X mesh for N = %d",N1));
    subplot(2,1,2) 
    plot(s,x_mat(:,1),'-r');
    xlabel('Computational domain [s]');
    ylabel('Physical domain [x]');
    SOLVE.display_mesh(x_mat(1,1), x_mat(end,1), -1, 1, x_mat(:,1), [-1, 1]);

    % wait for user input
    N_str = input('Continue? [y/n] : ','s');
    if N_str == 'y'
        % do nothing
    elseif N_str == 'n'
        fprintf("Trying a new mesh. . .\n");
        clear all
        close all
        continue;
    else
        fprintf("Invalid input. Try again. . .\n");
        clear all
        close all
        continue;
    end
    
    % Solve and plot Burger's equation with chosen nu value
    tol = 0.001;                % tolerance
    solve = @(N,nu)SOLVE.solve_burgers(N,nu,A,b,c);
    N2 = floor(1.5*N1);                % second set of points (larger than N1)
    fprintf("Searching for solution to Burgers equation with diffusivity nu = %6.2E. . .\n",nu);
    
    % compute target discretization using N1 points
    [x1, urt1] = solve(N1,nu);
    dx1_avg = mean(x1(2:end)-x1(1:end-1));

    % compute target discretization using N2 points
    [x2, urt2] = solve(N2,nu);
    dx2_avg = mean(x2(2:end)-x2(1:end-1));

    % interpolate N1 and N2 solutions to common domain
    x_common = linspace(x1(1),x1(end),(N1+N2)/2);
    u1c = interp1(x1,urt1,x_common);
    u2c = interp1(x2,urt2,x_common);

    % approximate error in N1 discretization
    n = dx1_avg/dx2_avg;
    diff = u2c - u1c;
    
    % Richardson's Extrapolation (assuming 2nd-order accuracy)
    error = (n^2/(n^2-1))*(diff);
    norminf = max(abs(error));
    fprintf("Error at N=%d is estimated to be %6.4f m/s.\n",N1,norminf);
    
    % Check if we are below tolerance
    if (norminf < tol)
        fprintf("Solution converged at N=%d. Error: %6.4e < %6.4e\n\n",N1,norminf,tol);
        fprintf("Parameters: A=%5.3f, b=%d, c=%5.3f, N=%d\n, nu=%6.2E\n",A,b,c,N1,nu);
        
        % wait for user input
        N_str = input('Save results? [y/n] : ','s');
        if N_str == 'y'
            % plot results
            figure()
            plot(x1,urt1,'or');
            xlabel('X coordinate [m]');
            ylabel("Function u [m/s]");
            title(sprintf('Solution to 1D Burgers equation; nu=%6.2E, N=%d',nu,N1));
            xlim([-1 1]);
            ylim([-1 1]);
            grid on
            save task2.mat x1 urt1 A b c N1 nu
            return;
        elseif N_str == 'n'
            fprintf("Closing. . .\n");
            close all
            clear all
            return;
        else
            fprintf("Invalid input. Closing. . .\n");
            return;
        end
    end
    
    % plot results
    figure()
    plot(x1,urt1,'or');
    hold on
    plot(x2,urt2,'-k');
    xlabel('X coordinate [m]');
    ylabel("Function u [m/s]");
    title(sprintf('Solution to 1D Burgers equation; nu=%6.2E',nu));
    legend(sprintf("N=%d",N1),sprintf("N=%d",N2));
    xlim([-1 1]);
    ylim([-1 1]);
    grid on
    
    % wait for user input
    N_str = input('Try again? [y/n] : ','s');
    if N_str == 'y'
        close all
        clear all
        continue;
    elseif N_str == 'n'
        fprintf("Closing. . .\n");
        return;
    else
        fprintf("Invalid input. Closing. . .\n");
        return;
    end
end