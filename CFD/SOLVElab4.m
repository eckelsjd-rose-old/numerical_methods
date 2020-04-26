%% SOLVE.m
% Class:       ME427 - Intro to CFD; functions for HW10
% Author:      Joshua Eckels
% Description: Class to provide helper functions for solving the 1D, steady
%              Burgers equation.
% Functions:   - solve_burgers
%              - residual_burgers
%              - generate_mesh
%              - state_residual_mesh
%              - state_deriv_mesh
%              - jacobian
%              - residual_deriv
%              - rk4 (runge-kutta 4)
%              - nr_1D (Newton-Raphson 1D)
%              - nr_MD (Newton-Raphson MD)
%              
% Calling:     SOLVE.function_name(arg1,arg2)
% Date:        4/11/2020

classdef SOLVE
    methods(Static)
        
        %% Estimate solution to 1D, steady, nonlinear, non-uniform mesh Burgers equation
        % using finite difference approximation
        % Parameters:   N   = number of discrete mesh points
        %               nu  = diffusivity
        % Return:       x   = domain mesh points in x-direction [m]
        %               urt = returned 'u' function values evaluated at x
        function [x, urt] = solve_burgers(N,nu,A,b,c)
            % Set boundary values
            umin = 1;
            umax = -1;
            % Generate the mesh for N points
            [~,x_mat,~] = SOLVE.generate_mesh(N+1,A,b,c);
            x = x_mat(:,1);
            % Define residual and Jacobian functions
            R = @(u)SOLVE.residual_burgers(u,x,nu,umin,umax);
            J = @(u)SOLVE.jacobian(u,R);

            % Perform Newton-Raphson method on initial guess of u0
            maxIter = 1000;
            toggle = 0;
            tol = 1*10^(-9);
            u0 = linspace(1,-1,N)';    % initial guess
            [urt, ~] = SOLVE.nr_MD(R,J,u0,tol,maxIter,toggle);
        end
        
        %% Computes a Residual vector for the nonlinear Burgers equation
        % using finite difference approximation
        % Parameters:
        %   u - vector of u values to compute residual at [u1;u2;u3;u4;u5...]
        %   x_mat - matrix of discrete x values and derivatives [x1...;x1'...;x1''...]
        %   nu- diffusivity
        %   umin - boundary condition 1
        %   umax - boundary condition 2
        % Calling:
        %   residual expected to be called as residual_burgers(u);
        %   x_mat,nu,umin,umax are provided
        function res = residual_burgers(u,x,nu,umin,umax)
            N = size(u,1);      % number of finite volume cells
            res = zeros(N,1);
            dx = x(2:end) - x(1:end-1);
            
            % residual evaluated at middle cells (variable mesh spacing)
            % interface values found with central interpolation
            % could replace these with higher order taylor approximations ->
            % boundary values
            res(1) = (1/2)*(umin^2-((dx(1)*u(2)+dx(2)*u(1))/(dx(2)+dx(1)))^2) ...
                -nu*((u(1)-umin)/(dx(1)/2) - (u(2)-u(1))/((dx(1)+dx(2))/2));
            for i = 2:N-1
                uin_half  = (dx(i)*u(i-1)+dx(i-1)*u(i))/(dx(i-1)+dx(i)); % ^
                uout_half = (dx(i)*u(i+1)+dx(i+1)*u(i))/(dx(i+1)+dx(i)); % ^
                duin_half = 2*(u(i)-u(i-1))/(dx(i-1)+dx(i));             % ^
                duout_half= 2*(u(i+1)-u(i))/(dx(i+1)+dx(i));             % ^
                res(i) = (1/2)*(uin_half^2-uout_half^2)-nu*(duin_half-duout_half);
            end
            res(N) = (1/2)*(((dx(N)*u(N-1)+dx(N-1)*u(N))/(dx(N)+dx(N-1)))^2 - umax^2) ...
                -nu*((u(N)-u(N-1))/((dx(N-1)+dx(N))/2)-(umax-u(N))/(dx(N)/2));
        end
        
        %% Compute a 1D non-uniform mesh using a grid-stretching ratio
        % Parameters:
        %   N       - number of mesh points
        % Outputs:
        %   s       - computational domain vector
        %   x_mat   - physical domain matrix
        %           - col 1: x(s), col 2: x'(s), col 3: x''(s)
        % Calling:
        %   @(N)generate_mesh(N);
        function [s,x_mat,fGSR] = generate_mesh(N,A,b,c)
            fprintf("Generating mesh with N=%d points . . .\n",N);
            % parameters
            xmin = -1;              % Lower x bound
            xmax = 1;               % Upper x bound
            param = [xmin,xmax,N,A,b,c];

            % Residual (vg=shooting parameter)
            R = @(vg)SOLVE.state_residual_mesh(vg,param);   
            dRdvg = @(vg)SOLVE.residual_deriv(vg,R);
            vg0 = 0.01;             % initial guess dx/ds(s=0)
            maxIter = 1000;         % Maximum iterations
            tol = 10^(-6);          % tolerance

            % Solve for shooting parameter with 1D Newton-Raphson
            [vgfin, ~] = SOLVE.nr_1D(R,dRdvg,vg0,tol,maxIter,0);

            % Use rk4 to compute the x-direction mesh from s=0 to s=1
            v0 = [xmin;vgfin];      % Initial conditions
            dvds = @(s,v_vec)SOLVE.state_deriv_mesh(s,v_vec,param);
            ds = 1/(N-1);
            [s, v_mat] = SOLVE.rk4(dvds, v0, ds, N-1);

            % Physical domain
            x_mat(:,1:2) = v_mat(:,1:2);
            
            % Grid-to-grid stretching ratio function (fGSR)
            f1 = (A/2)*(1-erf(b*(s-c)));
            f2 = (A/2)*(1-erf(b*(s-(1-c))));
            fGSR = f1+f2-A;
            
            % BVP definition: x''(s) - (1/ds)*fGSR*x'(s) = 0
            x_mat(:,3) = (1/ds)*fGSR.*x_mat(:,2);
            fprintf("Mesh complete.\n");
        end
        
        %% Computes the residual of grid-stretching second-order DE for vguess
        % Parameters:
        %   vg - independent variable for shooting parameter
        % Calling:
        %   residual expected to be called as R(vg); params are passed in
        function res = state_residual_mesh(vg,param)
            % parameters passed in
            xmin = param(1);        % initial condition for  x at s=0
            xmax = param(2);        % boundary condition for x at s=1
            N = param(3);           % number of points

            % Run rk4 to time-step state vector to s=1 using vguess
            v0 = [xmin;vg];         % Initial conditions
            dvds = @(s,v_vec)SOLVE.state_deriv_mesh(s,v_vec,param);
            ds = 1/(N-1);
            [s, v_mat] = SOLVE.rk4(dvds, v0, ds, N-1);

            % Get end results from rk4 (at edge of computational domain s=1)
            xend = v_mat(end,1);

            % Compute residual at boundary
            res = xend - xmax;
        end
        
        %% State derivative function for mesh BVP (fGSR mesh)
        % Parameters:
        %   s     - location value to evaluate dvds at
        %   v_vec - state vector to evaluate dvds at
        % Outputs:
        %   dvds  - the output vector of derivatives of state vector
        % Calling:
        %   @(s,v_vec)state_deriv_mesh(s, v_vec); params are passed in
        function dvds = state_deriv_mesh(s, v_vec,param)
            % parameters passed in
            N = param(3);   % number of points
            A = param(4);   % fGSR amplitude tuning            
            b = param(5);   % fGSR steepness tuning
            c = param(6);   % fGSR width tuning

            % Grid-to-grid stretching ratio function (fGSR)
            f1 = (A/2)*(1-erf(b*(s-c)));
            f2 = (A/2)*(1-erf(b*(s-(1-c))));
            fGSR = f1+f2-A;
            
            ds = 1/(N-1);
            n = size(v_vec,1); % v_vec comes in as column state vector
            dvds = zeros(n,1); % dvds will leave as column state vector
            dvds(1) = v_vec(2);
            dvds(2) = (1/ds)*fGSR*v_vec(2);
        end
        
        %% Computes an approximate Jacobian with 1% perturbation
        % Parameters:
        %   R - vector of residual functions; called as R(u)
        %   u - location to evaluate jacobian matrix at
        % Calling:
        %   jacobian expected to be called as jacobian(u); so R function handle must be
        %   passed in
        function ret = jacobian(u,R)
            N = size(u,1);
            ret = zeros(N);
            for i = 1:N
                u_new = u;
                u_new(i) = u(i)*1.01;
                ret(:,i) = (R(u_new) - R(u)) / (0.01*u(i));
            end
        end
        
        %% Computes an approximate derivative of a residual function
        % Parameters:
        %   R - function handle to residual function R(x)
        %   x - location to evaluate the derivative at
        % Calling:
        %   Expected to be called as dRdx(x); so R function handle must be
        %   passed in; 1% perturbation used
        function d = residual_deriv(x,R)
            d = (R(1.01*x) - R(x))/(0.01*x);
        end
        
        %% Runge-Kutta approximation with 4 slopes
        % Parameters:
        %   dv_fcn - function handle which evaluates derivative of state vector
        %   v0_fcn - vector with initial value of each state variable [nx1]
        %   dt     - time step 
        %   N      - number of time steps
        % Outputs:
        %   t_vec  - vector of time values [ (N+1) x 1 ]
        %   v_mat  - matrix with instantaneous value of each state variable 
        %            [ (N+1) x n ]
        % Calling:
        %   dv_fcn expected to be called as dv_fcn(t, v_vec), where v_vec is the
        %   current state vector
        function [t_vec, v_mat] = rk4(dv_fcn, v0_vec, dt, N)
            n = size(v0_vec,1);     % n = number of state variables
            t_vec = (0:dt:dt*N)';   % time column vector
            v_mat = zeros(N+1,n);   % columns for each state var; rows for time steps
            v_mat(1,:) = v0_vec';   % row=0; t=0; v = v0 (ICs)
            for i = 1:N
                k1 = dv_fcn(t_vec(i),v_mat(i,:)'); % d{v}/dt 1 (col vector)
                ystar1 = v_mat(i,:)' + (dt/2)*k1;

                k2 = dv_fcn(t_vec(i) + dt/2, ystar1); % d{v}/dt 2
                ystar2 = v_mat(i,:)' + (dt/2)*k2;

                k3 = dv_fcn(t_vec(i) + dt/2, ystar2); % d{v}/dt 3
                ystar3 = v_mat(i,:)' + dt*k3;

                k4 = dv_fcn(t_vec(i) + dt, ystar3); % d{v}/dt 4
                v_mat(i+1,:) = v_mat(i,:) + dt*(k1' + 2*k2' + 2*k3' + k4')/6;
            end
        end
        
        %% One-dimensional Newton-Raphson method
        % Parameters:
        %   R - residual function handle; R(x)
        %   dRdx - handle for derivative of residual function; dRdx(x)
        %   xi - initial guess
        %   tol - tolerance
        %   maxIter - number of iterations maximum
        %   toggle - 0; no print info. 1; print info
        function [soln, n_iter] = nr_1D(R, dRdx, xi, tol, maxIter, toggle)
        % correction factor
        corr = abs(R(xi)/dRdx(xi));
        xi_new = xi;
        n = 0; % number of iterations
        if (toggle == 1)
            fprintf('%5s %7s %9s %9s %9s\n','Count','xi','R(xi)','dRdx(xi)','corr');
        end

        while (abs(corr) > tol)
            n = n + 1;
            Rxi = R(xi);
            dRdxi = dRdx(xi);
            corr = -Rxi/dRdxi;
            xi_new = xi + corr;
            if (toggle == 1)
                fprintf('%5.0f %7.3f %9.3f %9.3f %9.3f\n', n, xi, Rxi, dRdxi, corr);
            end
            xi = xi_new;
            if (n >= maxIter)
                break;
            end
        end
        
        % found solution (or exceeded maxIter)
        soln = xi_new;
        n_iter = n;
        end
        
        %% Multi-dimensional Newton-Raphson method
        % Parameters:
        %   R - vector of residual functions; called as R(u)
        %   J - approximate Jacobian function; called as J(u)
        %   xi - initial unknown parameter vector
        %   tol - tolerance level to meet
        %   maxIter - maximum number of iterations to perform
        %   toggle - 0: hide print info, 1: display print info
        % Returns:
        %   xrt - column vector of final parameters
        %   n_iter - number of iterations it took to converge
        function [xrt, n_iter] = nr_MD(R, J, xi, tol, maxIter, toggle)
            N = size(xi,1);         % size of input residual vector
            corr = -J(xi)\R(xi);    % correction factor
            xi_new = xi;
            n = 0;                  % number of iterations
            if (toggle == 1)
                fprintf('%9s %12s   %1s\n','Iteration','Max Error','Current root approximation');
            end

            while (max(abs(corr)) > tol)
                n = n + 1;
                Rxi = R(xi);
                Jxi = J(xi);
                corr = -Jxi\Rxi;
                xi_new = xi + corr;
                if (toggle == 1)
                    g = sprintf('%9.3f', xi);
                    fprintf('%6i    %12.2e %s\n', n, max(abs(corr)), g);
                end
                xi = xi_new;
                if (n >= maxIter)
                    break;
                end
            end
            
            % found solution (or exceeded maxIter)
            xrt = xi_new;
            n_iter = n;
            if toggle == 1
                 fprintf("Error from %s at N=%d is estimated to be %6.4e m/s.\n\n","Newton-Raphson",N,max(abs(corr)));
            end
        end
        
        %% Plot an xy-plane structured mesh
        % Code provided by Dr. Lui
        % Parameters:
        %   x_min - minimum x boundary
        %   x_max - maximum y boundary
        %   y_min - minimum y boundary
        %   y_max - maximum y boundary
        %   x_vec - vector containing all x mesh points
        %   y_vec - vector containing all y mesh points
        % Outputs:
        %   Plots the resulting mesh on the xy-plane
        function [] = display_mesh(x_min, x_max, y_min, y_max, x_vec, y_vec)
            figure()
            hold on;
            Nx = length(x_vec);
            Ny = length(y_vec);

            for j = 1:Ny
                x_horizontal(j,:) = [x_min x_max];
                y_horizontal(j,:) = [y_vec(j) y_vec(j)];
                plot(x_horizontal(j,:), y_horizontal(j,:));
            end

            for i = 1:Nx
                x_vertical(i,:) = [x_vec(i) x_vec(i)];
                y_vertical(i,:) = [y_min y_max];
                plot(x_vertical(i,:), y_vertical(i,:));
            end

            hold off;
            axis([x_min x_max -5 5]);
            xlabel('x');
            ylabel('y');
            title(sprintf('Structured 1D Mesh (Nx = %3i)', Nx));
        end
        
    end
end