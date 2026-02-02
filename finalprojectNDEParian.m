clc; clear; close all;

% Define the domain
% Alternatively, Uncomment the pre-set values instead of getting them from user.

n_x = 20; 
%n_x = str2double(inputdlg('please enter n_x'));
n_y = 20; 
%n_y = str2double(inputdlg('please enter n_y'));

a = 0; 
%a = str2double(inputdlg('please enter a'));
b = 1; 
%b = str2double(inputdlg('please enter b'));
c = 0; 
%c = str2double(inputdlg('please enter c'));
d = 1; 
%d = str2double(inputdlg('please enter d'));

% Get the coefficients b and q from the user
b1 = 2; 
%b1 = str2double(inputdlg('Please enter a positive real number as b1:'));

b2 = 2; 
%b2 = str2double(inputdlg('Please enter a positive real number as b2:'));

q = 5; 
%q = str2double(inputdlg('Please enter positive real number q:'));

% define the grid
x = linspace(a, b, n_x);
y = linspace(c, d, n_y);
[X, Y] = meshgrid(x, y);

% define h_x and h_y
h_x = (b - a)/(n_x - 1);
h_y = (d - c)/(n_y - 1);

% Initialize the boundary condition arrays
g_bottom = zeros(size(x));
g_top = zeros(size(x));

% Approximation and iteration variables
tol = 1e-6;
error = tol+1;
U = zeros(n_y,n_x);

% Choose Mode
M = str2double(inputdlg('Please enter mode choice 1 or 2'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mode 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if M == 1

    % Get the function u(x, y) from the user
    %u = inputdlg('Please enter a function as u: ', 's');
    %u = str2func(['@(x,y)',u{1}]);

    % A few test functions
    u = @(x, y) sin(pi * x) .* sin(pi * y);
    %u = @(x, y) x.^2 + y.^2;
    %u = @(x, y) sin(pi * x) .* cos(pi * y);
    %u = @(x, y) exp(-(x.^2 + y.^2));
    %u = @(x, y) x + y;
    %u = @(x, y) exp(-(x.^2 + y.^2)) .* sin(pi * x) .* cos(pi * y);

    U_exact = u(X, Y);
    U_old = zeros(n_x ,n_y );
    grad = costumgradient(U_exact, h_x, h_y);
    lap = costumlaplacian(U_exact, h_x, h_y);

    % calculate f
    f = -lap + b1* grad(:,:,1) + b2* grad(:,:,2) + q* U_exact;

    %calculate partial derivatives of u in direction of x and y
    [du_dx, du_dy] = symbolicDerivatives(u);
    g_bottom = arrayfun(@(x_) du_dy(x_, c), x);
    g_top    = arrayfun(@(x_) du_dy(x_, d), x);

    iter = 0;

    % Iteration Loop
    while error>tol
        
        U_old = U;
        % Interior points
        C = 2/h_x^2 + 2/h_y^2 + q;
        for i = 2:n_y-1
            for j = 2:n_x-1
                U(i,j) = ( ...
                (U(i,j+1) + U(i,j-1)) / h_x^2 + ...
                (U(i+1,j) + U(i-1,j)) / h_y^2 - ...
                (b1 / (2*h_x)) * (U(i,j+1) - U(i,j-1)) - ...
                (b2 / (2*h_y)) * (U(i+1,j) - U(i-1,j)) + ...
                f(i,j) ...
                ) / C;
            end
        end

        % Boundary conditions
        % Neumann for top and bottom
        for j = 2:n_x-1
            U(1, j) = U(2, j) - h_y * g_bottom(j);
            U(end, j) = U(end-1, j) + h_y * g_top(j);
        end

        % Dirichle for left and right
        U(:, 1) = u(a, y);
        U(:, end) = u(b, y);
   
        % Calculate the error (Frobenius norm) to check for convergence
   
        error = norm(U - U_old, 'fro');
        iter = iter + 1;
    end

    % Plot the solution
    figure;
    surf(X, Y, U);
    title('Approximate Solution');
    xlabel('x');
    ylabel('y');

    % Plot the exact solution
    figure;
    surf(X, Y, U_exact);
    title('Exact Solution');
    xlabel('x');
    ylabel('y');


    %calculate the condition number over mesh refinement 20-200
    % you can do more but it gets a bit slow:(

    errors = zeros (10,1);
    condition = zeros(10,1);
    h_values = zeros(10,1);
    h_mesh = zeros(10,1);

    for m=1:10
        
        % 1. Use updated mesh size
        n_x = m * 20;
        n_y = m * 20;

        % 2. Define mesh
        x = linspace(a, b, n_x);
        y = linspace(c, d, n_y);
        [X, Y] = meshgrid(x, y);

        h_x = (b - a) / (n_x - 1);
        h_y = (d - c) / (n_y - 1);
        h_max = max(h_x, h_y);
        h_values(m) = 1 / h_max^2;
        h_mesh(m) = h_max;

        % 3. Define boundary data
        U = zeros(n_y, n_x);  % reset solution for new mesh
        U_exact = u(X, Y);
        grad = costumgradient(U_exact, h_x, h_y);
        lap = costumlaplacian(U_exact, h_x, h_y);
        f = -lap + b1 * grad(:,:,1) + b2 * grad(:,:,2) + q * U_exact;

        [du_dx, du_dy] = symbolicDerivatives(u);
        g_bottom = arrayfun(@(x_) du_dy(x_, c), x);
        g_top = arrayfun(@(x_) du_dy(x_, d), x);
        
        % Iteration Loop
        error_iter = tol + 1;  % initialize

        while error_iter > tol
            % Interior points
            U_old = U;
            C = 2/h_x^2 + 2/h_y^2 + q;
            for i = 2:n_y-1
                for j = 2:n_x-1
                    U(i,j) = ( ...
                    (U(i,j+1) + U(i,j-1)) / h_x^2 + ...
                    (U(i+1,j) + U(i-1,j)) / h_y^2 - ...
                    (b1 / (2*h_x)) * (U(i,j+1) - U(i,j-1)) - ...
                    (b2 / (2*h_y)) * (U(i+1,j) - U(i-1,j)) + ...
                    f(i,j) ...
                    ) / C;
                end
            end

            % Boundary conditions
            % Neumann for top and bottom
            for j = 2:n_x-1
              U(1, j) = U(2, j) - h_y * g_bottom(j);
              U(end, j) = U(end-1, j) + h_y * g_top(j);
            end

            % Dirichle for left and right
            U(:, 1) = u(a, y);
            U(:, end) = u(b, y);

            % Calculate the error (Frobenius norm) to check for convergence
            
            % Compute iteration error
            error_iter = norm(U - U_old, 'fro');
        end

        
        errors(m) = norm(U - u(X,Y), 'fro') * h_max;
        A= buildMatrixA(n_x, n_y, h_x, h_y, b1, b2, q);
        condition(m) =  condest(A);
       
    end

    % Plot the condition number
    %figure;
    %loglog(h_values, condition, '-o');
    %title('Condition Number vs mesh refinement');
    %xlabel('1/h^2');
    %ylabel('Condition Number');

    %plot the error convergence
    %figure;
    %loglog(h_values, errors, '-o');
    %title('Error Convergence');
    %xlabel('1/h^2');
    %ylabel('Error')


% Plot error convergence vs h
figure;
loglog(h_values, errors, 'b-o', 'LineWidth', 1.5); % Error vs h
hold on;
loglog(h_values, h_values.^2, 'r--', 'LineWidth', 1.5); % Reference h^2
xlabel('h');
ylabel('||u - U||_\infty');
title('Error Convergence');
legend('||u - U||_\infty', 'h^2', 'Location', 'southwest');
grid on;

% Plot condition number vs h
figure;
loglog(h_values, condition, 'c-o', 'LineWidth', 1.5); % Condition number vs h
hold on;
loglog(h_values, 1./(h_values.^2), 'r--', 'LineWidth', 1.5); % Reference 1/h^2
xlabel('h');
ylabel('K(A)');
title('Condition Number vs Mesh Size');
legend('K(A)', '1/h^2', 'Location', 'northwest');
grid on;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mode 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if M == 2

    % Get the function
    %prompt = {'Enter f(x,y):'};
    %response = inputdlg(prompt);
    %f = str2func(['@(x,y)', vectorize(response{1})]);
    %f = @(x,y) sin(pi*x).*sin(pi*y);
    %f = @(x,y) x.^2 + y.^2;
    f = @(x,y) exp(-x.^2 - y.^2);
    %f = @(x,y) cos(pi*x).*cos(pi*y);
    %f = @(x,y) x.*y;
    %f = @(x,y) x.^3 - y.^3;
    %f = @(x,y) log(x + 1) + log(y + 1);
    %f = @(x,y) x.^2 - y.^2;
    %f = @(x,y) sin(2*pi*x).*cos(2*pi*y);
    %f = @(x,y) ones(size(x));
    %f = @(x,y) 2*exp(x + y);
    %f  = @(x,y) 2*sin(x).*cos(y) + 2*cos(x).*cos(y) - 2*sin(x).*sin(y) + 5*sin(x).*cos(y);
    %f  = @(x,y) 5*ones(size(x));



    F = f(X,Y);
    figure;
    surf(X, Y, F);
    xlabel('x');
    ylabel('y');
    zlabel('f(x,y)');
    title('Plot of the entered function');


    % Iteration Loop
    error = tol + 1;
    while error > tol

        U_old = U; % Store the old solution
        % Interior points
        C = 2/h_x^2 + 2/h_y^2 + q;
        for i = 2:n_y-1
            for j = 2:n_x-1
                U(i,j) = ( ...
                (U(i,j+1) + U(i,j-1)) / h_x^2 + ...
                (U(i+1,j) + U(i-1,j)) / h_y^2 - ...
                (b1 / (2*h_x)) * (U(i,j+1) - U(i,j-1)) - ...
                (b2 / (2*h_y)) * (U(i+1,j) - U(i-1,j)) + ...
                F(i,j) ...
                ) / C;
            end
        end

        % Left/right Dirichlet (default u = 0)
        U(:,1)   = 0;
        U(:,end) = 0;

        % Top/bottom Neumann (zero normal derivative)
        for j = 2:n_x-1
            U(1,j)   = U(2,j);       
            U(end,j) = U(end-1,j);   
        end

        % Compute the error (norm of the difference)
        error = norm(U - U_old, 'fro');
    end


    % Plot the approximate solution
    figure;
    surf(X, Y, U);
    xlabel('x');
    ylabel('y');
    zlabel('U');
    title('Approximate Solution');

    % Vectorize solution
    z = U(:);
    x_plot = X(:);
    y_plot = Y(:);
   

    % Try a quadratic surface: 'poly22' (adjust degree as needed)
    U_poly = fit([x_plot, y_plot], z, 'poly22');

    % Display the fitted formula
    disp(U_poly);

    U_poly = U_poly(X,Y);

    % Plot the approximate solution
    figure;
    surf(X, Y, U_poly);
    xlabel('x');
    ylabel('y');
    zlabel('U_poly');
    title('Polynomial Solution');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lap = costumlaplacian(U, h_x, h_y)
     [n_y, n_x] = size(U);
     lap = zeros(n_y, n_x);
     for i = 2:n_y-1
         for j = 2:n_x-1
              lap(i,j) = (U(i+1,j) - 2*U(i,j) + U(i-1,j)) / h_y^2 + (U(i,j+1) - 2*U(i,j) + U(i,j-1)) / h_x^2;
         end
     end
end
   
function grad = costumgradient(U, h_x, h_y)
     [n_y, n_x] = size(U);
     grad_x = zeros(n_y, n_x);
     grad_y = zeros(n_y, n_x);
     for i = 2:n_y-1
         for j = 2:n_x-1
             grad_x(i,j) = (U(i,j+1) - U(i,j-1)) / (2*h_x);
             grad_y(i,j) = (U(i+1,j) - U(i-1,j)) / (2*h_y);
         end
     end
     grad = cat(3, grad_x, grad_y);
end

function A = buildMatrixA(n_x, n_y, h_x, h_y, b1, b2, q)
    N = (n_x - 2) * (n_y - 2); 
    A = sparse(N, N);          
    idx = @(i, j) (j - 1) * (n_x - 2) + i;  % convert (i,j) to linear index

    for j = 1:(n_y - 2)
        for i = 1:(n_x - 2)
            k = idx(i, j);

            % Diagonal (center)
            A(k, k) = 2/h_x^2 + 2/h_y^2 + q;

            % West (i-1)
            if i > 1
                A(k, idx(i-1, j)) = -1/h_x^2 - b1/(2*h_x);
            end

            % East (i+1)
            if i < (n_x - 2)
                A(k, idx(i+1, j)) = -1/h_x^2 + b1/(2*h_x);
            end

            % South (j-1)
            if j > 1
                A(k, idx(i, j-1)) = -1/h_y^2 - b2/(2*h_y);
            end

            % North (j+1)
            if j < (n_y - 2)
                A(k, idx(i, j+1)) = -1/h_y^2 + b2/(2*h_y);
            end
        end
    end
end

function [du_dx_func, du_dy_func] = symbolicDerivatives(u)
    syms x y
    u_sym = u(x, y);
    du_dx_func = matlabFunction(diff(u_sym, x), 'Vars', [x, y]);
    du_dy_func = matlabFunction(diff(u_sym, y), 'Vars', [x, y]);
end

