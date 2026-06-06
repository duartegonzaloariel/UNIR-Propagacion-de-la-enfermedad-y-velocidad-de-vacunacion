clear; clc; close all;

% =========================================================
% Código 1: Definición del campo vectorial
% =========================================================
beta_param = 0.0022;
mu         = 0.0137;
theta      = 0.4477;
gamma      = 0.04;

sir_vacunacion_rhs = @(t, y) [
    -beta_param * y(1) * y(2) - mu * y(1) + mu * sum(y) - gamma * y(1);
     beta_param * y(1) * y(2) - (theta + mu) * y(2);
     theta * y(2) - mu * y(3) + gamma * y(1)
];
disp('Sistema de primer orden definido como y'' = sir_vacunacion_rhs(t, y).');

% =========================================================
% Código 2: Condiciones iniciales
% =========================================================
N_total = 7000;
I0      = 1;
R0      = 0;
S0      = N_total - I0 - R0;
y0      = [S0; I0; R0];

T_condiciones_iniciales = table(S0, I0, R0, ...
    'VariableNames', {'S0', 'I0', 'R0'});
disp(T_condiciones_iniciales);

% =========================================================
% Código 3: Discretización del intervalo
% =========================================================
t0              = 0;
tf              = 10;
n_subintervalos = 500;
h               = (tf - t0) / n_subintervalos;

% =========================================================
% Código 4: Integración mediante Euler explícito
% =========================================================
[y_euler, t_euler] = euler_vector(t0, y0, sir_vacunacion_rhs, tf, ...
    n_subintervalos);

S_euler = y_euler(:, 1);
I_euler = y_euler(:, 2);
R_euler = y_euler(:, 3);

% =========================================================
% Código 5: Resultados y tabla de valores finales (Euler)
% =========================================================
[I_euler_max, idx_pico_euler] = max(I_euler);
t_pico_euler = t_euler(idx_pico_euler);

T_euler_final = table(S_euler(end), I_euler(end), R_euler(end), ...
    sum(y_euler(end, :)), I_euler_max, t_pico_euler, ...
    'VariableNames', {'S_10', 'I_10', 'R_10', 'N_10', ...
                      'I_max', 't_I_max'});
disp(T_euler_final);

% =========================================================
% Código 6: Figura de evolución temporal (Euler)
% =========================================================
figure('Color', 'w');
plot(t_euler, S_euler, 'LineWidth', 1.4);
hold on;
plot(t_euler, I_euler, 'LineWidth', 1.4);
plot(t_euler, R_euler, 'LineWidth', 1.4);
grid on;
xlabel('Tiempo');
ylabel('Número de individuos');
title('Solución aproximada mediante Euler explícito');
legend('Susceptibles S(t)', 'Infectados I(t)', 'Recuperados R(t)', ...
    'Location', 'best');

% =========================================================
% Código 7: Comprobaciones (Euler)
% =========================================================
error_conservacion = max(abs(sum(y_euler, 2) - N_total));
tol_conservacion   = 1e-8;

assert(error_conservacion <= tol_conservacion, ...
    'Euler no conserva la población total dentro de la tolerancia.');
assert(all(y_euler(:) >= 0), ...
    'Euler produce alguna población negativa.');

fprintf('Máximo error de conservación con Euler: %.3e\n', error_conservacion);

% =========================================================
% Código 8: Integración mediante RK4
% =========================================================
[y_rk4, t_rk4] = rk4(t0, y0, sir_vacunacion_rhs, tf, n_subintervalos);

S_rk4 = y_rk4(:, 1);
I_rk4 = y_rk4(:, 2);
R_rk4 = y_rk4(:, 3);

% =========================================================
% Código 9: Resultados y tabla de valores finales (RK4)
% =========================================================
[I_rk4_max, idx_pico_rk4] = max(I_rk4);
t_pico_rk4 = t_rk4(idx_pico_rk4);

T_rk4_final = table(S_rk4(end), I_rk4(end), R_rk4(end), ...
    sum(y_rk4(end, :)), I_rk4_max, t_pico_rk4, ...
    'VariableNames', {'S_10', 'I_10', 'R_10', 'N_10', ...
                      'I_max', 't_I_max'});
disp(T_rk4_final);

% =========================================================
% Código 10: Figura de evolución temporal (RK4)
% =========================================================
figure('Color', 'w');
plot(t_rk4, S_rk4, 'LineWidth', 1.4);
hold on;
plot(t_rk4, I_rk4, 'LineWidth', 1.4);
plot(t_rk4, R_rk4, 'LineWidth', 1.4);
grid on;
xlabel('Tiempo');
ylabel('Número de individuos');
title('Solución aproximada mediante Runge-Kutta de orden 4');
legend('Susceptibles S(t)', 'Infectados I(t)', 'Recuperados R(t)', ...
    'Location', 'best');

% =========================================================
% Código 11: Comparación con ode45
% =========================================================
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_ode45, y_ode45] = ode45(sir_vacunacion_rhs, [t0 tf], y0, opts);

S_ode45_interp = interp1(t_ode45, y_ode45(:,1), t_euler);
I_ode45_interp = interp1(t_ode45, y_ode45(:,2), t_euler);
R_ode45_interp = interp1(t_ode45, y_ode45(:,3), t_euler);

err_euler_S = abs(S_euler(end) - S_ode45_interp(end));
err_euler_I = abs(I_euler(end) - I_ode45_interp(end));
err_euler_R = abs(R_euler(end) - R_ode45_interp(end));
err_rk4_S   = abs(S_rk4(end)   - S_ode45_interp(end));
err_rk4_I   = abs(I_rk4(end)   - I_ode45_interp(end));
err_rk4_R   = abs(R_rk4(end)   - R_ode45_interp(end));

T_comparacion = table(...
    [S_euler(end); S_rk4(end); S_ode45_interp(end)], ...
    [I_euler(end); I_rk4(end); I_ode45_interp(end)], ...
    [R_euler(end); R_rk4(end); R_ode45_interp(end)], ...
    [err_euler_S;  err_rk4_S;  0], ...
    [err_euler_I;  err_rk4_I;  0], ...
    [err_euler_R;  err_rk4_R;  0], ...
    'VariableNames', {'S_10','I_10','R_10','err_S','err_I','err_R'}, ...
    'RowNames', {'Euler','RK4','ode45'});
disp(T_comparacion);

figure('Color', 'w');
plot(t_euler, S_euler,      'b--', 'LineWidth', 1.4); hold on;
plot(t_rk4,   S_rk4,        'r-',  'LineWidth', 1.4);
plot(t_ode45, y_ode45(:,1), 'k:',  'LineWidth', 1.4);
plot(t_euler, I_euler,      'b--', 'LineWidth', 1.4);
plot(t_rk4,   I_rk4,        'r-',  'LineWidth', 1.4);
plot(t_ode45, y_ode45(:,2), 'k:',  'LineWidth', 1.4);
plot(t_euler, R_euler,      'b--', 'LineWidth', 1.4);
plot(t_rk4,   R_rk4,        'r-',  'LineWidth', 1.4);
plot(t_ode45, y_ode45(:,3), 'k:',  'LineWidth', 1.4);
grid on;
xlabel('Tiempo');
ylabel('Número de individuos');
title('Comparación Euler, RK4 y ode45');
legend('S Euler','S RK4','S ode45', ...
       'I Euler','I RK4','I ode45', ...
       'R Euler','R RK4','R ode45','Location','best');

% =========================================================
% Código 12: Análisis efecto tasa de vacunación gamma = 0.8
% =========================================================
gamma_alt = 0.8;

sir_alt = @(t, y) [
    -beta_param * y(1) * y(2) - mu * y(1) + mu * sum(y) - gamma_alt * y(1);
     beta_param * y(1) * y(2) - (theta + mu) * y(2);
     theta * y(2) - mu * y(3) + gamma_alt * y(1)
];

[y_alt, t_alt] = rk4(t0, y0, sir_alt, tf, n_subintervalos);
I_alt = y_alt(:, 2);
[I_alt_max, idx_alt] = max(I_alt);
t_pico_alt = t_alt(idx_alt);

T_gamma = table(...
    [0.04; 0.8], ...
    [I_rk4_max;  I_alt_max], ...
    [t_pico_rk4; t_pico_alt], ...
    [y_rk4(end,1); y_alt(end,1)], ...
    [y_rk4(end,2); y_alt(end,2)], ...
    [y_rk4(end,3); y_alt(end,3)], ...
    'VariableNames', {'gamma','I_max','t_I_max','S_10','I_10','R_10'}, ...
    'RowNames', {'gamma_004','gamma_08'});
disp(T_gamma);

figure('Color', 'w');
plot(t_rk4, y_rk4(:,2), 'r-',  'LineWidth', 1.8); hold on;
plot(t_alt,  y_alt(:,2), 'b--', 'LineWidth', 1.8);
grid on;
xlabel('Tiempo');
ylabel('Número de individuos');
title('Infectados I(t): \gamma=0.04 vs \gamma=0.8');
legend('I(t) \gamma=0.04', 'I(t) \gamma=0.8', 'Location', 'best');

% =========================================================
% Funciones locales
% =========================================================
function [y, t] = euler_vector(t0, y0, f, T, n)
validateattributes(t0, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(T,  {'numeric'}, {'scalar', 'real', 'finite', '>', t0});
validateattributes(n,  {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(y0, {'numeric'}, {'vector', 'real', 'finite'});
if ~isa(f, 'function_handle')
    error('euler_vector:funcionInvalida', ...
        'f debe ser un identificador de función.');
end
y0 = y0(:).';
p  = numel(y0);
t  = linspace(t0, T, n + 1).';
h  = (T - t0) / n;
y  = zeros(n + 1, p);
y(1, :) = y0;
for k = 1:n
    pendiente = f(t(k), y(k, :).');
    pendiente = pendiente(:).';
    if numel(pendiente) ~= p
        error('euler_vector:dimensionInvalida', ...
            'f debe devolver un vector con %d componentes.', p);
    end
    y(k + 1, :) = y(k, :) + h * pendiente;
end
end

function [y, t] = rk4(t0, y0, f, T, n)
validateattributes(t0, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(T,  {'numeric'}, {'scalar', 'real', 'finite', '>', t0});
validateattributes(n,  {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(y0, {'numeric'}, {'vector', 'real', 'finite'});
if ~isa(f, 'function_handle')
    error('rk4:funcionInvalida', 'f debe ser un identificador de función.');
end
y0 = y0(:);
p  = numel(y0);
t  = linspace(t0, T, n + 1).';
h  = (T - t0) / n;
y  = zeros(n + 1, p);
y(1, :) = y0.';
for k = 1:n
    estado = y(k, :).';
    k1 = comprobar_dimension(f(t(k),           estado),           p);
    k2 = comprobar_dimension(f(t(k) + h/2, estado + h/2 * k1),   p);
    k3 = comprobar_dimension(f(t(k) + h/2, estado + h/2 * k2),   p);
    k4 = comprobar_dimension(f(t(k) + h,   estado + h   * k3),   p);
    incremento  = (k1 + 2*k2 + 2*k3 + k4) / 6;
    y(k+1, :)   = (estado + h * incremento).';
end
end

function valor = comprobar_dimension(valor, p)
valor = valor(:);
if numel(valor) ~= p
    error('rk4:dimensionInvalida', ...
        'f debe devolver un vector con %d componentes.', p);
end
end