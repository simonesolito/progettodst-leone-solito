clear; close all; clc;

%% --- Parametri Costanti ---
mu = 3.986e14;
a = 7000e3;
e = 0.1;
omega_0 = sqrt(mu/a^3);
T_orbit = 2 * pi / omega_0;

% Limiti di ricerca per il tempo di manovra (Outer Loop)
T_min = 0.5 * T_orbit;
T_max = 1.5 * T_orbit;

N = 5;
m0 = 200;
u_max = 1;
p = 800; % Risoluzione temporale bilanciata per l'ottimizzazione

% Pesi della funzione costo globale (Energia + Penalità Temporale)
beta = 0.01; 

%% --- Generazione Target (Sfera) ---
R = 100;
golden_angle = pi*(3-sqrt(5));
posizione = zeros(3, N);
for i = 1:N
    if N == 1
        z = 1;
    else
        z = 1 - ((i-1)/(N-1))*2;
    end
    phi = asin(z);
    theta = mod((i-1)*golden_angle, 2*pi);
    posizione(:,i) = R * [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
end
zeta_r = [posizione; zeros(3,N)];

%% --- Condizioni Iniziali ---
rng(42);
zeta_0 = [posizione + randn(3,N)*100; 0.02 * randn(3,N)];

% Matrice dei pesi per il controllo del terminale nell'inner loop
W = diag([10, 10, 10, 100, 100, 100]);

%% OUTER LOOP: OTTIMIZZAZIONE CONTINUA DEL TEMPO (fminbnd)
 
fprintf('=== Avvio Outer Loop: Ricerca del Tempo Ottimo (T_min = %.1f min, T_max = %.1f min) ===\n', ...
    T_min/60, T_max/60);

% Funzione obiettivo anonima per fminbnd basata sull'indice di costo del paper[cite: 7]
cost_function = @(T_cand) evaluate_maneuver_time(T_cand, zeta_0, zeta_r, a, e, mu, omega_0, m0, u_max, W, p, beta);

options = optimset('Display', 'iter', 'TolX', 1e-2, 'MaxIter', 30);
[T_opt, cost_opt] = fminbnd(cost_function, T_min, T_max, options);

fprintf('\n>>> TEMPO OTTIMO DI MANOVRA TROVATO: T* = %.2f s (%.2f minuti) <<<\n', T_opt, T_opt/60);

%% ESECUZIONE FINALE CON IL TEMPO OTTIMO T_OPT

fprintf('=== Esecuzione Simulazione con T* Ottimizzato ===\n');
[~, Traj_opt, u_opt, solver_info] = evaluate_maneuver_time(T_opt, zeta_0, zeta_r, a, e, mu, omega_0, m0, u_max, W, p, beta);

%% STAMPA ERRORE FINALE E STATISTICHE

fprintf('\n--- RISULTATI FINALI ---\n');
for i = 1:N
    e_pos = norm(Traj_opt(1:3, end, i) - zeta_r(1:3, i));
    e_vel = norm(Traj_opt(4:6, end, i) - zeta_r(4:6, i));
    fprintf('Sat %d [%s] -> Errore Posizione: %.6g m | Errore Velocita: %.6g m/s | Propellente: %.6g kg\n', ...
        i, solver_info(i).status, e_pos, e_vel, solver_info(i).fuel_kg);
end

%% PLOT 3D

figure('Color','w', 'Name', 'Traiettorie Ottimizzate a Tempo Libero');
hold on; grid on; view(3);
colors = lines(N);
[sx, sy, sz] = sphere(40);
surf(sx*R, sy*R, sz*R, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
for i = 1:N
    h1 = plot3(squeeze(Traj_opt(2,:,i)), squeeze(Traj_opt(1,:,i)), squeeze(Traj_opt(3,:,i)), 'LineWidth', 1.8, 'Color', colors(i,:));
    h2 = plot3(zeta_r(2,i), zeta_r(1,i), zeta_r(3,i), 'p', 'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:));
    h3 = plot3(zeta_0(2,i), zeta_0(1,i), zeta_0(3,i), 's', 'MarkerSize', 7, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:));
    h4 = plot3(Traj_opt(2,end,i), Traj_opt(1,end,i), Traj_opt(3,end,i), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
end
title(sprintf('Assemblaggio a Tempo Ottimo (T* = %.1f min)', T_opt/60));
xlabel('Y [m]'); ylabel('X [m]'); zlabel('Z [m]');
axis equal;
legend([h1, h3, h2, h4], {'Traiettoria', 'Inizio', 'Target (Desiderato)', 'Arrivo Reale'}, 'Location', 'best');


%% FUNZIONE DI VALUTAZIONE (INTERFACCIA OUTER-INNER LOOP)
 
function [total_cost, Traj_out, u_out, info] = evaluate_maneuver_time(T_cand, Z0, Zr, a, e, mu, omega_0, m0, u_max, W, p, beta)
    Ts = T_cand / p;
    T_vec = linspace(0, T_cand, p+1);
    
    M = omega_0 * (T_vec(1:p) + Ts/2);
    E_ecc = M;
    for i=1:10
        E_ecc = E_ecc - (E_ecc - e*sin(E_ecc) - M)./(1 - e*cos(E_ecc));
    end
    f = 2*atan2(sqrt(1+e)*sin(E_ecc/2), sqrt(1-e)*cos(E_ecc/2));
    df = omega_0*(1+e*cos(f)).^2 ./ (1-e^2)^(3/2);
    d2f = -2*omega_0^2*e*sin(f).*(1+e*cos(f)).^3 ./ (1-e^2)^3;
    rt = a*(1-e^2)./(1+e*cos(f));
    
    Adisc = zeros(6,6,p);
    Bdisc = zeros(6,3,p);
    Bcont = 1/m0 * [zeros(3,3); eye(3)];
    
    for i = 1:p
        A1 = [df(i)^2 + 2*mu/rt(i)^3, d2f(i), 0; -d2f(i), df(i)^2 - mu/rt(i)^3, 0; 0, 0, -mu/rt(i)^3];
        A2 = [0, 2*df(i), 0; -2*df(i), 0, 0; 0, 0, 0];
        Acont = [zeros(3,3), eye(3); A1, A2];
        Mexp = expm([Acont, Bcont; zeros(3,9)] * Ts);
        Adisc(:,:,i) = Mexp(1:6,1:6);
        Bdisc(:,:,i) = Mexp(1:6,7:9);
    end
    
    % Esecuzione dell'Inner Loop blindato
    [u_out, ~, Traj_out, info] = shooting_method(Z0, Zr, Adisc, Bdisc, u_max, T_cand, m0, W, p);
    
    % Calcolo dell'indice di costo globale del paper (Eq. 18):
    % J = (1/N) * sum_i sum_k ||u_i,k|| * Ts + beta * T[cite: 7]
    N = size(Z0, 2);
    sum_energy = 0;
    max_terminal_error = 0;
    
    for i = 1:N
        u_sat = u_out(:,:,i);
        norms_u = vecnorm(u_sat, 2, 1);
        sum_energy = sum_energy + sum(norms_u) * Ts;
        
        % Termine di penalità per eventuali fallimenti di convergenza nell'inner loop
        err_p = info(i).position_error;
        err_v = info(i).velocity_error;
        max_terminal_error = max(max_terminal_error, err_p + 100*err_v);
        
        if ~info(i).converged
            sum_energy = sum_energy + 1e6; % Penalità drastica se non converge
        end
    end
    
    energy_index = sum_energy / N;
    total_cost = energy_index + beta * T_cand + 1e5 * max_terminal_error;
end

%% INNER LOOP: SHOOTING CON REGOLARIZZAZIONE (INVARIATO E BLINDATO)

function [u_out, lambda_out, Traj_out, info] = shooting_method(Z0, Zr, Ak, Bk, usat, Tmax, m0, W, p)
    N = size(Z0,2);
    Ts = Tmax/p;
    I_sp = 300;
    g = 9.81;
    u_out = zeros(3,p,N);
    Traj_out = zeros(6,p+1,N);
    lambda_out = zeros(6,N);
    info = struct([]);
 
    eps_values = [0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001];
    max_iter = 100;
    tol_pos = 1e-5;
    tol_vel = 1e-8;
    InvAkT = zeros(6,6,p);
    for k = 1:p
        InvAkT(:,:,k) = Ak(:,:,k)' \ eye(6);
    end
 
    for i = 1:N
        z0 = Z0(:,i);
        zr = Zr(:,i);
        [lambda, lambda_scale, residual_scale] = get_gramian_lambda0(z0, zr, Ak, Bk, p);
 
        max_S = 0;
        lk_temp = lambda;
        for k = 1:p
            lk_temp = InvAkT(:,:,k)*lk_temp;
            max_S = max(max_S, norm(Bk(:,:,k)'*lk_temp));
        end
        if max_S > 0
            scale = fminbnd(@(alpha) norm(residual_scale \ ...
                (forward_sweep(z0, alpha*lambda, m0, p, Ts, Ak, Bk, usat, ...
                I_sp, g, InvAkT, eps_values(1))-zr)), ...
                0, 2/max_S, optimset('TolX',1e-10,'Display','off'));
            lambda = scale*lambda;
        end
        total_iter = 0;
 
        for eps_smooth = eps_values
            mu_reg = 1e-3;
            [zp, ~, u_seq, ms, z_hist] = forward_sweep(z0, lambda, m0, p, Ts, ...
                Ak, Bk, usat, I_sp, g, InvAkT, eps_smooth);
            ei = zp-zr;
            for iter = 1:max_iter
                total_iter = total_iter+1;
                if norm(ei(1:3)) < tol_pos && norm(ei(4:6)) < tol_vel, break; end
 
                q = lambda_scale \ lambda;
                J = zeros(6,6);
                for j = 1:6
                    delta_l = 1e-6*max(1,abs(q(j)));
                    lp = lambda+delta_l*lambda_scale(:,j);
                    lm = lambda-delta_l*lambda_scale(:,j);
                    zpp = forward_sweep(z0, lp, m0, p, Ts, Ak, Bk, usat, ...
                        I_sp, g, InvAkT, eps_smooth);
                    zpm = forward_sweep(z0, lm, m0, p, Ts, Ak, Bk, usat, ...
                        I_sp, g, InvAkT, eps_smooth);
                    J(:,j) = (zpp-zpm)/(2*delta_l);
                end
 
                JW = residual_scale \ J;
                residual = residual_scale \ ei;
                grad = JW'*residual;
                H = JW'*JW;
                A_lm = H+mu_reg*diag(max(diag(H),1e-12));
                if rcond(A_lm) < 1e-14
                    dlambda = -pinv(A_lm)*grad;
                else
                    dlambda = -(A_lm \ grad);
                end
                dlambda = lambda_scale*dlambda;
 
                accepted = false;
                if all(isfinite(dlambda))
                    for backtrack = 0:19
                        lambda_cand = lambda+(0.5^backtrack)*dlambda;
                        [zp_cand, ~, u_cand, ms_cand, z_cand] = forward_sweep( ...
                            z0, lambda_cand, m0, p, Ts, Ak, Bk, usat, ...
                            I_sp, g, InvAkT, eps_smooth);
                        ei_cand = zp_cand-zr;
                        if norm(residual_scale \ ei_cand) < norm(residual)
                            lambda = lambda_cand;
                            ei = ei_cand;
                            u_seq = u_cand;
                            ms = ms_cand;
                            z_hist = z_cand;
                            mu_reg = max(1e-10,mu_reg/3);
                            accepted = true;
                            break;
                        end
                    end
                end
                if ~accepted
                    mu_reg = min(1e12,mu_reg*10);
                    if mu_reg >= 1e12, break; end
                end
            end
            phase_converged = norm(ei(1:3)) < tol_pos && norm(ei(4:6)) < tol_vel;
            if ~phase_converged, break; end
        end
 
        converged = phase_converged && eps_smooth == eps_values(end) && norm(W*ei) < 1;
        u_out(:,:,i) = u_seq;
        Traj_out(:,:,i) = z_hist;
        lambda_out(:,i) = lambda;
        info(i).converged = converged;
        if converged
            info(i).status = 'CONVERGENTE';
        else
            info(i).status = 'NON CONVERGENTE';
        end
        info(i).position_error = norm(ei(1:3));
        info(i).velocity_error = norm(ei(4:6));
        info(i).weighted_error = norm(W*ei);
        info(i).fuel_kg = m0-ms;
        info(i).final_mass = ms;
        info(i).iterations = total_iter;
        info(i).epsilon = eps_smooth;
    end
end
 
function [z, lk, usq, ms, zhist] = forward_sweep(z0, l0, m0, p, Ts, Ak, Bk, usat, Isp, g, InvAkT, eps_smooth)
    z = z0; lk = l0; ms = m0;
    usq = zeros(3,p);
    zhist = zeros(6,p+1);
    zhist(:,1) = z;
    for k = 1:p
        lk = InvAkT(:,:,k)*lk;
        B_eff = Bk(:,:,k)*(m0/ms);
        S = B_eff'*lk;
        normS = norm(S);
        factor = 0.5*(1+tanh((normS-1)/eps_smooth));
        u = -usat*(S/max(normS,1e-12))*factor;
 
        z = Ak(:,:,k)*z+B_eff*u;
        ms = ms-(norm(u)/(Isp*g))*Ts;
        if ~isfinite(ms) || ms <= 0
            error('Shooting:InvalidMass','Massa non positiva durante il forward sweep.');
        end
        usq(:,k) = u;
        zhist(:,k+1) = z;
    end
end
 
function [lambda0, lambda_scale, residual_scale] = get_gramian_lambda0(z0, zr, Ak, Bk, p)
    Q = zeros(6,3,p);
    Phi_p0 = eye(6);
    for k = p:-1:1
        Q(:,:,k) = Phi_p0*Bk(:,:,k);
        Phi_p0 = Phi_p0*Ak(:,:,k);
    end
    Q = reshape(Q,6,3*p);
    G = Q*Q';
    [residual_scale, chol_status] = chol((G+G')/2,'lower');
    if chol_status ~= 0
        error('Shooting:Gramian','Gramiano non definito positivo: verificare il modello.');
    end
    lambda_scale = Phi_p0'/residual_scale';
    lambda0 = lambda_scale*(residual_scale \ (Phi_p0*z0-zr));
end