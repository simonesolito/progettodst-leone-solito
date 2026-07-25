clear; close all; clc;

% --- Parametri ---
mu = 3.986e14; a = 7000e3; e = 0.1; Ts = 10;
omega_0 = sqrt(mu/a^3); Tmax = 2*pi/omega_0;
T = 0 : Ts : Tmax; Nsamp = length(T); N = 8; m0 = 500; u_max = 1;

% --- Target ---
R = 100; golden_angle = pi*(3-sqrt(5)); posizione = zeros(3, N);
for i = 1:N
    z = 1 - ((i-1)/(N-1))*2; phi = asin(z); theta = mod((i-1)*golden_angle, 2*pi);
    posizione(:,i) = R * [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
end
zeta_r = [posizione; zeros(3,N)];

% --- Matrici T-H ---
M = omega_0*T; E_ecc = M; 
for i=1:10, E_ecc = E_ecc - (E_ecc - e*sin(E_ecc) - M)./(1 - e*cos(E_ecc)); end
f = 2*atan2(sqrt(1+e)*sin(E_ecc/2), sqrt(1-e)*cos(E_ecc/2));
df = omega_0*(1+e*cos(f)).^2 ./ (1-e^2)^(3/2);
d2f = -2*omega_0^2*e*sin(f).*(1+e*cos(f)).^3 ./ (1-e^2)^3;
rt = a*(1-e^2)./(1+e*cos(f));
Adisc = zeros(6,6,Nsamp); Bdisc = zeros(6,3,Nsamp); Bcont = 1/m0 * [zeros(3,3); eye(3)];

% --- DISCRETIZZAZIONE ---
for i = 1:Nsamp
    A1 = [df(i)^2 + 2*mu/rt(i)^3, d2f(i), 0; -d2f(i), df(i)^2 - mu/rt(i)^3, 0; 0, 0, -mu/rt(i)^3];
    A2 = [0, 2*df(i), 0; -2*df(i), 0, 0; 0, 0, 0];
    Acont = [zeros(3,3), eye(3); A1, A2];
    Mexp = expm([Acont, Bcont; zeros(3,9)]*Ts);
    Adisc(:,:,i) = Mexp(1:6,1:6); Bdisc(:,:,i) = Mexp(1:6,7:9);
end

% --- Esecuzione ---
zeta_0 = [posizione + randn(3,N)*5; zeros(3,N)];

disp('=== Avvio Risoluzione Controllo Ottimo PMP ===');
[u_opt, T_opt, m_fin, lambda_opt, Traj] = shooting_method(zeta_0, zeta_r, Adisc, Bdisc, u_max, Tmax, m0);

% --- Plot ---
figure('Color','w'); hold on; grid on; view(3); colors = lines(N);
for i=1:N
    plot3(squeeze(Traj(2,:,i)), squeeze(Traj(1,:,i)), squeeze(Traj(3,:,i)), 'LineWidth', 1.5, 'Color', colors(i,:));
    plot3(zeta_r(2,i), zeta_r(1,i), zeta_r(3,i), 'ko', 'MarkerFaceColor', colors(i,:));
end
title('Traiettorie di Assemblaggio'); xlabel('Y [m]'); ylabel('X [m]'); zlabel('Z [m]');
pos_fin_calc = [Traj(:,end,1), Traj(:,end,2), Traj(:,end,3), Traj(:,end,4), Traj(:,end,5), Traj(:,end,6), Traj(:,end,7), Traj(:,end,8)];



% --- FUNZIONI PER L'IMPOLEMENTAZIONE DELL'ALGORITMO DI SHOOTING METHOD ---
function [u_out, T_out, m_out, lambda_out, Traj_out] = shooting_method(Z0, Zr, Ak, Bk, usat, Tmax, m0)
N = size(Z0,2); p = size(Ak,3)-1; Ts = Tmax/p; I_sp = 300; g = 9.81;

u_out = zeros(3,p,N); Traj_out = zeros(6,p+1,N); m_out = zeros(1,N); lambda_out = zeros(6,N);

T_out = Tmax; 

% Normalmente le posizioni variano più delle velocità, controbilanciamo pesando maggiormente le velocità
W = diag([1,1,1,10,10,10]);

for i = 1:N
    z0 = Z0(:,i); zr = Zr(:,i); lambda = zeros(6,1);

    % Guess Iniziale Euristico
    lambda(4:6) = -(zr(1:3)-z0(1:3))/norm(zr(1:3)-z0(1:3)+1e-9)*1e-4; 
    mu_reg = 10;

    for iter = 1:50
        [zp, ~, u_seq, ~, z_hist] = forward_sweep(z0, lambda, m0, p, Ts, Ak, Bk, usat, I_sp, g, m0);
        ei = zp - zr;

        if norm(W*ei) < 1.0, break; end

        J = zeros(6,6);
        for j=1:6
            lp = lambda; lp(j) = lp(j) + 1e-6;
            [zpp, ~, ~, ~, ~] = forward_sweep(z0, lp, m0, p, Ts, Ak, Bk, usat, I_sp, g, m0);
            J(:,j) = (zpp - zp)/1e-6;
            
        end
        lambda = lambda - ((W*J)'*(W*J) + mu_reg*eye(6)) \ ((W*J)'*(W*ei));
    end

    u_out(:,:,i) = u_seq; Traj_out(:,:,i) = z_hist; lambda_out(:,i) = lambda;
    
    m_out(i) = m0 - sum(vecnorm(u_seq))/(I_sp*g)*Ts;
end
end

function [z, lk, usq, ms, zhist] = forward_sweep(z0, l0, m0, p, Ts, Ak, Bk, usat, Isp, g, minit)
z = z0; lk = l0; ms = m0; usq = zeros(3,p); zhist = zeros(6,p+1); zhist(:,1) = z;

for k=1:p
    lk = (Ak(:,:,k)') \ lk;
    u = -usat * (Bk(:,:,k)'*lk) / (norm(Bk(:,:,k)'*lk) + 1e-3);
    z = Ak(:,:,k)*z + Bk(:,:,k)*(minit/ms)*u;
    ms = ms - (norm(u)/(Isp*g))*Ts;
    usq(:,k) = u; zhist(:,k+1) = z;
end
end