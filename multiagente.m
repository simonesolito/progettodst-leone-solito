clear
close all

%% inizializzazione componenti del sistema

% Parametri orbitali
mu = 3.986e14;                                  % parametro gravitazionale standard (Terra)
a = 7000e3;                                     % lunghezza del semiasse maggiore (7000 km)
e = 0.1;                                        % Orbita ellittica
Ts = 10;                                        % Campionamento ogni 10 secondi
Tmax = 2 * pi * sqrt(a^3/mu);                   % durata dell'orbita

T = 0 : Ts : Tmax;                              % asse dei tempi
Nsamp = length(T);

omega_0 = sqrt(mu/a^3);                         % velocità angolare media del satellite
M = omega_0 * T;                                % anomalia media
m = 500;                                        % massa satellite

N = 8;                                          % numero di satelliti

u_max = 1;                                      % Limite di saturazione dell'input (N)

% Definizioone della matrice di adiacenza (collegamento ad anello ogni
% elemento comunica con i più vicini)

adjacency_matrix = zeros(N,N);
for i = 1:N-1
    adjacency_matrix(i, i+1) = 1;
    adjacency_matrix(i+1, i) = 1;
end
adjacency_matrix(N, 1) = 1; % Chiude l'anello
adjacency_matrix(1, N) = 1;

% 2. Calcolo del Grado (D) 
d = sum(adjacency_matrix, 2); 
D = diag(d);

% Calcolo della matrice Laplaciana
L = D - adjacency_matrix;

% Definizione della matrice di Pinning
G = eye(N);

% Calcolo della matrice di guadagno di accoppiamento
Gamma = (eye(N) + D + G) \ (L + G);

%% 3. Posizioni Desiderate e Parametri di Controllo  FUNZIONA SOLO PER 8 SATELLITI!!!!!!!!
R_target = 5; % Raggio della formazione sferica finale (metri)
Target_point = zeros(6, N);
pos_des = zeros(3, N);
 
golden_angle = pi * (3 - sqrt(5));
 
for i = 1:N
    % Calcolo di phi (Latitudine) e theta (Azimut)
    z_fib = 1 - ((i - 1) / (N - 1)) * 2;
    phi_fib = asin(z_fib);
    theta_fib = mod((i - 1) * golden_angle, 2*pi);
    % Assegnazione coordinate finali
    pos_des(1, i) = R_target * cos(phi_fib) * cos(theta_fib);   % X
    pos_des(2, i) = R_target * cos(phi_fib) * sin(theta_fib);   % Y
    pos_des(3, i) = R_target * sin(phi_fib);                    % Z
end
 
Target_point(1:3, :) = pos_des; % Le velocità di riferimento relative sono 0
 
% Matrici di peso e Guadagno di Accoppiamento
Q = diag([1, 1, 1, 1e3, 1e3, 1e3]);
R = diag([1e5, 1e5, 1e5]);
eta = 4;




%% calcolo dell'anomalia eccentrica E partendo dalla equazione M = E - e * sin(E) Newton-Raphson

E = M; 
for i = 1:10 
    E = E - (E - e.*sin(E) - M) ./ (1 - e.*cos(E));
end

% calcolo di posizione, velocità e accelerazione angolare
f = 2 * atan2(sqrt(1 + e) * sin(E/2) , sqrt(1 - e) * cos(E/2));
df = omega_0 * (1 + e * cos(f)).^2/(1 - e^2)^(3/2);
d2f = -2* omega_0^2 * e * sin(f) .* (1 + e * cos(f)).^3/(1 - e^2)^3;

rt = a * (1 - e^2)./(1 + e * cos(f));               % distanza satellite terra 

A1 = zeros(3,3,Nsamp);                          % LA PARTE CONINUA LA POTREI MEMORIZZARE MOMENTANEAMENTE
A2 = A1;                                            % RIDUCENDO LE MATRICI A DUE DIMENSIONI
Acont = zeros(6,6,Nsamp);
Bcont = 1/m * [zeros(3,3) ; eye(3)];

Adisc = Acont;
Bdisc = Bcont;


for i = 1:Nsamp 
    % calcolo delle matrici del sistema continuo
    A1(:,:,i) = [df(i)^2 + 2*mu/rt(i)^3, d2f(i), 0;
                 -d2f(i), df(i)^2 - mu/rt(i)^3, 0;
                 0,     0,      -mu/rt(i)^3         ];

    A2(1,2,i) = 2*df(i);
    A2(2,1,i) = -A2(1,2,i);

    Acont(:,:,i) = [zeros(3,3) , eye(3) ; A1(:,:,i) , A2(:,:,i)];

    % discretizzazione del sistema
    Maux = [Acont(:,:,i) , Bcont ; zeros(3,9)] * Ts;
    Mexp = expm(Maux);

    Adisc(:,:,i) = Mexp(1:6,1:6);
    Bdisc(:,:,i) = Mexp(1:6,7:9);

end

%% Calcolo della traiettoria
P = zeros(6,6,Nsamp);
P(:,:,Nsamp) = Q;
NiterRiccati = 5;

for i = 1 : NiterRiccati
    for k = Nsamp-1:-1:1
        Ak = Adisc(:,:,k);
        Bk = Bdisc(:,:,k);
        Pk_next = P(:,:,k+1);
            
        % Termini intermedi per l'Eq. 32 (per leggibilità e velocità)
        temp_inv = (R + Bk' * Pk_next * Bk);
            
        % Equazione di Riccati Discreta (Eq. 32 del paper)
        P(:,:,k) = Q + Ak' * Pk_next * Ak - ...
                (Ak' * Pk_next * Bk) * (temp_inv \ (Bk' * Pk_next * Ak));
    end
    % Per la periodicità: l'inizio del giro precedente diventa la fine del nuovo
    P(:,:,Nsamp) = P(:,:,1);
end

flag = zeros(1,3);
Kk = zeros(3,6,Nsamp);

for k = 1:Nsamp-1
    Ak = Adisc(:,:,k);
    Bk = Bdisc(:,:,k);
    Pk = P(:,:,k);
    Pk_next = P(:,:,k+1);
    
    
    Kk(:,:,k) = (R + Bk' * Pk_next * Bk) \ (Bk' * Pk_next * Ak);

    % Analisi di stabiliità
    Acl = Ak - Bk * Kk(:,:,k);
    [~, flag(1)] = chol(Pk);
    Q_test = Pk - transpose(Acl)*Pk_next*Acl;
    [~, flag(2)] = chol(Q_test);

    if flag(1) > 0 || flag(2) > 0
        fprintf('La matrice del campione numero %d non è definita positiva\n', k);
        flag(3)=1;
    end
end

if flag(3) == 0
    disp('Il sistema è asintoticamente stabile');
end


Kk(:,:,Nsamp) = Kk(:,:,1); % Continuità periodica


%% condizioni iniziali (casuali)
X = zeros(6, N, Nsamp);
rng(42);                                                    % Seed per riproducibilità
for i = 1:N
    X(1:3, i, 1) = 100 * (2*rand(3,1) - 1);
    X(4:6, i, 1) = 1 * (2*rand(3,1) - 1);
end


%% Definizione del sistema
for k = 1 : Nsamp-1
    for i = 1 : N
        Ak = Adisc(:,:,k);
        Bk = Bdisc(:,:,k);

        u_ir = pinv(Bk) * (eye(6) - Ak) * Target_point(:, i);

        % Variazioni rispetto al riferimento
        zeta_bar_i = X(:, i, k) - Target_point(:, i);
        % Termine di consenso
        sum_consensus = zeros(6, 1);
        for j = 1:N
            if adjacency_matrix(i,j) == 1
                zeta_bar_j = X(:, j, k) - Target_point(:, j);
                sum_consensus = sum_consensus + adjacency_matrix(i,j) * (zeta_bar_j - zeta_bar_i);
            end
        end
        % Legge di controllo con consenso locale 
        zeta_bar_0 = zeros(6,1); % Target è all'origine
        u_tilde = eta * (1 + D(i,i) + G(i,i))^-1 * Kk(:,:,k) * (G(i,i) * (zeta_bar_0 - zeta_bar_i) + sum_consensus);
        % Input totale e Saturazione
        u_tot = u_ir + u_tilde;
        u_sat = max(min(u_tot, u_max), -u_max);
        % Propagazione degli stati
        X(:, i, k+1) = Ak * X(:, i, k) + Bk * u_sat;
    end
end



%% 6. Visualizzazione dei Risultati
figure('Name', 'Traiettorie 3D');
hold on; grid on; view(3);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
 
colors = lines(N);
for i = 1:N
    plot3(squeeze(X(1,i,:)), squeeze(X(2,i,:)), squeeze(X(3,i,:)), 'Color', colors(i,:), 'LineWidth', 1.5);
    scatter3(X(1,i,1), X(2,i,1), X(3,i,1), 50, colors(i,:), 'o', 'filled');
    scatter3(Target_point(1,i), Target_point(2,i), Target_point(3,i), 100, colors(i,:), 'x', 'LineWidth', 2);
end
title('Assemblaggio Multi-Satellite: Traiettorie in Orbita Ellittica');
legend('Traiettoria', 'Inizio', 'Target (Vertice Cubo)');





%Scelta posizione satelliti DA CORREGGERE
%{
% definizione target per ogni satellite
N = 5;     % da togliere
R = 100; % Raggio della sfera
golden_angle = pi * (3 - sqrt(5)); % Circa 137.5 gradi

theta = zeros(N, 1);
phi = zeros(N, 1);
posizione = zeros(N, 3);

for i = 1:N
    % 1. Calcolo di phi (Latitudine) - Corretto con parentesi
    % z varia linearmente da 1 a -1
    z = 1 - ((i - 1) / (N - 1)) * 2; 
    phi(i) = asin(z); 
    
    % 2. Calcolo di theta (Azimut) - Corretto con parentesi
    theta(i) = (i - 1) * golden_angle;
    theta(i) = mod(theta(i), 2*pi); % Normalizzazione
    
    % 3. Conversione corretta in coordinate Cartesiane
    posizione(i, 1) = R * cos(phi(i)) * cos(theta(i)); % X
    posizione(i, 2) = R * cos(phi(i)) * sin(theta(i)); % Y
    posizione(i, 3) = R * sin(phi(i));                 % Z (altezza)
end

% Visualizzazione tabella
risultati = table((1:N)', rad2deg(theta), rad2deg(phi), ...
    'VariableNames', {'Satellite', 'Theta_deg', 'Phi_deg'});
disp(risultati);

% --- Plot ---
figure('Color', 'w');
[sx, sy, sz] = sphere(30); 
mesh(sx*R, sy*R, sz*R, 'EdgeColor', [0.8 0.8 0.8], 'FaceAlpha', 0.05); 
hold on;

scatter3(posizione(:,1), posizione(:,2), posizione(:,3), 100, 'filled', 'MarkerFaceColor', '#0072BD');

for i = 1:N
    text(posizione(i,1)+5, posizione(i,2)+5, posizione(i,3)+5, ['Sat ', num2str(i)], 'FontSize', 10);
end
axis equal; grid on; view(3);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Distribuzione Ottimale Satelliti (Fibonacci Spiral)');
%}