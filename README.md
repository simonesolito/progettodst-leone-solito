# progettodst-leone-solito
programma progetto dst
% Parametri orbitali
mu = 3.986e14;                    % parametro gravitazionale standard (Terra)
a = 7000e3;                       % lunghezza del semiasse maggiore (7000 km)
e = 0.1;                          % Orbita ellittica
Ts = 10;                          % Campionamento ogni 10 secondi
Tmax = 2 * pi * sqrt(a^3/mu);     % durata dell'orbita

T = 0 : Ts : Tmax;                % asse dei tempi

omega_0 = sqrt(mu/a^3);           % velocità angolare media del satellite
M = omega_0 * T;                  % anomalia media
m = 500e3;                        % massa satellite

% calcolo dell'anomalia eccentrica E partendo dalla equazione M = E - e * sin(E)
% Newton-Raphson

E = M;
for i = 1:10
    E = E - (E - e*sin(E) - M) / (1 - e*cos(E));
end

% calcolo di posizione, velocità e accelerazione angolare
f = 2 * atan2(sqrt(1 + e) * sin(E/2) , sqrt(1 - e) * cos(E/2));
df = omega_0 * (1 + e * cos(f)).^2/(1 - e^2)^(3/2);
d2f = -2* omega_0^2 * e * sin(f) .* (1 + e * cos(f)).^3/(1 - e^2)^3;

rt = a * (1 - e^2)./(1 + e * cos(f));            % distanza satellite terra

A1 = zeros(3,3,length(T));                       % LA PARTE CONTINUA LA POTREI MEMORIZZARE MOMENTANEAMENTE
A2 = A1;                                         % RIDUCENDO LE MATRICI A DUE DIMENSIONI
Acont = zeros(6,6,length(T));
Bcont = 1/m * [zeros(3,3) ; eye(3)];

Adisc = Acont;
Bdisc = Bcont;
