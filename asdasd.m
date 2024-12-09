clc;
clear;
close all;

% Steady state analysis
% Develop an steady-state model of a Type 1 wind turbine connected to the grid

%% Squirrel cage induction motor equivalent scheme

%         Xs      Rs              Xr    Rr/s
% ------|||||--/\/\/\-----------|||||--/\/\/\----.
% +        --->         |   |    --->            |
%           Is          X   /     Ir             |
% V                 Xm  X   \ Rm                 |
%                       X   /                    |
% -                     |   |                    |
% -----------------------------------------------·

% Induction Machine parameters (from Assignment 2 tables) 

Vnom=960;               % Nominal grid voltage
Pnom=2;
fs = 50;                 % Synchronous freaquency [Hz]
pols=2;                 % Pole pairs
Rs=0.005;               % Stator resistance
Xs=2 * pi * fs * 4e-4;  % Leakage stator inductance (impedance)
Rr=0.009;               % Rotor resistance
Xr=2*pi*50*3e-4 ;       % Leakage rotor inductance (impedance)
Xm=2 * pi * fs * 15e-3; % Magnetizing branch inductance
Rm=140;                 % Iron branch resistance
V=Vnom/sqrt(3);         % .......
ws=2*pi*fs/pols;        % Synchronous speed [rad s^-1]

% Equations workflow to generate the torque-speed plot

s = -1 : 1e-4 : 1; 
wg = sort(( 1 - s ) .* ws); % [rad s^-1]  % rps of the rotor
Zm = ( Rm * 1j * Xm) / (Rm + 1j * Xm); 
paralel = (((Rr ./ s + 1j * Xr) * Zm) ./ ((Rr ./ s + 1j * Xr) + Zm )); 
% the parallel reactance changes depending on the freaquency fr (s)
imp = paralel + 1j*Xs + Rs;
% the total impedance has values that changes from negative to positive
Is = V ./ imp;  %stator current
Vr = Is .* paralel; 
Ir = Vr ./ (Rr ./ s + 1j.*Xr);  
% The rotor current is positive when it's in motor behaviour
Telec1 = (3 .* Rr .* abs(Ir).^2) ./ (s .* ws);
Pelec1 = Telec1 .* ( 1 - s ) .* ws; 
% Calcolo della potenza complessa iniettata nella rete
S = 3 * V * conj(Is);  % Potenza complessa trifase

% Calcolo della potenza attiva (parte reale)
P = real(S);

% Calcolo della potenza reattiva (parte immaginaria)
Q = imag(S);

figure (1); 
plot(wg * 30 / pi, -Telec1,'LineWidth',2);grid on;
title('Electric Torque - fast shaft')
xlabel('\omega_g [rpm]','FontSize',14);
ylabel('T_{elec} [Nm]','FontSize',14);

figure (2); 
plot(wg * 30 / pi, -Pelec1 ./ 1e6,'LineWidth',2);grid on;
title('Electric Power - fast shaft')
xlabel('\omega_g [rpm]','FontSize',14);
ylabel('P_{elec} [MW]','FontSize',14);

%% Wind Turbine Aereodynamics

% Create a script to obtain the Speed-Power and Speedtorque curves 
% for a wind turbine using the Assignment 2 data

% Cp parameters for fixed speed turbine
c1=0.44;
c2=125;
c3=0;c4=0;c5=0;
c6=6.94;
c7=16.5;
c8=0;
c9=-0.002;

% operations limit 
v_in = 3; % cut-in speed [m/s]
v_out = 20; % cut-in speed [m/s]

% WT parameters
R_turbina=38;          % WT diameter is 76 [m] 
A=pi * R_turbina^2;    % Area
rho=1.225;             % Air density
angle_pitch=0;         % Pitch angle
n_multiplicador = 80;  % Transmission ratio

% Calculation for different shaft speed and for different wind speeds         
wg_rpm = 60 .* wg ./ (2*pi);  % rotation per minute of the fast shaft
for ii=1:1:12       % Creation of the 'for' loop. It will execute it 6 times, as defined (1:1:6)
    vw=1+ii;     % Wind speed (for simplicity, we calculate it based on the iteration number)
    tsr= wg * R_turbina / (n_multiplicador * vw); % Calculation of the Tip speed ratio
    k1=( tsr+c8*angle_pitch).^(-1) - c9/(1+angle_pitch^3); % aux variable for the cp calculation
    cp = max(0,c1*(c2*k1-c3*angle_pitch-c4*angle_pitch^c5-c6).*exp(-c7*k1)); % Calculation of the Cp
    %constant for the current iteration
    T_turbina(ii,:) = (1/n_multiplicador)*0.5*rho*A*cp*vw^3./(wg/n_multiplicador); % Turbine torque (fast shaft)
    txt{ii}=['Wind=' num2str(vw) ' m/s']; % Creation of a 'txt' vector for the legend
end

P_turbina = T_turbina .* wg;

% Torque - Speed plot
figure(3)
plot(wg_rpm,T_turbina,'LineWidth',2)
xlabel('\omega_g fast shaft [rpm]','FontSize',14);
ylabel('T_{mech} [Nm]','FontSize',14);
legend(txt) % legend
grid on
title('Mechanical Torque - fast shaft') % title

% Power - Speed plot
figure(4)
plot(wg_rpm,P_turbina ./ 1e6,'LineWidth',2)
xlabel('\omega_g fast shaft [rpm]','FontSize',14);
ylabel('P_{mech} [MW]','FontSize',14);
legend(txt) % legend
grid on
title('Mechanical Power - fast shaft') % title

% With the mouse pointer check on the graphs that the maximum power
% for a certain wind speed doesn't lay on the same max torque wg_rpm

%% Superimpose the two graphs
% Trova i valori di T_elec negativi

% Trova il massimo di T_elec_neg e il suo indice
[T_elec_max, idx_max] = max(Telec1);

% Trova tutti gli indici con T_elec_neg > 0 prima del picco massimo
idx_positive = find(Telec1(1:idx_max) > 0); % Solo fino al picco

% Crea un vettore ridotto mantenendo NaN per i valori indesiderati
T_elec_reduce = Telec1;                % Copia il vettore originale
T_elec_reduce(setdiff(1:length(Telec1), idx_positive)) = NaN; % Imposta NaN dove Telec1 < 0

% Inizializza le intersezioni
wg_rpm_intersections = [];  % Velocità alle intersezioni
T_intersections = [];       % Coppia alle intersezioni
vw_intersections = [];      % Velocità del vento associata

% Loop sulle curve di coppia per ogni velocità del vento
for ii = 1:size(T_turbina, 1) % Per ogni curva colorata
    % Differenza tra T_turbina (curva corrente) e T_elec_reduce
    DeltaT = T_turbina(ii, :) - T_elec_reduce;
    
    % Trova i punti di cambio segno (escludi i NaN)
    idx_intersection = find(~isnan(DeltaT(1:end-1)) & ~isnan(DeltaT(2:end)) & ...
                            (DeltaT(1:end-1) .* DeltaT(2:end) < 0));
    
    % Registra i valori validi (velocità e coppia)
    wg_rpm_intersections = [wg_rpm_intersections, wg_rpm(idx_intersection)];
    T_intersections = [T_intersections, T_turbina(ii, idx_intersection)];
    vw_intersections = [vw_intersections, ones(1, length(idx_intersection)) * (1 + ii)]; % Velocità del vento
end

% Crea una tabella con i dati
IntersectionTable = table(vw_intersections', T_intersections', wg_rpm_intersections', ...
    'VariableNames', {'WindSpeed_vw', 'Torque_T', 'ShaftSpeed_wg'});

% Salva la tabella in un file CSV
writetable(IntersectionTable, 'IntersectionData.csv');

% Mostra i risultati
disp('Tabella delle intersezioni salvata:');
disp(IntersectionTable);

% Visualizza le intersezioni sul grafico
figure(5);
plot(wg_rpm, T_turbina, 'LineWidth', 2); grid on; % Coppia turbina
hold on;
plot(wg_rpm, T_elec_reduce, 'k', 'LineWidth', 2); % Coppia elettrica ridotta
scatter(wg_rpm_intersections, T_intersections, 'ro', 'filled'); % Intersezioni
title('Superimpose generator and turbine characteristics with intersections');
xlabel('\omega_g [rpm]', 'FontSize', 14);
ylabel('Torque [Nm]', 'FontSize', 14);
legend([txt, "Electric Torque", "Intersections"]);

%% Power curve 

P_t = zeros(size(T_intersections)); % Potenza calcolata

% Calcola la potenza per ogni intersezione
for i = 1:length(T_intersections)
    P_t(i) = T_intersections(i) * (wg_rpm_intersections(i) * 2 * pi / 60); % Calcola P = T * ω
    if P_t(i) > Pnom * 1e6 % Se supera la potenza nominale
        P_t(i) = Pnom * 1e6; % Limita la potenza
    end
end

% Assicurati che vw_intersections e P_t siano colonne
vw_intersections = vw_intersections(:); % Converte in colonna se non lo è
P_t = P_t(:); % Converte in colonna se non lo è

% Aggiungi i valori iniziali e finali
v_in = vw_intersections(1); % Velocità iniziale del vento
vw_intersections = [v_in; vw_intersections; 20; 20]; % Aggiungi velocità iniziale e finale
P_t = [0; P_t; Pnom * 1e6; 0]; % Aggiungi potenza iniziale (0) e finale (2 MW, poi 0)

% Limiti per il plot
x_lim = max(vw_intersections) + 1; % Limite asse x (aggiungi margine)
y_lim = Pnom + 1; % Limite asse y (aggiungi margine)

% Grafico della potenza
figure(6);
plot(vw_intersections, P_t / 1e6, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b'); % Potenza in MW
grid on;
xlabel('Wind Speed [m/s]', 'FontSize', 14);
ylabel('Power [MW]', 'FontSize', 14);
title('Turbine Power vs Wind Speed', 'FontSize', 14);
legend('Turbine Power', 'Location', 'NorthWest');
xlim([0 x_lim])
ylim([0 y_lim])

%% Power coeffficient 

% Calcolo della potenza ideale
P_ideal = 0.5 * rho * A * (vw_intersections.^3); % Potenza teorica
C_p = P_t ./ P_ideal; % Fattore di potenza (efficienza)

% Limiti per il plot
x_lim = max(vw_intersections) + 1; % Limite asse x (aggiungi margine)
y_lim = max(C_p) + 0.1; % Limite asse y (aggiungi margine)

% Grafico del fattore di potenza
figure(7);
plot(vw_intersections, C_p, 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
xlabel('Wind Speed [m/s]', 'FontSize', 14);
ylabel('Power Factor', 'FontSize', 14);
title('Power Factor vs Wind Speed', 'FontSize', 14);
xlim([0 x_lim])
ylim([0 y_lim])

%% Tip speed ratio 

tsr = R_turbina .* (mean(wg_rpm_intersections)./60) ./ (vw_intersections .* n_multiplicador);

% Limiti per il plot
x_lim = max(vw_intersections) + 1; % Limite asse x (aggiungi margine)
y_lim = max(tsr) + 0.1; % Limite asse y (aggiungi margine)

figure(8);
plot(vw_intersections, tsr, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
xlabel('Wind Speed [m/s]', 'FontSize', 14);
ylabel('TRS (-)', 'FontSize', 14);
title('Tip speed ratio vs Wind Speed', 'FontSize', 14);
xlim([0 x_lim ]);
ylim([0 y_lim]); % Il coefficiente di potenza massimo teorico è Betz limit (~0.59)

%% Grid active and reactive power 
%hola hola 
%modifica un acazzata 