%% import
%%% sroots
load("sroots.mat");
% programmatic matfile access init goes here

%% initial params
%%% from lab 1
alpha = 25;
beta = 2334;

%%% given model
A = [0, 1;
    0, -alpha];
B = [0;
    beta];
C = [1, 0];
D = 0;

%%%% derived / calculated characteristics
nOrder = min(size(A));
eigA = eig(A);

%% begin design
%%% initial selected design params
Ts = 0.25; % init val 0.3

%%%% addn dynamics
phia = 1;
gammaa = 1;

%%% calc design params
TbetaMax = max(imag(eigA));
T = min(Ts / (20 * (nOrder + 1)), ...
    pi / (5 * TbetaMax));

[phi, gamma] = c2d(A, B, T);

% programmatic bessel selection goes here
sPoles = s3 / Ts;
zPoles = exp(T * sPoles);

[K1, K2, delta1, delta2] = tsd(phi, gamma, C, phia, gammaa, zPoles, T, 'place')

%% settling time and overshoot
max_y = max(y(:, 2))
PO = (max_y - 60) / 60 * 100 % percent overshoot
max_u = max(u(:, 2))
ind = find(y(:, 2) > 0.99 * 60);
ST = y(min(ind), 1) % settling time

% %% design optimization goes here

% %% emulate gain discrepancy q for B
% 
% UGM = 0.15;
% q = 10 ^ (UGM / 20);
