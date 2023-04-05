load("sroots.mat");

%% initial params

alpha = 25;
beta = 2334;

constA = 41.5;
constC = alpha;
constD = beta;
n = 495;
g = 9.81;

A = [0, 1, 0, 0;
    constA, 0, 0, -(constA .* constC)./(n .* g);
    0, 0, 0, 1;
    0, 0, 0, -constC];
B = [0; (constA .* constD)./(n .* g); 0; constD];
C = [1, 0, 0, 0;
    0, 0, 1, 0];
D = zeros(2,1);

x0 = [0.17; 0; 0; 0];

T = 0.01; % carried over from lab 2
[phi, gamma] = c2d(A, B, T);
%% regulator vars
% (from lab 2)

Ts = .9; % default 1, between .1 and 2

%% regulator poles

sPoles = [-25 s3/Ts];

zpoles = exp(T * sPoles);

K = place(phi, gamma, zpoles);

disp('Regulator DSM')
dsm(phi, gamma, K)

%% observer vars

Tso = Ts / 5; % 5-times faster observer speed factor

%% observer poles

soPoles = s4/Tso; % 4th-order bessel poles, 5-times faster

zoPoles = exp(T * soPoles);

L = place(phi', C', zoPoles)';

disp('Observer DSM')
dsm_regob(phi, gamma, C, K, L)

%% observer coefficients and conditions

Ao = phi - L * C - gamma * K;
Bo = L;
Co = -K;
Do = zeros(1,2);

xo0 = C \ (C * x0);

%% emulate gain discrepancy q for B

UGM = 0.15;
q = 10 ^ (UGM / 20);
