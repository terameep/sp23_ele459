load("sroots.mat");

%% initial params
alpha = 25;
beta = 2334;

constA = 41.5;
constC = alpha;
constD = beta;
n = 495;
g = 9.81;

x0 = [0.17; 0; 0; 0];
noisePower = 0.002;
T = 0.01;

A = [0, 1, 0, 0;
    constA, 0, 0, -(constA .* constC)./(n .* g);
    0, 0, 0, 1;
    0, 0, 0, -constC];
B = [0; (constA .* constD)./(n .* g); 0; constD];
C = eye(4);
D = zeros(4,1);

%% slider variables
Ts = .66; % default 1, between .1 and 2
f = 5; % greater than or equal default 5

%% filter
Tsf = Ts/f; % filter settling time; choose f to be greater than or equal to 5
bwfp = (s1+1j*s1)/Tsf; % 2nd-order Butterworth filter pole scaled to Tsf-sec settling time
den = poly([bwfp conj(bwfp)]);
num = [den(end) 0]; % numerator is the constant term in the denominator polynomial
vfilter_a = tf(num,den); % analog differentiator with lowpass filter
vfilter_d = c2d(vfilter_a,T,'tustin'); % digital differentiator with lowpass filter
[Af,Bf,Cf,Df]=ssdata(vfilter_d); % convert vfilter_d to state-space model

%% poles
sPoles = s4/Ts;
% sPoles = [-25 s3/Ts]

zpoles = exp(T * sPoles);

[phi, gamma] = c2d(A, B, T);

K = place(phi, gamma, zpoles);

dsm(phi, gamma, K)
