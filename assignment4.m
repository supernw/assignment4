set(0, 'defaultFigureWindowStyle', 'docked')
clear
clf
close all

addpath('./code/');

%Get R3

%This section takes a long time to run
%{
volt = linspace(0.1, 10, 5);
curr = zeros(1, length(volt));
for i= 1:length(volt)
    curr(i) = PART23_Func(volt(i), 40e-9);
end

plot(volt, curr)
title('Voltage vs Current - Varrying Voltage')
xlabel('Voltage (V)')
ylabel('Current (A/m)')
c = polyfit(volt, curr, 1);
invR = c(1);
R3 = invR^-1/200e-9;

%}
R3 = 23.57; %Resistance from Assignment 3 (Ohm)
%Constants
G1 = 1/1;
C_const = 0.25;
G2 = 1/2;
L = 0.2;
G3 = 1/R3;
alpha = 100;
G4 = 1/0.1;
G0 = 1/1000;

Vs_var = linspace(-10,10,1000);
%G MATRIX
    %V1 V2 V3 I3 V4 Vo
G = [1 0 0 0 0 0 ; %1
    -G1 G1+G2 0 1 0 0 ; %2
    0 1 -1 0 0 0; %3
    0 0 G3 -1 0 0 ; %4
    0 0 0 -alpha 1 0 ; %5
    0 0 0 0 -G4 G4+G0]; %6

%C MATRIX
C = [0 0 0 0 0 0;
     -C_const C_const 0 0 0 0;
     0 0 0 -L 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];

%F MATRIX
F = [1;0;0;0;0;0];

%PLOT
figure(1)
V3 = zeros(1,length(Vs_var));
Vo = zeros(1,length(Vs_var));
for i = 1:length(Vs_var)
    F(1,1) = Vs_var(i);
    V = G\F;
    V3(i) = V(3);
    Vo(i) = V(6);
end

figure(1)
plot(Vs_var, V3);
hold on;
plot(Vs_var, Vo);
hold off;
xlabel('Input Voltage (V)')
ylabel('Voltage (V)')
title('V3 and Vo Plot')
legend('V3','Vo');
%AC
w = linspace(0,1000,10000);
F(1) = 1;
V1AC = zeros(1,length(w));
VoAC = zeros(1,length(w));
db = zeros(1,length(w));
for i = 1:length(w)
    V = (G+1i*w(i)*C)\F;
    V1AC(i) = V(1);
    VoAC(i) = V(6);
    db(i) = 20*log(abs(VoAC(i)/V1AC(i)));
end

figure(2)
subplot(1,2,1)
plot(w, abs(VoAC));
xlabel('Omega (Deg)')
ylabel('Voltage (V)')
title('Output Voltage (AC)')

subplot(1,2,2)
semilogx(w, db);
ylabel('Gain (dB)')
xlabel('Omega (Freq)')
title('Gain plot')

%Random Perturbations on C
x = linspace(0,5,1000);
 cDist = C_const + 0.05*randn(5000,1);
 
V1ACR = zeros(1,length(x));
VoACR = zeros(1,length(x));
dbACR = zeros(1,length(x));
for i = 1:length(x)
    
    C_rand = [0 0 0 0 0 0;
     -cDist(i) cDist(i) 0 0 0 0;
     0 0 0 -L 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];
    V = (G+1i*pi*C_rand)\F;
    V1ACR(i) = V(1);
    VoACR(i) = V(6);
    dbACR(i) = abs(V1ACR(i)/VoACR(i));
end
figure(4)
histogram(dbACR,50);
ylabel('Gain (dB)')
xlabel('Omega (freq)')
title('Gain plot')

%Part 4
num_step = 1000; %Number of steps
dt = 0.001; %Time step
time_test = linspace(0,1,1000);
%Step
timeSum = 0;
F(1) = 0;
V = zeros(6,1);
VoStep = zeros(1,num_step);
VinStep = zeros(1,num_step);
time = zeros(1,num_step);
freq = zeros(1,num_step);
for i=1:num_step
    if round(timeSum,4) == 0.03
           F(1) = 1;
    end
    time(i) = timeSum;
    freq(i) = 1/timeSum;
    timeSum = timeSum + dt;
    V_old = V;
    V = (C./dt + G)\(C*(V_old./dt) + F);
    VoStep(i) = V(6);
    VinStep(i) = V(1);
end
figure(5)
subplot(1,2,1);
plot(time, VoStep);
hold on;
plot(time, VinStep);
hold off;
xlabel('Time (seconds)')
ylabel('Voltage (Volts)')
title('Transient Simulation - Step Function')

subplot(1,2,2);
FVinStep = fftshift(fft(VinStep));
FVoStep = fftshift(fft(VoStep));
plot(20*log(abs(FVinStep)))
hold on;
plot(20*log(abs(FVoStep)))
hold off;
xlabel('Frequency')
ylabel('Gain')
title('Frequency Content of Input and Output Signals')
%Sine Wave
time_old = 0;
figure(6)
Vsin = zeros(6,1);
VoSin = zeros(1,num_step);
VinSin = zeros(1,num_step);
time = zeros(1,num_step);
freq = zeros(1,num_step);
timeSum = 0;
for i=1:num_step
    F(1) = sin(2*pi*(1/0.03)*timeSum);
    time(i) = timeSum;
    freq(i) = 1/timeSum;
    timeSum = timeSum + dt;
    V_oldSin = Vsin;
    Vsin = (C/dt + G)\(C*(V_oldSin/dt) + F);
    VoSin(i) = Vsin(6);
    VinSin(i) = Vsin(1);
end
subplot(1,2,1);
plot(time, VinSin);
hold on;
plot(time, VoSin);
xlabel('Time (seconds)')
ylabel('Voltage (Volts)')
title('Sine Wave Input')

subplot(1,2,2);
FVinSin = fftshift(fft(VinSin));
FVoSin = fftshift(fft(VoSin));
plot(20*log(abs(FVinSin)))
hold on;
plot(20*log(abs(FVoSin)))
hold off;
xlabel('Frequency')
ylabel('Gain (dB)')
title('Transient Simulation - Sine Function')
%Gaussian Pulse
gsig = 0.03;
gmu = 0.06;
gp = @(t) 1*exp( (-0.5*(t - gmu).^2)/(gsig^2) );

tSum = 0;
time_vec = zeros(1,num_step);
VoldG = zeros(6,1);
VGin = zeros(1,num_step);
VGo = zeros(1,num_step);
for i=1:num_step
    F(1) = gp(tSum);
    time_vec(i) = tSum;
    V = (C/dt + G)\(C*(VoldG/dt) + F);
    VoldG = V;
    tSum = tSum + dt;
    VGin(i) = V(1);
    VGo(i) = V(6);
end
figure(7)
subplot(1,2,1);
plot(time, VGin);
hold on;
plot(time, VGo);
xlabel('Time (seconds)')
ylabel('Gain (dB)')
title('Gaussian Pulse Input')

subplot(1,2,2);
FVGin = fftshift(fft(VGin));
FVGo = fftshift(fft(VGo));
plot(20*log(abs(FVGin)))
hold on;
plot(20*log(abs(FVGo)))
hold off;
xlabel('Frequency')
ylabel('Gain (dB)')
title('Transient Simulation - Gaussian Pulse Function')
%G MATRIX
    %V1 V2 V3 I3 V4 Vo
GN = [1 0 0 0 0 0 ; %1
    -G1 G1+G2 0 1 0 0 ; %2
    0 1 -1 0 0 0; %3
    0 0 G3 -1 0 0 ; %4
    0 0 0 -alpha 1 0 ; %5
    0 0 0 0 -G4 G4+G0]; %6

%C MATRIX
Cn_const = 0.00001;
CN = [0 0 0 0 0 0;
     -C_const C_const 0 0 0 0;
     0 0 0 -L 0 0;
     0 0 Cn_const 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];

%F MATRIX
FN = [1;0;-1;0;0;0];
tSum = 0;
time_vecN = zeros(1,num_step);
VoldGN = zeros(6,1);
VGinN = zeros(1,num_step);
VGoN = zeros(1,num_step);
for i=1:num_step
    FN(1) = gp(tSum);
    FN(3) = -0.001*randn;
    time_vecN(i) = tSum;
    V = (CN/dt + GN)\(CN*(VoldGN/dt) + FN);
    VoldGN = V;
    tSum = tSum + dt;
    VGinN(i) = V(1);
    VGoN(i) = V(6);
end
figure(8)
subplot(1,2,1);
plot(time, VGinN);
hold on;
plot(time, VGoN);
hold off;
xlabel('Time (seconds)')
ylabel('Voltage (Volts)')
title('Gaussian Pulse Input - Time Response')

subplot(1,2,2);
FVGinN = fftshift(fft(VGinN));
FVGoN = fftshift(fft(VGoN));
plot(20*log(abs(FVGinN)))
hold on;
plot(20*log(abs(FVGoN)))
hold off;
xlabel('Frequency')
ylabel('Gain (dB)')
title('Gaussian Pulse Input - Frequency Response')
