clc;
clear;
N = 3; %Total no of events of the discrete random process
syms tau; %Time difference
x1 = -1;
x2 = 0;
x3 = 1;

X = [x1 x2 x3];

%Conditional probabilities of event occurring

p_i_eq_j = (1+2*exp(-abs(tau)))/3;
p_i_not_eq_j = (1-exp(-abs(tau)))/3;

%Probabilities of each events

p1 = 1/3;
p2 = 1/3;
p3 = 1/3;

P = [p1 p2 p3];

%ACF of discrete random process

acf = 0;
for i = 1:N
    for j = 1:N
        if i==j
            acf = acf + X(i)*X(j)*p_i_eq_j*P(j);
        else
            acf = acf + X(i)*X(j)*p_i_not_eq_j*P(j);
        end
    end
end
acf;
figure('Name','Autocorrelation function of discrete random process');
fplot(tau,acf)
grid on
title('Autocorrelation function of discrete random process');
ylim([0 0.8])
xlabel('Time difference in second');
ylabel('Amplitude');

