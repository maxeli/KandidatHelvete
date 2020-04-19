clc
clf
Amp = 10;
f = 0.1;%:2:80;
t = linspace (0,100);
A = out.get('phase_A');
B = out.get('phase_B');
C = out.get('phase_C');
hold on, grid on
axis([0 10 -11 11])
plot(A)
plot(B)
plot(C)