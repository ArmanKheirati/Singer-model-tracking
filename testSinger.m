clearvars; clc; close all;
randn('seed',5);
N = 2000;
snr = 10;
Ts = 0.01;
alpha = 1.6;
%%%%%%%%%%%%%%  3e-3 < \alpha < 1e3 %%%%%%%%%%%%%%%%%%%%%%%%%
Q = 0.001^2; % Process noise
R = [0.01^2 0;0 0.005^2]; % Measurement noise
beta = exp(-alpha*Ts);
B = [Ts^3/6 Ts^2/2 Ts]';
f = inline('[1 Ts (alpha*Ts - 1 + exp(-alpha*Ts))/alpha^2;0 1 (1-exp(-alpha*Ts))/alpha;0 0 exp(-alpha*Ts)]*x','x', 'alpha', 'Ts' ,'t');
h = inline('[1 0 0;0 0 1]*x','x');
pe = inline('exp(-(1/2)*(sum(e.*(inv(U*[10^2 0;0 0.5^2])*e))))', 'e', 'U');
x(:,1) = [1 1 0.1]'; % Initial state used in simulation
xd(:,1) = x(:,1);
A = [1 Ts (alpha*Ts- 1 + beta)/alpha^2;0 1 (1-beta)/alpha;0 0 beta];
H = [1 0 0]; 
for t=1:N % Simulate the model
    xd(:,t+1) = A*xd(:,t);
    xdR(:,t) = H*xd(:,t) ;
end
for t=1:N % Simulate the model
    x(:,t+1) = A*x(:,t) + B*sqrtm(Q)*randn(1);
    xn(:,t) = H*x(:,t) + snr*sqrtm(R)*randn(2,1);
end
% figure, plot(xn(1,:))

y = xn(1,:)';
% Calculate coefficients m1, m2, and m3
m1 = sqrt(1/(2*alpha^5)*(1 - exp(-2*alpha*Ts) + 2*alpha*Ts + (2*alpha^3*Ts^3)/3 - 2*alpha^2*Ts^2 - 4*alpha*Ts*exp(-alpha*Ts)));
m2 = sqrt(1/(2*alpha^3)*(4*exp(-alpha*Ts) -3 - exp(-2*alpha*Ts) + 2*alpha*Ts));
m3 = sqrt(1/(2*alpha)*(1 - exp(-2*alpha*Ts)));

phi1 = -(beta+2);
phi2 = 1+2*beta;
phi3 = -beta;
zeta1 = m1;
zeta2 = -m1*(beta + 1) +m2*Ts +m3*(alpha*Ts - 1 + beta)/alpha^2;
zeta3 = m1*beta -m2*Ts*beta +m3*(Ts*(1-beta)/alpha - (alpha*Ts -1 + beta)/alpha^2);
a = [phi3 phi2 phi1 1]
b = [zeta3 zeta2 zeta1]

lam = 1e-4;

[xhat1, xhat2, xhat3, xhat4] = l2ARsmoothing(y, lam, a, b);
%[xhat1, xhat2, xhat3, xhat4] = l2ARsmoothing2(y, lam, a, b);
figure, hax=axes;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 6], 'PaperUnits', 'Inches', 'PaperSize', [10, 6])
hold on
% Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(x(1,:)), hold all, plot(xhat1), plot(xhat3)
xlabel('Samples','FontSize',14,'FontName','Times New Roman','interpreter','latex');
ylabel('Amplitude','FontSize',14,'FontName','Times New Roman','interpreter','latex');
legend('Interpreter', 'latex', 'FontSize', 18);
h = legend('True position', 'Matrix Smoothing', 'Forward-backward filtering')
figure, hax=axes;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 6], 'PaperUnits', 'Inches', 'PaperSize', [10, 6])
hold on
% Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(x(1,:)), hold all, plot(xhat2, '--'), plot(xhat4)
xlabel('Samples','FontSize',14,'FontName','Times New Roman','interpreter','latex');
ylabel('Amplitude','FontSize',14,'FontName','Times New Roman','interpreter','latex');
legend('Interpreter', 'latex', 'FontSize', 18);
h = legend('True position', 'Matrix FIR filter', 'Forward filtering')