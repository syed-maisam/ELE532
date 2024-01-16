clear;
% Name: Syed Maisam Abbas, Student number: 501103255, Section: 01
%%
% ELE532_LAB2: System Properties and Convolution

%% NOTE: If an error is found in any given section: you must run the code in a seperate temp file
%% Part A: Impulse Response

%Problem A.1
% Set component values:
R = [1e4, 1e4, 1e4]; 
C = [1e-6, 1e-6];
% Determine coefficients for characteristic equation:
A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
% Determine characteristic roots:
lambda = roots(A);
%root of the matrix 
poly(lambda);
disp ("The lambda values for the poly function are:");
disp(lambda);

%%
%Problem A.2: 
t = (0:0.0005:0.1);
R = [1e4, 1e4, 1e4]; 
C = [1e-6, 1e-6];
A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
lambda = roots(A);
poly(lambda);
u = @(t) 1 * (t>=0);
h = @(t) (C(1).* exp(lambda(1).*t) + C(2).* exp(lambda(2).*t)).*(u(t));
figure; 
plot (t, h(t));
xlabel('time'); ylabel('Amplitude'); title('Problem A1 plot'); grid;
hold off;

%%
%Problem A.3: 
function [lambda] = CH2MP2(R,C)
% Function M-file finds characteristic roots of op-amp circuit.
% INPUTS:   R = length-3 vector of resistances
%           C = length-2 vector of capacitances
% OUTPUTS:  lambda = characteristic roots

% Determine coefficients for characteristic equation:
A = [1, (1/R(1)+1/R(2)+1/R(3))/C(2), 1/(R(1)*R(2)*C(1)*C(2))];
% Determine characteristic roots:
lambda = roots(A);
end

%% Part B: Convolution

%Problem B.1
% Script M-file graphically demonstrates the convolution process.
figure(1); % Create figure window and make visible on screen
u = @(t) 1.0*(t>=0);
x = @(t) 1.5*sin(pi*t).*(u(t)-u(t-1));
h = @(t) 1.5*(u(t)-u(t-1.5))-u(t-2)+u(t-2.5);
dtau = 0.005; tau = -1:dtau:4;
ti = 0; tvec = -.25:.1:3.75;
y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
for t = tvec,
    ti = ti+1; % Time index
    xh = x(t-tau).*h(tau); lxh = length(xh);
    y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
    subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
        [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
        [.8 .8 .8],'edgecolor','none');
    xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
    c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
    subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1.0 2.0]); grid;
    pause; %drawnow changed to pause
end

%Problem B.2
figure(2);
u = @(t) 1.0*(t>=0);
x = @(t) u(t)-u(t-2);
h = @(t) (t+1).*(u(t+1)-u(t));
dtau = 0.005; tau = -2:dtau:2;
ti = 0; 
tvec = -5:.1:5; y = NaN*zeros(1,length(tvec)); 
% Pre-allocate memory
for t = tvec
ti = ti+1; % Time index
xh = x(t-tau).*h(tau); lxh = length(xh);
y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
axis([tau(1) tau(end) -2.0 2.5]);
patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
[zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
[.8 .8 .8],'edgecolor','none');
xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
axis([tau(1) tau(end) -1.0 2.0]); grid;
pause;
end

%Problem B.3

%PART A (Figure 2.4-27, part a)
figure(3);
u = @(t) 1.0*(t>=0);
%Assume the values of A and B: 
A = 1; B = 2;
x = @(t) A*(u(t-4)-u(t-6));
h = @(t) B*(u(t+5)-u(t+4));
dtau = 0.005; tau = -6:dtau:3;
ti = 0; tvec = -5:.1:5;
y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
for t = tvec
ti = ti+1; % Time index
xh = x(t-tau).*h(tau); lxh = length(xh);
y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
axis([tau(1) tau(end) -2.0 2.5]);
patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
[zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
[.8 .8 .8],'edgecolor','none');
xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
axis([tau(1) tau(end) -1.0 2.0]); grid;
drawnow; %drawnow used for simplicity
end

%PART B (Figure 2.4-27, part b)
figure(4);
u = @(t) 1.0*(t>=0);
%Assume the values of A and B: 
A = 1; B = 0.5;
x = @(t) A*(u(t-3)-u(t-5));
h = @(t) B*(u(t+5)-u(t+3));
dtau = 0.005; tau = -6:dtau:3;
ti = 0; tvec = -5:.1:5;
y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
for t = tvec
ti = ti+1; % Time index
xh = x(t-tau).*h(tau); lxh = length(xh);
y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
axis([tau(1) tau(end) -2.0 2.5]);
patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
[zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
[.8 .8 .8],'edgecolor','none');
xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
axis([tau(1) tau(end) -1.0 2.0]); grid;
drawnow; %drawnow used for simplicity
end

%PART C (Figure 2.4-27, part h)
figure(5);
u = @(t) 1.0*(t>=0);
x = @(t) exp(t).*(u(t+2)-u(t));
h = @(t) exp(-2*t) .*(u(t)-u(t-1));
dtau = 0.005; tau = -6:dtau:3;
ti = 0; tvec = -5:.1:5;
y = NaN*zeros (1,length(tvec)); % Pre-allocate memory
for t = tvec
ti = ti+1; % Time index
xh = x(t-tau).*h(tau); lxh = length(xh);
y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
axis([tau(1) tau(end) -2.0 2.5]);
patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
[zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
[.8 .8 .8],'edgecolor','none');
xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
axis([tau(1) tau(end) -1.0 2.0]); grid;
drawnow; %drawnow used for simplicity
end

%% Part C: System Behavior and Stability
%Problem C.1

t = (-1:0.001:5);
u = @(t) 1.0.*(t>=0);
h1 = @(t) exp(t/5).*u(t);
h2 = @(t) 4.*exp(-t/5).*u(t);
h3 = @(t) 4.*exp(-t).*u(t);
h4 = @(t) 4.*(exp(-t/5)-exp(-t)).*u(t);
figure;
plot (t, h1(t), t, h2(t), t, h3(t), t, h4(t));
xlabel('time'); ylabel('hn(t)'); title('LDI Systems'); grid;
legend('S1 = h1(t)', 'S2 = h2(t)', 'S3 = h3(t)', 'S4 = h4(t)');
hold on;

%Problem C.2 

%The characteristic values (eigenvalues) of the 
%systems S1-S4 are found using the exponential value 
%and are provided in the lab report

%Problem C.3

%PART A (Implementation of S1)
u = @(t) 1.0*(t>=0);
x = @(t) sin(5*t).*(u(t) - u(t-3));
%Impulse responce function
h = @(t) exp(t/5).*(u(t) - u(t-20));
%Modified CH2MP2 section
dtau = 0.005; tau = 0:dtau:20;
ti = 0; tvec = 0:0.1:20;
y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
for t = tvec
ti = ti+1; % Time index
xh = x(t-tau).*h(tau); lxh = length(xh);
y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
axis([tau(1) tau(end) -2.0 2.5]);
patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
[zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
[.8 .8 .8],'edgecolor','none');
xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
axis([tau(1) tau(end) -1.0 2.0]); grid;
drawnow;
end

%PART B (Implementation of S2)
u = @(t) 1.0*(t>=0);
x = @(t) sin(5*t).*(u(t) - u(t-3));
%Impulse responce function h2(t)
h = @(t) 4*exp(-t/5).*(u(t) - u(t-20));
%Modified CH2MP2 section
dtau = 0.005; tau = 0:dtau:20;
ti = 0; tvec = 0:0.1:20;
y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
for t = tvec
ti = ti+1; % Time index
xh = x(t-tau).*h(tau); lxh = length(xh);
y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
axis([tau(1) tau(end) -2.0 2.5]);
patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
[zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
[.8 .8 .8],'edgecolor','none');
xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
axis([tau(1) tau(end) -1.0 2.0]); grid;
drawnow;
end

%PART C (Implementation of S3)
u = @(t) 1.0*(t>=0);
x = @(t) sin(5*t).*(u(t) - u(t-3));
%Impulse responce function h3(t)
h = @(t) 4*exp(-t).*(u(t) - u(t-20));
%Modified CH2MP2 section
dtau = 0.005; tau = 0:dtau:20;
ti = 0; tvec = 0:0.1:20;
y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
for t = tvec
ti = ti+1; % Time index
xh = x(t-tau).*h(tau); lxh = length(xh);
y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
axis([tau(1) tau(end) -2.0 2.5]);
patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
[zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
[.8 .8 .8],'edgecolor','none');
xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
axis([tau(1) tau(end) -1.0 2.0]); grid;
drawnow;
end

%PART D (Implementation of S4)
u = @(t) 1.0*(t>=0);
x = @(t) sin(5*t).*(u(t) - u(t-3));
%Impulse responce function h4(t)
h = @(t) 4*(exp(-t/5)-exp(-t)).*(u(t) - u(t-20));
%Modified CH2MP2 section
dtau = 0.005; tau = 0:dtau:20;
ti = 0; tvec = 0:0.1:20;
y = NaN*zeros(1,length(tvec)); % Pre-allocate memory
for t = tvec
ti = ti+1; % Time index
xh = x(t-tau).*h(tau); lxh = length(xh);
y(ti) = sum(xh.*dtau); % Trapezoidal approximation of convolution integral
subplot(2,1,1),plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
axis([tau(1) tau(end) -2.0 2.5]);
patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)],...
[zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)],...
[.8 .8 .8],'edgecolor','none');
xlabel('\tau'); title('h(\tau) [solid], x(t-\tau) [dashed], h(\tau)x(t-\tau) [gray]');
c = get(gca,'children'); set(gca,'children',[c(2);c(3);c(4);c(1)]);
subplot(2,1,2),plot(tvec,y,'k',tvec(ti),y(ti),'ok');
xlabel('t'); ylabel('y(t) = \int h(\tau)x(t-\tau) d\tau');
axis([tau(1) tau(end) -1.0 2.0]); grid;
drawnow;
end

%The plots for S2, S3, and S4 all produce similar waveforms 
%in the convolution with h(t). These waveforms all have a 
%similar structure as it starts off with some sort of sin wave and 
%saturates into a straight line.


%% Part D: Discussion

%This section is discusssed in the lab report




