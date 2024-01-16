clear;
% Name: Syed Maisam Abbas, Student number: 501103255, Section: 01
%%
% ELE532_LAB4: The Fourirer Trasnfrom: Properties and Applications 

%% PREPERATION

% Generated Signal from Figure 1

N = 100; yhrrr
PulseWidth = 10;
t = [0:1:(N-1)];
x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
f = [-(N/2):1:(N/2)-1]*(1/N);

% Plot x(t) 

figure('name','A.0');
stairs(t,x); 
title('Generated Signal from Figure 1');
xlabel('t');
ylabel('x(t)');
grid on; axis([-10,110,-0.1,1.1])

%% Problem A: The Fourier Transform and its Properties 

%% Problem A.1

figure(2);
u = @(t) 1.0*(t>=0);
x = @(t) 1.*(u(t)-u(t-10));
h = @(t) 1.*(u(t)-u(t-10));
dtau = 0.05;
tau = -1:dtau:25;
ti = 0;
tvec = -0.25:0.1:25;
y = NaN*zeros(1,length(tvec));
for t = tvec
    ti = ti+1;
    xh = x(t-tau).*h(tau);
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], [0.8 0.8 0.8], 'edgecolor', 'none');
    xlabel('\tau');
    title('x(\tau) [solid], x(t-\tau) [dashed], x(\tau)x(t-\tau) [gray]');
    c = get(gca, 'children');
    set(gca, 'children', [c(2);c(3);c(4);c(1)]);
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    title('z(t) Time Domain');
    xlabel('t');
    ylabel('z(t) = \int x(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1 15]);
    grid;
    drawnow;
end

%% Problem A.2

% Compute X(w) 
Xf = fft(x);

% Calculate Z(w) by multiplying in the frequency domain
zf = Xf.*Xf;

%% Problem A.3

figure('name','Plot of z(w) = X(w)*X(w)');

subplot(211); plot(f,fftshift( abs(zf))); 
grid on;
title('Magnitude-Spectra Plot')
xlabel('w')
ylabel('|z(w)|')

subplot(212); plot(f,fftshift(angle(zf))); 
grid on;
title('Phase-Spectra Plot');
xlabel('w');
ylabel('∠z(w)');

%% Problem A.4

%% Time Domain z(t)
figure(2);
u = @(t) 1.0*(t>=0);
x = @(t) 1.*(u(t)-u(t-10));
h = @(t) 1.*(u(t)-u(t-10));
dtau = 0.05;
tau = -1:dtau:25;
ti = 0;
tvec = -0.25:0.1:25;
y = NaN*zeros(1,length(tvec));
for t = tvec
    ti = ti+1;
    xh = x(t-tau).*h(tau);
    lxh = length(xh);
    y(ti) = sum(xh.*dtau);
    subplot(2,1,1), plot(tau,h(tau),'k-',tau,x(t-tau),'k--',t,0,'ok');
    axis([tau(1) tau(end) -2.0 2.5]);
    patch([tau(1:end-1);tau(1:end-1);tau(2:end);tau(2:end)], [zeros(1,lxh-1);xh(1:end-1);xh(2:end);zeros(1,lxh-1)], [0.8 0.8 0.8], 'edgecolor', 'none');
    xlabel('\tau');
    title('x(\tau) [solid], x(t-\tau) [dashed], x(\tau)x(t-\tau) [gray]');
    c = get(gca, 'children');
    set(gca, 'children', [c(2);c(3);c(4);c(1)]);
    subplot(2,1,2), plot(tvec,y,'k',tvec(ti),y(ti),'ok');
    title('z(t) Time Domain');
    xlabel('t');
    ylabel('z(t) = \int x(\tau)x(t-\tau) d\tau');
    axis([tau(1) tau(end) -1 15]);
    grid;
    drawnow;
end

%% Frequency Domain Z(w) 🚨🚨🚨

zf = ifft(zt);

figure('name','Plotting z(w) derived from z(t)'); 
plot(t,zf); 

title('z(t) from z(w)');
xlabel('t');
ylabel('z(t)');
grid on; 
axis([-10,110,-0.1,10.1]);

%% Problem A.5

%% PART 1: For pulse width of 5

% signal generation, x(t)
N = 100; 
PulseWidth = 5;
t = [0:1:(N-1)];
xt = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
f = [-(N/2):1:(N/2)-1]*(1/N);

zt = conv(xt,xt);

% Compute X(w) 
xf = fft(x);

% Calculate Z(w) by multiplying in the frequency domain
zf = xf.*xf;

% Plot magnitude- and phase- spectra for z(w)
figure('name','Plot of z(w) = X(w)*X(w)');

subplot(211); plot(f,fftshift( abs(zf))); 
grid on;
title('Magnitude-Spectra Plot')
xlabel('w')
ylabel('|z(w)|')

subplot(212); plot(f,fftshift(angle(zf))); 
grid on;
title('Phase-Spectra Plot');
xlabel('w');
ylabel('∠z(w)');

%% Part 2: For pulse width of 25

% signal generation, x(t)
PulseWidth = 25;
N = 100; 
t = [0:1:(N-1)];
xt = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
f = [-(N/2):1:(N/2)-1]*(1/N);

zt = conv(xt,xt);

% Compute X(w) 
Xf = fft(x);

% Calculate Z(w) by multiplying in the frequency domain
zf = Xf.*Xf;

% Plot magnitude- and phase- spectra for z(w)
figure('name','Plot of z(w) = X(w)*X(w)');

subplot(211); plot(f,fftshift( abs(zf))); 
grid on;
title('Magnitude-Spectra Plot')
xlabel('w')
ylabel('|z(w)|')

subplot(212); plot(f,fftshift(angle(zf))); 
grid on;
title('Phase-Spectra Plot');
xlabel('w');
ylabel('∠z(w)');

%% Problem A.6

%% PART 1: Magnitude and Phase Plot for w+(t)

N = 100;
PulseWidth = 10;
t = [0:1:(N-1)];
x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
wplus = x.*exp(1j.*(pi/3).*t)
Xf = fft(wplus);
f = [-(N/2):1:(N/2)-1]*(1/N);
figure();
subplot(311); 
plot(f,fftshift(Xf));
title('X(w)');
xlabel('w');
subplot(312); 
plot(f,fftshift(abs(Xf)));
title('|X(w)|');
xlabel('w');
subplot(313); 
plot(f,fftshift(angle(Xf)));
title('angle X(w)');
xlabel('w');

%% PART 2: Magnitude and Phase Plot for w-(t)

N = 100;
PulseWidth = 10;
t = [0:1:(N-1)];
x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
wminus = x.*exp(-1j.*(pi/3).*t);
Xf = fft(wminus);
f = [-(N/2):1:(N/2)-1]*(1/N);
figure();
subplot(311); 
plot(f,fftshift(Xf));
title('X(w)');
xlabel('w');
subplot(312); 
plot(f,fftshift(abs(Xf)));
title('|X(w)|');
xlabel('w');
subplot(313); 
plot(f,fftshift(angle(Xf)));
title('angle X(w)');
xlabel('w');

%% PART 3: Magnitude and Phase Plot for wc(t)

N = 100;
PulseWidth = 10;
t = [0:1:(N-1)];
x = [ones(1,PulseWidth), zeros(1,N-PulseWidth)];
wc = x.*cos((pi/3).*t);
Xf = fft(wc);
f = [-(N/2):1:(N/2)-1]*(1/N);
figure();
subplot(311); 
plot(f,fftshift(Xf));
title('X(w)');
xlabel('w');
subplot(312); 
plot(f,fftshift(abs(Xf)));
title('|X(w)|');
xlabel('w');
subplot(313); 
plot(f,fftshift(angle(Xf)));
title('angle X(w)');
xlabel('w');



