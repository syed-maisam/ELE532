clear;
% Name: Syed Maisam Abbas, Student number: 501103255, Section: 01
%%
% ELE532_LAB2: Fourier Series Analysis 

%% Problem A.1

%Given the period signal, x1(t) in the lab manual, the following 
%expression for Dn can be derived:
D1 = (1/2) * [sin((3-n).*pi) + sin((3+n).*pi) + (1/2)*sin((1+n).*pi) + (1/2)*sin((1-n).*pi)];

%% Problem A.2

%Given the period signals, x2(t) and x3(t) in the lab manual, 
%the following expressions for Dn can be derived respectively:
D2 = (1/(n.*pi)) * sin ((n*pi)/2);
D3 = (1/(n.*pi)) * sin ((n*pi)/4);

%NOTE: The derivations for problems A1 and A2 are provided in the lab report

%% Problem A.3

function [D] = ELE532_LAB3_(D,n)
D1 = (1/2) * [sin((3-n).*pi) + sin((3+n).*pi) + (1/2)*sin((1+n).*pi) + (1/2)*sin((1-n).*pi)];
D2 = (1/(n.*pi)) * sin ((n*pi)/2);
D3 = (1/(n.*pi)) * sin ((n*pi)/4);

if (d == 1)
    D = D1;
end 

if (d == 2)
    D = D2;
end 

if (d == 3)
    D = D3;
end 

%% Problem A.4

%PART A 

%% x_1(t)
clf;
n = (-5:5);
D_n = 1./2.*((1./(pi.*n)).*sin((3-n).*pi )) + (1./pi.*n).*sin((3+n).*pi)  ...
    + (1./(2.*n.*pi).*sin((1+n).*pi)) + (1./(2.*n.*pi).*sin((1-n).*pi)) ;
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle Dn [rad]');

%% x_2(t)
clf;
n = (-5:5);
D_n = (1./(n.*pi).*sin((n.*pi)./2)); 
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle D_n [rad]');

%% x_3(t)
clf;
n = (-5:5);
D_n = (1./(n.*pi).*sin((n.*pi)./4)); 
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle D_n [rad]');

%PART B

%% x_1(t)
clf;
n = (-20:20);
D_n = 1./2.*((1./(pi.*n)).*sin((3-n).*pi )) + (1./pi.*n).*sin((3+n).*pi)  ...
    + (1./(2.*n.*pi).*sin((1+n).*pi)) + (1./(2.*n.*pi).*sin((1-n).*pi)) ;
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle Dn [rad]');

%% x_2(t)
clf;
n = (-20:20);
D_n = (1./(n.*pi).*sin((n.*pi)./2)); 
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle D_n [rad]');

%% x_3(t)
clf;
n = (-20:20);
D_n = (1./(n.*pi).*sin((n.*pi)./4)); 
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle D_n [rad]');

%PART C

%% x_1(t)
clf;
n = (-50:50);
D_n = 1./2.*((1./(pi.*n)).*sin((3-n).*pi )) + (1./pi.*n).*sin((3+n).*pi)  ...
    + (1./(2.*n.*pi).*sin((1+n).*pi)) + (1./(2.*n.*pi).*sin((1-n).*pi)) ;
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle Dn [rad]');

%% x_2(t)
clf;
n = (-50:50);
D_n = (1./(n.*pi).*sin((n.*pi)./2)); 
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle D_n [rad]');

%% x_3(t)
clf;
n = (-50:50);
D_n = (1./(n.*pi).*sin((n.*pi)./4)); 
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle D_n [rad]');

%PART D

%% x_1(t)
clf;
n = (-500:500);
D_n = 1./2.*((1./(pi.*n)).*sin((3-n).*pi )) + (1./pi.*n).*sin((3+n).*pi)  ...
    + (1./(2.*n.*pi).*sin((1+n).*pi)) + (1./(2.*n.*pi).*sin((1-n).*pi)) ;
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle Dn [rad]');

%% x_2(t)
clf;
n = (-500:500);
D_n = (1./(n.*pi).*sin((n.*pi)./2)); 
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle D_n [rad]');

%% x_3(t) 
clf;
n = (-500:500);
D_n = (1./(n.*pi).*sin((n.*pi)./4)); 
subplot(1,2,1); stem(n,abs(D_n),'.k'); 
xlabel('n'); ylabel('|D_n|'); 
subplot(1,2,2); stem(n,angle(D_n),'.k'); 
xlabel('n'); ylabel('\angle D_n [rad]');

end

%% Problem A.5

function x = a5(d ,Dn)
if(d == 1)
    w = pi/10;
elseif (d == 2)
    w = pi/10;
elseif (d == 3)
    w = pi/20;
end

t = -300:1:300;
x = zeros(size(t));
for i = 1:length(x)
    total = 0;
    j = 1;
    for n = -500:500
        total = total + Dn(j) * exp(1i* n * w * t(i));
        j = j+1;
    end
    x(i) = total;
end

figure(1);
plot(t, x, 'b')
xlabel('t (s)');
ylabel('x(t)');
if(d ~= 1)
    axis([-300 300 -1 2]);
end

title('Reconstructed Fourier Coefficients');
grid;

%% Problem A.6

%for x1(t)

%% PART A
clf
t = -300:1:300;
x = 0;
for n = -5:5
    D_n = 0;
if (n==3 || n==-3)
    D_n=(1/2);
end
if (n==1||n==-1)
    D_n=(1/4);
end 
x=x+D_n.*(exp(sqrt(-1)*n*(pi/10)*t)); 
end
x = real(x);
plot(t,x);
xlabel('t');
ylabel('x1(t)');
title('A6 1(a) from D_n and n=-5:5');

%% PART B

clf
t = -300:1:300;
x = 0;
for n = -20:20
    D_n = 0;
if (n==3 || n==-3)
    D_n=(1/2);
end
if (n==1||n==-1)
    D_n=(1/4);
end 
x=x+D_n.*(exp(sqrt(-1)*n*(pi/10)*t)); 
end
x = real(x);
plot(t,x);
xlabel('t'); 
ylabel('x1(t)');
title('A6 1(b) from D_n and n=-20:20');

%% PART C

clf
t = -300:1:300;
x = 0;
for n = -50:50
    D_n = 0;
if (n==3 || n==-3)
    D_n=(1/2);
end
if (n==1||n==-1)
    D_n=(1/4);
end 
x=x+D_n.*(exp(sqrt(-1)*n*(pi/10)*t)); 
end
x = real(x);
plot(t,x);
xlabel('t');
ylabel('x1(t)');
title('A6 1(c) from D_n and n=-50:50');

%% PART D

clf
t = -300:1:300;
x = 0;
for n = -500:500
    D_n = 0;
if (n==3 || n==-3)
    D_n=(1/2);
end
if (n==1||n==-1)
    D_n=(1/4);
end 
x=x+D_n.*(exp(sqrt(-1)*n*(pi/10)*t)); 
end
x = real(x);
plot(t,x);
xlabel('t');
ylabel('x1(t)');
title('A6 1(d) from D_n and n=-500:500');

%%

%for x2(t)

%% PART A

D_n=[-5:5];
nleftlim = -5;
nrightlim = 5;
x = 5+1;
for n = [nleftlim:nrightlim];
if n == 0,
    D_n(x) = 0.05;
else
    D_n(n-nleftlim+1) = (sin(n.*pi*0.5)./(n.*pi));
end
end
n = [nleftlim:nrightlim];
W0 = pi/10;
t = -300:300;
s = 300+1;
b=length(t);
x = zeros(1,b); %Pre-allocate Memory for recreated signal
for t=-300:300
for n=nleftlim:nrightlim
    x(t+s) = x(t+s) + real(D_n(n-nleftlim+1).*exp(n.*1i*W0*t));
end
end
t=-300:300;
plot(t,real(x),'b');
ylabel('x(t)');
xlabel('t');
title('Regenerated A6 2(a) from D_n for -5:5');
legend('Signal x_2(t)');
grid;

%% PART B

D_n=[-20:20]; 
nleftlim = -20;
nrightlim = 20;
x = 20+1;
for n = [nleftlim:nrightlim];
if n == 0,
    D_n(x) = 0.05;
else
    D_n(n-nleftlim+1) = (sin(n.*pi*0.5)./(n.*pi));
end
end
n = [nleftlim:nrightlim];
W0 = pi/10;
t = -300:300;
s = 300+1;
b=length(t);
x = zeros(1,b);
for t=-300:300
for n=nleftlim:nrightlim
    x(t+s) = x(t+s) + real(D_n(n-nleftlim+1).*exp(n.*1i*W0*t));
end
end
t=-300:300;
plot(t,real(x),'b');
ylabel('x(t)');
xlabel('t');
title('Regenerated A6 2(b) from D_n for -20:20');
legend('Signal x_2(t)');
grid;

%% PART C

D_n=[-50:50];
nleftlim = -50;
nrightlim = 50;
x = 50+1;
for n = [nleftlim:nrightlim];
if n == 0,
    D_n(x) = 0.05;
else
    D_n(n-nleftlim+1) = (sin(n.*pi*0.5)./(n.*pi));
end
end
n = [nleftlim:nrightlim];
W0 = pi/10;
t = -300:300;
s = 300+1;
b=length(t);
x = zeros(1,b);
for t=-300:300
for n=nleftlim:nrightlim
    x(t+s) = x(t+s) + real(D_n(n-nleftlim+1).*exp(n.*1i*W0*t));
end
end
t=-300:300;
plot(t,real(x),'b');
ylabel('x(t)');
xlabel('t');
title('Regenerated A6 2(c) from D_n for -50:50');
legend('Signal x_2(t)');
grid;

%% PART D

D_n=[-500:500];
nleftlim = -500;
nrightlim = 500;
x = 500+1;
for n = [nleftlim:nrightlim];
if n == 0,
    D_n(x) = 0.05;
else
    D_n(n-nleftlim+1) = (sin(n.*pi*0.5)./(n.*pi));
end
end
n = [nleftlim:nrightlim];
W0 = pi/10;
t = -300:300;
s = 300+1;
b=length(t);
x = zeros(1,b); %Pre-allocate Memory for recreated signal
for t=-300:300
for n=nleftlim:nrightlim
    x(t+s) = x(t+s) + real(D_n(n-nleftlim+1).*exp(n.*1i*W0*t));
end
end
t=-300:300;
plot(t,real(x),'b');
ylabel('x(t)');
xlabel('t');
title('Regenerated A6 2(d) from D_n for -500:500');
legend('Signal x_2(t)');
grid;

%%

%for x3(t)

%% PART A

D_n=[-5:5];
nleftlim = -5;
nrightlim = 5;
x = 5+1;
for n = [nleftlim:nrightlim];
if n == 0,
    D_n(x) = 0.025;
else
    D_n(n-nleftlim+1) = (sin(n.*pi*0.25)./(n.*pi));
end
end
n = [nleftlim:nrightlim];
W0 = pi/20;
t = -300:300;
s = 300+1;
b=length(t);
x = zeros(1,b);
for t=-300:300
for n=nleftlim:nrightlim
    x(t+s) = x(t+s) + real(D_n(n-nleftlim+1).*exp(n.*1i*W0*t));
end
end
t=-300:300;
plot(t,real(x),'b');
ylabel('x(t)');
xlabel('t');
title('Regenerated A6 3(a) from D_n for -5:5');
legend('Signal x_3(t)');
grid;

%% PART B
D_n=[-20:20];
nleftlim = -20;
nrightlim = 20;
x = 20+1;
for n = [nleftlim:nrightlim];
if n == 0,
    D_n(x) = 0.025;
else
    D_n(n-nleftlim+1) = (sin(n.*pi*0.25)./(n.*pi));
end
end
n = [nleftlim:nrightlim];
W0 = pi/20;
t = -300:300;
s = 300+1;
b=length(t);
x = zeros(1,b);
for t=-300:300
for n=nleftlim:nrightlim
    x(t+s) = x(t+s) + real(D_n(n-nleftlim+1).*exp(n.*1i*W0*t));
end
end
t=-300:300;
plot(t,real(x),'b');
ylabel('x(t)');
xlabel('t');
title('Regenerated A6 3(b) from D_n for -20:20');
legend('Signal x_3(t)');
grid;

%% PART C
D_n=[-50:50]; 
nleftlim = -50; 
nrightlim = 50; 
x = 50+1;
for n = [nleftlim:nrightlim];
if n == 0,
    D_n(x) = 0.025;
else
    D_n(n-nleftlim+1) = (sin(n.*pi*0.25)./(n.*pi));
end
end
n = [nleftlim:nrightlim]; 
W0 = pi/20; 
t = -300:300;
s = 300+1;
b=length(t);
x = zeros(1,b);
for t=-300:300
for n=nleftlim:nrightlim
    x(t+s) = x(t+s) + real(D_n(n-nleftlim+1).*exp(n.*1i*W0*t));
end
end
t=-300:300;
plot(t,real(x),'b');
ylabel('x(t)');
xlabel('t');
title('Regenerated A6 3(c) from D_n for -50:50');
legend('Signal x_3(t)');
grid;

%% PART D
D_n=[-500:500];
nleftlim = -500;
nrightlim = 500;
x = 500+1;
for n = [nleftlim:nrightlim];
if n == 0,
    D_n(x) = 0.025;
else
    D_n(n-nleftlim+1) = (sin(n.*pi*0.25)./(n.*pi));
end
end
n = [nleftlim:nrightlim];
W0 = pi/20;
t = -300:300;
s = 300+1;
b=length(t);
x = zeros(1,b);
for t=-300:300
for n=nleftlim:nrightlim
    x(t+s) = x(t+s) + real(D_n(n-nleftlim+1).*exp(n.*1i*W0*t));
end
end
t=-300:300;
plot(t,real(x),'b');
ylabel('x(t)');
xlabel('t');
title('Regenerated A6 3(d) from D_n for -500:500');
legend('Signal x_3(t)');
grid;
end

%% Problem B

%The listed problems B.1-B.7 are answered in the provided lab report
%as a PDF file 





