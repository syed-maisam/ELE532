clear;
% Name: Syed Maisam Abbas, Student number: 501103255, Section: 01
%%
% ELE532_LAB4: The Fourirer Trasnfrom: Properties and Applications 

%% Problem B: Application of the Fourier Transform 

%% Problem B.1

% Display Signal
 figure('name','MagSpect Signal');
 MagSpect(xspeech);

% Display Filter
%  figure('name','MagSpect 2kHz Passband');
%  MagSpect(hLPF2000);

% Display Filter
%  figure('name','MagSpect 2.5kHz Passband');
%  MagSpect(hLPF2500);

% Display Channel
% figure('name','MagSpect hChannel');
% MagSpect(hChannel);

t=(1:80000);

%% ENCODE: FREQUENCY SHIFT RIGHT XSPEECH

w0=.4;
xspeech_right=zeros(80000,1);
for t=1:80000
    xspeech_right(t,1)=xspeech(t,1)*exp(1i*w0*t);
end
figure('name','MagSpect xspeech_right');
MagSpect(xspeech_right);
movegui('northeast');

%% ENCODE: FREQUENCY SHIFT LEFT XSPEECH

w0=-.4;
xspeech_left=zeros(80000,1);
for t=1:80000
    xspeech_left(t,1)=xspeech(t,1)*exp(1i*w0*t);
end
figure('name','MagSpect xspeech_left');
MagSpect(xspeech_left);
t=(1:80000);
movegui('northwest');
pause;
close all;

%% ENCODE: APPLY FILTERS TO REPECTIVE SIDE OF THE SIGNAL

xspeech_rightcut=conv(xspeech_right,hLPF2000_right);
xspeech_leftcut=conv(xspeech_left,hLPF2000_left);
figure('name','MagSpect xspeech_rightcut');
MagSpect(xspeech_rightcut);
movegui('northeast');
figure('name','MagSpect xspeech_leftcut');
MagSpect(xspeech_leftcut);
movegui('northwest');

figure('name','MagSpect hLPF2000_right');
MagSpect(hLPF2000_right);
movegui('southeast');
figure('name','MagSpect hLPF2000_left');
MagSpect(hLPF2000_left);
movegui('southwest');

%% ENCODE: ADD THE TWO HALF SINGALS TO CREATE ENCODED SINGAL
xspeech_new=xspeech_leftcut+xspeech_rightcut;
figure('name','MagSpect new');
MagSpect(xspeech_new);

pause;
close all;

%% TRNASMISSION: CONVOLUDE THE SIGNAL WITH PCHANNEL TO TRANSMIT

xspeech_encoded=conv(xspeech_new,hChannel);
figure('name','MagSpect Encoded');
MagSpect(xspeech_encoded);

pause;
close all;

%% DECODE: ISOLATE RESPECTIVE HALFS TO DECODE

xspeech_rightcutdec=conv(xspeech_encoded,hLPF2000_right);
xspeech_leftcutdec=conv(xspeech_encoded,hLPF2000_left);

%% DECODE: FREQUENCY SHIFT RIGHT XSPEECH

w0=-.4;
xspeech_rightdec=zeros(82230,1);
for t=1:82230
    xspeech_rightdec(t,1)=5*xspeech_rightcutdec(t,1)*exp(1i*w0*t);
end
figure('name','MagSpect decoded right');
MagSpect(xspeech_rightdec);
movegui('northeast');

%% DECODE: FREQUENCY SHIFT LEFT XSPEECH

w0=.4;
xspeech_leftdec=zeros(82230,1);
for t=1:82230
    xspeech_leftdec(t,1)=5*xspeech_leftcutdec(t,1)*exp(1i*w0*t);
end
figure('name','MagSpect decoded left');
MagSpect(xspeech_leftdec);
movegui('northwest');

%% DECODE: ADD HALF SINGALS TO RECREATE ORIGNAL SIGNAL

xspeech_reconstruct=xspeech_leftdec+xspeech_rightdec;
figure('name','MagSpect reconstruct');
MagSpect(xspeech_reconstruct);


% sound(xspeech);
% pause;
% sound(xspeech_reconstruct);

pause;
close all;











