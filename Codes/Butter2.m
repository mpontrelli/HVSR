%% Butter2
  % a Bandpass butterworth filter to detrend groundmotion microtremor data
    %INPUTS
    %x - waveform
       
    %OUTPUTS
    %y - butterworth filtered waveform

function [y] = Butter2(x)
% Here are the filter parameters
LowCorner=0.5;
%HighCorner=9;
HighCorner=49;
Npoles=4;  % Corner for 1 pass of the two-pass filter

fN=100/2;
Lowcut=LowCorner/fN;
Highcut=HighCorner/fN;
[bb1 aa1]=butter(Npoles,[Lowcut Highcut]);
%freqz(bb1,aa1)
y=filtfilt(bb1,aa1,(x-mean(x)));

freqz(bb1,aa1)
% plot filter
fs = 100;
[h, f] = freqz(bb1,aa1, 256, fs);

%compute magnitude squared of the frequency response
mag = abs(h);
mag = mag.^2;
ydb = 20*log10(mag);
figure
subplot(2,1,1)
plot(f,ydb)
title('Magnitude response')
xlabel('frequency (Hz)')
ylabel('|H(\omega)|^{2}')
ylim([-100 20])
xlim([0.2 49])
set(gca, 'XScale', 'log')
grid on

%compute phase
ang = angle(h);
ang = (ang*180)/pi;
xlim([0.2 49])
subplot(2,1,2)
Fig3 = plot(f, ang);
title('Phase')
xlabel('frequency (Hz)')
ylabel('Phase')
set(gca, 'XScale', 'log')
xlim([0.2 49])
grid on
%saveas(Fig3, 'Fig3.jpg')

end