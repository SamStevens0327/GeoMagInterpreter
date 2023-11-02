% HER SQUID from slowdata - natural oscillations, xcorr with HBK/KMH
% ver20230728 1349
clear
close all
clc

%%%%% SET AND DOUBLE CHECK ALL following PARAMS %%%%%%
% ---------------start----------------------%%
% ------------------------------------------%%

%cd 'I:\refill-1+1';  %working main directory. the files for -1 0 +1 days
%after refill should reside in a subdirectory i.e. in 2023-12-31 if refill
% was on 31st.

%### HER, HBK, and KMH files download ###%
% https://imag-data.bgs.ac.uk/GIN_V1/GINForms2
%- choose date, period:second data, then view, then download IAGA2002
% !! caution in the vsec files, replace : semicolons with empty space in
% a editor% otherwise loading will fail !!

year = 2022; % REFILLDAY
month = 04; % date when refill happened
day = 10; %

whichsite = 'kmh';  % correlate with hbk or kmh - will load the vsec data.
whichday = 1; %% 0, +1, or -1 days relative to refill.
secondsdifference = -41; % difference from true 00 00 00 file start see
% the header of LVM file!! -102 secs is for a file from 23:58:18

% ---------------end----------------------%%
% ----------------------------------------%%


fs = 2; % 2Hz sampling rate

t0 = datenum(year, month, day);

switch(whichday)
   case 0
      t = t0;
   case 1
      t = t0 + 1;
   case -1
      t = t0 - 1;
end

filename = 'slowDataUTC_23-07-03_0000.lvm'; % hersquid
filename2 = 'hbk20230703.txt'; % hbk
filename3 = 'her20230703.txt'; % her

opts = delimitedTextImportOptions('VariableNamesLine', 0, 'DataLines', 24);
T = readtable(filename, opts);
SLOW = table2array(T(:, 1:10));

SQUID = circshift(SLOW(:, 2:4), secondsdifference); % SQUID, Hhz
CTUMAG = circshift(SLOW(:, 5:7), secondsdifference); % CTUMAG, HhZ

opts = delimitedTextImportOptions('VariableNamesLine', 0, 'DataLines', 19);
T = readtable(filename2, opts);
hbkk = table2array(T(:, 4:6)); % 1 second HBK HDZ
hbkH = hbkk(:, 1);
hbkh = tan(pi .* (hbkk(:, 2) / 60) / 180) .* hbkH;  % mins to degs, degs to rads, D to h.
hbkZ = hbkk(:, 3);

opts = delimitedTextImportOptions('VariableNamesLine', 0, 'DataLines', 19);
T = readtable(filename3, opts);
herr = table2array(T(:, 4:6)); % 1 second HER HDZ
herH = herr(:, 1);
herh = tan(pi .* (herr(:, 2) / 60) / 180) .* herH;  % mins to degs, degs to rads, D to h.
herZ = herr(:, 3);

hbkH(hbkH == 0) = NaN; hbkH(hbkH == 9999) = NaN;
hbkh(hbkh == 0) = NaN; hbkh(hbkh == 9999) = NaN;
hbkH = fillmissing(hbkH, 'linear'); hbkh = fillmissing(hbkH, 'linear');

SQUIDZ = fillmissing(SQUID(:, 3), 'linear');
SQUIDH = fillmissing(SQUID(:, 1), 'linear');
sqH = 0.685 * decimate(SQUIDH, fs, 'fir'); % HEYA HEYA 0.685 factor for SQH!
sqZ = decimate(SQUIDZ, fs, 'fir');

herH(herH == 0) = NaN; herH(herH == 9999) = NaN;
herh(herh == 0) = NaN; herh(herh == 9999) = NaN;
herH = fillmissing(herH, 'linear'); herh = fillmissing(herH, 'linear');

hbk = [hbkH hbkh hbkZ];
her = [herH herh herZ];

for n = 1:length(hbk)
    timehour(n) = n / 3600;
end

bpFilt = designfilt('bandpassiir', 'FilterOrder', 20, ...
    'HalfPowerFrequency1', 3e-3, 'HalfPowerFrequency2', 0.15, ...
    'SampleRate', 1); % avoid 5mHz and its 2nd harmonic

fsqH = filtfilt(bpFilt, sqH - mean(sqH));
fsqH(1:2000) = 0; fsqH(end - 550:end) = 0;
fhbkH = filtfilt(bpFilt, hbkH - mean(hbkH));
fhbkH(1:2000) = 0; fhbkH(end - 550:end) = 0;
fherH = filtfilt(bpFilt, herH - mean(herH));
fherH(1:2000) = 0; fherH(end - 550:end) = 0;

% Plotting code continues...

% figure 3
figure;
plot(fhbkH);
hold;
plot(fsqH);
title([whichsite 'H and SquidH hipassed']);
legend(whichsite, 'SQUID');
xlabel('Seconds');
ylabel('Nano Tesla');

% figure 4
figure;
plot(fhbkH);
hold;
plot(fherH);
title([whichsite 'H and HerH hipassed']);
legend(whichsite, 'her');
xlabel('Seconds');
ylabel('Nano Tesla');

Nfft = 1 * 2048;

% figure 5
[NNXX1, F] = pwelch(fhbkH(1000:end) - mean(fhbkH(1000:end)), Nfft, [], Nfft, 1);
figure;
loglog(F, sqrt(NNXX1));
hold;
grid;
[NNXX1, F] = pwelch(fherH(1000:end) - mean(fherH(1000:end)), Nfft, [], Nfft, 1);
loglog(F, sqrt(NNXX1));
[NNXX1, F] = pwelch(fsqH(1000:end) - mean(fsqH(1000:end)), Nfft, [], Nfft, 1);
loglog(F, sqrt(NNXX1));
legend(whichsite, 'herH', 'sqH');
xlabel('freq');
ylabel('amplitude asd');

decfact = 5; % 2Hz/10=200mHz
nfft = 2048;
Fs = 1;

% figure 6
figure;
plot(hbkH - mean(hbkH));
hold;
plot(herH - mean(herH));
plot(sqH - mean(sqH));
grid;
hold;
legend('HBK', 'HER', 'SQUID');
xlabel('Seconds');
ylabel('NanoTesla');

% figure 7
figure;
plot(timehour, fhbkH);
hold;
plot(timehour, fherH);
plot(timehour, fsqH);
grid;
hold;

[X, Z, ~] = alignsignals(fhbkH, fsqH, 10, 'truncate');
xspectrogram(decimate(X, decfact, 'fir'), decimate(Z, decfact, 'fir'), 128, [], Nfft, Fs / decfact, 'yaxis');
clim([-50 30]);
colormap('jet');

% figure 8
[X, Z, Y] = alignsignals(fhbkH, fsqH, 10, 'truncate');
figure;
plot(X .* Z);

% figure 9
% figure;
% plot(moving((X .* Z), 500));
% ylim([-0.01 0.25]);
% title('cross product lowpass filtered -quazisyncdet');

% MGsquared coherence
clear CXY30;
Nfft = 512;
length = 3600;  % slice len in seconds.
iterations = 24;
decfact = 1;
Fs = 1 / decfact;

for i = 0:(iterations - 1)  % iteration no
    zac = length * i + 1;
    kon = length * i + length;
    if kon > 86384
        kon = 86384;
    end

    XX = fhbkH(zac:kon);
    ZZ = fherH(zac:kon);
    [X, Z, Y] = alignsignals(XX, ZZ, 7, 'truncate'); % align signal, max lag 10s.

    % figure 10
    figure;
    plot(X);
    hold;
    plot(Z);
    title('signals in the window');

    Fs = 1;
    % coherence
    [Cxy, F] = mscohere(X, Z, 300, 280, Nfft, Fs);

    if i == 0
        CXY30 = Cxy;
    else
        CXY30 = [CXY30 Cxy];  % append
    end
end

% figure 11
figure;
surf(CXY30(1:200, 2:(iterations - 1)));
view(0, 90);
title('Coherence');
ylabel('Frequency in mHz');
durationmins = sprintf(' %d', length / 60);
xlabel([strcat('Time in', durationmins, '-min slices')]);

% MANUALLY MGsquared coherence if required

% close all
i = 10;            % iteration no
length = 3600;  % slice len in seconds.
zac = length * i + 1;
kon = length * i + 1 + length;

XX = fhbkH(zac:kon);
ZZ = (fherH(zac:kon));
[X, Z, D] = alignsignals(XX, ZZ, 10, 'truncate'); % alighn signal, max lag 10s.
figure;
plot(X);
hold;
plot(Z);
title('signals in the window');

Nfft = 1024;

Fs = 1;
% coherence
[Cxy, F] = mscohere(X, Z, 300, 200, Nfft, Fs);

% figure 12
figure;
plot(F, Cxy);
grid;
title('mag squared coh hbk and her squid H axis');
xlim([2e-3 0.18]); % limit X to bpf corner.
ylim([0 1]);
