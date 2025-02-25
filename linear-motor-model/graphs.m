%% Simulering 10-100 Hz
clear all;clf;
fig = figure(1);
set(fig, 'Position', [0 0 500 200]);
hold on;

T = readtable('10-100 cgf 3.csv','Format','%f%f%f%f%f%f%f%f', 'HeaderLines',5);

f = T.Var2;
t = T.Var3;
% 3 2 1 eco12 eco14
y = [T.Var4,T.Var6, T.Var7,T.Var5,T.Var8];

markers = ["o","x","+","*","s","d","^","v",">","<","p","h" ];
color = ["blue", "green", "red", "cyan", "magenta"];

y_mean = zeros(size(unique(f),1),size(y,2));

f_i = 1;
for freq=unique(f)'
    x = f(f==freq);
    y_freq = y(f==freq,:);
    for i=1:size(y_freq,2)
        o = (i-2);
        y_wind = -y_freq(:,i);
        if max(abs(y_wind)) > 20
            y_mean(f_i,i) = NaN;
            continue;
        end
        plot(x(11:end)+o,y_wind(11:end), color(i),'HandleVisibility','off');
        y_mean(f_i,i) = mean(y_wind(11:end));
    end
    f_i = f_i+1;
end

for i=1:size(y_mean,2)
    o = (i-2);
    plot(unique(f)'+o,y_mean(:,i), append(color(i),"-", markers(i)));
end

%title("Vågform på spänning över spolarna");
%axis([15/sample 75/sample -0.8 0.8]);
lgd = legend(["3 Lager","2 Lager","1 Lager","Ekonomisk 12","Ekonomisk 14"]);
lgd.Location = 'southeast';
xlabel("Frekvens (Hz)");
ylabel("Kraft (N)");
axis([0 100 0 10]);
print("./img/sim-frekvens-kraft", '-dpng');

%% Simulering 10-100 Hz eco
clear all;clf;
fig = figure(1);
set(fig, 'Position', [0 0 500 200]);
hold on;

T = readtable('10-100 cgf 3.csv','Format','%f%f%f%f%f%f%f%f', 'HeaderLines',5);

f = T.Var2;
t = T.Var3;
y = -T.Var8;

leg = [];
for freq=fliplr(unique(f)')
    plot(t(f==freq).*f(f==freq)*2,y(f==freq));
    leg = [leg sprintf("%d Hz",freq)];
end


legend(leg);
xlabel("Perioder");
ylabel("Kraft (N)");
axis([0 4 0 10]);
print("./img/sim-frekvens-kraft-eco14", '-dpng');

%% Simulering instabilitet

clear all;clf;
fig = figure(1);
set(fig, 'Position', [0 0 250 150]);
hold on;

T = readtable('10-100 cgf 3.csv','Format','%f%f%f%f%f%f%f%f', 'HeaderLines',5);

f = T.Var2;
t = T.Var3;
% 3 2 1 eco12 eco14
y = [T.Var4,T.Var8];

markers = ["o","s"];
color = ["blue", "magenta"];

y_mean = zeros(size(unique(f),1),size(y,2));

f_i = 1;
for freq=unique(f)'
    x = f(f==freq);
    y_freq = y(f==freq,:);
    for i=1:size(y_freq,2)
        o = (i-2);
        y_wind = -y_freq(:,i);
        plot(x(11:end)+o,y_wind(11:end), color(i),'HandleVisibility','off');
        y_mean(f_i,i) = mean(y_wind(11:end));
    end
    f_i = f_i+1;
end

for i=1:size(y_mean,2)
    o = (i-2);
    plot(unique(f)'+o,y_mean(:,i), append(color(i),"-", markers(i)));
end

%title("Vågform på spänning över spolarna");
%axis([15/sample 75/sample -0.8 0.8]);
lgd = legend(["3 Lager","Ekonomisk 14"]);
lgd.Location = 'northwest';
xlabel("Frekvens (Hz)");
ylabel("Kraft (N)");
print("./img/sim-unstable", '-dpng');

%% Simulering instabilitet 1 frekvens

clear all;clf;
fig = figure(1);
set(fig, 'Position', [0 0 250 150]);
hold on;

T = readtable('10-100 cgf 3.csv','Format','%f%f%f%f%f%f%f%f', 'HeaderLines',5);

f = T.Var2;
t = T.Var3;
y = -T.Var4;

markers = ["o","s"];
color = ["blue", "magenta"];

freq = 60;
tf = t(f==freq);
yf = y(f==freq);


yyaxis left; hold on;
plot(tf, yf, append("-", markers(1)));
ylabel("Kraft (N)");
axis([0 0.034 4 6]);

yyaxis right; hold on;
plot(tf, yf, append("-", markers(2)));
ylabel("Kraft (N)");
axis([0 0.034 -1500 2500]);


lgd = legend(["3 Lager, 60 Hz, inzoomad","3 Lager, 60 Hz, utzoomad"]);
lgd.Location = 'southeast';
xlabel("Frekvens (Hz)");
ylabel("Kraft (ln N)");
print("./img/sim-unstable+single", '-dpng');

%% Signal Spänning
clear all;clf;hold on;
dsample = 10;
sample = 10000/dsample;

fig = figure(1);
set(fig, 'Position', [0 0 500 100]);

T = readtable('signal-1v.csv', 'HeaderLines',29);
t = downsample(T.Var2,1);
y = downsample(T.Var3,1)-1.6725;
plot(t,y);
plot(t,sin((t+0.002)*sample/pi)*0.7);
%title("Vågform på spänning över spolarna");
axis([15/sample 75/sample -0.8 0.8]);
lgd = legend(["Genererad signal 50 Hz", "Referenssignal sinus 50Hz"]);
lgd.Location = 'northwest';
ylabel("Spänning (V)");
xlabel("Tid (s)");
commas();
print("./img/mes-sig-spanning", '-dpng');

%% Signal Ström???
% clear all;clf;hold on;
% dsample = 10;
% sample = 10000/dsample;
% T = readtable('ström-signal-10v.csv', 'HeaderLines',29);
% t = downsample(T.Var2,1);
% y = downsample(T.Var3,1);
% y = bandpass(y,[15 17], sample);
% plot(t,y);

%% Kalibrering lastcell
clear all;clf;
fig = figure(1);
set(fig, 'Position', [0 0 500*0.65 200*0.65]);
hold on;
T = readtable('Mätresultat - Lastcell.csv', 'HeaderLines',2);
y = [0; T.Var12(1:21)]/1000;
x = [0; T.Var16(1:21)];

P = polyfit(x,y,1);
yfit = P(1)*x+P(2);
hold on;
plot(x/1000,y,'ro');
plot(x/1000,yfit);
ylabel("Relativ Kraft (N)");
xlabel("Relativ Spänning (V)");
lgd = legend(["Mätvärden","Linjär approximation"]);
lgd.Location = 'northwest';
commas();
print("./img/mes-kalibrering", '-dpng');

%% 48 V - Frekvens kraft
clear all;clf;
fig = figure(1); hold on;
set(fig, 'Position', [0 0 500 200]);

T = readtable('Mätresultat - 48 V.csv', 'HeaderLines',1);
f = T.Var2;
yFrc = T.Var9/1000;

plot(f,yFrc, '-o');

axis([15 50 0 10]);
ylabel("Kraft (N)");
xlabel("Frekvens (Hz)");
legend(["Uppmätt kraft"]);
print("./img/mes-48-kraft-frekvens", '-dpng');

%% 48 V - Ström spänning
clear all;clf;
fig = figure(1); hold on;
set(fig, 'Position', [0 0 500 200]);

T = readtable('Mätresultat - 48 V.csv', 'HeaderLines',1);
f = T.Var2;
yVDC = T.Var3;
yARMS = str2double(strrep(T.Var4,",","."));
yVRMS = str2double(strrep(T.Var5,",","."));

yyaxis left; hold on;
plot(f,yARMS, '-x');
ylabel("Ström (A)");
axis([15 50 0 11]);

yyaxis right; hold on;
plot(f,yVDC, '-o');
plot(f,yVRMS, '-+');
ylabel("Spänning (V)");
axis([15 50 40 51]);

xlabel("Frekvens (Hz)");
lgd = legend([ "Ström Spole (A RMS)", "Matningsspänning (V)", "Spänning Spole (V RMS)"]);
lgd.Location = 'southeast';

print("./img/mes-48-ström-spänning", '-dpng');

%% Temperaturberoende
clear all;clf;
fig = figure(1);
set(fig, 'Position', [0 0 500 200]);
hold on;

T = readtable('Mätresultat - Temperatur.csv','HeaderLines',1);
u = str2double(strrep(T.Var2(47+10:68),",","."));
i = str2double(strrep(T.Var3(47+10:68),",","."));
temp = str2double(strrep(T.Var4(47+10:68),",","."));
t = T.Var5(47+10:68);
f = str2double(strrep(T.Var12(47+10:68),",","."));

yyaxis left; hold on;
plot(temp,f,"-o");
ylabel("Kraft (N)");
commas();

yyaxis right; hold on;
plot(temp,i,"-x");
ylabel("Ström (A)");
axis([35 70 8.5 9]);

xlabel("Temperatur (ºC)");
lgd = legend(["Uppmätt kraft", "Uppmätt Ström (A RMS)"]);
lgd.Location = 'north';
print("./img/mes-temp-kraft", '-dpng');

%% Frekvensanalys Upmätt
clear all;clf;
fig = figure(1);
set(fig, 'Position', [0 0 500 200]);
hold on;
T = readtable('Mätresultat - Frekvensanalys.csv', 'HeaderLines',1);
freq = T.Var2;
a = T.Var4;
forc = T.Var9/1000;

prevF = freq(1);
start = 1;
leg = [];
markers = ["o","x","+","*","s","d","^","v",">","<","p","h" ];
n = 1;
for i=1:length(freq)
    f = freq(i);
    if (f ~= prevF) | (i == length(freq))
        plot(a(start:i-1),forc(start:i-1),append("-", markers(n)));
        leg = [leg sprintf("%d Hz",prevF)];
        start = i;
        n = n + 1;
    end
    prevF = f;
end
lgd = legend(leg);
lgd.Location = 'northwest';

xlabel("Ström (A)");
ylabel("Kraft (N)");
axis([0 7 0 12]);
print("./img/mes-frekvens-kraft", '-dpng');

%% Frekvensanalys Simulerad
clear all;clf;
fig = figure(1);
set(fig, 'Position', [0 0 500 200]);
hold on;

T = readtable('ström cfg 3 10-60 Hz.csv','Format','%f%f%f%f%f%f%f%f%f', 'HeaderLines',5);
f = T.Var2;
a = T.Var3;
t = T.Var4;
y = -T.Var9;

markers = ["o","x","+","*","s","d","^","v",">","<","p","h" ];
color = ["blue", "green", "red", "cyan", "magenta"];

f_i = 2;
leg = [];
plot(0, 0,append("-", markers(1)));
leg = [leg sprintf("5 Hz")];
for freq=unique(f)'
    if freq > 60
        break;
    end
    y_freq = y(f==freq,:);
    a_freq = a(f==freq);
    a_means = zeros(size(unique(a_freq),1),1);
    a_i = 1;
    for amp=unique(a_freq)'
        i = logical((f==freq) .* a==amp);
        y_a = y(i,:);
        if max(abs(y(i))) > 4
            a_means(a_i) = NaN;
        else 
            a_means(a_i) = mean(y_a(11:end));
        end
        a_i = a_i + 1;
    end
    plot(unique(a_freq)',a_means,append("-", markers(f_i)));
    leg = [leg sprintf("%d Hz",freq)];
    f_i = f_i+1;
end

lgd = legend(leg);
lgd.Location = 'northwest';
xlabel("Ström (A)");
ylabel("Kraft (N)");
axis([0 7 0 4]);
print("./img/sim-frekvens-strom-kraft", '-dpng');

%% Mätosäkerhet
clear all;clf;hold on;
sample = 10000;

fig = figure(1);
set(fig, 'Position', [0 0 500 200]);

T = readtable('temp-force.csv', 'HeaderLines',29);
y = [T.Var3(1*sample:3*sample-1); T.Var3(20*sample:22*sample)];
t = linspace(0,4,length(y));
plot(t,y);
%title("Vågform på spänning över spolarna");
axis([0 4 1 1.13]);
lgd = legend(["Lastcell"]);
lgd.Location = 'southeast';
ylabel("Spänning (V)");
xlabel("Tid (s)");
commas();
print("./img/mes-osäkerhet", '-dpng');