%% Load force data
T = readtable('temp-force.csv', 'HeaderLines',29);
t = downsample(T.Var2(20*sample*dsample:end),dsample)-11;
y = downsample(T.Var3(20*sample*dsample:end),dsample)-1.11;

%% Load 20 Hz 6 A data
T = readtable('6a_20hz.csv', 'HeaderLines',29);
t = downsample(T.Var2,dsample);
y = downsample(T.Var3,dsample)-1.0519;

%% Load 20 Hz data
T = readtable('4a_20hz.csv', 'HeaderLines',29);
t = downsample(T.Var2,dsample);
y = downsample(T.Var3,dsample)-1.03;

%% Filter
Y=y;
for i=1:10
    Wo = i*20/(sample/2);  BW = Wo/10;
    [b,a] = iirnotch(Wo,BW);  
    Y = filter(b,a,Y);
end
plot(t,Y);
y=Y;


%% FFT
clf;
L = length(y);
Y = fft(y);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = sample*(0:floor(L/2))/L;

plot(f(1:floor(length(y)/5)),P1(1:floor(length(y)/5)));

%% Difference between 
clf; hold on;
plot(t,y);
step = floor(sample/1000);
steps = floor(length(y)/sample*1000)-1;
err = 0;
t_2 = zeros(steps, 1);
val_rms = zeros(steps, 1);
val_v = zeros(steps, 1);
for i=1:steps
    s = floor(i*step);
    e = floor(s+step);
    d = y(s:e);
    v = abs(max(d)+min(d))/2;
    exp = rms(d);
    diff = v-exp;
    if abs(diff) > abs(err)
        err = diff;
    end
    
    t_2(i) = t(s);
    val_v(i) = v;
    val_rms(i) = exp;
end

plot(t_2, val_v-val_rms);
%plot(t_2, val_rms);


%% AVG Per second
for i=5:14
    %disp(round(t((i)*sample))+28)
    %fprintf(strrep(sprintf("%f\n",y((i)*sample)),'.',','));
    s = floor(i*sample-sample/4);
    e = floor(s+sample/4);
    v = max(y(s:e));
    fprintf(strrep(sprintf("%f\n",v),'.',','));
end
for i=15:5:60
    %disp(round(t((i)*sample))+28)
    %fprintf(strrep(sprintf("%f\n",y((i)*sample)),'.',','));
    s = floor(i*sample-sample/4);
    e = floor(s+sample/4);
    v = max(y(s:e));
    fprintf(strrep(sprintf("%f\n",v),'.',','));
end


%% Bandpass filter
%b = (1/sample)*ones(1,sample);
%a = 1;
%yfilt = filter(b,a,y);
%plot(t(sample:end),yfilt(sample:end));
y = bandpass(y,[20,50],sample);
plot(t,y);
