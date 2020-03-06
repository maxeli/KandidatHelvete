% kraft beroende av L_s
clf
clc
L_s=linspace(0,1,100);
K_ts = [2 1/2];%[1 0.7 0.5 0.3 0.1]; %tand/slot förhållande

f1=10;% frekvens fixerat vid 10Hz
u_0=4*pi*10^-7;
w=2*pi*f1; % elektrisk vinkelfrekvens?

Wse = 0.03; %statorbredd.

K_ta=1; % anta att transverse effects ?r f?rsumbara. sida 86 LMES
s=1; %slip
p=4;
K_w=0.6;
sigma_i=0;
delta_i=0.03;
K_ti=1;
d_a=0.003;
sigma=37.7*10^6;
g_e=0.005; % nu har vi ansatt ekvivalent luftgap till bara luftgap, se sida 89 i LMES
g_0=g_e;
g=g_e;
tau=@(L_s)L_s./p; % pole pitch, kolla s? den ?r konsistent med boken
sigma_e=sigma*((1/K_ta)+((sigma_i*delta_i)/(sigma*K_ti*d_a)));
I_1=16;
x=@(L_s)(L_s./(13+14*K_ts(i)));
for i=1: length(K_ts)
    N=@(L_s) (L_s./(13+14*K_ts(i))).*(0.05/(pi*(0.001^2)));
    G_e=@(L_s)(u_0.*w.*(tau(L_s).^2).*sigma_e.*d_a)/((pi^2).*g_e); % goodness factor
    X_m=@(L_s)(6.*u_0.*w.*Wse.*tau(L_s).*L_s.*(K_w.*N(L_s)).^2.)/(p.*g.*pi^2); % magnetiseringsreakktans
    % Försvinner inte statorbredden i Xm? Jmf. andra ekv...
    
    R_2=@(L_s) X_m(L_s)./G_e(L_s); % rotorresistans
    F_x=@(L_s) 3.*(I_1.^2).*R_2(L_s)./((s.*2.*tau(L_s).*f1).*((1./s.*G_e(L_s).^2+1)));
    hold on
        if i==1
            plot(L_s,F_x(L_s),'g');
        elseif i == 2
            plot(L_s,F_x(L_s),'-.b');
        end
    %plot(0.4,F_x(0.4),'-x','linewidth',1)

end

axis([0,1,0,15])
ylabel('Thrust [N]');
xlabel('statorlängd [m]');
legend('K_{ts} 2:1','K_{ts} 1:2');
grid on
F_x(0.4);

%% Trådlängd
clf
clc
L_s = linspace(0,1,100);
K_ts = 3/4; %tand:slot faktor

Wse = 0.03; %statorbredd.



K_w = 0.6; %lindningsfaktor

x =@(L_s)(L_s./(13+14*K_ts));

N =@(L_s) (x(L_s).*0.05.*K_w)./(2.*pi.*(0.0005^2));
L_t =@(L_s) N(L_s).*2.*(Wse+x(L_s)).*12.*2;
plot(L_s,L_t(L_s))
grid on, hold on
% plot(L_s,N(L_s))
axis([0 0.5 0 700])
















%% test testsson
clear all
m=3; %antal faser
I1=10;% fasstr?m, tror det ?r RMS. godtyckligt v?rde just nu. enhet: A
S=1; %slip
u0=4*pi*10^-7; % konstant permeabilitet f?r omagnetiska material
f=20; % frekvens, justera enligt behov.
Wse=0.03; % ekvivalent statorbredd. oklart vad skillnaden ?r mot statorbredd. begr?nsas iaf av r?lsens h?jd
Kw=0.9; %Winding factor, oklart hur den ber?knas f?r linj?rmotor. Ansatt till 0.9
N1=200; % Lindningar per slot
r=0.1; % wtf ?r den h?r egentligen. h?rlett till meter. radie p? lindningar?
p=4; %antal poler, totalt
ge= 0.005; % ekvivalent luftgap. samma sak som luftgap? ansatt till 5mm.
Ls=0.5; %statorl?ngd
tau=Ls/p; %statorl?ngd per fas
pr= 26.5*10^-9; %volymresistivitet f?r aluminium. konstant
d= 0.003; % aluminiumets tjocklek
Vs=2*f*tau;

G=2*u0*f*tau^2/(pi*(pr/d)*ge);
Xm=24*u0*pi*f*Wse*Kw*(N1^2)*r/((pi^2)*p*ge);
R2=Xm/G; % rotorresistans
F=(m*((I1)^2)*R2)./(((1/(S*G).^2)+1)*Vs*S);

%% 

m=3; %antal faser
I1=10;% fasstr?m, tror det ?r RMS. godtyckligt v?rde just nu. enhet: A
S=1; %slip
u0=4*pi*10^-7; % konstant permeabilitet f?r omagnetiska material

Wse=0.03; % ekvivalent statorbredd. oklart vad skillnaden ?r mot statorbredd. begr?nsas iaf av r?lsens h?jd
Kw=0.9; %Winding factor, oklart hur den ber?knas f?r linj?rmotor. Ansatt till 0.9
N1=250; % Lindningar per slot
r=0.1; % wtf ?r den h?r egentligen. h?rlett till meter. radie p? lindningar?
p=4; %antal poler, per fas. 2, 4, 6, 8, 10 osv
ge= 0.004; % ekvivalent luftgap. samma sak som luftgap? ansatt till 5mm.
Ls=0.5; %statorl?ngd
f1=1;% frekvens, justera enligt behov.
n=150;
pr= 26.5*10^-9; %volymresistivitet f?r aluminium. konstant
d= 0.003; % aluminiumets tjocklek
F=linspace(0,0,n);
fplot=linspace(0,0,n);
for i=1:n
f=f1*i;
fplot(i)=f;
tau=Ls/p; %statorl?ngd per fas
Vs=2*f*tau;
G=2*u0*f*tau^2/(pi*(pr/d)*ge);
Xm=24*u0*pi*f*Wse*Kw*(N1^2)*r/((pi^2)*p*ge);
R2=Xm./G; % rotorresistans
F(i)=(m*((I1)^2).*R2)./(((1/(S*G).^2)+1).*Vs*S);
end
[Fmax,i1]=max(F)
Vs_m=2*tau*fplot(i1)
fplot(i1)
plot(fplot,F./Fmax);
hold on
grid on








%%  reluktans

Ls=0.5; %
p=4; % poler per fas

t_yoke=0.025;
t_stator=0.03;
A_rals=0.01*t_stator; 
u0=4*pi*10^-7;
u_rals= 200*u0; 
tau=Ls/p;
u_core=4000*u0;
A_yoke=t_stator*t_yoke;
t_tand=0.01;
t_air=0.01;
A_air=(2*t_stator*t_air); % ansatt tand till 1 cm bred
l_air=0.005;
A_tand=2*t_stator*t_tand;

R_tand=t_yoke/(u_core*(A_tand));
R_air=l_air/(u0*A_air);
R_yoke=tau/(u_core*(A_yoke));
R_rals=tau/(u_rals*A_rals);
R_tot=R_rals+R_yoke+2*R_air+2*R_tand;

N=200;
I=10;

phi=2*N*I/R_tot;
B_tand=phi/A_tand;
B_rals=phi/A_rals;
B_air=phi/A_air;

%% DIsponibel area
Da=((Ls/12)-(t_tand))*t_yoke/2;
dk=1*10^-3;
Ak=(dk/2)^2*pi;
Fill=0.5;
N1=Fill*Da/Ak;

%% kopparresistans
% 
p_cu=1.68*10^-8;
l_wire=4*N*(2*0.06+2*t_stator);
A_cu=0.0005^2*pi;
R_cu=p_cu*l_wire/A_cu;
P_cu=I_1^2*R_cu;


%% Kraft beroende av frekvensen
L_s=0.33; %statorlängd
N=((36/(49*12))*(25/pi)).*(L_s)*1000;
f1=linspace(0,10,100);% frekvens
u_0=4*pi*10^-7;
w=@(f1)2*pi*f1; % elektrisk vinkelfrekvens?
K_ta=1; % anta att transverse effects ?r f?rsumbara. sida 86 LMES
s=1; %slip
p=4;
K_w=0.9;
sigma_i=0;
delta_i=0.03;
K_ti=1;
d_a=0.003;
sigma=37.7*10^6;
g_e=0.005; % nu har vi ansatt ekvivalent luftgap till bara luftgap, se sida 89 i LMES
g_0=g_e;
g=g_e;
tau=L_s/p; % pole pitch, kolla s? den ?r konsistent med boken
sigma_e=sigma*((1/K_ta)+((sigma_i*delta_i)/(sigma*K_ti*d_a)));
I_1=10;

G_e=@(f1)(u_0.*w(f1).*(tau^2).*sigma_e.*d_a)/((pi^2).*g_e); % goodness factor
X_m=@(f1)(6.*u_0.*w(f1).*tau.*L_s./(p.*g.*pi^2)).*(K_w*N)^2; % magnetiseringsreakktans
R_2=@(f1) X_m(f1)./G_e(f1); % rotorresistans

F_x=@(f1) 3.*(I_1.^2).*R_2(f1)./((s.*2.*tau.*f1).*((1./s.*G_e(f1).^2+1)));

plot(f1,F_x(f1));
grid on
% % Vs=2*f1*tau;
%fsurf(F_x(L_s,f1),[0,1,0,100]);
