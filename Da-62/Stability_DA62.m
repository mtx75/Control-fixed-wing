clear all ;clc;
S=17.28;      %wing area
Cbar=1.24;  %mean chord of the wing
b=14.53;      %span
AR=11.94;      %aspect ratio for wing
clalpha=6.8755;     %app. gradient of cl - alpha 2D
%CLa=clalpha./(1+clalpha./(pi*AR));   %converting cl-a from 2d to 3d using prandl for high AR
CLa=5.608; 
alpha0=-7.5; %alpha for 0 lift taken from airfoil tools website
CL0=0.76547;
%CL0=CLa*-alpha0*pi./180; %coeff. of lift at 0 angle of attack for 3D , pi/180 converting from rad to deg
Xacw=0.25.*Cbar; %assumed aerodynamic center of wing
Xacb=0;
%Xacb=-0.0942.*Cbar; %app. aerodynamic center shift due to fuselage deduced from marcello ref.
Xac=(Xacw+Xacb);   %total aerodyanmic center wing + fuselage
Xcg=0.344.*Cbar;  %app. value exact cant be found
Cmacw=-0.18;       %assumed coeff. of pitching moment at aerodynamic center 
alpha= -10 :2:10;
CLw=CL0 +CLa.*alpha*pi./180; %Coeff. of lift eq. for wing
Cmow=Cmacw+CL0*(Xcg-Xac)./Cbar; %const. term in Coeff. of moment eq. for wing
Cmaw=CLa*(Xcg-Xac)./Cbar; %slope of moment coeff. with alpha
figure (1)
plot(alpha,CLw)
grid minor
xlabel('angle of attack')
ylabel('coeff. of lift due to wing')
alpha2= 0 :15;
Cmw=Cmow+Cmaw.*alpha2.*pi./180;
figure (2)
plot(alpha2,Cmw,'LineWidth',3)
xlabel('angle of attak')
ylabel('coeff. of wing pitching moment ')
ax = gca;
ax.FontSize = 15;
grid minor
%% horizontal tail
ettah=1.07;    %tail eff. range (0.8-1.2) nelson ref.
Sth=2.91;     %area of horizontal tail
lth=4.6;  %app. distance from mean chord of tail to the CG range(3.407-4.654)
vh=Sth.*lth./(S*Cbar); %volume ratio
ebsy0=2*CL0./(pi*AR);   %downwash angle at zero angle of attack
ebsya=2*CLa./(pi*AR);   %slope of downwash with angle of attack
ith=-2 *pi./180;    %incidence angle of horizontal tail relative to longitudinal axis of airplane in rad
clat=5.7295;     %assumed gradient of lift tail coeff. with alpha for low Aspect ratio for 2d naca 0018
bht=3.75;     %app. horizontal wing span
ARht=bht.^2./Sth;   %app. aspect ratio for horizontal wing
%Clat=(pi.*ARht)./(1+sqrt(1+(pi.*ARht./clat)^2)); %conversion from 2D to 3D valid for low Aspect ratio
%Clat=clat./(1+clat./(pi*ARht));
Clat=4.22;
iw=0*pi/180;    %assumed incidence angle of wing in rad
Cmot=ettah.*vh.*Clat.*(ebsy0+iw-ith); %coeff. of tail moment at zero alpha
Cmat=-ettah.*vh.*Clat.*(1-ebsya);     %gradient of tail moment with alpha
a=0:15;
Cmt=Cmot+Cmat.*a.*pi/180;     %coeff. of tail moment equation 
figure(3)
plot(a,Cmt)
xlabel('angle of attack')
ylabel('coeff. of tail moment')
grid minor
Cmwt=Cmot+Cmow +(Cmat+Cmaw).*a.*pi/180;   %total coeff. of moment wing +tail
figure (4)
plot (a,Cmwt,a,Cmw,a,Cmt,'LineWidth',2)
xlabel('angle of attack')
ylabel('coeff. of pitching moment')
legend('wing+tail contribution','wing only','tail only')
ax = gca;
ax.FontSize = 25;
grid minor
%%
CLot=-ettah.*Clat.*(ebsy0+iw-ith).*Sth./S;  %coeff. of tail lift at zero alpha
Clatt=ettah.*Clat.*(1-ebsya).*Sth./S;   %gradient of tail lift due to downwash
CLaT=CLa+Clatt;
CMAT=Cmat+Cmaw;
CLA=CLaT.*a*pi/180+CLot+CL0;
%plot(a,CLA)
%% inaccurate values
Cmafus=0.395509; % slope of pitching moment due to fuselage deduced by nelson ref.
Cmafusr=0.42516; %by using imperical formula from raymer ref
Cmof=-0.045132; %coeff. of pitching moment at zero angle of attack due to fuselage
Cmf=Cmof+Cmafus.*a.*pi./180; % coeff.of pitching moment due to fuselage 
Cmwt=Cmot+Cmow+Cmof +(Cmafus+Cmat+Cmaw).*a.*pi/180;   %total coeff. of moment wing +tail
figure (4)
plot (a,Cmwt,a,Cmw,a,Cmt,a,Cmf,'LineWidth',2)
xlabel('angle of attack')
ylabel('coeff. of pitching moment')
legend('Airplane contribution','wing only','tail only','fuselage only','Location','Best')
ax = gca;
ax.FontSize = 25;
grid minor
%%
Clfus=0.2;   % assumed lift slope due to fuselage 
SM=-(Cmaw+Cmat+Cmafusr)./(CLa+Clatt+Clfus)    %static margin range [0.15- 0.25] for commercial aviation and military transport  
Xn=(SM+Xcg./Cbar).*Cbar       %neutral point where moment is 0 range (0.4988-0.7482)
%% trial 1 
NP=@(XN)[XN-(Cbar./CLa).*(ettah.*Sth.*(5.3-XN).*Clat.*(1-ebsya))./(S.*Cbar)-Xac+(Cmafus*Cbar)./CLa]; 
intialguess= 0;
XNP=fsolve(NP,intialguess);
display(XNP)
Sm=(XNP-Xcg)./Cbar
%% trial 2
Xnp=(((Xac./Cbar)+Clat.*ettah.*Sth*(1-ebsya).*5.334./(Cbar.*CLa.*S))./(1+Clat.*ettah.*Sth.*(1-ebsya)./(CLa.*S)))*Cbar %eq. from marcello ref
sm=(Xnp-Xcg)./Cbar
%% trial 3
XnP=Xac +ettah*vh.*Clat.*(1-ebsya).*Cbar./CLa -(Cmafus.*Cbar)./CLa  %by approximating Xact -Xnp to lt
SM2=(XnP-Xcg)./Cbar
%% elevator 
Tawh=0.44;   %app. value for elevator effectivness from taw - area ratio curve
CLdeltae=ettah.*Sth.*Clat.*Tawh./S;  %slope of coeff. of tail lift with elevator deflection 
%deltaelevator= linspace(-10,10,16);    %variation for different deflection for elevator deflecting downward is +ve deflection
al=-10:10;
for deltaelevator=-20:10:20
   CL= CL0+CLot +(CLa+Clatt+Clfus).*al*pi./180+CLdeltae.*deltaelevator.*pi./180; 
   figure (5)
plot(al,CL,'linewidth',2)
hold on 
legend('\delta=-20','\delta=-10','\delta=0','\delta=10','\delta=20')
xlabel('Angle of attak')
ylabel('coeff. of Lift for Airplane')
ax = gca;
ax.FontSize = 25;
grid minor
end
%%
CMdeltae=-ettah.*vh.*Clat.*Tawh;    %slope of coeff. of tail moment with elevator deflection known as elevator control power
%deltaelevator=1:16;
for deltaelevator=-20:10:20
Cm=Cmot+Cmow+Cmof +(Cmat+Cmaw+Cmafus).*a.*pi/180+CMdeltae.*deltaelevator.*pi./180; 
figure (6)
plot(a,Cm,'linewidth',2)
hold on 
legend('\delta=-20','\delta=-10','\delta=0','\delta=10','\delta=20')
xlabel('Angle of attak')
ylabel('coeff. of pitching for Airplane')
grid minor
ax = gca;
ax.FontSize = 25;
end
%% for trimming during cruise 
V=160*1.852.*5./18;    %appropiate assumption speed is 160knots provided 
W=2279.6695227.*9.81;       %by assuming that the weight doesnt change in cruise
Rhoh=0.7965;              %density for 14000 ft=4267.2m
CLtrim=(2.*W)./(Rhoh.*S.*V.^2);  %lift coeff. for cruising required for triming at specific speed
CMcgtrim=0;            %coeff. of moment at CG for trimming
syms alphatrim deltatrim
eq1=(Cmat+Cmaw+Cmafus).*alphatrim +CMdeltae.*deltatrim==(CMcgtrim-Cmot-Cmow-Cmof).*180./pi;
eq2=(CLa+Clatt+Clfus).*alphatrim+CLdeltae.*deltatrim==(CLtrim-CL0-CLot)*180./pi;
[z,y]= equationsToMatrix([eq1,eq2],[alphatrim,deltatrim]);
x=linsolve(z,y);
v=double(x);
alphatrim=v(1)
deltatrim=v(2)

deltaelevators=-20:10:20;
alphatrimd=(CMcgtrim-Cmot-Cmow-Cmof-CMdeltae.*deltaelevators.*pi./180)./((Cmat+Cmaw+Cmafus).*pi./180) % alpha trim when deflecting elevator with various deflections 
%% to trim a/c at max. angle of attack for landing max forward limit to trim
alphamax=15.*pi./180; %by assuming max. angle of attack is 17 deg
deltamax=-15.*pi./180; %assuming max deflection upward
CMAlphacg=(CMcgtrim-Cmot-Cmow-Cmof-CMdeltae.*deltamax)./alphamax;
Sm2=-CMAlphacg./(CLa+Clatt+Clfus)
XCgmax=Xn-Sm2.*Cbar %max forward limit 
%% another solution for trim
alphamax2=15;
CMOT=Cmot+Cmof+Cmow;
CMATT=CMAT+Cmafus;
deltamax2=-((CMOT+CMATT.*alphamax2*pi./180)./CMdeltae)*180./pi
%%
dlta=-12;
alp=-((CMOT+CMdeltae*dlta*pi/180)./CMATT)*180/pi
%% directional static stability 
lvt=4.9; %app. distance from mean chord of vertical tail to CG
Svt=2.31; %area of vertical tail provided by the manual
vv=(Svt.*lvt)./(S.*b); %volume ratio for vertical tail
sweep=1; %app. sweep angle at quarter chord 
Dmax=1.3136; %app max. diameter of fuselage 
zw=0.3; %app. distance from root of the wing to the FRL - indicates that the wing is low wing
et=0.724+3.06.*(Svt./S)./(1+cosd(sweep))+(0.4.*zw)./Dmax+0.009.*AR; %imperical formula to cal. ettah(1+dsigma/dbeta)
clav=2.*pi;%appropiate assumption as the vertical tail is assumed symmetric
bvt=1.79; % app. span for vT
ARvt=bvt.^2./Svt;   %app. aspect ratio for horizontal wing
k=(pi.*ARvt)./clav;
Clav=(pi.*ARvt)./(1+sqrt(1+k.^2)); %conversion from 2D to 3D valid for low Aspect ratio
CNbeta=vv.*Clav.*et; %slope of yaw moment with side slip +ve for stability
beta1=-10:10;
CN=CNbeta.*beta1.*pi./180;
figure (7)
plot(beta1,CN,'linewidth',2)
xlabel('side slip')
ylabel('coeff. of yaw moment')
ax = gca;
ax.FontSize = 25;
grid minor
%% Rudder control +ve deflection left direction
tawv=0.52; %rudder effectivness parameter by assuming rudder area by 0.74
ettavh=1;% vertical tail efficiency
CNdr=-ettavh.*Clav.*vv.*tawv; %slope of yaw moment with rudder deflection -ve due to +ve def. of rudder resu;t in -ve yaw moment
%CNdr known as rudder control power
vin=46.3; %speed for one engine inoperative in m/sec 90knots
P=0.9*132.*10.^3; %power 90% for inoperative engine in 132kw
T=P./vin; %thrust in newtons
yt=1.68812; %app. moment arm (m) between the operative engine and FRl
Nengine=T.*yt; %yaw moment due to inoperative engine
rhoh2=0.9048; %density at 10000ft =3048 m in kg/m^3
CNd=Nengine./(0.5.*rhoh2.*vin.^2.*S.*b); %coeff. of yaw moment due to engine inoperative to get max rudder deflection
deltarmax=(CNd./CNdr).*180./pi %max calc. rudder deflection 29 left, 29 right
%% Roll static Stability 
%right aileron up for +ve delta aileron
cr=1.86213;      %app.root chord
ct=0.89419;     %app. tip chord
lamda=ct./cr; %app. tip to root ratio 
diangle=5.2; % dihedral angle wrt to y axis 
syms y
c=int(cr.*(1+((lamda-1)./(b.*0.5)).*y).*y,0,b./2); %chord function of y for trapizodial wing 
C=double(c);
CLb=(-2.*CLa.*diangle.*C.*pi)./(180.*S.*b); %slope of Roll moment with side slip -ve for stability
beta2=-10:10;
CLL=CLb.*beta2.*pi./180; %coeff. of roll moment with beta
figure (8)
plot(beta2,CLL,'linewidth',2)
xlabel('side slip')
ylabel('coeff. of Roll moment')
ax = gca;
ax.FontSize = 25;
grid minor
%% aileron 
tawa=0.35; %aileron effectivness parameter for app. avg aileron chord to wing chord is 0.1902
syms Y
c2=int(cr.*(1+((lamda-1)./(b.*0.5)).*Y).*Y,4.802,6.44); %chord function of y for trapizodial wing 
C2=double(c2);
% NOTE 4.5 , 6.47 are the distances between the root chord and the start ,
% end of aileron all dimensions are app. 
CLDA=(2.*CLa.*tawa.*C2)./(S.*b); %Slope of coeff. of roll moment with aileron deflection
% known as aileron control power 
DA=-15:25; %for right aileron max up 25 and 15 down 
CLDAa=CLDA.*DA.*pi./180;
%figure(9)
%plot(DA,CLDAa)

DA2=-25:15; % for left aileron
CLDAa2=CLDA.*DA2.*pi./180;
%figure(10)
%plot(DA2,CLDAa2)
%% dynamic
U1 = 90; %m/s
M=0.3; %mach no.
e=0.777; %oswald eff. assumed
m=2279.66922; %kg
g=9.81;
% mass moment of inertia calulated by assuming the body as cylinder and
% wing as flat plate method deduced from M.Sedery ref.by assuming mass of
% fuselage 10% of total mass and mass of wing of 400kg
IXXB = 7167.02; % (kg.m^2) 
IyyB = 2685.4636; % (kg.m^2) 
IZZB = 9718.5467; % (kg.m^2) 

%Longitudinal Stability Derivatives 
CD0 = 0.041; 
CDu =0; %can be ignored for low subsonic 
%CDalpha = 0.22;
CTxu = -CD0;
q1 = 0.5*U1.^2.*Rhoh; %dynamic pressure at 14000ft in pa
CXu=-(CDu+2*CD0)+CTxu;
CXalpha= CL0*(1-(2.*CLa)./(e.*pi.*AR)); %/rad
CZu=-(M.^2./(1-M.^2)).*CL0 -2.*CL0;
CZalpha=-(CLa+CD0); %/rad
CZalphadot=-2.*Clat.*ettah.*vh.*ebsya; %coeff. of z forces with rate of alpha /rad
CZq=-2*Clat.*ettah.*vh; %change of coeff. of Z-forces WRT non-dimensonal pitch /rad
CZdeltae=-Clat.*Tawh.*ettah.*Sth./S%coeff. of Z- forces WRT elevator angle /rad
Cmalpha=CMATT%change of pitching moment wrt angle of attack
Cmu=0; %change of pitching moment wrt speed  for low speeds assumed 0
Cmalphadot=CZalphadot.*lth./Cbar;
Cmq=CZq.*lth./Cbar;
Cmdeltae=CZdeltae.*lth./Cbar
CDdeltaE=0;
Cmtu=0;
CmTalpha=0;
%longitudinal derivatives 
Xu = ( (  q1) *S* (CXu)) /(m*U1); % /sec
XTu = ( (q1*S* (CTxu))) / (m*U1) ; 
Xalpha = ( (  q1) *S* (CXalpha)) /m; % m/sec^2
Zu = ( (q1) *S* CZu) / (m*U1); % /sec 
Zalpha = (( q1)*S*CZalpha)/m; %(m/sec^2) 
Zalphadot =  (q1*S*Cbar*CZalphadot) /(2*m*U1); % ( m/sec^2) 
Zq =  (q1*S*Cbar*CZq)/(2*m*U1); % ( m/sec) 
Mu = (q1*S*Cbar*(Cmu))/(IyyB*U1); %(rad/sec^2)/(m/sec) 
MTu = (q1*S*Cbar* (Cmtu)) / (IyyB*U1); % (rad/sec^2) / (ft/sec) 
Malpha = (q1*S*Cbar*Cmalpha) /IyyB;  % (/sec^2)
MTalpha = (q1*S*Cbar*CmTalpha)/ IyyB; % {rad/sec^2) /rad 
Malphadot = (q1*S*Cbar^2*Cmalphadot) / (2*IyyB*U1); % (/sec) 
Mq = (q1*S*Cbar^2*Cmq)/(2*IyyB*U1); % (/sec) 
% Longitudinal Dimensional Control Derivatives 
XdeltaE = (q1*S*CDdeltaE) /m; 
ZdeltaE =  (q1*S*CZdeltae) /m; % (m/sec^2 ) /rad 
MdeltaE = (q1*S*Cbar*Cmdeltae) /IyyB; % (/sec^2)  
% conversion of stability dervatives between alpha to w 
Xw=Xalpha./U1;
Zw=Zalpha./U1;
Mw=Malpha./U1;
Mwdot=Malphadot./U1;

% State Variable Representation of the Equations of Motion
% A is stability matrix
A=[Xu Xw 0 -g; Zu Zw U1 0; Mu+Mwdot*Zu Mw+Mwdot*Zw Mq+Mwdot*U1 0; 0 0 1 0] 
Dyn=eig(A)
% B matrix cal. assuming all control derviative in x-force and Mw dot
% neglected

B=[XdeltaE ;ZdeltaE; MdeltaE+Mwdot*ZdeltaE; 0]; %control matrix
% C-matrix output (desired) matrix
C3=[1 0 0 0];
%D-matrix feedforward usually null matrixhe
D=[0];
sys1=ss(A,B,C3,D)
uss=tf(sys1)
C4=[0 1 0 0];
sys2=ss(A,B,C4,D);
wss=tf(sys2)
C5=[0 0 0 1];
sys3=ss(A,B,C5,D);
theta_ss=tf(sys3)
uuss=tf([-10/57.3],[1 1])
T1=tf([-10],[1 10])
TT=T1*theta_ss
uuss=tf([-10/57.3],[1 1])
TT1=feedback(TT,1)
TT2=TT/(1+TT);
TT3=minreal(TT2)
step(TT,TT3)
figure
rlocus(TT)
[RR,KK]=rlocus(TT3);
figure
plot(KK,RR)
figure 
bode(TT3)
figure
nyquist(TT3)
%controlSystemDesigner(TT)
%%
C6=[0 0 1 0];
sys4=ss(A,B,C6,D)
Qss=tf(sys4)
figure
rlocus(uss)
t=0:0.01:300;
for i=1:300
    u(i)=0;
end
for i=301:500
    u(i)=10/57.3;
end
for i=501:2000
    u(i)=0;
end
for i=2001:2400
    u(i)=-10/57.3;
end
for i=2401:30001
    u(i)=0;
end
%step(sys1)
long_input=[t' u'];
%%
xo=[0 0 0 0];
VV=lsim(sys1,u,t,xo);
%%
subplot (211) ; plot(t,VV,'linewidth',2)
title("U vs time ")
xlabel("time")
ylabel("U m/s")
ax = gca;
ax.FontSize = 20;
grid minor

u2=u.*57.3;
subplot (212);
plot(t,u2,'linewidth',2)%elevator deflection wrt time 
title ("elevator deflection vs time ")
xlabel("time")
ylabel("elevator deg")
ax = gca;
ax.FontSize = 20;
grid minor
%%
AL=lsim(sys2,u,t,xo);
subplot (223);plot(t,AL,'linewidth',2) %W velocity wrt elevator deflection
title("W vs time ")
xlabel("time")
ylabel("W m/s")
ax = gca;
ax.FontSize = 20;
grid minor
TH=lsim(sys3,u,t,xo)
TH2=TH.*57.3;
subplot(224); plot(t,TH2,'linewidth',2)
title ("Theta vs time ")
xlabel("time")
ylabel("Theta deg")
ax = gca;
ax.FontSize = 20;
grid minor
%% Coefficients of the NUM(s) of u-Transfer FUnction 
alpha1 = 1.5/57.3; % <rad) 
theta1 = alpha1; 
Au = XdeltaE*(U1- Zalphadot);
Bu = - XdeltaE*(((U1-Zalphadot)*Mq)+Zalpha+(Malphadot*(U1+Zq)) +(ZdeltaE*Xalpha)); 
Cu = (XdeltaE*((Mq*Zalpha)+(Malphadot*g*sin(theta1))- ((Malpha+MTalpha)*(U1+Zq)))+ (ZdeltaE*((-Malphadot*g*cos(theta1))-(Xalpha*Mq))+(MdeltaE*((Xalpha*(U1+Zq))- ((U1 - Zalphadot)*g*cos(theta1))))));
Du = (XdeltaE*(Malpha+MTalpha)*g*sin(theta1))-(ZdeltaE*Malpha*g*cos(theta1))+(MdeltaE*((Zalpha*g*cos(theta1))- (Xalpha*g*sin(theta1))));
NU = [Au Bu Cu Du]; 
% Coefficients of t he NUM(s) of alpha-Transfer FWlction 
Aalpha = ZdeltaE; 
Balpha = (XdeltaE*Zu)+(ZdeltaE*(-Mq-(Xu+XTu))) + (MdeltaE*(U1+Zq)); 
Calpha = (XdeltaE*(((U1+Zq)*(Mu+MTu))-(Mq*Zu)))+(ZdeltaE*Mq*(Xu+XTu))+(MdeltaE*((-g*sin(theta1))-((U1+Zq)*(Xu+XTu)))); 
Dalpha = -(XdeltaE*(Mu+MTu)*g*sin(theta1))+(ZdeltaE*(Mu+MTu)*g*cos(theta1))+(MdeltaE*(((Xu+XTu)*g*sin(theta1))-(Zu*g*cos(theta1)))); 
Nalpha = [Aalpha Balpha Calpha Dalpha] ; 
% Coefficients of t he NUM(s) of theta-Transfer Function 
Atheta = (ZdeltaE*Malphadot) + (MdeltaE* (U1 - Zalphadot) ); 
Btheta = (XdeltaE * ((Zu * Malphadot) + ((U1 - Zalphadot) * (Mu + MTu)))) + (ZdeltaE * ((Malpha + MTalpha) - (Malphadot * (Xu + XTu)))) + (MdeltaE * (-Zalpha - ((U1 - Zalphadot)*(Xu + XTu))));
Ctheta = (XdeltaE * (((Malpha + MTu) * Zu) - (Zalpha * (Mu+ MTu)))) + (ZdeltaE * ((- (Malpha + MTu) * (Xu + XTu)) + (Xalpha * (Mu + MTu)))) + (MdeltaE * ((Zalpha * (Xu + XTu)) - (Xalpha * Zu)));
Ntheta = [Atheta Btheta Ctheta];
% Coefficients of the Longitudinal Characteristic Equation (DBN(s)))
A1 = U1 - Zalphadot;
B1 = - (U1 - Zalphadot) * (Xu + XTu + Mq) - Zalpha - (Malphadot * (U1 + Zq));
C1 = Xu + (XTu * (Mq * (U1 - Zalphadot) + Zalpha + (Malphadot * (U1 + Zq)))) + (Mq * Zalpha) - (Zu * Xalpha) + (Malphadot * g * sin(theta1)) - ((Malpha + MTalpha) * (U1 + Zq));
D1 = g * sin(theta1) * (Malpha + MTalpha - (Malphadot * (Xu + XTu))) + (g * cos(theta1) * ((Zu * Malphadot) + ((Mu + MTu) * (U1 - Zalphadot)))) + ((Mu + MTu) * (-Xalpha * (U1 + Zq))) + (Zu * Xalpha * Mq) + ((Xu + XTu) * (((Malpha + MTalpha) * (U1 + Zq)) - (Mq * Zalpha)));
E1 = (g * cosd(theta1) * ((Malpha + MTalpha) * Zu - (Zalpha * (Mu + MTu)))) + (g * sin(theta1) * ((Mu + MTu) * Xalpha) - ((Xu + XTu) * (Malpha + MTalpha)));
Dbar1= [A1 B1 C1 D1 E1];
%% Simulation 
sys_1=tf (Nalpha, Dbar1); 
alpha=lsim(sys_1,u,t); 
sys_2=tf(NU, Dbar1); 
UU=lsim(sys_2,u,t);
sys_3=tf(Ntheta,Dbar1); 
theta=lsim(sys_3,u,t); 
figure
plot (t, UU, 'g'); 
  title ('Small Perturbation velocity vs.Time'); 
  xlabel ('Time (s)'); 
  ylabel ('Velocity(ft/s)'); 
  grid on; 
  % plot the poles on the s-domain 
poles = roots (Dbar1) ; 
figure, 
plot(poles, '*' ); 
grid on 
title('LOngitudinal Poles in the s- domain'); 
xlabel ('Real AXis'); 
ylabel ('Imaginary AXis' ) ; 
omega_sp=sqrt(abs(poles(1,1)^2)) 
omega_ph=sqrt(abs(poles(3,1)^2))
damp_sp=abs(real(poles(1,1)))/omega_sp 
damp_ph=abs(real(poles(3,1)))/omega_ph 
%%
[N,d]=ss2tf(A,B,C3,D)
G=tf(N,d)
G2=tf(NU,Dbar1)
rlocus(G2)

%% compare between transfer function and State space
Ad=1;
Bd=-Mq-U1*Mwdot-Zw-Xu
Cd=Zw*Mq-U1*Mw-Xw*Zu+Xu*(Mq+U1*Mwdot+Zw)
Dd=-Xu*(Zw*Mq-U1*Mw)+Zu*(Xw*Mq+g*Mwdot)-Mu*(U1*Xw-g)
Ed=g*(Zu*Mw-Mu*Zw)
AU=XdeltaE
BU=-XdeltaE*(Zw+Mq+U1*Mwdot)+ZdeltaE*Xw
CU=XdeltaE*(Zw*Mq-U1*Mw)-ZdeltaE*(Xw*Mq+g*Mwdot)+MdeltaE*(U1*Xw-g)
DU=g*(MdeltaE*Zw-ZdeltaE*Mw)
G3=tf([AU BU CU DU],[Ad Bd Cd Dd Ed])
lsim(G3,u,t)
poles2=roots([Ad Bd Cd Dd Ed])
omega_sp2=sqrt(abs(poles2(1,1)^2)) 
omega_ph2=sqrt(abs(poles2(3,1)^2))
damp_sp2=abs(real(poles2(1,1)))/omega_sp 
damp_ph2=abs(real(poles2(3,1)))/omega_ph 
%%
%step(G3)
u1=1/57.3.*ones(size(t));
lsim(G3,u,t)
%%
rlocus(G3)
rlocus(G2)
%% subplot (121)
%plot(t,velocity)
%plot(t,u)
%% lateral dynamcis (lateral+ directional)
% coeff. due to Y-forces (side force)
ettav=1; %asumed eff. for vertical tail
CYbeta_tail=-Clav*ettav*Svt./S %slope of side force with side slip angle acts on V.tail
C_ybeta=-et*Clav*Svt./S %due to side slip
C_yp=0; %assumed due to very small sweep angle coeff. due to roll rate(p)
C_yr=-2*CYbeta_tail*lvt./b %coeff. due to yawing rate(r)
C_ydeltaa=0; %coeff. due to aileron deflection
C_ydeltar=tawv*Clav*Svt./S %coeff. due to rudder deflection
%% coeff. due to rolling moment (L)
VMAX=90; %max true speed at 14000ft =4267m
RHOH=0.79721036; %density at 4267m
CLift_h=W./(0.5*RHOH*VMAX^2*S) %ref. lift at certain altitude
clbeta_gammaratio=-0.0002; %approximated from curves 
deltaclbeta=-0.0002; %assume max. ordinates at upper surface
C_lbeta=clbeta_gammaratio*diangle*180./pi+deltaclbeta*180./pi %roll due to side slip
C_lp=(-CLa*(1+3*lamda))./(12*(1+lamda)) %roll moment due to roll rate needs to recheck
%C_lp2=(-4*CLa./(S*b^2))*C
Zv=1; %assumed distance between centerline of pressure of vertical tail to FRL
C_lr=CLift_h/4-2*(lvt*Zv*CYbeta_tail)./b^2 %roll due to yaw rate
C_ldeltaa=(2*CLa*tawa*C2)./(S*b)%roll due to aileron deflection
C_ldeltar=(Svt*Zv*tawa*CLa)./(S*b) %roll due to rudder deflection

%% coeff. due to Yawing moment (N)
Length=9.17; %fuse length
Sfs=6.13347; % projected side area of fuselage
XM=2.6; %approximated distance from front nose to CG
h1=1.06; h2=0.5;
KN=0.0005; %wing interference body factor
neu=2.076945*10^-5; %kinamatic viscosity @ 4267m unit m^2/s
ReLfus=VMAX*Length/neu; %fuselage reyonlds 
KRl=1.8; %correction factor
CNBetawingfuse=(-KN*KRl*Sfs*Length./(S*b))*180/pi %wing +fuse contribution to directional stability
C_nbeta=CNBetawingfuse+et*vv*Clav %due to side slip
C_np=-CLift_h/8 %due to roll rate
C_nr=-2*ettav*vv*Clav*lvt./b % due to yaw rate
K=-0.1; %assumed empirical factor 
C_ndeltaa=2*K*CLift_h*CLDA %due to aileron deflection
C_ndeltar=-vv*ettav*Clav*tawv %due to rudder deflection
%Lateral Derivatives
Y_beta = (q1*S*C_ybeta)/m
Y_p = (q1*S*b*C_yp) / (2*m*U1)
Y_r = (q1*S*b*C_yr) / (2*m*U1)
L_beta = (q1*S*b*C_lbeta) / IXXB
L_p = (q1*S*b.^2.*C_lp) / (2*IXXB*U1)
L_r = (q1*S*b^2*C_lr) / (2*IXXB*U1)
N_beta = (q1*S*b*C_nbeta) / IZZB
N_p = (q1*S*b^2*C_np) / (2*IZZB*U1)
N_r = (q1*S*b^2*C_nr) / (2*IZZB*U1)
Ydelta_a = (q1*S*b^2*C_ydeltaa) / (m)
Ydelta_r = (q1*S*C_ydeltar) / (m)
Ndelta_a=(q1*S*b*C_ndeltaa)/IZZB
Ndelta_r=(q1*S*b*C_ndeltar)/IZZB
Ldelta_a=(q1*S*b*C_ldeltaa)/IXXB
Ldelta_r=(q1*S*b*C_ldeltar)/IXXB
A_Lateral=[ Y_beta/U1 Y_p/U1 -(1-(Y_r/U1)) g*cosd(2)/U1 0; L_beta L_p L_r 0 0; N_beta N_p N_r 0 0; 0 1 0 0 0;0 0 1 0 0]
Dyn2 = eig(A_Lateral)
omega_n_DR = sqrt((Y_beta*N_r - N_beta*Y_r + U1*N_beta) / U1)
zeta_DR = (-1/(2*omega_n_DR)) * ((Y_beta + U1*N_r) / U1)

B_lateral=[Ydelta_a/U1 Ydelta_r/U1;Ldelta_a Ldelta_r; Ndelta_a Ndelta_r; 0 0;0 0 ]
C_lateral=eye(5)
Q =[100 0 0 0 0;0 400 0 0 0; 0 0 100 0 0;0 0 0 100 0;0 0 0 0 100]
R=100
[Klqrgains,s,E]=lqr(A_Lateral,B_lateral,Q,R)
eig(A_Lateral-B_lateral*Klqrgains)
Kx=Klqrgains(:,1:3)
Ki=Klqrgains(4)
