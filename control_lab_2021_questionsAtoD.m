clear all
close all

%parametroi sustimatos
var_1=20.6;
var_2 = -0.5;
A=[0 1 0 0; var_1 0 0 0;0 0 0 1; var_2 0 0 0];
B=[0;-1;0;0.5];
C=[1 0 0 0;0 0 1 0];
D=[0;0];
x0=[-0.2 -0.06 0.01 0.3];
t=0:0.01:8;
r=ones(size(t)); %eisodos tou neou systimatos 


states={'th','th_der','x','x_der'};
input={'u'};
output={'y1','y2'};

sys = ss(A,B,C,D,'statename',states,'inputname',input,'outputname',output); %arxiko systima

%controllability
C_array=ctrb(A,B);
C_rank=rank(C_array); %dinei 4, ara systima controllable

%observability
Q_array=obsv(A,C);
Q_rank=rank(Q_array); %dinei 4, ara systima observable

%Erwtima A-compute
poles_before=pole(sys); %enas astathis polos

z=0.5;
ts=1.2; %epilogi diki mas, prepei na einai <2sec
a=4; %epilogi gia na exoume 2% tis telikis timis
wn=a/z*ts;
wd=wn*sqrt(1-z^2);
p1=-z*wn+wd*1i; %epithimitoi poloi
p2=-z*wn-wd*1i;
p3=-50; %oi p3,p4 epilegontai wste na min epireazoun tin sxediasi(diladi me Re<0)
p4=-100;
poles_after=[p1 p2 p3 p4]; %oloi eustatheis

K_acker=acker(A,B,poles_after); %vriskoume kerdos mesw typou ackermann

A_new=A-B*K_acker;
sys_A=ss(A_new,B,C,D,'statename',states,'inputname',input,'outputname',output);
[y_A,t,x_A]=lsim(sys_A,r,t,x0);

%Erwtima A-display
figure;
plot(t,x_A);
legend('theta','theta_{der}','x','x_{der}');
title('State Response for pole selection');
xlabel('Time(sec)');
ylabel('Amplitude');

figure;
u_A = - K_acker* x_A'; %vriskoume antistrofo gia na ginei o pollaplasiasmos
plot(t, u_A);
title('Input u(t) for pole selection');
xlabel('Time (sec)');
ylabel('Amplitude');

%erwtima B-compute
Q=eye(4);%xol_T*xol=x1^2+x2^2+x3^2+x4^2
R=1;%logw u^2
N=0;
[K_B,S,e] = lqr(sys,Q,R,N);
A_new2=A-B*K_B;
sys_B=ss(A_new2,B,C,D,'statename',states,'inputname',input,'outputname',output);


%erwtima B-display
figure;
[y_B,t,x_B]=lsim(sys_B,r,t);
plot(t, x_B);
legend('theta','theta_{der}','x','x_{der}');
xlabel('Time (sec)');
ylabel('Amplitude');
title('State Response for LQR');

u_B = - K_B * x_B';
figure;
plot(t, u_B);
title('Input u(t) for LQR');
xlabel('Time (sec)');
ylabel('Amplitude');

%erwtima C - compute 
xr = [0 0 1 0]';
xf = inv(A - B * K_acker) * B * K_acker * xr;
u_C = -K_acker * x_A' - K_acker * xr;

x_C = x_A;
[rows, cols] = size(x_C);
for i = 1 : rows
 x_C(i, :) = x_C(i, :) + xf';
end

p = u_C'.*x_C(:,4); %4hgrammi einai taxutita x_der, kai Power=Force*Velocity

%erwtima C- display
figure;
plot(t, u_C);
title('Input u(t) to send system to desired position');
xlabel('Time (sec)');
ylabel('Amplitude');

figure;
plot(t, x_C);
title('State Response to send system to desired position');
legend('theta', 'theta_{der}', 'x', 'x_{der}');
xlabel('Time (sec)');
ylabel('Amplitude');

figure;
plot(t, p);
title('Power P(t) to send system to desired position');
xlabel('Time (sec)');
ylabel('Watts');

%erwtima D 

%DA -compute
x0_DA=[x0 [0,0,0,0]]; %Ypothetoume arxiko lathos mideniko
L=place(A',C',poles_after)'; %epilgoi polwn wste to e_der=(A-LC)e na nai xei Re(si)<0
A_DA = [(A-B*K_acker) (B*K_acker); zeros(size(A)) (A-L*C)];
B_DA = [B; zeros(size(B))];
C_DA = [C zeros(size(C))];
D_DA = [0; 0];
state_DA = {'theta' 'theta_dot' 'x' 'x_dot' 'e1' 'e2' 'e3' 'e4'};
input_DA = {'u'};
sys_DA = ss(A_DA,B_DA,C_DA,D_DA,'statename',state_DA,'inputname',input_DA,'outputname',output);
[y_DA,t,x_DA]=lsim(sys_DA,r,t,x0_DA);


%Erwtima DA-display
figure;
plot(t,y_DA);
legend('theta','x');
title('Output Response for pole selection on Luenberg Observer');
xlabel('Time(sec)');
ylabel('Amplitude');

%erwtima DB -compute
x0_DB=x0_DA;
A_DB = [(A-B*K_B) (B*K_B); zeros(size(A)) (A-L*C)]; %Xrisimopoioume to kerdos tou erwtimatos B
sys_B=ss(A_new2,B,C,D,'statename',states,'inputname',input,'outputname',output);
state_DB = {'theta' 'theta_dot' 'x' 'x_dot' 'e1' 'e2' 'e3' 'e4'};
input_DB = {'u'};
sys_DB = ss(A_DB,B_DA,C_DA,D_DA,'statename',state_DB,'inputname',input_DB,'outputname',output);
[y_DB,t,x_DB]=lsim(sys_DB,r,t,x0_DB);


%erwtima DB-display
figure;
plot(t, y_DB);
legend('theta','x');
xlabel('Time (sec)');
ylabel('Amplitude');
title('Output Response for LQR on Luenberg Observer');












