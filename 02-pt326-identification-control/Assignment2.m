%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Diagnosis and Self Regulation Algorithms - Assignment 2               %
%   MIEB 2019/2020                                                        %
%   Giulio Steyde
%   Fernando Catania
%   Hugo Sebastiao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all, clc, clear all
%% import data

load 'adar_data_S7_v1' 
data_1=y;
load 'adar_data_S7_v2' 
data_2=y;
load 'adar_data_S7_v3'
data_3=y;
clear y

%%  we removed the first two samples from all the signals, given that transients, related to the starting of the acquisitions may occur 

data_1(2:3,1:2)=0;
data_2(2:3,1:2)=0;
data_3(2:3,1:2)=0;


%% Plot data
% figure (1)
% subplot(3,1,1)
% plot(data_1(1,:),data_1(2,:),data_1(1,:),data_1(3,:))
% legend('Input','Output')
% xlabel('time[s]')
% ylabel('Voltage [V]')
% title('data1');
% subplot(3,1,2)
% plot(data_2(1,:),data_2(2,:),data_2(1,:),data_2(3,:))
% title('data2');
% xlabel('time[s]')
% ylabel('Voltage [V]')
% subplot(3,1,3)
% plot(data_3(1,:),data_3(2,:),data_1(1,:),data_3(3,:))
% title('data3');
% xlabel('time[s]')
% ylabel('Voltage [V]')

%% Pre-process your data by removing 1st and 2nd order trends

data_1=[data_1(1,:);detrend(data_1(2,:),'linear',2);detrend(data_1(3,:),'linear',2)];
data_2=[data_2(1,:);detrend(data_2(2,:),'linear',2);detrend(data_2(3,:),'linear',2)];
data_3=[data_3(1,:);detrend(data_3(2,:),'linear',2);detrend(data_3(3,:),'linear',2)];


%% we used iddata to create objects for each data package
Ts=0.08;
obj_1=iddata(data_1(3,:)',data_1(2,:)',Ts);
obj_2=iddata(data_2(3,:)',data_2(2,:)',Ts);
obj_3=iddata(data_3(3,:)',data_3(2,:)',Ts);


%% ARX model
 %Consider a variation between [1,20] for the number of poles and zeros. 
 %For the delay consider values between [1,10].


%V = arxstruc(ze,zv,NN) ze:estimation zv:validation NN: number of different ARX-model structures nn = [na nb nk]
%NN=[na,nb,nk]
na=[1:20]; nb=[1:20]; nk=[1:10];
NN = struc(na,nb,nk); 

figure

% ZE = 1; ZV = 2; ZT = 3
arxmodel123 = arxstruc(obj_1, obj_2,NN);
bestmodel123 = selstruc(arxmodel123,'aic');
model123 = arx([obj_1 ; obj_2],bestmodel123);
subplot(6,1,1)
compare(obj_3, model123)

% ZE = 2; ZV = 1; ZT = 3
arxmodel213 = arxstruc(obj_2, obj_1,NN);
bestmodel213 = selstruc(arxmodel213,'aic');
model213 = arx([obj_2 ; obj_1],bestmodel213);
subplot(6,1,2)
compare(obj_3, model213)

% ZE = 1; ZV = 3; ZT = 2
arxmodel132 = arxstruc(obj_1, obj_3,NN);
bestmodel132 = selstruc(arxmodel132,'aic');
model132 = arx([obj_1 ; obj_3],bestmodel132);
subplot(6,1,3)
compare(obj_2, model132)

% ZE = 3; ZV = 1; ZT = 2
arxmodel312 = arxstruc(obj_3, obj_1,NN);
bestmodel312 = selstruc(arxmodel312,'aic');
model312 = arx([obj_3 ; obj_1],bestmodel312);
subplot(6,1,4)
compare(obj_2, model312)

% ZE = 2; ZV = 3; ZT = 1
arxmodel231 = arxstruc(obj_2, obj_3,NN);
bestmodel231 = selstruc(arxmodel231,'aic');
model231 = arx([obj_2 ; obj_3],bestmodel231);
subplot(6,1,5)
compare(obj_1, model231)

% ZE = 3; ZV = 2; ZT = 1
arxmodel321 = arxstruc(obj_3, obj_2,NN);
bestmodel321 = selstruc(arxmodel321,'aic');
model321 = arx([obj_3 ; obj_2],bestmodel321);
subplot(6,1,6)
compare(obj_1, model321)


%% Second Part

%create new model using all data
new_model=arx([obj_1;obj_2;obj_3],[2 2 3]);
% A(z)y(t) = B(z)u(t) + e(t)
%copy from workspace
%use block Discrete Filter
%G num-> new_model.B [0 0 0 0.057939582979414 0.193667573950829]
%  den-> new_model.A. [1 -0.835566853971648 -0.046530787036599]
%
%D num-> new_model.A [1 -0.835566853971648 -0.046530787036599]
%  den -> subs(B,z,1)- new_model.B [6291/25000 0 0 -0.057939582979414 -0.193667573950829]

%% solve on paper and inplement on Simulink Deadbeat

%% Three Blocks controller: Diophantine equation
syms p0 p1 p2 p3 l0 l1 z

syms p0 p1 p2 p3 l0 l1
eq1 = p0 == 1;
eq2 = -0.835566853971648*p0+p1==-0.72 - 0.42;
eq3 = -0.046530787036599*p0 -0.835566853971648*p1 +p2==0.42*0.72;
eq4 = -0.046530787036599*p1 +p3+0.057939582979414*l0-0.83556685397164*p2==0;
eq5 = -0.046530787036599*p2-0.835566853971648*p3+0.057939582979414*l1+0.193667573950829*l0==0;
eq6 = -0.046530787036599*p3+0.193667573950829*l1==0;

solution = solve([eq1, eq2, eq3, eq4, eq5, eq6], [p0, p1, p2, p3, l0, l1]);
p0=double(vpa(solution.p0));
p1=double(vpa(solution.p1));
p2=double(vpa(solution.p2));
p3=double(vpa(solution.p3));
l0=double(vpa(solution.l0));
l1=double(vpa(solution.l1));

H=((1-0.42)*(1-0.72))/(0.05794+0.1937);

%% RMSE. APPLY THIS CODE FOR EVERY CASE
ref=ScopeData{1}.Values.Data;
sig=ScopeData{2}.Values.Data;
RMSE = sqrt(mean((sig - ref).^2));

%without the transient
ref_notr=ScopeData{1}.Values.Data(round(end/10):end);
sig_notr=ScopeData{2}.Values.Data(round(end/10):end);
RMSE_notr = sqrt(mean((sig - ref).^2));



