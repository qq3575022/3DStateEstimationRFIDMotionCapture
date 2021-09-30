clc, clear, close all
data=readtable('2.csv','Delimiter', ',');  g=9.7953;
load('RF.mat');
% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

time   = table2array(data(:,1));               
data_x = table2array(data(:,3));  data_y = table2array(data(:,4));  data_z = table2array(data(:,5));

acc_time   = time(find_acc_data);  acc_time   = (acc_time - acc_time(1))/1000000000;                
gyro_time  = time(find_gyro_data); gyro_time  = (gyro_time - gyro_time(1))/1000000000;
mag_time   = time(find_mag_data);  mag_time   = (mag_time - mag_time(1))/1000000000;

% start time
xS = find(abs(acc_time-47.9531)<0.002); xS = xS(1);
xE = find(abs(acc_time-51.4569)<0.001); xE = xE(1);

xxS = find(abs(gyro_time-47.9531)<0.002); xxS = xxS(1);
xxE = find(abs(gyro_time-51.4569)<0.001); xxE = xxE(1);

xSS = find(abs(mag_time-47.9531)<0.007); xSS = xSS(1);
xEE = find(abs(mag_time-51.4569)<0.009); xEE = xEE(1);

yS = find(abs(acc_time-46)<0.002); yS = yS(1);
yE = find(abs(acc_time-48.5)<0.002); yE = yE(1);

yyS = find(abs(gyro_time-46.1784)<0.002);  yyS = yyS(1);
yyE = find(abs(gyro_time-48.5)<0.002); yyE = yyE(1);

ySS = find(abs(mag_time-46.1784)<0.009);   ySS = ySS(1);
yEE = find(abs(mag_time-48.5)<0.009); yEE = yEE(1);

zS = find(abs(acc_time-95.5)<0.002); zS = zS(1);
zE = find(abs(acc_time-98)<0.002); zE = zE(1);

zzS = find(abs(gyro_time-95.5027)<0.002); zzS = zzS(1);
zzE = find(abs(gyro_time-98)<0.002); zzE = zzE(1);

zSS = find(abs(mag_time-95.5027)<0.009); zSS = zSS(1);
zEE = find(abs(mag_time-98)<0.009); zEE = zEE(1);

%time = unique(sort([accT; gyroT; magT]),'rows');

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = data_x(find_acc_data);
acc_data_y = data_y(find_acc_data)-g;
acc_data_z = data_z(find_acc_data);

%%
gyro_data_x = data_x(find_gyro_data); gyro_data_y = data_y(find_gyro_data); gyro_data_z = data_z(find_gyro_data);
mag_data_x  = data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);

accx = acc_data_x(xS:xE) - mean(acc_data_x(xS:xE));
accy = acc_data_y(xS:xE) - mean(acc_data_y(xS:xE));
accz = acc_data_z(xS:xE) - mean(acc_data_z(xS:xE));
accT = acc_time(xS:xE);

gyrox = gyro_data_x(xxS:xxE) - mean(gyro_data_x(xxS:xxE));
gyroy = gyro_data_y(xxS:xxE) - mean(gyro_data_y(xxS:xxE));
gyroz = gyro_data_z(xxS:xxE) - mean(gyro_data_z(xxS:xxE));
gyroT = gyro_time(xxS:xxE);

magx = mag_data_x(xSS:xEE) - mean(mag_data_x(xSS:xEE));
magy = mag_data_y(xSS:xEE) - mean(mag_data_y(xSS:xEE));
magz = mag_data_z(xSS:xEE) - mean(mag_data_z(xSS:xEE));
magT = mag_time(xSS:xEE);
% ========== Start End Time Acceleration ========== 


time = unique(sort([accT; gyroT; magT]),'rows');

%
[SimPP1, SimVV1, SimAA1] = groundtruth1Dx(accT-accT(1));

[PP1, VV1, AA1] = groundtruth1Dx2(time - time(1));

% Start from 0
% SimAA1 = SimAA1 + 0.14273*ones(1, length(acc_time));
% SimAA1 = awgn(SimAA1,35,'measured');
% 
% PP1 = PP1 - PP1(1) + 1.03;

% 1D - 3D coordinates
coord3 = zeros(3, length(time));
%
% x
coord3(1,:) = PP1;

%y;
coord3(2,:) = 1.31;
coord3(3,:) = 1.09;
%
% Parameters of tag
Gt = 14.62;    % tag's antenna gain
X = 0.85;      % polarization mismatch
M = 4;         % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

x1 = [0,    0,  0.865];%[2.6256, 0.0889,0.858];%[-0.2,-0.2, 0.4];
x2 = [2.29, 0,  1.27];%[ 0.9, 0.8, 0.2];
x3 = [2.29,2.52, 0.865];%[ 0.8,-0.3,-0.2];
x4 = [0, 2.52,  1.27];

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.00012; 

% phase cconcatenation
% global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 

% Get RSS and phase from each reader observing the moving tag
z = NaN(3,length(coord3)-1); z_prev = NaN(3,length(coord3)-1);

z_prev(1,:) = coord3(1,1:end-1); z(1,:) = coord3(1,2:end);% x coordinate
z_prev(2,:) = coord3(2,1:end-1); z(2,:) = coord3(2,2:end);% y coordinate
z_prev(3,:) = coord3(3,1:end-1); z(3,:) = coord3(3,2:end);% z coordinate

zdiff1 = z(1,:)-z_prev(1,:);
zdiff2 = z(2,:)-z_prev(2,:);
zdiff3 = z(3,:)-z_prev(3,:);

%
H1 = NaN(1,length(PP1));     phi1 = NaN(1,length(PP1));       phi_mu1 = NaN(1,length(PP1));
H2 = NaN(1,length(PP1));     phi2 = NaN(1,length(PP1));       phi_mu2 = NaN(1,length(PP1));
H3 = NaN(1,length(PP1));     phi3 = NaN(1,length(PP1));       phi_mu3 = NaN(1,length(PP1));
H4 = NaN(1,length(PP1));     phi4 = NaN(1,length(PP1));       phi_mu4 = NaN(1,length(PP1));

r_sim1 = NaN(1,length(PP1)); r_sim1_ = NaN(1,length(PP1)); rdot_sim1 = NaN(1,length(PP1)); rdot_sim1_ = NaN(1,length(PP1));  diff1 = NaN(1,length(PP1));
r_sim2 = NaN(1,length(PP1)); r_sim2_ = NaN(1,length(PP1)); rdot_sim2 = NaN(1,length(PP1)); rdot_sim2_ = NaN(1,length(PP1));  diff2 = NaN(1,length(PP1));
r_sim3 = NaN(1,length(PP1)); r_sim3_ = NaN(1,length(PP1)); rdot_sim3 = NaN(1,length(PP1)); rdot_sim3_ = NaN(1,length(PP1));  diff3 = NaN(1,length(PP1));
r_sim4 = NaN(1,length(PP1)); r_sim4_ = NaN(1,length(PP1)); rdot_sim4 = NaN(1,length(PP1)); rdot_sim4_ = NaN(1,length(PP1)); diff4 = NaN(1,length(PP1));

for k = 1:1:length(time)-1  
[H1(k+1),phi_mu1(k+1),r_sim1(k+1),r_sim1_(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),diff1(k+1)] = noisysim(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,time(k+1)-time(k));
[H2(k+1),phi_mu2(k+1),r_sim2(k+1),r_sim2_(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),diff2(k+1)] = noisysim(x2,f2,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,time(k+1)-time(k));
[H3(k+1),phi_mu3(k+1),r_sim3(k+1),r_sim3_(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),diff3(k+1)] = noisysim(x3,f3,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,time(k+1)-time(k));   
[H4(k+1),phi_mu4(k+1),r_sim4(k+1),r_sim4_(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),diff4(k+1)] = noisysim(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,time(k+1)-time(k));  
end

% figure
% subplot(4,1,1),plot(time, H1,'LineWidth',2);title('Simulated H from Reader $\#1$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
% subplot(4,1,2),plot(time, H2,'LineWidth',2);title('Simulated H from Reader $\#2$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
% subplot(4,1,3),plot(time, H3,'LineWidth',2);title('Simulated H from Reader $\#3$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
% subplot(4,1,4),plot(time, H4,'LineWidth',2);title('Simulated H from Reader $\#4$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')

%r_sim1_ = r_sim1+0.1*rand(1, length(r_sim1))-0.5*0.1; r_sim2_ = r_sim2+rand(1, length(r_sim2)); r_sim3_ = r_sim3+rand(1, length(r_sim3)); r_sim4_ = r_sim4+rand(1, length(r_sim4));
r_sim1(1) = r_sim1(2); r_sim2(1) = r_sim2(2); r_sim3(1) = r_sim3(2); r_sim4(1) = r_sim4(2); 
rdot_sim1(1) = rdot_sim1(2); rdot_sim2(1) = rdot_sim2(2); rdot_sim3(1) = rdot_sim3(2); rdot_sim4(1) = rdot_sim4(2); 
figure
subplot(8,1,1),plot(time, r_sim1_,'LineWidth',2);hold on;plot(time, r_sim1,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_1$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,2),plot(time, rdot_sim1_,'-','LineWidth',2);hold on;plot(time, rdot_sim1,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_1$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,3),plot(time, r_sim2_,'LineWidth',2);hold on;plot(time, r_sim2,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_2$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,4),plot(time, rdot_sim2_,'-','LineWidth',2);hold on;plot(time, rdot_sim2,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_2$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,5),plot(time, r_sim3_,'LineWidth',2);hold on;plot(time, r_sim3,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_3$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,6),plot(time, rdot_sim3_,'-','LineWidth',2);hold on;plot(time, rdot_sim3,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_3$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,7),plot(time, r_sim4_,'LineWidth',2);hold on;plot(time, r_sim4,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_4$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,8),plot(time, rdot_sim4_,'-','LineWidth',2);hold on;plot(time, rdot_sim4,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_4$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
% Nonlinear - Extended kalman filter
% Load ground truth orientation, angular velocity, acc_x, acc_y, acc_z
%[phi,phi_gt,AXY,XX,Tmag,tdmag,Tacc,tdacc,T3,td3] = groundtruth3D();


mag1 = medfilt1(mag1,120);
mag2 = medfilt1(mag2,120);
mag3 = medfilt1(mag3,120);
mag4 = medfilt1(mag4,120);

time1 = t1(21511:38779) +13.5167; time3 = t3(21511*2:38779*2)+13.5167;
mag11 = mag1(21511:38779); mag22 = mag2(21511:38779); mag33 = mag3(21511*2:38779*2); mag44 = mag4(21511:38779); 

%
figure
subplot(411), plot(t1, mag1)
subplot(412), plot(t1, mag2)
subplot(413), plot(t3, mag3)
subplot(414), plot(t1, mag4)
%%
phase11 = phase1(21511:38779); phase22 = phase2(21511:38779); phase33 = phase3(21511*2:38779*2); phase44 = phase4(21511:38779); 


phaseDiff1 = diff(phase11); phaseDiff2 = diff(phase22);  phaseDiff3 = diff(phase33);  phaseDiff4 = diff(phase44); 

phaseDiff1 = -medfilt1(phaseDiff1,800)*1/(42e-03);
phaseDiff2 = -medfilt1(phaseDiff2,800)*1/(42e-03);
phaseDiff3 = -medfilt1(phaseDiff3,800)*1/(42e-03);
phaseDiff4 = -medfilt1(phaseDiff4,800)*1/(42e-03);

%%

measr1 = 0.016*mag11.^(-0.5)-0.1; measr2 = 0.006*mag22.^(-0.5)+1.19; 
measr3 = 0.007*mag33.^(-0.5)+1.2; measr4 = 0.01*mag44.^(-0.5)+0.5;

% measr1 = medfilt1(measr1,10);
% measr2 = medfilt1(measr2,10);
% measr3 = medfilt1(measr3,10);
% measr4 = medfilt1(measr4,10);



measrdot1 = phaseDiff1; measrdot2 = phaseDiff2; measrdot3 = phaseDiff3; measrdot4 = phaseDiff4; 

dmeasr1 = zeros(1, length(time));dmeasr2 = zeros(1, length(time));
dmeasr3 = zeros(1, length(time));dmeasr4 = zeros(1, length(time));

dmeasrdot1 = zeros(1, length(time));  dmeasrdot2 = zeros(1, length(time)); 
dmeasrdot3 = zeros(1, length(time));  dmeasrdot4 = zeros(1, length(time)); 

for i = 3:1:length(time)
    difT = time(i) - time(1);
    index = round(difT/(2.0632*0.0001));
    index2 = min(round(difT/(1.0141*0.0001)), length(measrdot3));
    dmeasr1(i) = measr1(index);
    dmeasr2(i) = measr2(index);
    dmeasr3(i) = measr3(index2);
    dmeasr4(i) = measr4(index);
    
    dmeasrdot1(i) = measrdot1(index);
    dmeasrdot2(i) = measrdot2(index);
    dmeasrdot3(i) = measrdot3(index2);
    dmeasrdot4(i) = measrdot4(index);
end

dmeasr1(1) = measr1(1);dmeasr2(1) = measr2(1);dmeasr3(1) = measr3(1);dmeasr4(1) = measr4(1);
dmeasrdot1(1) = measrdot1(1);dmeasrdot2(1) = measrdot2(1);dmeasrdot3(1) = measrdot3(1);dmeasrdot4(1) = measrdot4(1);
dmeasr1(2) = measr1(2);dmeasr2(2) = measr2(2);dmeasr3(2) = measr3(2);dmeasr4(2) = measr4(2);
dmeasrdot1(2) = measrdot1(2);dmeasrdot2(2) = measrdot2(2);dmeasrdot3(2) = measrdot3(2);dmeasrdot4(2) = measrdot4(2);


figure
subplot(411), plot(time1,  measr1, 'LineWidth', 1), hold on, plot(time, r_sim1, 'LineWidth', 1), title('Radial Distance $r1$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([47.9531, 51.4569]);
subplot(412), plot(time1,  measr2, 'LineWidth', 1), hold on, plot(time, r_sim2, 'LineWidth', 1),title('Radial Distance $r2$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([47.9531, 51.4569]);
subplot(413), plot(time3,  measr3, 'LineWidth', 1), hold on, plot(time, r_sim3, 'LineWidth', 1),title('Radial Distance $r3$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([47.9531, 51.4569]);
subplot(414), plot(time1,  measr4, 'LineWidth', 1), hold on, plot(time, r_sim4, 'LineWidth', 1),title('Radial Distance $r4$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([47.9531, 51.4569]);



figure
subplot(411), plot(time1(2:end), phaseDiff1); hold on; plot(time, rdot_sim1, 'LineWidth', 1), title('Radial Velocity $\dot r_1$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([47.9531, 51.4569]);
subplot(412), plot(time1(2:end), phaseDiff2); hold on; plot(time, rdot_sim2, 'LineWidth', 1), title('Radial Velocity $\dot r_2$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([47.9531, 51.4569]);
subplot(413), plot(time3(2:end), phaseDiff3); hold on; plot(time, rdot_sim3, 'LineWidth', 1), title('Radial Velocity $\dot r_3$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([47.9531, 51.4569]);
subplot(414), plot(time1(2:end), phaseDiff4); hold on; plot(time, rdot_sim4, 'LineWidth', 1), title('Radial Velocity $\dot r_4$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([47.9531, 51.4569]);

%
%%
% Measurement Vector [acc_x, acc_y, acc_z, orientation, angular velocity]

% State vector [xdotdot, ydotdot, orientation, angular velocity]
%x = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]'*ones(1,length(accT));
%o1 a1 o2 a2 o3 a3
x = zeros(15, length(time)); x(1,1) = 1.03;
% Intial P matrix
P = eye(length(x(:,1)));
% covariance for w
Q = diag([0.1, 0.1, 10, 0.1, 0.1, 10, 0.1, 0.1, 10, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,]);


% covariance for v
RR = [300, 300, 300, 300, 300, 300, 300, 300, 0.03*var(magx), 0.03*var(gyrox), 0.03*var(magy), 0.03*var(gyroy), 0.03*var(magz), 0.03*var(gyroz), 0.03*var(accx), 0.03*var(accy), 0.03*var(accz)]; 

i = 1;
j = 1;
k = 1;
for m = 1:1:length(time)
    %[index, i, j, k] = getIndex(i, j, k, mag_time, gyro_time, acc_time, m, time);
 
    [y, i, j, k, index] = getyNPVA(dmeasr1(m), dmeasrdot1(m), dmeasr2(m), dmeasrdot2(m), dmeasr3(m), dmeasrdot3(m), dmeasr4(m), dmeasrdot4(m), magx, gyrox, magy, gyroy, magz, gyroz, accx, accy, accz, accT, gyroT, magT, time, i, j, k, m);
    
    if m == 1
        x(:,m) = x(:,m);
    else
        F = getF(time, m);
        %x(:,m-1)
        x(:,m) = F*x(:,m-1);
        %xx = x(:,m)
        
        P = F*P*F' + Q;
        H = getJacoN(y,x(:,m),index);

        R = getR(RR, index);
        % back tracking
        G = P*H'/(H*P*H' + R);%lambda
        % e
        x(:,m) = x(:,m) + G*(y - getH(x(:,m), index));
        P = (eye(length(x(:,1))) - G*H)*P;
    end
    
    
    
    
end
%%
% figure
% subplot(311), plot(time, x(1,:), 'LineWidth', 3); hold on; plot(accT(2:end), posX, 'LineWidth', 2); hold on; plot(time, PP1, 'LineWidth', 2); legend('estimated','integration','ground truth'); title('Position Along x Axis [m]'); xlabel('t [s]'); ylabel('position [m]'); grid on; grid minor;%hold on; plot(time, XX(1,:))
% subplot(312), plot(time, x(7,:), 'LineWidth', 2); hold on; plot(accT(2:end), posY, 'LineWidth', 2); legend('estimated','integration'); title('Position Along y Axis [m]'); xlabel('t [s]'); ylabel('position [m]'); grid on; grid minor;%hold on; plot(time, XX(4,:))
% subplot(313), plot(time, x(4,:), 'LineWidth', 2); hold on; plot(accT(2:end), posZ, 'LineWidth', 2); legend('estimated','integration'); title('Position Along z Axis [m]'); xlabel('t [s]'); ylabel('position [m]'); grid on; grid minor;%hold on; plot(time, XX(7,:))
% csvwrite('data/x.csv',x);   csvwrite('data/XX.csv',XX);   csvwrite('data/phi_gt.csv', phi_gt); 
% csvwrite('data/T3.csv',x);   csvwrite('data/td3.csv',XX); 

phase11 = phase1(21511:38779); 
time1 = t1(21511:38779);

phaseDiff = diff(phase11);
phaseDiff = medfilt1(phaseDiff,1000);

figure
subplot(311), plot(time, x(1,:), 'LineWidth', 2); hold on; plot(time, PP1, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','SouthEast'); title('Position Along x Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([47.9408, 51.456])
subplot(312), plot(time, x(2,:), 'LineWidth', 2); hold on; plot(time, VV1, 'LineWidth', 2); legend('Estimated Velocity','Ground Truth'); title('Velocity Along x Axis [m/s]'); xlabel('t [s]'); ylabel('Velocity [m/s]'); grid on; grid minor;xlim([47.9408, 51.456])%xlim([3.5, 7.5])
subplot(313), plot(time, x(3,:), 'LineWidth', 2); hold on; plot(time, AA1, 'LineWidth', 2); legend('Estimated Acceleration','Ground Truth'); title('Acceleration Along x Axis [m/s^2]'); xlabel('t [s]'); ylabel('Acceleration [m/s^2]'); grid on; grid minor;xlim([47.9408, 51.456])

rmsErrorX = rms(x(1,15:end)-PP1(15:end))
%%
figure
subplot(9,1,1), plot(time, x(1,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along x axis'); grid on;
subplot(9,1,2), plot(time, x(2,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along y axis'); grid on;
subplot(9,1,3), plot(time, x(3,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(9,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along z axis'); grid on;
subplot(9,1,4), plot(time, x(4,:), 'b', 'LineWidth', 2); %hold on; plot(time,, phi_gt(1,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along x axis'); grid on;
subplot(9,1,5), plot(time, x(5,:), 'b', 'LineWidth', 2); %hold on; plot(time, phi_gt(2,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along x axis'); grid on;
subplot(9,1,6), plot(time, x(6,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along y axis'); grid on;
subplot(9,1,7), plot(time, x(7,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(4,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along y axis'); grid on;
subplot(9,1,8), plot(time, x(8,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(5,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along z axis'); grid on;
subplot(9,1,9), plot(time, x(9,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along z axis'); grid on;
