clc; clear; close all;

%% =========================
immune_capacity = 1; %% 2 是正常状态
self_constant = 1;
% 参数定义
% =========================
par.p1  = 1e2;          % 这里改成 0
par.p2  = 0.5;
par.p3  = 0.05;
par.p4  = 10;
par.p5  = 1e20;% 
par.p6  = 1e-9*self_constant;% 5e-9
par.p7  = 1e8; %  1e10
par.p8  = 1e-12*self_constant;
par.p9  = 0.1;
par.p10 = 0.9;
par.p11 = 1e9;
par.p12 = 9*1.22e10;
par.p13 = 0.1;
par.p14 = 1e-9*self_constant;
par.p15 = 1.0;
par.p16 = 0.01;
par.p17 = 1e-10*self_constant;
par.p18 = 1;
par.p19 = 1e8;
par.p20 = 0.01;
par.p21 = 1e1*immune_capacity;
par.p22 = 1e1*immune_capacity;
par.p23 = 0.01;
par.p24 = 0.01;
par.p25 = 1e7;
par.p26 = 1e8;
par.p27 = 0.01;
par.p28 = 0.5;
par.p29 = 1e1*immune_capacity;
par.p30 = 1;% 8e-8
par.p31 = 1;
par.p32 = 1;
par.p33 = 0.1;
par.p34 = 1e4;%% 1e4
par.p35 = 1e9;
par.p36 = 9*1.22e10;
par.p37 = 1e-11*self_constant;
par.p38 = 1.1;
par.p39 = 1e-10*self_constant;
par.p40 = 2;
par.p41 = 0.01;
par.p42 = 1e12;
par.p43 = 0.01;
par.p44 = 0.01;
par.p45 = 0.9;
par.p46 = 0.01;
par.p47 = 1;
par.p48 = 0.01;
par.p49 = 1;
par.p50 = 1;
par.p51 = 1e5*immune_capacity;
par.p52 = 50*immune_capacity;
par.p53 = 0.01;
par.p54 = 0;
par.p55 = 0;
par.p56 = 10; %% 转录系数
par.p57 = 1e2; %% 翻译系数

par.p58 = 1e2;
par.p59 = 1e2;
par.p60 = 0;

par.p61 = 1e8;
par.p62 = 1;
par.p63 = 1e4;


% 肿瘤清零阈值
par.Tu_cutoff = 1e-4;

%% =========================
% 初值
% 变量顺序：
% 1 Tu, 2 N, 3 C, 4 F1, 5 F2, 6 F3, 7 W, 8 Wp, 9 X, 10 A,
% 11 Tc, 12 Tc_non, 13 Tcp, 14 Z, 15 Zp, 16 R1 MHC-2, 17 R2 MHC-1,
% 18 F1p, 19 F2p, 20 F3p, 21 Y, 22 Yp, 23 B, 24 Bp, 25 P,
% 26 C1, 27 M1, 28 M2, 29 M3, 30 Xp 31 RNA  32 N_inside  33 complex_1
% 34 complex_2  35 tumor membrane protein 36 tumor MHC-1 complex
%% =========================
y0 = zeros(36,1);

y0(1)  = 100;
y0(2)  = 0;
y0(3)  = 0;
y0(4)  = 0;
y0(5)  = 1e10;
y0(6)  = 0;
y0(7)  = 1e6;
y0(8)  = 9999*1e6;
y0(9)  = 0;
y0(10) = 1e10;
y0(11) = 1e6;
y0(12) = 9999*1e6;
y0(13) = 0;
y0(14) = 0;
y0(15) = 0;
y0(16) = 0;
y0(17) = 0;
y0(18) = 0;
y0(19) = 1e10;
y0(20) = 0;
y0(21) = 0;
y0(22) = 0;
y0(23) = 1e10;
y0(24) = 9999*1e10;
y0(25) = 1e6;
y0(26) = 0;
y0(27) = 0;
y0(28) = 1e10;
y0(29) = 0;
y0(30) = 0;

y0(31) = 0;
y0(32) = 0;
y0(33) = 0;
y0(34) = 0;


y0(35) = 0; % tumor membrane protein
y0(36) = 0; % tumor MHC-1 complex


%% 时间区间
tspan = [0 5000];

%% 求解器设置
opts = odeset('NonNegative',1:36);

%% 求解
[t1,y1] = ode15s(@(t,y) tumor_immune_interaction_herv(t,y,par), tspan, y0,opts);




figure
plot(t1,(y1(:,1)),'r','LineWidth',2); 
figure
plot(t1,(y1(:,10)),'r','LineWidth',2); 
figure
plot(t1,(y1(:,32)),'r','LineWidth',2); 
% hold on
% 
% plot(t1,log10(y1(:,32)./y1(:,1)),'r','LineWidth',2); 






