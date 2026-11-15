function dydt = tumor_immune_interaction_herv(t,y,par)

% =========================
% 变量拆包
% =========================
Tu     = y(1);
N      = y(2);
C      = y(3);
F1     = y(4);
F2     = y(5);
F3     = y(6);
W      = y(7);
Wp     = y(8);
X      = y(9);
A      = y(10);
Tc     = y(11);
Tc_non = y(12);
Tcp    = y(13);
Z      = y(14);
Zp     = y(15);
R1     = y(16);
R2     = y(17);
F1p    = y(18);
F2p    = y(19);
F3p    = y(20);
Y      = y(21);
Yp     = y(22);
B      = y(23);
Bp     = y(24);
P      = y(25);
C1     = y(26);
M1     = y(27);
M2     = y(28);
M3     = y(29);
Xp     = y(30);
R      = y(31);
N_in   = y(32);
complex_1  = y(33);
complex_2  = y(34);
N_sur =  y(35);
N_sur_2 = y(36);

new_capacity = 1;


% =========================
% 参数拆包
% =========================
p1=par.p1;   p2=par.p2;   p3=par.p3;   p4=par.p4;   p5=par.p5;
p6=par.p6;   p7=par.p7;   p8=par.p8;   p9=par.p9;   p10=par.p10;
p11=par.p11; p12=par.p12; p13=par.p13; p14=par.p14; p15=par.p15;
p16=par.p16; p17=par.p17; p18=par.p18; p19=par.p19; p20=par.p20;
p21=par.p21*new_capacity; p22=par.p22*new_capacity; p23=par.p23; p24=par.p24; p25=par.p25;
p26=par.p26; p27=par.p27; p28=par.p28; p29=par.p29*new_capacity; p30=par.p30;
p31=par.p31; p32=par.p32; p33=par.p33; p34=par.p34; p35=par.p35;
p36=par.p36; p37=par.p37; p38=par.p38; p39=par.p39; p40=par.p40;
p41=par.p41; p42=par.p42; p43=par.p43; p44=par.p44; p45=par.p45;
p46=par.p46; p47=par.p47; p48=par.p48; p49=par.p49; p50=par.p50;
p51=par.p51*new_capacity; p52=par.p52*new_capacity; p53=par.p53; p54=par.p54; p55=par.p55;
p56=par.p56; p57=par.p57; p58=par.p58; p59=par.p59; p60=par.p60; p61=par.p61; p62 = par.p62;p63 = par.p63;
Tu_cutoff = par.Tu_cutoff;

eps_val = 1e-30;
dydt = zeros(36,1);

% =========================
% 肿瘤清零规则：
% 若 Tu < 1e-10，则认为Tu=0，且以后保持为0
% =========================
tumor_dead = (Tu < Tu_cutoff);

if tumor_dead
    Tu_eff = 0;
else
    Tu_eff = Tu;
end

%% 1. Tumor
if tumor_dead
    dydt(1) = 0;
else
    dydt(1) = p1 + p2*Tu_eff - p3*Tu_eff ...
              - p4*((complex_1/Tu_eff)^4/((complex_1/Tu_eff)^4+p5))*Tu_eff ...
              -p4*((complex_2/Tu_eff)^4/((complex_2/Tu_eff)^4+p7))*Tu_eff;
end

%% 2. Neoantigen
if tumor_dead
    dydt(2) = - p8*N*A + p9*C - p8*N*B + p9*C1;
else

    dydt(2) =  p3*N_in +p4*(((complex_1/Tu_eff)^4/((complex_1/Tu_eff)^4+p5))+((complex_2/Tu_eff)^4/((complex_2/Tu_eff)^4+p7)))*N_in ...
        - p8*N*A + p9*C - p8*N*B + p9*C1;
end
%% 3. Immune complex C
dydt(3) = p8*N*A - p9*C - p10*C;

%% 4. F1
dydt(4) = p11*C/(C+p12) - p13*F1 - p14*F1*F2;

%% 5. F2
dydt(5) = -p14*F1*F2 + p15*F3 - p14*R1*F2 + p61 - p53*F2;

%% 6. F3
dydt(6) = p14*F1*F2 - p15*F3 + p14*R1*F2 ...
          - p16*F3 - p17*1e4*F3*W + p18*X;

%% 7. W
dydt(7) = p19*W/(W+Wp+eps_val) - p20*W ...
          + (-p17*F3*1e4*W + p18*X)/1000 ...
          + p21*X ...
          - (p17*M3*1e4*W + p18*Xp)/10000 ...
          + p22*Xp;

%% 8. W'
dydt(8) = p19*Wp/(W+Wp+eps_val) - p20*Wp;

%% 9. X
dydt(9) = p17*F3*1e4*W - p18*X - p23*X ...
          - p30*X*(Y+Yp) + p31*(Z+Zp);

%% 10. A
dydt(10) = -p24*A - p8*N*A + p9*C + p25*P+p60;

%% 11. Tc
if tumor_dead
    dydt(11) = p26*Tc/(Tc+Tc_non+eps_val) - p27*Tc ...
           + p29*Z ...
           + (-p39*F3p*1e4*Tc + p40*Y)/10000;
else
    dydt(11) =p26*Tc/(Tc+Tc_non+eps_val) - p27*Tc ...
           - p28*(p4*((complex_2/Tu_eff)^4/((complex_2/Tu_eff)^4+p7))*Tu_eff) + p29*Z ...
           + (-p39*F3p*1e4*Tc + p40*Y)/10000;
end

%% 12. Tc_non
dydt(12) = p26*Tc_non/(Tc+Tc_non+eps_val) - p27*Tc_non;

%% 13. Tc'
if tumor_dead
    dydt(13) = -p27*Tcp ...
           + (-p39*F3p*1e4*Tcp + p40*Yp)/10000;
else
   dydt(13) = -p27*Tcp + p28*((p4*((complex_2/Tu_eff)^4/((complex_2/Tu_eff)^4+p7))*Tu_eff)) ...
           + (-p39*F3p*1e4*Tcp + p40*Yp)/10000;
end

%% 14. Z
dydt(14) = p30*X*Y - p31*Z - p32*Z;

%% 15. Z'
dydt(15) = p30*X*Yp - p31*Zp - p32*Zp;

%% 16. R1
dydt(16) = p54-p33*R1 - p14*R1*F2 + p15*F3; %% ...
           %% - p14*R1*M2 + p15*M3;

%% 17. R2
dydt(17) = p55-p33*R2 - p37*R2*F2p + p38*F3p;

%% 18. F1'
dydt(18) = p35*C/(C+p36) - p13*F1p - p37*F1p*F2p;

%% 19. F2'
dydt(19) = -p37*F1p*F2p + p38*F3p - p37*R2*F2p + p61 - p53*F2p;

%% 20. F3'
dydt(20) = p37*F1p*F2p - p38*F3p + p37*R2*F2p ...
           - p16*F3p - p39*1e4*F3p*(Tc+Tcp) + p40*(Y+Yp);

%% 21. Y
dydt(21) = p39*1e4*F3p*Tc - p40*Y - p41*Y - p30*X*Y + p31*Z;

%% 22. Y'
dydt(22) = p39*1e4*F3p*Tcp - p40*Yp - p41*Yp - p30*X*Yp + p31*Zp;

%% 23. B
dydt(23) =p42*B/(B+Bp+eps_val) - p43*B ...
           - p8*N*B + p9*C1 ...
           + p51*(1 - p50*C1/(C1+B+eps_val))*Xp;

%% 24. B'
dydt(24) = p42*Bp/(B+Bp+eps_val) - p43*Bp;

%% 25. P
dydt(25) = p52*(p50*C1/(C1+B+eps_val))*Xp - p44*P + p34;

%% 26. C1
dydt(26) = p8*N*B - p9*C1 - p45*C1;

%% 27. M1
dydt(27) = p46*C1 - p13*M1 - p14*M1*M2;

%% 28. M2
dydt(28) = -p14*M1*M2 ... %% -p14*M1*M2 - p14*R1*M2 + p15*M3 ...
           + p47*p51*(1 - p50*C1/(C1+B+eps_val))*Xp ...
           + p61 - p53*M2;

%% 29. M3
dydt(29) = p14*M1*M2-p48*M3... %% + p14*R1*M2 - p15*M3 - p48*M3 ...
           - p17*1e4*M3*W + p18*Xp;

%% 30. X'
dydt(30) = p17*1e4*M3*W - p18*Xp - p49*Xp;

%% 31. RNA
if tumor_dead 
dydt(31) = 0;
else
dydt(31) = p56*Tu- p3*R - p4*(((complex_1/Tu_eff)^4/((complex_1/Tu_eff)^4+p5))+((complex_2/Tu_eff)^4/((complex_2/Tu_eff)^4+p7)))*R + p62*R*(N_in/Tu_eff)/(N_in/Tu_eff+p63);
end

%% 32 N_in

if tumor_dead 
dydt(32) = 0;
else
dydt(32) = p57*R- p3*N_in - p4*(((complex_1/Tu_eff)^4/((complex_1/Tu_eff)^4+p5))+((complex_2/Tu_eff)^4/((complex_2/Tu_eff)^4+p7)))*N_in;
end


%% complex_1
if tumor_dead 
dydt(33) = 0;
else
dydt(33) =  p8*A*N_sur - p9*complex_1 - p4*((complex_1/Tu_eff)^4/((complex_1/Tu_eff)^4+p5))*complex_1-0.01*N_sur-p3*complex_1;
end

%% complex_2
if tumor_dead 
dydt(34) = 0;
else
dydt(34) = p39*Tc*N_sur_2 - p40*complex_2 - p4*((complex_2/Tu_eff)^4/((complex_2/Tu_eff)^4+p7))*complex_2-0.01*N_sur_2-p3*complex_2;

end

%% tumor free neoantigen on membrane
if tumor_dead 
dydt(35) = 0;
else
dydt(35) =  N_in/p58 - p8*A*N_sur + p9*complex_1 - p4*((complex_1/Tu)^4/((complex_1/Tu)^4+p5))*N_sur-0.01*N_sur-p3*N_sur;
end


%% tumor free MHC-1-neoantigen-peptiede complex on membrane
if tumor_dead 
dydt(36) = 0;
else
dydt(36) = N_in/p59 - p39*Tc*N_sur_2 + p40*complex_2 - p4*((complex_2/Tu)^4/((complex_2/Tu)^4+p7))*N_sur_2 -0.01*N_sur_2-p3*N_sur_2;
end

end