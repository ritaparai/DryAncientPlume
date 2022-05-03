%% This is a fresh page from which to do the plume paper calculations
%% RParai Oct 20 2020

clear; close all;

%% Initialize rand and load components
s = RandStream('mt19937ar','Seed', sum(1e4*clock));
RandStream.setGlobalStream(s);
load components2
load Xemasses
load o
Ms = 0.4e27;
run.late_t = 1e6*[0:0.1:200 201:1:1302 1307:5:4567]';
run.BestQ = findBestQ_plumes(o,run.late_t,Ms);
whatCCx = 1;
eval(strcat('dCCX = o.dCC', num2str(whatCCx),';'));
tBestQ = run.BestQ(1,1,sub2ind([3,3],1,whatCCx));

mcnum = 5e2;

mm = Xe_masses;
run.whatM = 1;          %% 50%, 75%, 90%

adjustI = 1;
lv = 1;

%% decay constants
lambda129 = 4.4150e-8;
lambda232 = 4.9334e-11;
lambda235 = 9.8486e-10;
lambda238 = 1.5514e-10;
lambda244 = 8.6643e-9;

%% fission yields
yPu136 = 7e-5;
yPu132 = yPu136./components(5,4,1);     %% use endmember fission compositions to compute partial yields from 136 yields
yPu131 = yPu132.*components(3,4,1);
yPu134 = yPu132.*components(4,4,1);
yU136 = 3.43e-8;
yU132 = yU136./components(5,3,1);
yU131 = yU132.*components(3,3,1);
yU134 = yU132.*components(4,3,1);

Ne21He4prod = 4.5e-8;                   %% Yatsevich and Honda 1997

%% Parent concentrations in the BULK UPPER MANTLE-CONTINENT (= BSE)
BSEU = 21e-09;                      %% 21 ppb U [g/g]
Ui = (BSEU*0.992745/238)/exp(-lambda238*4.567e9)*6.02e23;   %% initial U conc [atoms/g]
PuU = 0.0068;                       %% 244Pu/238U solar from Jacobsen and Harper 1996
Pui = (PuU*Ui);                     %% initial Pu conc [atoms/g]
I129127 = 1.1e-4;                   %% initial ratio
Xe129130C = 6.436;                  %% Q Xe (Pepin and Porcelli, 2006)
T = 4.567e9;
t = 1e6*(0:0.01:4567)';
U235U238 = (0.0072*exp(lambda235*(4.567e9-t)))./(0.992745*exp(lambda238*(4.567e9-t)));  %% 235U/238U ratio as a function of time
Th232U238 = ((80e-09/232)*exp(lambda232*(4.567e9-t)))./((BSEU*0.992745/238)*exp(lambda238*(4.567e9-t))); %% 232Th/238U ratio as a function of time
Pu244U238 = (Pui*exp(-lambda244*t))./(Ui*exp(-lambda238*t));

%% Set up the parameter space
%% Accretion-related parameters
ACC.t = [3 30];                     %% e-folding timescale for accretion
ACC.fr = 0.99;                      %% mass fraction accreted before closure

%% Pre-closure parameters
preC.retI = [0 1];
preC.delI = [0 1];
preC.retXe = [0 15];
preC.delXe = [0 1];
preC.IC = 7./adjustI;               %% 7 ppb 127I

%% Single instantaneous loss event (Giant Impact)
SSL.retI = [0 1];
SSL.retXe = [0 1];

%% Pre-closure parameters
postC.retI = [0 1];                
postC.delI = [0 1];
postC.retXe = [0 10];               %% want to weight this towards 1, so use 1-log10(-N)
postC.delXe = [0 1];
postC.IC = 450./adjustI;            %% just prescribe 450 ppb 127I

run.accfr = ACC.fr;
run.IC = preC.IC;
run.IC_pc = postC.IC;

%% ==========================================================
%% EEXe-blank
%% ==========================================================
%% Accretion-related parameters
run.acct = 11;

%% Pre-closure parameters
run.retI = 1;
run.delI = 1;
run.retXe = 0;
run.delXe = 1;
% run.IC = (max(preC.IC) - min(preC.IC))*rand + min(preC.IC);

%% Single instantaneous loss event (Giant Impact)
run.retI_SSL = 1;
run.retXe_SSL = 0;

%% Post-closure parameters
run.retI_pc = 1;
run.delI_pc = 1;
run.retXe_pc = 1;
run.delXe_pc = 1;
%% ==========================================================

LVs = [0.001 0.005 0.01];   %% late veneer fractions


%% Downwelling Xe concentration time series
%% ========================================================================
form = 'sigmoidal';           %% Choose the form of the downwelling Xe function
highearly = 0;                  %% Choose 0 or 1 depending on whether you want to account for high early Xe atmospheric concentrations and Henry's Law

chavrit = 9e7;                  %% max atoms per gram 130Xe in subducting slab; Chavrit et al., 2016
ours = 5e8/4;                  %% spend more time searching in around the right ballpark
maximum = ours;

switch form
    case 'sigmoidal'
%         r1 = [4.3e5/maximum 1]; %% K always 0 to 1, scaling factor for the max
        r1 = [0 1];
        r2 = [8e7 1e10];         %% M inflection point
        r3 = [1e-8 1e-10];       %% Z rise rate
    case 'exponential'
        res3 = 1;
%         r1 = [4.3e5/maximum 1];
        r1 = [0 1];        
        r2 = [1e-11 5e-8];
        r3 = NaN;
    case 'parabolic'
        r1 = [1e9 5e9];
        r2 = [1e3 1e6];
        r3 = [1e-16 1e-12];
end

%% Boundaries of the solution space
%% ========================================================================
CS = [4.3e5 9.2e5];             %% concentration success boundaries
% CS = Xe130M_T.*[0.5 2];       %% concentration success boundaries
RS = zeros(7,2);                %% based solely on DICE data (M12)
RS(1,:) = [0.4725 0.474];        
RS(2,:) = [0.0703 0.0706];
RS(3,:) = [1.024 1.040];
RS(4,:) = [0.1488 0.1490];
RS(5,:) = [0.7813 0.7833];
RS(6,:) = [0.4004 0.4020];
RS(7,:) = [0.3450 0.3470];


%% Modern and initial atmospheric Xe compositions
%% ========================================================================
Xe128130air = components(1,1,1)./components(2,1,1);

%% Atmospheric Xe isotopic time series (x-axis)
ff = (-0.0777/2.5e9).*run.late_t + 1.0777; ff(ff<1)=1;                 %% fractionations in a vector based on run.late_t; matched 128/130 to U-Xe
ff(:)=1;    %% TURN OFF Xe ATMOSPHERIC EVOLUTION
Xe128130St = Xe128130air.*ff;
a = (1-( (mm(3)./mm([2 4:7])).^(1/2) ))./(1-sqrt(mm(3)/mm(1)));
b = ((mm(3)./mm([2 4:7])).^(1/2)) ./ (sqrt(mm(3)/mm(1)).^a);
b = b./b;
Xe129130St = 6.496.*b(1).* ((ff).^a(1));
Xe131130St = components(3,1,1)./components(2,1,1).*b(2).* ((ff).^a(2));
Xe132130St = 1./components(2,1,1).*b(3).* ((ff).^a(3));
Xe134130St = components(4,1,1)./components(2,1,1).*b(4).* ((ff).^a(4));
Xe136130St = components(5,1,1)./components(2,1,1).*b(5).* ((ff).^a(5));



%% Start with the parameter ranges
Nrange = linspace(1,4,30);
initrange = linspace(1e8,3e9,30);
mapvx = zeros(numel(Nrange),numel(initrange));
Xetotals = zeros(numel(Nrange),numel(initrange));
Netotals = zeros(numel(Nrange),numel(initrange));
Hetotals = zeros(numel(Nrange),numel(initrange));
elemtotals = zeros(numel(Nrange),numel(initrange));
f130ares =  zeros(numel(Nrange),numel(initrange));
% grid1 = 1; grid2 = 14;
Nres = 3.3;
Xe130i = 2e9;
for grid1 = 1:numel(Nrange)
    Nres = Nrange(grid1);
    for grid2 = 1:numel(initrange)
        Xe130i = initrange(grid2);

%% choose a 130Xe/22Ne value for the Iceland mantle
% test13022M = linspace(0.00276,0.012,100);
test_res = 0.02;
%% fraction 130Xe atm
Icef130a = 0.93;
adjustment = 1;
% while test_res > 0.01

% %% Accreted volatile parameters
% ACC.Ne2022sol_l = [13.36-0.09 13.36+0.09];  %% limits; sticking with Heber
% ACC.Ne2022cc_l = [9.03-2.46 9.03+2.46];
% 
% ACC.Xe130Ne22sol_l = [1.03e-6 1.97e-5];     %% ott 2014; Heber and Meshik recalculated WM19
% ACC.Xe130Ne22cc_l = [0.015 0.055];

Ne2022atm = 9.8;
Xe130Ne22atm = 0;
Xe130Ne22atm_nom = 0.0021;

%% Order is SOLAR CC ATM
Ne2022 = [13.36 8.9 9.8]; 
Xe130Ne22 = [1.03e-6 0.033 0];
A = [Ne2022; Xe130Ne22; [1 1 1]];


%% Iceland WILLIAMS AND MUKHOPADHYAY 2019
%% 20Ne/22Ne 
Ice2022 = 13.17;
%% 130Xe/22Ne
Ice13022 = 0.003;
% %% fraction 130Xe atm
% Icef130a = 0.97;
b = [Ice2022; Ice13022.*(1-Icef130a); 1];
%% 3He/22Ne
He3Ne22sol = 1.46;  %% Tucker and Mukhopadhyay 2014
He3Ne22cc = 0.89;
He3Ne22atm = 4.3e-6;

xPLUME = A\b;

check1 = sum(xPLUME .* Ne2022') == Ice2022;
Ne2022acc = (sum(xPLUME(1:2) .* Ne2022(1:2)'))./sum(xPLUME(1:2));
Ne22M = 1;
Ne22sol = Ne22M.*xPLUME(1);
Ne22cc = Ne22M.*xPLUME(2);
Ne22atm = Ne22M.*xPLUME(3);

Xe130M = Ne22M.*Ice13022; %% Total 130Xe in the mantle today assuming 22Ne = 1
Xe130sol = Ne22sol.*Xe130Ne22(1);
Xe130cc = Ne22cc.*Xe130Ne22(2);
Xe130atm_part = Ne22atm.*Xe130Ne22atm_nom;
XewithNe = Xe130sol + Xe130cc + Xe130atm_part;
Xe130excess = Xe130M - XewithNe;
Xe130atm = Xe130atm_part + Xe130excess;
Xe130Ne22_R = Xe130atm ./ Ne22atm;
check2 = abs(Xe130atm./Xe130M - Icef130a) < 1e-15;

Ice3130 = 900;
He3Xe130 = Ice3130./(1-Icef130a);   %% the 3He/130Xe ratio at accretion
% He3M = Xe130M.*He3Xe130;
He3sol = Ne22sol.*He3Ne22sol;
He3cc = Ne22cc.*He3Ne22cc;
He3atm = Ne22atm.*He3Ne22atm;
% He3Ne22acc = (He3sol + He3cc + He3atm) ./ Ne22M < Ice322

% He3Ne22_afterMO = Ice322 * (1 ./ (1-xPLUME(3)));
% FrFr = 1.8;
% He3Ne22_afterMO = FrFr*(He3sol + He3cc + He3atm) ./ Ne22M;

%% So the mantle 3/22 had to be raised by magma ocean degassing after

% He3Xe130acc = FrFr*(He3sol + He3cc) ./ (Xe130sol + Xe130cc);
He3Xe130acc = He3Xe130;
Xe130Ne22acc = (Xe130sol + Xe130cc) ./ (Ne22sol + Ne22cc);
% Xe130Ne22acc = (1./He3Ne22_afterMO) .* He3Xe130acc;

Ice322 = Ice3130 .* Ice13022;


%% ICELAND
Neslope = (13.8-9.8)/(0.0375-0.029);
He3Xe130i = He3Xe130acc;
Xe129Xe130 = 10.5;
Xe129136star = 2.8;
Xe136PuXe130 = (Xe129Xe130 - Xe129130C)./ Xe129136star;
NeXeslabs = 1./Xe130Ne22_R;
targetHe = [40000; 44000];
targetNe = [Ice2022; 0.029+(Ice2022-9.8)/Neslope];
targetHeXe = Ice3130 * [0.9 1.1];
plotx = [0.475 0.001; 1.032 8e-3; 0.1489 1e-4; 0.7823 1e-3; 0.4012 8e-4; 0.3460 1e-3];



%% ==EARLY===========================================================================================================
%% TIME EVOLUTION MODEL====================================================
%% Continuous accretion

a = 1/(run.acct*1e6);               %% accretion time constant
Mt = (1-exp(-a.*t));                %% M(t) mass of growing planet
dx = find(Mt>(1-LVs(3)),1);         %% dx is start index for tclosure (0.99 accreted)
dm = find(Mt>(1-LVs(2)),1);         %% dm is end index for tclosure (0.995 accreted)
dy = find(Mt>(1-LVs(1)),1);         %% dy is end index for tclosure (0.999 accreted)
dend = find(Mt>(1-1e-4),1);         %% dend is index of time at end of accretion (0.9999 accreted)
if t(dend) < 200e6
    dend = floor(dend/10)*10 + 1; 	%% round down to nearest 0.1 Myr to hand off to long-term loop
else
    dend = floor(dend/100)*100 + 1; %% round down to nearest 1 Myr to hand off to long-term loop
end
dc = 1e6;                           %% resolution on candidate tclosure
tc = ([round(t(dx)/dc) round(t(dm)/dc) round(t(dy)/dc)]).*dc;	%% vector of candidate tclosure 1Myr resolution
tc = tc(lv);
X = ones(size(tc));                 %% tc-sized dummy of ones

%% Initialize last values at zero
Ilast = 0; Xe129last = 0; Pulast = 0; Xe136last = 0; Ulast = 0; XeUlast = 0;
Xe130last = 0; Xe129Clast = 0; Xe130FDlast = 0; Xe129CFDlast = 0;
INOGlast = 0; Xe129NOGlast = 0; Xe136NOGlast = 0; Xe130NOGlast = 0;
Xe131last = 0; Xe132last = 0; Xe134last = 0;
XeU131last = 0; XeU132last = 0; XeU134last = 0;

%% Initialize pre-closure parameters for all tc (make a vector the length of tc)
retfrI = run.retI.*X; retfrXe = run.retXe.*X; delfrI = run.delI.*X; delfrXe = run.delXe.*X; Ii = (run.IC*1e-09/127)*I129127.*X.*6.02e23; %% 129I [atoms/g]
istc = false(1,numel(tc));

%% EARLY DEGASSING TIME LOOP
for k = 2:dend                      %% t vector starts at zero
    tnow = t(k);
    tlast = t(k-1);
    dt = tnow-tlast;
    Mlast = Mt(k-1);
    Mnow = Mt(k);
    dM = Mnow-Mlast;
    %% Loop through and for each tc, replace with SSL or post-closure parameters once tnow == tc or tnow > tc
    e2 = tnow == tc;
        retfrI(e2) = run.retI_SSL; retfrXe(e2) = run.retXe_SSL; Ii(e2) = run.IC_pc*1e-09/127*I129127;
        istc(e2) = 1;
    e3 = tnow > tc;
        retfrI(e3) = run.retI_pc; retfrXe(e3) = run.retXe_pc; delfrI(e3) = run.delI_pc; delfrXe(e3) = run.delXe_pc;
        istc(e3) = 0;
  
    Inow = (Ilast.*exp(-lambda129*dt)*Mlast.*retfrI + Ii.*exp(-lambda129*tnow)*dM.*delfrI)/Mnow;
    INOGnow = (INOGlast.*exp(-lambda129*dt)*Mlast + Ii*exp(-lambda129*tnow)*dM.*delfrI)/Mnow;
    
    Xe129now = ((Xe129last+Ilast.*(1-exp(-lambda129*dt)))*Mlast.*retfrXe)/Mnow;
    Xe129NOGnow = ((Xe129NOGlast+Ilast.*(1-exp(-lambda129*dt)))*Mlast)/Mnow;
    
    Punow = (Pulast.*exp(-lambda244*dt)*Mlast + Pui*exp(-lambda244*tnow)*dM)/Mnow;
    Xe136now = ((Xe136last+Pulast.*(1-exp(-lambda244*dt))*yPu136)*Mlast.*retfrXe)/Mnow;
    Xe136NOGnow = ((Xe136NOGlast+Pulast.*(1-exp(-lambda244*dt))*yPu136)*Mlast)/Mnow;
    
    Xe131now = ((Xe131last+Pulast.*(1-exp(-lambda244*dt))*yPu131)*Mlast.*retfrXe)/Mnow;
    Xe132now = ((Xe132last+Pulast.*(1-exp(-lambda244*dt))*yPu132)*Mlast.*retfrXe)/Mnow;
    Xe134now = ((Xe134last+Pulast.*(1-exp(-lambda244*dt))*yPu134)*Mlast.*retfrXe)/Mnow;
    
    Unow = (Ulast.*exp(-lambda238*dt)*Mlast + Ui*exp(-lambda238*tnow)*dM)/Mnow;
    XeUnow = ((XeUlast+Ulast.*(1-exp(-lambda238*dt))*yU136)*Mlast.*retfrXe)/Mnow;
    XeU131now = ((XeU131last+Ulast.*(1-exp(-lambda238*dt))*yU131)*Mlast.*retfrXe)/Mnow;
    XeU132now = ((XeU132last+Ulast.*(1-exp(-lambda238*dt))*yU132)*Mlast.*retfrXe)/Mnow;
    XeU134now = ((XeU134last+Ulast.*(1-exp(-lambda238*dt))*yU134)*Mlast.*retfrXe)/Mnow;
    
    
    Xe130now =  (Xe130last.*Mlast + Xe130i*dM.*delfrXe).*retfrXe/Mnow;
    Xe129Cnow = (Xe129Clast.*Mlast + Xe129130C*Xe130i*dM.*delfrXe).*retfrXe/Mnow;
    Xe130NOGnow =  (Xe130NOGlast.*Mlast + Xe130i*dM.*delfrXe)/Mnow;
    Xe129CFDnow = (Xe129CFDlast.*Mlast + Xe129130C*Xe130i*dM).*retfrXe/Mnow;
    Xe130FDnow =  (Xe130FDlast.*Mlast + Xe130i*dM).*retfrXe/Mnow;
    
    Ilast = Inow; Xe129last = Xe129now; Pulast = Punow; Xe136last = Xe136now;  Ulast = Unow; XeUlast = XeUnow;
    Xe131last = Xe131now; Xe132last = Xe132now; Xe134last = Xe134now; 
    XeU131last = XeU131now; XeU132last = XeU132now; XeU134last = XeU134now;
    Xe130last = Xe130now; Xe129Clast = Xe129Cnow; Xe130FDlast = Xe130FDnow; Xe129CFDlast = Xe129CFDnow;
    INOGlast = INOGnow; Xe129NOGlast = Xe129NOGnow; Xe136NOGlast = Xe136NOGnow; Xe130NOGlast = Xe130NOGnow;
    
end

% Xe130last

%% INTERMEDIATE AND LONG-TERM DEGASSING TIME LOOPS
tlast = t(dend);
t = run.late_t;
startidx = find(t>tlast,1);
Mproc = findMproc_plumes_spec(run.late_t,startidx,Nres,Ms);

tblock = t;                                                     %% tblock (and therefore Cs130) based on run.late_t
z = numel(t(startidx:end));

%% Puff everything down and out and back (Nres x tc x Mres:dCC)
Mres = repmat(Ms,[1 numel(tc)]);       %% z1 (1) rows by tc columns
Pulast = repmat(Pulast,[1 numel(tc)]); Ulast = repmat(Ulast,[1 numel(tc)]);
Z = ones(size(Pulast));

%% Initialize the ballpark sanity check values for 3He/4He degassing
U235last = Ulast.*U235U238(dend);       %% [moles/g]
Thlast = Ulast.*Th232U238(dend);        %% [moles/g]

He3last = Xe130last.*He3Xe130i;     %% [atoms/g] THIS IS THE PLUME SOURCE 3He concentration after accretion; Harper and Jacobsen (1995) suggest 7.3e10 as a lower limit
He4last = He3last./(120*1.39e-06);      %% [atoms/g]

Ne22last = Xe130last./Xe130Ne22acc;
Ne20last = Ne22last.*Ne2022acc;      %% 20Ne/22Ne at the end of accretion, prior to regassing
Ne21last = Ne22last.*0.0324;         %% solar, largely; Williams and Mukhopadhyay (2019) recalc of Heber (2012)
Ne22Xe130S = NeXeslabs;


U235U238 = (0.0072*exp(lambda235*(4.567e9-t)))./(0.992745*exp(lambda238*(4.567e9-t)));  %% 235U/238U ratio as a function of time
Th232U238 = ((80e-09/232)*exp(lambda232*(4.567e9-t)))./((BSEU*0.992745/238)*exp(lambda238*(4.567e9-t))); %% 232Th/238U ratio as a function of time
Pu244U238 = (Pui*exp(-lambda244*t))./(Ui*exp(-lambda238*t));

%% LONG-TERM DEGASSING TIME LOOPS

%% Regassing budgets start at zero
XeR130last = 0;
XeR128last = 0;
XeR129last = 0;
XeR131last = 0;
XeR132last = 0;
XeR134last = 0;
XeR136last = 0;

NeR20last = 0;
NeR21last = 0;
NeR22last = 0;

Xe130Rlast = 0;
Xe130FDlast = 0;
Xe129Clast = 0;
Xe129CFDlast = 0;
SaveEarly_outputs = num2cell([Ulast U235last Thlast Pulast XeUlast XeU131last XeU132last XeU134last Xe136last Xe131last Xe132last Xe134last XeR130last XeR128last XeR129last XeR131last XeR132last XeR134last XeR136last Xe130Rlast Xe129last Xe130last Xe130FDlast Xe129Clast Xe129CFDlast He3last He4last Ne20last Ne21last Ne22last NeR20last NeR21last NeR22last]);

    
%% Allocate memory for time series or reset to zeros
tr = rand(3,mcnum);     %% call rand once per loop
storeparams = zeros(3,mcnum);
store128130M = zeros(z,mcnum);
store128130R = zeros(z,mcnum);
store128132M = zeros(z,mcnum);
store129130M = zeros(z,mcnum);
store129132M = zeros(z,mcnum);
store130132M = zeros(z,mcnum);
store131130M = zeros(z,mcnum);
store131132M = zeros(z,mcnum);
store132130M = zeros(z,mcnum);
store134130M = zeros(z,mcnum);
store134132M = zeros(z,mcnum);
store136130M = zeros(z,mcnum);
store136132M = zeros(z,mcnum);
store4He3HeM = zeros(z,mcnum);
storePuUM = zeros(mcnum,1);
store130M = zeros(z,mcnum);
store3HeM = zeros(z,mcnum);
storeCs130 = zeros(z,mcnum);
store130Rsteps = zeros(z,mcnum);
store2022M = zeros(z,mcnum);
store2122M = zeros(z,mcnum);
store43M = zeros(z,mcnum);
storeHe3M = zeros(z,mcnum);
storeNe22M = zeros(z,mcnum);
store3130M = zeros(z,mcnum);
store13022M = zeros(z,mcnum);
store322M = zeros(z,mcnum);
store13022R = zeros(z,mcnum);
store_err = zeros(10,mcnum);
Sidx  = zeros(1,mcnum);
store130fracs = zeros(2,mcnum);
store132fracs = zeros(4,mcnum);
store22fracs = zeros(2,mcnum);
store13022x = zeros(2,mcnum);

for mc = 1:mcnum
    switch form
        case 'sigmoidal'
            K = (max(r1) - min(r1))*tr(1,mc) + min(r1);
            M = (max(r2) - min(r2))*tr(2,mc) + min(r2);
            Z = (max(r3) - min(r3))*tr(3,mc) + min(r3);
            storeparams(:,mc) = [K;M;Z];
            Cs130 = (maximum).*(0 + (K-0) ./ (1 + exp(-Z.*(tblock-M))).^1);
        case 'exponential'
            mu = (max(r1) - min(r1))*tr(1,mc) + min(r1);
            beta = exp((log(max(r2)) - log(min(r2)))*tr(2,mc) + log(min(r2)));
%                         beta = (max(r2) - min(r2))*tr(2,mc) + min(r2);
            storeparams(:,mc) = [mu;beta;NaN];
            Cs0 = maximum.*exp(-beta.*T)';         %% calc initial Cs based on present-day; apply scaling factor in next step
            Cs130 = mu.* repmat(Cs0,[size(t,1) 1]).*exp(beta.*tblock);
        case 'parabolic'
            a = (max(r1) - min(r1))*tr(1,mc) + min(r1);
            b = (max(r2) - min(r2))*tr(2,mc) + min(r2);
            c = (max(r3) - min(r3))*tr(3,mc) + min(r3);
            storeparams(:,mc) = [a;b;c];
            Cs130 = c.*(tblock-a).^2+b;

    end
    
    %% NEED TO RE-INITIALIZE WITHIN THE MC LOOP
    %% Add second part to time series
    IT2 = zeros(z,numel(tc));
    Xe129T2 = zeros(z,numel(tc));
    PuT2 = zeros(z,numel(tc));
    Xe136T2 = zeros(z,numel(tc));
    Xe131T2 = zeros(z,numel(tc));
    Xe132T2 = zeros(z,numel(tc));
    Xe134T2 = zeros(z,numel(tc));
    UT2 = zeros(z,numel(tc));
    XeUT2 = zeros(z,numel(tc));
    Xe131UT2 = zeros(z,numel(tc));
    Xe132UT2 = zeros(z,numel(tc));
    Xe134UT2 = zeros(z,numel(tc));
    Xe128RT2 = zeros(z,numel(tc));
    Xe129RT2 = zeros(z,numel(tc));
    Xe130RT2 = zeros(z,numel(tc));
    Xe131RT2 = zeros(z,numel(tc));
    Xe132RT2 = zeros(z,numel(tc));
    Xe134RT2 = zeros(z,numel(tc));
    Xe136RT2 = zeros(z,numel(tc));
    Xe130T2 = zeros(z,numel(tc));
    Xe129CT2 = zeros(z,numel(tc));
    Xe130RaT2 = zeros(z,numel(tc));
    Xe130Rsteps = zeros(z,numel(tc));
    He3T2 = zeros(z,numel(tc)); 
    He4T2 = zeros(z,numel(tc));
    Ne20T2 = zeros(z,numel(tc));
    Ne21T2 = zeros(z,numel(tc));
    Ne22T2 = zeros(z,numel(tc));
    NeR20T2 = zeros(z,numel(tc));
    NeR21T2 = zeros(z,numel(tc));
    NeR22T2 = zeros(z,numel(tc));
    
    %% Reinitialize last values with outputs from early degassing loop using deal function
    [Ulast U235last Thlast Pulast XeUlast XeU131last XeU132last XeU134last Xe136last Xe131last Xe132last Xe134last XeR130last XeR128last XeR129last XeR131last XeR132last XeR134last XeR136last Xe130Rlast Xe129last Xe130last Xe130FDlast Xe129Clast Xe129CFDlast He3last He4last Ne20last Ne21last Ne22last NeR20last NeR21last NeR22last] = deal(SaveEarly_outputs{:});

    for j = startidx:numel(t)
        tnow = t(j);
        tlast = t(j-1);
        dt = tnow-tlast;

        dM = repmat(Mproc(:,:,(j-startidx+1)),[1 numel(tc)]);     %% block of Mass processed between tlast and tnow; dim1 = #res, dim2 = tclosure, dim3 = Mres
%         dM = dt.*Qp*10;
        Q = repmat(tBestQ,[1 numel(tc)]);   %% block of Q's Nres x tc x Mres:dCC 
        dCC = repmat(dCCX(j),[1 numel(tc)]);%% block of dCC Nres x tc x Mres:dCC 
        dUCC = Q.*dCC.*6.02e23;                      %% atoms/g U lost from mantle in this step

        %% Radioactive parent isotope evolution (decay and sequestration in CC)
        Unow = Ulast.*exp(-lambda238*dt) - dUCC;
        U235now = U235last.*exp(-lambda235*dt) - U235U238(j)*dUCC;
        Thnow = Thlast.*exp(-lambda232*dt) - Th232U238(j)*dUCC;
        Punow = Pulast.*exp(-lambda244*dt) - Pu244U238(j)*dUCC;
        
        %% U fission Xe isotope evolution (production and degassing)
        XeUnow = (XeUlast + Ulast.*(1-exp(-lambda238*dt)).*yU136) .* (1 - dM./Mres);        %Xe136U
        XeU131now = (XeU131last + Ulast.*(1-exp(-lambda238*dt)).*yU131) .* (1 - dM./Mres);  %Xe131U
        XeU132now = (XeU132last + Ulast.*(1-exp(-lambda238*dt)).*yU132) .* (1 - dM./Mres);  %Xe132U
        XeU134now = (XeU134last + Ulast.*(1-exp(-lambda238*dt)).*yU134) .* (1 - dM./Mres);  %Xe134U

        %% Pu fission Xe isotope evolution (production and degassing)
        Xe136now = (Xe136last + Pulast.*(1-exp(-lambda244*dt)).*yPu136) .* (1 - dM./Mres);
        Xe131now = (Xe131last + Pulast.*(1-exp(-lambda244*dt)).*yPu131) .* (1 - dM./Mres);
        Xe132now = (Xe132last + Pulast.*(1-exp(-lambda244*dt)).*yPu132) .* (1 - dM./Mres);
        Xe134now = (Xe134last + Pulast.*(1-exp(-lambda244*dt)).*yPu134) .* (1 - dM./Mres);
        
        %% Regassed isotope evolution (regassing and degassing)
        XeR130now = XeR130last + (dM./Mres).*(Cs130(j)-XeR130last);
        XeR128now = XeR128last + (dM./Mres).*((Cs130(j).*Xe128130St(j)) - XeR128last);
        XeR129now = XeR129last + (dM./Mres).*((Cs130(j).*Xe129130St(j)) - XeR129last);
        XeR131now = XeR131last + (dM./Mres).*((Cs130(j).*Xe131130St(j)) - XeR131last);
        XeR132now = XeR132last + (dM./Mres).*((Cs130(j).*Xe132130St(j)) - XeR132last);
        XeR134now = XeR134last + (dM./Mres).*((Cs130(j).*Xe134130St(j)) - XeR134last);
        XeR136now = XeR136last + (dM./Mres).*((Cs130(j).*Xe136130St(j)) - XeR136last);
        
        Xe130Rnow = Xe130Rlast + dM.*Cs130(j);      %% ATOMS COUNT -- atoms/g * non-normalized mass processed (should end up not mattering)         
        
        Xe129now = (Xe129last) .* (1 - dM./Mres);   %% treating this as a primordial isotope for plumes
        
        %% Primordial 130Xe budget
        Xe130now = Xe130last.*(1 - dM./Mres); Xe130FDnow = Xe130FDlast.*(1 - dM./Mres);
        Xe129Cnow = Xe129Clast.*(1 - dM./Mres); Xe129CFDnow = Xe129CFDlast.*(1 - dM./Mres);

        He3now = He3last .* (1 - dM./Mres);
        Herad = 8*Ulast.*(1-exp(-lambda238*dt)) + 7*U235last.*(1-exp(-lambda235*dt)) + 6*Thlast.*(1-exp(-lambda232*dt));
        He4now = (He4last + Herad) .* (1 - dM./Mres);
        
        Ne21now = (Ne21last + Herad.*Ne21He4prod) .* (1 - dM./Mres);
        Ne20now = Ne20last .* (1 - dM./Mres);
        Ne22now = Ne22last .* (1 - dM./Mres);
        
        NeR22now = NeR22last + (dM./Mres).*((Cs130(j).*Ne22Xe130S)-NeR22last);
        NeR21now = NeR21last + (dM./Mres).*((Cs130(j).*Ne22Xe130S.*0.029)-NeR21last);       %% assume atmospheric Ne isotopic evolution is negligible
        NeR20now = NeR20last + (dM./Mres).*((Cs130(j).*Ne22Xe130S.*9.8)-NeR20last);
        
        Ulast = Unow; Pulast = Punow; %Ilast = Inow;
        XeUlast = XeUnow; Xe136last = Xe136now; Xe129last = Xe129now;

        Xe130last = Xe130now; Xe129Clast = Xe129Cnow; Xe130FDlast = Xe130FDnow; Xe129CFDlast = Xe129CFDnow;
        Xe131last = Xe131now; Xe132last = Xe132now; Xe134last = Xe134now; 
        XeU131last = XeU131now; XeU132last = XeU132now; XeU134last = XeU134now;
        XeR130last = XeR130now; XeR128last = XeR128now; XeR129last = XeR129now; XeR131last = XeR131now; XeR132last = XeR132now; XeR134last = XeR134now; XeR136last = XeR136now; 
        Xe130Rlast = Xe130Rnow; %% regassing atom count
        U235last = U235now; Thlast = Thnow;
        He3last = He3now; He4last = He4now;
        Ne20last = Ne20now; Ne21last = Ne21now; Ne22last = Ne22now;
        NeR20last = NeR20now; NeR21last = NeR21now; NeR22last = NeR22now;
 
        %% Fill in time series
        %IT2((j-startidx+1),:) = Inow;
        Xe129T2((j-startidx+1),:) = Xe129now;
        PuT2((j-startidx+1),:) = Punow;
        Xe136T2((j-startidx+1),:) = Xe136now;   %% Pu-fission 136Xe
        Xe131T2((j-startidx+1),:) = Xe131now;
        Xe132T2((j-startidx+1),:) = Xe132now;
        Xe134T2((j-startidx+1),:) = Xe134now;
        UT2((j-startidx+1),:) = Unow;
        XeUT2((j-startidx+1),:) = XeUnow;       %% U-fission 136Xe
        Xe131UT2((j-startidx+1),:) = XeU131now;
        Xe132UT2((j-startidx+1),:) = XeU132now;
        Xe134UT2((j-startidx+1),:) = XeU134now;
        Xe128RT2((j-startidx+1),:) = XeR128now;
        Xe129RT2((j-startidx+1),:) = XeR129now;
        Xe130RT2((j-startidx+1),:) = XeR130now;
        Xe131RT2((j-startidx+1),:) = XeR131now;
        Xe132RT2((j-startidx+1),:) = XeR132now;
        Xe134RT2((j-startidx+1),:) = XeR134now;
        Xe136RT2((j-startidx+1),:) = XeR136now;
        Xe130T2((j-startidx+1),:) = Xe130now;
        Xe130RaT2((j-startidx+1),:) = Xe130Rnow;    %% TOTAL
        Xe130Rsteps((j-startidx+1),:) = dM.*Cs130(j);	%% how many atoms of 130Xe are recycled per step?
        Xe129CT2((j-startidx+1),:) = Xe129Cnow;
        He3T2(j-startidx+1,:) = He3now;
        He4T2(j-startidx+1,:) = He4now;
        Ne20T2(j-startidx+1,:) = Ne20now;
        Ne21T2(j-startidx+1,:) = Ne21now;
        Ne22T2(j-startidx+1,:) = Ne22now;
        NeR20T2(j-startidx+1,:) = NeR20now;
        NeR21T2(j-startidx+1,:) = NeR21now;
        NeR22T2(j-startidx+1,:) = NeR22now;
        
    end
% toc
% %% Traditional stuff
% output.tc = tc;
% output.Xe129Xe136 = (Xe129now + Inow)./(Xe136now + Punow.*yPu136);
% output.XePuU = Xe136now./(XeUnow+Xe136now);
% output.Xe130 = Xe130now;
% output.Xe129130 = (Xe129Cnow+Xe129now)./Xe130now; 
% output.Xe130DG = Xe130DG;           %% Xe130DG = (Xe130now./Xe130NOGnow)';
% output.Xe129DG = Xe129DG;
% output.Xe129DGLGI = Xe129DGsnapshot;
% output.Xe130DGLGI = Xe130DGsnapshot;
% output.Xe128130 = (XeR128now+(Xe130now.*components(1,2,1)./components(2,2,1)))./(XeR130now+Xe130now);


%% Make primordial evolution time series
XeP132_2 = Xe130T2./components(2,2,1);
XeP128_2 = XeP132_2.*components(1,2,1);
XeP131_2 = XeP132_2.*components(3,2,1);
XeP134_2 = XeP132_2.*components(4,2,1);
XeP136_2 = XeP132_2.*components(5,2,1);
He3Xe130P = He3T2./Xe130T2;
Ne22Xe130P = Ne22T2./Xe130T2;

Xe128M = (Xe128RT2+(Xe130T2.*components(1,2,1)./components(2,2,1)));
Xe129M = (Xe129RT2+(Xe130T2.*Xe129130C)+Xe129T2);
Xe130M = (Xe130RT2+Xe130T2);
Xe131M = (Xe131RT2+(Xe130T2.*components(3,2,1)./components(2,2,1))+Xe131T2+Xe131UT2);
Xe132M = (Xe132RT2+(Xe130T2./components(2,2,1))+Xe132T2+Xe132UT2);
Xe134M = (Xe134RT2+(Xe130T2.*components(4,2,1)./components(2,2,1))+Xe134T2+Xe134UT2);
Xe136M = (Xe136RT2+(Xe130T2.*components(5,2,1)./components(2,2,1))+Xe136T2+XeUT2);
Ne20M = Ne20T2 + NeR20T2;
Ne21M = Ne21T2 + NeR21T2;
Ne22M = Ne22T2 + NeR22T2;

Xe136Putotal = Xe136now./(XeUnow+Xe136now);
%Xe136T2(end)./(XeUT2(end)+Xe136T2(end));

Xe128130M = Xe128M./Xe130M;
Xe129130M = Xe129M./Xe130M;
Xe131130M = Xe131M./Xe130M;
Xe132130M = Xe132M./Xe130M;
Xe134130M = Xe134M./Xe130M;
Xe136130M = Xe136M./Xe130M;

Xe128132M = Xe128M./Xe132M;
Xe129132M = Xe129M./Xe132M;
Xe130132M = Xe130M./Xe132M;
Xe131132M = Xe131M./Xe132M;
Xe134132M = Xe134M./Xe132M;
Xe136132M = Xe136M./Xe132M;

Ne2022M = Ne20M./Ne22M;
Ne2122M = Ne21M./Ne22M;

He43M = He4T2./He3T2;

He3Xe130M = He3T2./Xe130M;
Xe13022NeM = Xe130M./Ne22M;
He3Ne22 = He3T2./Ne22M;

%% Store results
store128130M(:,mc) = Xe128130M; store128132M(:,mc) = Xe128132M;
store128130R(:,mc) = Xe128RT2./Xe130RT2;
store129130M(:,mc) = Xe129130M; store129132M(:,mc) = Xe129132M;
store130132M(:,mc) = Xe130132M; store132130M(:,mc) = Xe132130M;
store131130M(:,mc) = Xe131130M; store131132M(:,mc) = Xe131132M;
store134130M(:,mc) = Xe134130M; store134132M(:,mc) = Xe134132M;
store136130M(:,mc) = Xe136130M; store136132M(:,mc) = Xe136132M;
storePuUM(mc) = Xe136Putotal;
store130M(:,mc) = Xe130M;
store130fracs(:,mc) = [Xe130T2(end) Xe130RT2(end)]./(Xe130RT2(end)+Xe130T2(end));
store132fracs(:,mc) = [Xe130T2(end)./components(2,2,1) Xe132T2(end) Xe132UT2(end) Xe132RT2(end)]./Xe132M(end);
store130Rsteps(:,mc) = Xe130Rsteps;
storeCs130(:,mc) = Cs130(startidx:end);
store2022M(:,mc) = Ne2022M;
storeNe22M(:,mc) = Ne22M;
store22fracs(:,mc) = [Ne22T2(end) NeR22T2(end)]./Ne22M(end);
store13022x(:,mc) = [Xe130T2(end)/Ne22T2(end) Xe130RT2(end)/NeR22T2(end)];
store2122M(:,mc) = Ne2122M;
store43M(:,mc) = He43M;
storeHe3M(:,mc) = He3T2;
store3130M(:,mc) = He3Xe130M;
store13022M(:,mc) = Xe13022NeM;
store322M(:,mc) = He3Ne22;
store13022R(:,mc) = Xe130RT2./NeR22T2;
uuu = t(abs(Cs130-(Cs130(1)+0.5*(Cs130(end)-Cs130(1)))) == min(abs(Cs130-(Cs130(1)+0.5*(Cs130(end)-Cs130(1))))));
                   
end


itscloseidx = abs(store130132M(end,:) - plotx(3,1)) < 2*plotx(3,2);
closetwice = itscloseidx & (abs(store134132M(end,:) - plotx(5,1)) < 2*plotx(5,2));

Xescore = ((store130132M(end,:) - plotx(3,1))./ (plotx(3,2))).^2 + ...
    ((store131132M(end,:) - plotx(4,1))./ (plotx(4,2))).^2 + ...
    ((store134132M(end,:) - plotx(5,1))./ (plotx(5,2))).^2 + ...
    ((store136132M(end,:) - plotx(6,1))./ (plotx(6,2))).^2 + ...%;
    ((store128130M(end,:) - plotx(1,1))./ (plotx(1,2))).^2;

% BestXeidx = Xescore == min(Xescore);
BestXeidx = abs(store130fracs(2,:) - Icef130a) == min(abs(store130fracs(2,:) - Icef130a));

% test_res = abs(mean(store130fracs(2,BestXeidx)) - Icef130a);
% Icef130a = mean(store130fracs(2,BestXeidx));
% end

score = zeros(10,size(store130132M(end,BestXeidx),2));
score(1,:) = 2*((store130132M(end,BestXeidx) - plotx(3,1)) ./ plotx(3,2)).^2;
score(2,:) = 2*((store134132M(end,BestXeidx) - plotx(5,1)) ./ plotx(5,2)).^2;
score(3,:) = 2*((store136132M(end,BestXeidx) - plotx(6,1)) ./ plotx(6,2)).^2;
score(4,:) = 2*((store128130M(end,BestXeidx) - plotx(1,1)) ./ plotx(1,2)).^2;
score(5,:) = ((store2022M(end,BestXeidx) - Ice2022)./0.02).^2 + ((store2122M(end,BestXeidx) - targetNe(2))./0.0003).^2;

% score(5,:) = abs(((store2022M(end,BestXeidx) - 9.8)./(store2122M(end,BestXeidx) - 0.029)) - ((12.6-9.8)/(0.0590-0.029)))./8; %%ERROR MORB NUMBERS?
score(6,:) = ((store43M(end,BestXeidx) - mean(targetHe))./3000).^2;
% if mean(store43M(end,closetwice)) > targetHe(1) && mean(store43M(end,closetwice)) < targetHe(2)
%     score(6,:) = 0;
% else
%     score(6,:) = 100;
% end
score(7,:) = ((store3130M(end,BestXeidx) - Ice3130)./100).^2;
score(8,:) = ((store13022M(end,BestXeidx) - Ice13022)./0.0005).^2;
score(9,:) = ((store322M(end,BestXeidx) - Ice322)./0.1).^2;
score(10,:) = ((store130fracs(2,BestXeidx) - Icef130a)./0.02).^2;

totalscores = sum(score,1);
bestscore = totalscores == min(totalscores);
mapvx(grid1,grid2) = min(totalscores);
Xetotals(grid1,grid2) = sum(score(1:4,bestscore));
Netotals(grid1,grid2) = sum(score(5,bestscore));
Hetotals(grid1,grid2) = sum(score(6,bestscore));
elemtotals(grid1,grid2) = sum(score(7:9,bestscore));
f130ares(grid1,grid2) = sum(score(10,bestscore));


    end
end

% 
save OIBgrid_01012021_acct11_Mres3_Nres2to45_w128_nowhile_005_098.mat

% % 
% % 
% % figure;
% % tblock = repmat(t(startidx:end),[1 mcnum]);
% % plot(tblock,store3130M);
% % plot(store3130M(end,:),store130fracs(2,:),'k.')
% % 
% % itscloseidx = abs(store130132M(end,:) - plotx(3,1)) < 2*plotx(3,2);
% % closetwice = itscloseidx & (abs(store134132M(end,:) - plotx(5,1)) < 2*plotx(5,2));
% % 
% % store130fracs(2,closetwice)
% % 

itscloseidx = BestXeidx;
closetwice = false(size(closetwice));%abs(store130fracs(2,:) - Icef130a) < 0.01;
figure; 
subplot(2,3,1)
plot(store3130M(end,:),store43M(end,:),'k.'); hold on;
rectangle('Position',[targetHeXe(1) targetHe(1) targetHeXe(2)-targetHeXe(1) targetHe(2)-targetHe(1)],'EdgeColor','g');
plot(store3130M(end,itscloseidx),store43M(end,itscloseidx),'c.');
plot(store3130M(end,closetwice),store43M(end,closetwice),'m.');
xlabel('3He/130Xe'); ylabel('4He/3He');
xlim([0 1500]); axis square; box on;

subplot(2,3,2)
plot(store2122M(end,:),store2022M(end,:),'k.'); hold on;
plot(0.0356,12.9,'r*');
plot(0.0362,13.17,'g*');
plot(0.029,9.8,'ks');
plot(store2122M(end,itscloseidx),store2022M(end,itscloseidx),'c.');
plot(store2122M(end,closetwice),store2022M(end,closetwice),'m.');
xlabel('21Ne/22Ne'); ylabel('20Ne/22Ne');
cz = 'k'; axis square; box on;

subplot(2,3,3)
plot(Ice13022,Ice322,'gs'); xlabel('130Xe/22Ne'); ylabel('3He/22Ne'); axis square; box on; hold on;
plot(store13022M(end,:),store322M(end,:),'k.');
plot(store13022M(end,itscloseidx),store322M(end,itscloseidx),'c.');
plot(store13022M(end,closetwice),store322M(end,closetwice),'m.');

subplot(2,3,4)
errorbar(plotx(3,1),plotx(4,1),plotx(4,2),plotx(4,2),plotx(3,2),plotx(3,2)); xlabel('130Xe/132Xe'); ylabel('131Xe/132Xe'); axis square; box on;
hold on; plot(plotx(3,1),plotx(4,1),'yd');
plot(components(2,1,1),components(3,1,1),'bs','Color',[0 182/255 1],'MarkerFaceColor',[0 182/255 1],'MarkerSize',7); plot(gca,components(2,2,1),components(3,2,1),'kd','Color',[107 61 6]./255,'MarkerFaceColor',[107 61 6]./255,'MarkerSize',7); plot(components(2,2,3),components(3,2,3),'rs', 'Color',[255 124 0]./255,'MarkerFaceColor',[255 124 0]./255,'MarkerSize',7);
plot([components(2,1,1) components(2,3,1)],[components(3,1,1) components(3,3,1)],'b-');
plot([components(2,1,1) components(2,4,1)],[components(3,1,1) components(3,4,1)],'r-');
scatter(store130132M(end,:),store131132M(end,:),4,cz,'filled');
scatter(store130132M(end,itscloseidx),store131132M(end,itscloseidx),4,'c','filled');
scatter(store130132M(end,closetwice),store131132M(end,closetwice),4,'m','filled');
xlim([0.145 0.165]); ylim([0.76 0.83]);

subplot(2,3,5)
errorbar(plotx(3,1),plotx(5,1),plotx(5,2),plotx(5,2),plotx(3,2),plotx(3,2)); xlabel('130Xe/132Xe'); ylabel('134Xe/132Xe'); axis square; box on
hold on; plot(plotx(3,1),plotx(5,1),'yd');
plot(components(2,1,1),components(4,1,1),'bs','Color',[0 182/255 1],'MarkerFaceColor',[0 182/255 1],'MarkerSize',7); plot(gca,components(2,2,1),components(4,2,1),'kd','Color',[107 61 6]./255,'MarkerFaceColor',[107 61 6]./255,'MarkerSize',7); plot(components(2,2,3),components(4,2,3),'rs', 'Color',[255 124 0]./255,'MarkerFaceColor',[255 124 0]./255,'MarkerSize',7);
plot([components(2,1,1) components(2,3,1)],[components(4,1,1) components(4,3,1)],'b-');
plot([components(2,1,1) components(2,4,1)],[components(4,1,1) components(4,4,1)],'r-');
scatter(store130132M(end,:),store134132M(end,:),4,cz,'filled');
scatter(store130132M(end,itscloseidx),store134132M(end,itscloseidx),4,'c','filled');
scatter(store130132M(end,closetwice),store134132M(end,closetwice),4,'m','filled');
xlim([0.145 0.165]); ylim([0.375 0.415]);

subplot(2,3,6)
errorbar(plotx(3,1),plotx(6,1),plotx(6,2),plotx(6,2),plotx(3,2),plotx(3,2)); xlabel('130Xe/132Xe'); ylabel('136Xe/132Xe'); axis square; box on
hold on; plot(plotx(3,1),plotx(6,1),'yd');
plot(components(2,1,1),components(5,1,1),'bs','Color',[0 182/255 1],'MarkerFaceColor',[0 182/255 1],'MarkerSize',7); plot(gca,components(2,2,1),components(5,2,1),'kd','Color',[107 61 6]./255,'MarkerFaceColor',[107 61 6]./255,'MarkerSize',7); plot(components(2,2,3),components(5,2,3),'rs', 'Color',[255 124 0]./255,'MarkerFaceColor',[255 124 0]./255,'MarkerSize',7);
plot([components(2,1,1) components(2,3,1)],[components(5,1,1) components(5,3,1)],'b-');
plot([components(2,1,1) components(2,4,1)],[components(5,1,1) components(5,4,1)],'r-');
scatter(store130132M(end,:),store136132M(end,:),4,cz,'filled');
scatter(store130132M(end,itscloseidx),store136132M(end,itscloseidx),4,'c','filled');
scatter(store130132M(end,closetwice),store136132M(end,closetwice),4,'m','filled');
xlim([0.145 0.165]); ylim([0.3 0.38]);

extable = [Nres Xe130i Icef130a Ice13022 Ice3130 Xe130Ne22 xPLUME']';

figure;
subplot(1,3,1); hold on;
plot(t(startidx:end),storeHe3M(:,BestXeidx),'r-');
xlabel('time'); ylabel('3He atoms/g'); axis square; box on; xlim([0 T]);
subplot(1,3,2); hold on;
plot(t(startidx:end),storeNe22M(:,BestXeidx),'r-');
xlabel('time'); ylabel('22Ne atoms/g'); axis square; box on; xlim([0 T]);
subplot(1,3,3); hold on;
plot(t(startidx:end),store130M(:,BestXeidx),'r-');
xlabel('time'); ylabel('130Xe atoms/g'); axis square; box on; xlim([0 T]);

storeHe3M(end,BestXeidx)
storeNe22M(end,BestXeidx)
store130M(end,BestXeidx)
abs(Ice3130 - storeHe3M(end,BestXeidx)./store130M(end,BestXeidx))

store132fracs(:,BestXeidx)
(ans(2)*components(5,4,1))/((ans(2)*components(5,4,1))+(ans(3)*components(5,3,1)))
tc/1e6

% Spreidx = zeros(4,mcnum);
% Spreidx(1,:) = abs(store130132M(end,:) - plotx(3,1)) < 2*plotx(3,2);
% Spreidx(2,:) = abs(store134132M(end,:) - plotx(5,1)) < 2*plotx(5,2);
% Spreidx(3,:) = abs(store136132M(end,:) - plotx(6,1)) < 2*plotx(6,2);
% Spreidx(4,:) = store2022M(end,:) > 12.9;
% 
% Sidx = sum(Spreidx,1) == 4; %% all criteria are true
% figure;
% plot(t(startidx:end),storeCs130(Sidx),'b-')

% sndit.Xe128130M = store128130M(end,:);
% sndit.Xe128132M = store128132M(end,:);
% sndit.Xe129130M = store129130M(end,:);
% sndit.Xe129132M = store129132M(end,:);
% sndit.Xe130132M = store130132M(end,:);
% sndit.Xe132130M = store132130M(end,:);
% sndit.Xe131130M = store131130M(end,:);
% sndit.Xe131132M = store131132M(end,:);
% sndit.Xe134130M = store134130M(end,:);
% sndit.Xe134132M = store134132M(end,:);
% sndit.Xe136130M = store136130M(end,:);
% sndit.Xe136132M = store136132M(end,:);
% sndit.Xe130M = store130M(end,:);
% sndit.Ne2022M = store2022M(end,:);
% sndit.Ne2122M = store2122M(end,:);
% sndit.Neslope = (store2022M(end,:) - 9.8)./(store2122M(end,:) - 0.029);
% sndit.He43M = store43M(end,:);
% sndit.He3130M = store3130M(end,:);
% sndit.Xe13022M = store13022M(end,:);
% sndit.He322M = store322M(end,:);
% 
% output = assess_residuals_plume_target(sndit,plotx,RS);

% mean(store132fracs(4,itscloseidx))
% store13022x(1,itscloseidx)
% store13022x(2,itscloseidx)
% store130fracs(2,itscloseidx)
% store22fracs(2,itscloseidx)

% xPLUME
% end
