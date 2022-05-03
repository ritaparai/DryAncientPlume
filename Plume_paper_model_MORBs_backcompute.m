%% This is a fresh page from which to do the plume paper calculations where
%% I remove the MC and back-compute regassing to make the grid
%% RParai June 24 2021

clear; close all;

%% Initialize rand and load components
s = RandStream('mt19937ar','Seed', sum(1e4*clock));
RandStream.setGlobalStream(s);
load components2
load Xemasses
load o

Ms = 3.6e27;
run.acct = 18;
PRf130a = 0.75;
PR3130 = 1000;

run.late_t = 1e6*[0:0.1:200 201:1:1302 1307:5:4567]';
run.BestQ = findBestQ_MORBs(o,run.late_t,Ms);
whatCCx = 1;
eval(strcat('dCCX = o.dCC', num2str(whatCCx),';'));
tBestQ = run.BestQ(1,1,sub2ind([3,3],1,whatCCx));

mm = Xe_masses;
run.whatM = 1;          %% 50%, 75%, 90%

adjustI = 1;
lv = 1;         %% 0.01 0.005 0.001 (they're in chronological order, and the big number happens first)

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
% run.acct = 18;

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

%% Boundaries of the solution space
%% ========================================================================
Xe130i = 5.38e-14*6.02e23;      %% 130Xe [atoms/g] (Marty, 2012) Murchison and Orgueil average
Xe130M_T = 9.45e-19*6.02e23;    %% 130Xe present-day upper mantle [atoms/g]; computed based on Bianchi and 3He/130Xe
CS = [4.3e5 9.2e5];             %% concentration success boundaries
% CS = Xe130M_T.*[0.5 2];         %% concentration success boundaries
RS = zeros(7,2);
RS(1,:) = [0.475 0.478];        %% 2s ratio success boundaries based on wellgasbd20bsm, caffeeharding, eqdm, swirwest and sqwireast
RS(2,:) = [0.0690 0.0710];
RS(3,:) = [1.058 1.134];
RS(4,:) = [0.1445 0.1493];
RS(5,:) = [0.7608 0.7786];
RS(6,:) = [0.4082 0.4302];
RS(7,:) = [0.3559 0.3835];


%% Modern and initial atmospheric Xe compositions
%% ========================================================================
Xe128130air = components(1,1,1)./components(2,1,1);

%% Atmospheric Xe isotopic time series (x-axis)
ff = (-0.0777/2.5e9).*run.late_t + 1.0777; ff(ff<1)=1;                 %% fractionations in a vector based on run.late_t; matched 128/130 to U-Xe
% ff(:)=1;    %% TURN OFF Xe ATMOSPHERIC EVOLUTION
Xe128130St = Xe128130air.*ff;
a = (1-( (mm(3)./mm([2 4:7])).^(1/2) ))./(1-sqrt(mm(3)/mm(1)));
b = ((mm(3)./mm([2 4:7])).^(1/2)) ./ (sqrt(mm(3)/mm(1)).^a);
% b = b./b;   %% TURN OFF Xe ATMOSPHERIC EVOLUTION
Xe129130St = 6.496.*b(1).* ((ff).^a(1));
Xe131130St = components(3,1,1)./components(2,1,1).*b(2).* ((ff).^a(2));
Xe132130St = 1./components(2,1,1).*b(3).* ((ff).^a(3));
Xe134130St = components(4,1,1)./components(2,1,1).*b(4).* ((ff).^a(4));
Xe136130St = components(5,1,1)./components(2,1,1).*b(5).* ((ff).^a(5));

%% Primordial isotope calculations ========================================

Ne2022atm = 9.8;
Xe130Ne22atm = 0;
Xe130Ne22atm_nom = 0.0021;

%% Order is SOLAR CC ATM
Ne2022 = [13.36 8.9 9.8]; 
Xe130Ne22 = [1.03e-6 0.033 0];  %% CHANGE TO 0.015 IF NEEDED
A = [Ne2022; Xe130Ne22; [1 1 1]];


%% POPPING ROCK || (PARAI AND MUKHOPADHYAY, 2020)
%% 20Ne/22Ne 
PR2022 = 12.6;
%% 130Xe/22Ne
PR13022 = 0.0078;
% %% fraction 130Xe atm
% PRf130a = 0.97;
b = [PR2022; PR13022.*(1-PRf130a); 1];
%% 3He/22Ne
He3Ne22sol = 1.46;  %% Tucker and Mukhopadhyay 2014
He3Ne22cc = 0.89;
He3Ne22atm = 4.3e-6;

xMORB = A\b;

check1 = sum(xMORB .* Ne2022') == PR2022;
Ne2022acc = (sum(xMORB(1:2) .* Ne2022(1:2)'))./sum(xMORB(1:2));
Ne2122acc = (sum(xMORB(1:2) .* [0.0324; 0.0268]))./sum(xMORB(1:2));     %% Heber et al 2012 adjusted following Williams and Mukhopadhyay 2019
Ne22M = 1;
Ne22sol = Ne22M.*xMORB(1);
Ne22cc = Ne22M.*xMORB(2);
Ne22atm = Ne22M.*xMORB(3);

Xe130M = Ne22M.*PR13022; %% Total 130Xe in the mantle today assuming 22Ne = 1
Xe130sol = Ne22sol.*Xe130Ne22(1);
Xe130cc = Ne22cc.*Xe130Ne22(2);
Xe130atm_part = Ne22atm.*Xe130Ne22atm_nom;
XewithNe = Xe130sol + Xe130cc + Xe130atm_part;
Xe130excess = Xe130M - XewithNe;
Xe130atm = Xe130atm_part + Xe130excess;
Xe130Ne22_R = Xe130atm ./ Ne22atm;
check2 = abs(Xe130atm./Xe130M - PRf130a) < 1e-15;

He3Xe130 = PR3130./(1-PRf130a);   %% the 3He/130Xe ratio at accretion
He3sol = Ne22sol.*He3Ne22sol;
He3cc = Ne22cc.*He3Ne22cc;
He3atm = Ne22atm.*He3Ne22atm;

He3Xe130acc = He3Xe130;
Xe130Ne22acc = (Xe130sol + Xe130cc) ./ (Ne22sol + Ne22cc);

PR322 = PR3130 .* PR13022;


%% POPPING ROCK
He3Xe130i = He3Xe130acc;
Xe129Xe130 = 10.5;
Xe129136star = 2.8;
Xe136PuXe130 = (Xe129Xe130 - Xe129130C)./ Xe129136star;
NeXeslabs = 1./Xe130Ne22_R;
targetHe = [80000; 90000];
targetNe = [PR2022; 0.0590];
targetHeXe = PR3130 * [0.9 1.1];
plotx = [0.4775 0.0025; 1.104 3e-3; 0.1489 1e-4; 0.7763 7e-4; 0.4145 4e-4; 0.3643 4e-4];  %% PR 2PD43 using 40Ar/36Ar of 27000
% plotx = [NaN NaN; 1.101 3e-3; 0.1489 1e-4; 0.7766 7e-4; 0.4138 4e-4; 0.3634 4e-4];  %% PR 2PD43 using 40Ar/36Ar of 23900
% plotx = [NaN NaN; 1.118 8e-3; 0.1449 4e-4; 0.7626 17e-4; 0.4306 11e-4; 0.3851 12e-4];  %% Eq Atlantic DM Tucker et al 2012

%% ===============================

%% Start with the parameter ranges
% Nrange = linspace(4.2,5.2,10);
% initrange = linspace(4.2e9,6.2e9,10);
Nrange = linspace(3.8,4.8,100);
initrange = linspace(4e9,7e9,100);

mapvx = zeros(numel(Nrange),numel(initrange));
Xetotals = nan(numel(Nrange),numel(initrange));
Netotals = nan(numel(Nrange),numel(initrange));
Hetotals = nan(numel(Nrange),numel(initrange));
elemtotals = nan(numel(Nrange),numel(initrange));
% grid1 = 50; grid2 = 43;
% Nres = 4.5;
% Xe130i = 7.8e9; 
for grid1 = 1:numel(Nrange)
    Nres = Nrange(grid1);
    for grid2 = 1:numel(initrange)
        Xe130i = initrange(grid2);


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

% Xe130last;

%% INTERMEDIATE AND LONG-TERM DEGASSING TIME LOOPS
tlast = t(dend);
t = run.late_t;
startidx = find(t>tlast,1);
Mproc = findMproc_MORBs_spec(run.late_t,startidx,Nres,Ms);

tblock = t;                                                     %% tblock (and therefore Cs130) based on run.late_t
z = numel(t(startidx:end));

% decoarse = [[0:0.1:200].*0+0.1 [201:1:1302].*0+1 [1307:5:4567].*0+5]';
% plot(tblock(startidx:end),reshape(Mproc,[numel(Mproc) 1 1])./decoarse(startidx:end)./1e6,'b-')
% sum(Mproc)./Ms

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
Ne21last = Ne22last.*Ne2122acc;         %% solar, largely; Williams and Mukhopadhyay (2019) recalc of Heber et al. (2012)
Ne22Xe130S = NeXeslabs;


U235U238 = (0.0072*exp(lambda235*(4.567e9-t)))./(0.992745*exp(lambda238*(4.567e9-t)));  %% 235U/238U ratio as a function of time
Th232U238 = ((80e-09/232)*exp(lambda232*(4.567e9-t)))./((BSEU*0.992745/238)*exp(lambda238*(4.567e9-t))); %% 232Th/238U ratio as a function of time
Pu244U238 = (Pui*exp(-lambda244*t))./(Ui*exp(-lambda238*t));

%% LONG-TERM DEGASSING TIME LOOPS

Xe130FDlast = 0;
Xe129Clast = 0;
Xe129CFDlast = 0;
SaveEarly_outputs = num2cell([Ulast U235last Thlast Pulast XeUlast XeU131last XeU132last XeU134last Xe136last Xe131last Xe132last Xe134last Xe129last Xe130last Xe130FDlast Xe129Clast Xe129CFDlast He3last He4last Ne20last Ne21last Ne22last]);
    
%% Reinitialize last values with outputs from early degassing loop using deal function
[Ulast U235last Thlast Pulast XeUlast XeU131last XeU132last XeU134last Xe136last Xe131last Xe132last Xe134last Xe129last Xe130last Xe130FDlast Xe129Clast Xe129CFDlast He3last He4last Ne20last Ne21last Ne22last] = deal(SaveEarly_outputs{:});

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


    Ulast = Unow; Pulast = Punow; %Ilast = Inow;
    XeUlast = XeUnow; Xe136last = Xe136now; Xe129last = Xe129now;

    Xe130last = Xe130now; Xe129Clast = Xe129Cnow; Xe130FDlast = Xe130FDnow; Xe129CFDlast = Xe129CFDnow;
    Xe131last = Xe131now; Xe132last = Xe132now; Xe134last = Xe134now; 
    XeU131last = XeU131now; XeU132last = XeU132now; XeU134last = XeU134now;

    U235last = U235now; Thlast = Thnow;
    He3last = He3now; He4last = He4now;
    Ne20last = Ne20now; Ne21last = Ne21now; Ne22last = Ne22now;

end

%% Compute primordial budgets
XeP132 = Xe130now./components(2,2,1);
XeP128 = XeP132.*components(1,2,1);
XeP131 = XeP132.*components(3,2,1);
XeP134 = XeP132.*components(4,2,1);
XeP136 = XeP132.*components(5,2,1);
He3Xe130P = He3now./Xe130now;
Ne22Xe130P = Ne22now./Xe130now;

%% Compute regassed budgets
XeR130 = (PRf130a * Xe130now)./(1-PRf130a);
XeR128 = XeR130 * components(1,1,1)./components(2,1,1);
XeR131 = XeR130 * components(3,1,1)./components(2,1,1);
XeR132 = XeR130 / components(2,1,1);
XeR134 = XeR130 * components(4,1,1)./components(2,1,1);
XeR136 = XeR130 * components(5,1,1)./components(2,1,1);

NeR22 = XeR130 * NeXeslabs;
NeR20 = NeR22 * 9.8;
NeR21 = NeR22 * 0.029;

%% Compute total budgets
Xe130M = Xe130now + XeR130;
Xe128M = XeP128 + XeR128;
Xe131M = XeP131 + XeR131 + Xe131now + XeU131now;
Xe132M = XeP132 + XeR132 + Xe132now + XeU132now;
Xe134M = XeP134 + XeR134 + Xe134now + XeU134now;
Xe136M = XeP136 + XeR136 + Xe136now + XeUnow;

Ne20M = Ne20now + NeR20;
Ne21M = Ne21now + NeR21;
Ne22M = Ne22now + NeR22;

%% Compute ratios
store130132M = Xe130M./Xe132M;
store131132M = Xe131M./Xe132M;
store134132M = Xe134M./Xe132M;
store136132M = Xe136M./Xe132M;
store128130M = Xe128M./Xe130M;
store2022M = Ne20M./Ne22M;
store2122M = Ne21M./Ne22M;
store43M = He4now./He3now;
store3130M = He3now./Xe130M;
store13022M = Xe130M./Ne22M;
store322M = He3now./Ne22M;

itscloseidx = abs(store130132M(end,:) - plotx(3,1)) < 2*plotx(3,2);
closetwice = itscloseidx & (abs(store134132M(end,:) - plotx(5,1)) < 2*plotx(5,2));

Xescore = ((store130132M(end,:) - plotx(3,1))./ (plotx(3,2))).^2 + ...
    ((store131132M(end,:) - plotx(4,1))./ (plotx(4,2))).^2 + ...
    ((store134132M(end,:) - plotx(5,1))./ (plotx(5,2))).^2 + ...
    ((store136132M(end,:) - plotx(6,1))./ (plotx(6,2))).^2 + ...%;
    ((store128130M(end,:) - plotx(1,1))./ (plotx(1,2))).^2;

score = zeros(9,1);
score(1,:) = 2*((store130132M - plotx(3,1)) ./ plotx(3,2)).^2;
score(2,:) = 2*((store134132M - plotx(5,1)) ./ plotx(5,2)).^2;
score(3,:) = 2*((store136132M - plotx(6,1)) ./ plotx(6,2)).^2;
score(4,:) = 2*((store131132M - plotx(4,1)) ./ plotx(4,2)).^2;
score(5,:) = ((store2022M - PR2022)./0.02).^2 + ((store2122M - targetNe(2))./0.0003).^2;
score(6,:) = ((store43M - mean(targetHe))./3000).^2;
score(7,:) = ((store3130M - PR3130)./100).^2;
score(8,:) = ((store13022M - PR13022)./0.0005).^2;
score(9,:) = ((store322M - PR322)./0.1).^2;

totalscores = sum(score,1);
mapvx(grid1,grid2) = min(totalscores);
Xetotals(grid1,grid2) = sum(score(1:4));
Netotals(grid1,grid2) = sum(score(5));
Hetotals(grid1,grid2) = sum(score(6));
elemtotals(grid1,grid2) = sum(score(7:9));

    end
end

save MORBgrid_10212021_acct18_f130a0p75_3He130Xe1000_Mres3p6e27_no128w131_fine

figure
colormap hot
pcolor(mapvx)
axis square; shading flat; box on;
caxis([0 50])


figure
colormap hot
pcolor(mapvx)
axis square; shading flat; box on;
caxis([0 100])

bidx = find(mapvx == min(min(mapvx)));
[grid1 grid2] = ind2sub([numel(Nrange) numel(initrange)], bidx)
Nres = Nrange(grid1)
Xe130i = initrange(grid2)