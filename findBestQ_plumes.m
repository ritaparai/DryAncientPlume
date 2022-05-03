function output = findBestQ_plumes(o,t,Mres)
%% Q is the scaling factor between dCC and dUCC that gives you the right U content at the end of the model run
%% BestQ depends only on Mres and CCmodel
%% It is independent of Nres because the mass processed over a time step doesn't matter -- we are setting the net drop in U content given partial melting and subduction of seds

%% Assume that 10% of the CC U comes from plume source

% Qp = 2.33 * 2.9e15 * 10;       %% g/yr
lambda238 = 1.5514e-10;
BSEU = 21e-09;                      %% 21 ppb U [g/g]
Ui = (BSEU*0.992745/238)/exp(-lambda238*4.567e9);	%% initial U conc [moles/g]

qn = 1e4;
Qs = linspace(2e-12,6e-11,qn);
n = numel(t);
% 
% ch = figure; subplot(3,3,1);
% cj = figure; subplot(3,3,1);
BestQ = zeros(1,1,9);
for H = 1%:9
    [midx didx] = ind2sub([3,3],H);
    eval(strcat('dCC = o.dCC', num2str(didx),';'));
    
    UPM = 21e-09 - (2.2e25*(1.3e-06)*0.1)./Mres;
    tols = zeros(qn,1);
    for q = 1:qn
        Q = Qs(q);
        Ulast = Ui;
        for j = 2:n
            tnow = t(j);
            tlast = t(j-1);
            dt = tnow-tlast;
            dUCC = Q*dCC(j);
            Unow = Ulast.*exp(-lambda238*dt) - dUCC;
            Ulast = Unow;
%             Us(j) = Unow;
%             CCU(j) = CCU(j-1).*exp(-lambda238*dt) + dUCC*Mres(m)/2.2e25;
        end
        tols(q,:) = Unow*238 - UPM;
    end
%     min(abs(tols))
    idx = find(abs(tols)==min(abs(tols)));
    BestQ(1,1,H) = Qs(idx);
    figure; plot(Qs,tols)
%     Q = Qs(idx);
%     Ulast = Ui;
%     Us = ones(n,1)*Ui;
%     CCU = zeros(n,1);
%     for j = 2:n
%         tnow = t(j);
%         tlast = t(j-1);
%         dt = tnow-tlast;
%         dUCC = Q*dCC(j);
%         Unow = Ulast.*exp(-lambda238*dt) - dUCC;
%         Ulast = Unow;
%         Us(j) = Unow;
%         CCU(j) = CCU(j-1).*exp(-lambda238*dt) + dUCC*Mres/2.2e25;
%     end
%     figure(ch); subplot(3,3,H);
%     plot(t,CCU*238,'r-')
%     xlabel('time (Myr)'); ylabel('U ppb'); title('CC evolution');
%     figure(cj); subplot(3,3,H);
%     plot(t,Us*238,'b-');
%     hold on;
%     plot(t, 238*Ui.*(exp(-lambda238.*t)),'k-');
%     xlabel('time (Myr)'); ylabel('U ppb'); legend({'decay and extraction','decay only'}); 
%     title('mantle evolution');
end

output = BestQ;