function dydt = DLC_human_ode(y, t)

global g
global L %number of repeating variables
global M %number of compounds in the mixture
global N %starting index position of repeating variables




%% ------------------------- STATE NAME MAPPING----------------------------%%
%repeating variables
MST               = y(N : L : N+(M-1)*L);
LYMLUM            = y(N+1 : L : N+(M-1)*L+1);
LIMLUM            = y(N+2 : L : N+(M-1)*L+2);
AURI              = y(N+3 : L : N+(M-1)*L+3);
AFB               = y(N+4 : L : N+(M-1)*L+4);
AF                = y(N+5 : L : N+(M-1)*L+5);
ARBB              = y(N+6 : L : N+(M-1)*L+6);
ARB               = y(N+7 : L : N+(M-1)*L+7);
ALB               = y(N+8 : L : N+(M-1)*L+8);
AL                = y(N+9 : L : N+(M-1)*L+9);
A_Dioxin_AHR      = y(N+10 : L : N+(M-1)*L+10);
A_Dioxin_CYP1A2   = y(N+11 : L : N+(M-1)*L+11);

%Non-repeating variables
CYP1A2_mRNA       = y(N+M*L);   #Concentration of CYP1A2 mRNA
CYP1A2            = y(N+1+M*L); #Concentration of free CYP1A2


%% --------------------------- Derived Parameters ------------------------ %%
day         = t/24;
year        = day/365;
ageInYears  = g.init_age + g.growswitch*year;

BWmale = g.male*(g.RN1*0.00058*ageInYears.^3 - g.RN2*0.0948*ageInYears.^2 + g.RN3*4.8434*ageInYears + g.RN4*2.2785);      %male body weight in kg as a function of ageInYears
BWfemale = (1 - g.male)* (g.RN1*0.0006*ageInYears.^3 - g.RN2*0.0912*ageInYears.^2 + g.RN3*4.3200*ageInYears + g.RN4*3.6520);  %female body weight in kg as a function of ageInYears
BW        = BWmale + BWfemale;
BWgram    = BW*1.0e3;    %body weight in grams

HEIGHT      = -2e-5*ageInYears^4+4.2e-3*ageInYears^3.0-0.315*ageInYears^2.0+9.7465*ageInYears+72.098;   % Height in cm
HEIGHTmeter = (HEIGHT/100.0);   %Height in meter
BMI         = BW/(HEIGHTmeter^2.0);    %body mass index

MSTOT_NM  = g.MSTOT./g.MW;            %AMOUNT IN NMOL/KG
MSTT		  = MSTOT_NM * BW;     %AMOUNT IN NMOL
ORAL_DOSE_RATE = MSTT/g.ORAL_DURATION;


VLC = (3.59e-2 -(4.76e-7*BWgram)+(8.50e-12*BWgram^2.0)-(5.45e-17*BWgram^3.0));        %Liver volume fraction of BW (Wang et al (1997)), luecke paper (2007)
VFCmale = g.male*(-5.26e-20*BWgram^4.0 +1.09e-14*BWgram^3.0 -6.99e-10*BWgram^2.0 +1.59e-5*BWgram+3.95e-2);   %male fat volume fraction of BW
VFCfemale = (1-g.male)*(-6.36e-20*BWgram^4.0 +1.12e-14*BWgram^3.0 -5.8e-10*BWgram^2.0 +1.2e-5*BWgram+5.91e-2);   %female fat volume fraction of BW
VFC = VFCfemale + VFCmale;    %fat fraction of BW


g.VRBC = (0.91 - (g.VLBC*VLC + g.VFBC*VFC + VLC + VFC))/(1+g.VRBBC);      %NO BODY BLOOD VOLUME CONSIDERED HERE


g.kfAHR = g.kbAHR./g.KDAHR;    %Association rate constant for AHR binding (1/nM/hour)
g.kf1A2 = g.kb1A2./g.KD1A2;    %Association rate constant for 1A2 binding (1/nM/hour)
g.b = 1.0-g.a;         %portal fraction absoption of dioxin in liver
VL = VLC * BW;   %Liver volume
VF = VFC * BW;     %Fat volume
VRB = g.VRBC * BW;   %RB volume
VLB = g.VLBC * VL; %Liver blood volume
VFB = g.VFBC * VF;    %Fat blood volume
VRBB = g.VRBBC * VRB; %RB blood volume
QC = g.QCC*(BW)^0.75;   %Cardiac output (L/HR)
QF = g.QFC*QC;    %Fat blood flow rate (L/HR)
QL = g.QLC*QC;    %Liver blood flow rate (L/HR)
QRB = g.QRBC*QC;    %RB blood flow rate (L/HR)
QTTQ = QF+QRB+QL; %Flow balance check
PAF = g.PAFC.*QF;      %ADIPOSE (L/HR)
PARB = g.PARBC.*QRB;   %REST OF THE BODY (L/HR)
PAL = g.PALC.*QL;   %LIVER TISSUE (L/HR)


%SYSTEMIC Blood compartment========
CB = (QF*(AFB/VFB) + QRB*(ARBB/VRBB) + QL*(ALB/VLB) + g.KABS.*MST.*g.a) ./ (QC+g.CLURI);     %By only dealing with steady state situation assuming fast equilirium, the blood volume is not needed here. It only matters when the amount in the blood needs to be known and using ODE to describe the change. 
CA = CB;



%% ------------------------------ ODEs-------------------------------------%%

dydt = zeros(length(y),1); %make dydt as a column vector as required by MatLab ode function

%1
%AMOUNT IN THE GI TRACK
%MST
dydt(N : L : N+(M-1)*L) = - (g.KST+g.KABS).*MST + ORAL_DOSE_RATE;	 

%2
%AMOUNT ABSORBED INTO LYMPH CIRCULATION
%LYMLUM
dydt(N+1 : L : N+(M-1)*L+1) = g.KABS.*MST.*g.a;

%3
%AMOUNT ABSORBED INTO PORTAL CIRCULATION
%LIMLUM
dydt(N+2 : L : N+(M-1)*L+2) = g.KABS.*MST.*g.b;

%4
%URINARY EXCRETION BY KIDNEY
%AURI
dydt(N+3 : L : N+(M-1)*L+3) = g.CLURI .* CB;



%=========Adipose tissue compartment========
%5
%AFB
dydt(N+4 : L : N+(M-1)*L+4) = QF*(CA - AFB/VFB) - PAF.*(AFB/VFB.*g.Fu - (AF/VF)./g.PF.*g.Fu); 		%(nmol/hr)

%6
%AF
dydt(N+5 : L : N+(M-1)*L+5) = PAF.*(AFB/VFB.*g.Fu - (AF/VF)./g.PF.*g.Fu);  %(nmol/hr)



%=============Rest of the body compartment========
%7
%ARBB
dydt(N+6 : L : N+(M-1)*L+6)= QRB*(CA-ARBB/VRBB)  - PARB.*(ARBB/VRBB.*g.Fu  - ARB/VRB./g.PRB.*g.Fu);   %(nmol/hr)

%8
%ARB
dydt(N+7 : L : N+(M-1)*L+7) = PARB.*(ARBB/VRBB.*g.Fu - ARB/VRB./g.PRB.*g.Fu);    %(nmol/hr)



%=================Liver compartment===================

%Amount to concentration conversion
CFLL = AL/VL./g.PL.*g.Fu;
Dioxin_CYP1A2 = A_Dioxin_CYP1A2/VL;
Dioxin_AHR = A_Dioxin_AHR/VL;

%9
%ALB (Amount in liver blood)
dydt(N+8 : L : N+(M-1)*L+8) = QL*(CA-ALB/VLB) - PAL.*(ALB/VLB.*g.Fu - CFLL) + g.KABS.*MST.*g.b;    %(nmol/hr);

%10
%AL (Amount of free PCDD and nonspecifc bound PCDD in liver tissue proper)
dydt(N+9 : L : N+(M-1)*L+9) = PAL.*(ALB/VLB.*g.Fu - CFLL) + (-g.kelim.*CFLL.*(CYP1A2 + sum(Dioxin_CYP1A2)-g.CYP1A2_1BASAL)/g.CYP1A2_1BASAL - g.kfAHR .* CFLL * (g.AHRtot - sum(Dioxin_AHR)) + g.kbAHR .* Dioxin_AHR - g.kf1A2 .* CFLL * CYP1A2 + g.kb1A2 .* Dioxin_CYP1A2) * VL;
%dydt(N+9 : L : N+(M-1)*L+9) = PAL.*(ALB/VLB.*g.Fu - CFLL) + (-g.kelim.*CFLL.*(CYP1A2 + sum(Dioxin_CYP1A2)-g.CYP1A2_1BASAL)/g.CYP1A2_1BASAL - g.kfAHR .* CFLL * (g.AHRtot - sum(Dioxin_AHR)) + g.kbAHR .* Dioxin_AHR - g.kf1A2 .* CFLL * CYP1A2 + g.kb1A2 .* Dioxin_CYP1A2 + g.kdegCYP1A2 * Dioxin_CYP1A2) * VL; %For version where CYP1A2 within Dioxin_CYP1A2 can be degraded and dioxin recycled

%11
%A_Dioxin_AHR (Amount of Dioxin_AHR complex)
dydt(N+10 : L : N+(M-1)*L+10) = (g.kfAHR .* CFLL * (g.AHRtot - sum(Dioxin_AHR))  - g.kbAHR .* Dioxin_AHR)*VL;

%12
%A_Dioxin_CYP1A2 (Amount of Dioxin_CYP1A2) Note: this step requires using amount not concentration, otherwise it creates discrepancy from the BM model
dydt(N+11 : L : N+(M-1)*L+11) = (g.kf1A2 .* CFLL * CYP1A2 - g.kb1A2 .* Dioxin_CYP1A2)*VL;
%dydt(N+11 : L : N+(M-1)*L+11) = (g.kf1A2 .* CFLL * CYP1A2 - g.kb1A2 .* Dioxin_CYP1A2 - g.kdegCYP1A2 * Dioxin_CYP1A2)*VL; %For version where CYP1A2 within Dioxin_CYP1A2 can be degraded and dioxin recycled

%----------------Non-repeating variables-------------------%
%13
%CYP1A2_mRNA (Concentration of CYP1A2 mRNA)
dydt(N+M*L) = g.ktranscription_1A2 * (1 + g.CYP1A2_1EMAX * sum(Dioxin_AHR)^0.6 / (g.CYP1A2_1EC50^0.6 + sum(Dioxin_AHR)^0.6)) - g.kdegCYP1A2_mRNA * CYP1A2_mRNA;

%14
%CYP1A2 (Concentration of free CYP1A2)
dydt(N+1+M*L) = g.ktranslation_1A2 * CYP1A2_mRNA - g.kdegCYP1A2 * CYP1A2 - sum(g.kf1A2 .* CFLL * CYP1A2) + sum(g.kb1A2 .* Dioxin_CYP1A2);




