clear all
clc
close all

tic


%Unit: time = hour, chemical mass = ng (nmol), body weight = kg, volume = L 

pkg load io
pkg load statistics
pkg load optim

global g;

% =======================USER INPUT: SELECT GENDER, CONGENERS, EXPOSURE, AND SIMULATION PARAMETERS================ %
%Gender
g.male = 0;           %Set it to 1 for male; set it to 0 for female

%Select specific congeners from the congener listed in "DLC_human_parameters.xlsx" for mixture simulations
selected_congeners = [1:11]; %{1,5]. In the form of [a, b , c...] or [a:b, c...]
num_compounds = length(selected_congeners); %number of selected compounds

%Exposure parameters
g.ORAL_DURATION = 24;      %Oral exposure duration per day (hours)
g.MSTOT = []; %reset
g.MSTOT(1:num_compounds,1) = 0.02;  %%low dose: 7E-4, high dose: 0.02; 		%Oral exposure dose (ng/kg bw/day)
%exposure_frequency = 24*4; %in hours

%Simulation parameters
simulation_hours = 24*365*70;
output_interval_years = 1/2;
output_interval_days = 365/2;
output_interval_hours = 24*365/2;
tspan = [0:output_interval_hours:simulation_hours];
day = tspan/24;
week = day/7;
year = day/365;


% =======================PHYSIOLOGICAL PARAMETERS================ %
g.init_age = 0;       %age (years) at beginning of simulation
g.growswitch = 1.0;   %1 for growth, 0 for no growth

%CARDIAC OUTPUT CONSTANT USED IN (QC) EQUATION
g.QCC = 15.36;    %CONSTANT (L/HOUR/KG)

%TISSUE BLOOD FLOW FRACTION OF CARDIAC OUTPUT
g.QLC = 0.26;            %LIVER (Wang et al (1997))
g.QFC = 0.05;            %ADIPOSE (Wang et al (1997))
g.QRBC = 1 - g.QFC - g.QLC;

%FRACTION OF TISSUE BLOOD VOLUME
g.VLBC = 0.266;    %LIVER Wang et al (1997)
g.VFBC = 0.050;     %ADIPOSE TISSUE Wang et al (1997)
g.VRBBC = 0.030;    %REST OF THE BODY Wang et al (1997)

%AHR and CYP1A2 related parameters
g.AHRtot = 0.35;
g.CYP1A2_mRNA_1BASAL = 160;
g.CYP1A2_1BASAL = 1600;
g.CYP1A2_1EC50 = 130;
g.CYP1A2_1EMAX = 9300;
g.ktranscription_1A2 = 16; %Transcription rate constant
g.kdegCYP1A2_mRNA = 0.1;  %0.1 is the value for CYP1A2 degradation before the two-step delay, claimed as "FIRST ORDER RATE CONSTANT OF DEGRADATION (1/hour) Wang et al (1997)"
g.ktranslation_1A2 = 1.25;  %Translation rate constant
g.kdegCYP1A2 = 0.125;     %0.125 is half of the value 0.25 used in Emond for the two-step time delay.

parameters_default = g;


%-------Growth curve--------%
ageInYears = g.init_age + g.growswitch*year';    %time in term of year

%Body weight
BWmale = g.male*(0.00058*ageInYears.^3 - 0.0948*ageInYears.^2 + 4.8434*ageInYears + 2.2785);      %male body weight in kg as a function of ageInYears
BWfemale = (1 - g.male)* (0.0006*ageInYears.^3 - 0.0912*ageInYears.^2 + 4.3200*ageInYears + 3.6520);  %female body weight in kg as a function of ageInYears
BW =  BWmale + BWfemale;
BWgram = BW*1.0e3;    %body weight in grams

%Liver volume fraction
VLC = (3.59e-2 -(4.76e-7*BWgram)+(8.50e-12*BWgram.^2.0)-(5.45e-17*BWgram.^3.0));        %Liver fraction of BW (Wang et al (1997)), luecke paper (2007)

%Fat volume fraction
VFCmale = g.male*(-5.26e-20*BWgram.^4.0 +1.09e-14*BWgram.^3.0 -6.99e-10*BWgram.^2.0 +1.59e-5*BWgram+3.95e-2);   %male fat fraction of BW
VFCfemale = (1-g.male)*(-6.36e-20*BWgram.^4.0 +1.12e-14*BWgram.^3.0 -5.8e-10*BWgram.^2.0 +1.2e-5*BWgram+5.91e-2);   %female fat fraction of BW
VFC = VFCfemale + VFCmale;    %fat fraction of BW

%Rest of Body (RB) fraction
VRBC = (0.91 - (g.VLBC*VLC + g.VFBC*VFC + VLC + VFC))/(1+g.VRBBC);      %CHECK LATER SHOULD THERE BE BODY BLOOD VOLUME

VL = VLC .* BW;       %Liver volume
VF = VFC .* BW;       %Fat volume
VRB = VRBC .* BW;     %RB volume

VLB = g.VLBC * VL;    %Liver blood volume
VFB = g.VFBC * VF;    %Fat blood volume
VRBB = g.VRBBC * VRB; %RB blood volume

QC = g.QCC*(BW).^0.75;   %Cardiac output (L/HR)
QL = g.QLC*QC;    %Liver blood flow rate (L/HR)
QF = g.QFC*QC;    %Fat blood flow rate (L/HR)
QRB = g.QRBC*QC;    %RB blood flow rate (L/HR)



% =======================LOAD CHEMICAL SPECIFIC PARAMETERS================ %
[num,txt] = xlsread('DLC_human_parameters.xlsx');

param_values = num;
compound_names = txt(2:end,1);
param_names = txt(1,2:end);
num_param = length(param_names);

%assign parameter values to local parameter names respectively
for j=1:1:num_param
  %Select specific congeners for mixture simulations
  temp = param_values(selected_congeners,j);
  eval(strcat('g.',param_names{j}, '=temp;'));  %param_name{j} returns string; param_name(j) returns cell
end

normalized_exposure_TEQ = sum(g.MSTOT.*g.TEF)/mean(g.MSTOT)/num_compounds;

% =======================OTHER PARAMETERS================ %
global L %number of repeating variables
global M %number of compounds in the mixture
global N %starting index position of repeating variables

L = 12; %number of repeating variables
M = num_compounds; %number of compounds in the mixture
N = 1; %starting index position of repeating variables
% =======================COLOR SCHEME FOR PLOTTING================ %
color_scheme = [
                [0, 0.4470, 0.7410];
                [0.8500, 0.3250, 0.0980];	          	
                [0.9290, 0.6940, 0.1250];
                [0.4940, 0.1840, 0.5560];	          	
                [0.4660, 0.6740, 0.1880];	          	
                [0.3010, 0.7450, 0.9330];	          	
                [0.6350, 0.0780, 0.1840];
                [0, 0.5, 0];
                [0, 0.75, 0.75];
                [0.75, 0, 0.75];
                [0.75, 0.75, 0]
                ];
color_scheme = color_scheme';               

%Color preview
figure(222)
for i=1:11
  line([1,10],[i,i], 'color', color_scheme(:,i))
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------Simulation of mixture of arbitrary number of compounds---------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ======================= INITIAL CONDITION ======================= %
init = struct();
init.MST = zeros(M, 1);
init.LYMLUM = zeros(M, 1);
init.LIMLUM = zeros(M, 1);
init.AURI = zeros(M, 1);
init.AFB = zeros(M, 1);
init.AF = zeros(M, 1);
init.ARBB = zeros(M, 1);
init.ARB = zeros(M, 1);
init.ALB = zeros(M, 1);
init.AL = zeros(M, 1);
init.A_Dioxin_AHR = zeros(M, 1);
init.A_Dioxin_CYP1A2 = zeros(M, 1);
init.CYP1A2_mRNA = g.CYP1A2_mRNA_1BASAL; %initial CYP1A2_mRNA concentration
init.CYP1A2 = g.CYP1A2_1BASAL;     %initial CYP1A2 concentration

y0 = [];
for i = 1 : 1 : M
  y0(N+(i-1)*L) = init.MST(i);
  
  y0(N+1+(i-1)*L) = init.LYMLUM(i);
  
  y0(N+2+(i-1)*L) = init.LIMLUM(i);
  
  y0(N+3+(i-1)*L) = init.AURI(i);
  
  y0(N+4+(i-1)*L) = init.AFB(i);
  
  y0(N+5+(i-1)*L) = init.AF(i);
  
  y0(N+6+(i-1)*L) = init.ARBB(i);
  
  y0(N+7+(i-1)*L) = init.ARB(i);
  
  y0(N+8+(i-1)*L) = init.ALB(i);
  
  y0(N+9+(i-1)*L) = init.AL(i);
  
  y0(N+10+(i-1)*L) = init.A_Dioxin_AHR(i);
  
  y0(N+11+(i-1)*L) = init.A_Dioxin_CYP1A2(i);
end

%non-repeating variables
  y0(N+M*L) = init.CYP1A2_mRNA;
  
  y0(N+1+M*L) = init.CYP1A2;
  


% ======================= RUN SIMULATION ======================= %
[y, istate, msg] = lsode("DLC_human_ode", y0, tspan);


% ======================= Display simulation result  ======================= %
for i = 1 : 1 : M
      %plotting compound amount(nmol) in MST
      V = 1;
      figure(V);
      plot(year, y(:,(N-1)+V+(i-1)*L), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      hold on
      xlabel("year");
      ylabel("Amount (nmol)");
      title("MST");
      %legend(txt{2:12,1})

      %plotting compound concentration(nM) in Fat tissue
      V = 6;
      figure(V);
      CF(:,i) = y(:,(N-1)+V+(i-1)*L)./VF;
      plot(year, CF(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      pbaspect([1.25,1,1])
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("Fat");

      %plotting compound concentration(nM) in RB tissue
      V = 8;
      figure(V);
      CRB(:,i) = y(:,(N-1)+V+(i-1)*L)./VRB;
      plot(year, CRB(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      pbaspect([1.25,1,1])
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("Rest of Body");  
      
      %plotting compound concentration (nM) in Liver blood
      V = 9;
      figure(V);
      plot(year, y(:,(N-1)+V+(i-1)*L)./VLB, '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("Liver blood");

      %plotting free DLC + nonspecific bound DLC concentration in liver tissue proper
      V = 10;
      figure(V);
      plot(year, y(:,(N-1)+V+(i-1)*L)./VL, '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("Liver free and nonspecific bound");
      
      %plotting free DLC in liver tissue proper
      V = 10;
      figure(V*10);
      CL_free(:,i) = y(:,(N-1)+V+(i-1)*L)./VL/g.PL(i);
      plot(year, CL_free(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("Liver free PCDD");
         
      %plotting Dioxin_AHR complex concentration (nM) in liver tissue proper
      V = 11;
      figure(V);
      plot(year, y(:,(N-1)+V+(i-1)*L)./VL, '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("Dioxin-AhR");
          
      %plotting Dioxin_CYP1A2 in liver tissue proper
      V = 12;
      figure(V);
      plot(year, y(:,(N-1)+V+(i-1)*L)./VL, '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("Dioxin-CYP1A2");
          
      %plotting total compound concentration (nM) in liver tissue proper
      V1 = 10;  %AL (Amount of free PCDD and nonspecific bound PCDD in liver tissue proper)
      V2 = 11;  %A_Dioxin_AHR (Amount of Dioxin_AHR complex)
      V3 = 12;  %A_Dioxin_CYP1A2 (Amount of Dioxin_CYP1A2)
      CL_total(:,i) = (y(:,(N-1)+V1+(i-1)*L) + y(:,(N-1)+V2+(i-1)*L) + y(:,(N-1)+V3+(i-1)*L))./VL;  
      figure(V1+V2+V3);
      plot(year, CL_total(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      pbaspect([1.25,1,1])
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("Liver Total");
      
      %plotting CA
      CA(:,i) = (QF.*(y(:, N+4+(i-1)*L)./VFB) + QRB.*(y(:, N+6+(i-1)*L)./VRBB) + QL.*(y(:, N+8+(i-1)*L)./VLB) + y(:,N+(i-1)*L).*g.KABS(i)'.*g.a(i)) ./ (QC+g.CLURI(i)); 
      figure(1000);
      plot(year, CA(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 2);
      set(gca, 'fontsize', 20)
      pbaspect([1.25,1,1])
      hold on
      xlabel("year");
      ylabel("Concentration (nM)");
      title("CA");
            
endfor

figure(101);
plot(year, BW, '-', 'color', color_scheme(:,1), 'linewidth', 2);
hold on
set(gca, 'fontsize', 20)
xlabel("year");
ylabel("Weight (Kg)");
title("BW");


%plotting CYP1A2 mRNA concentration
V = 13;
figure(V);
plot(year, y(:,(N+M*L)), '-', 'color', color_scheme(:,1), 'linewidth', 2);
hold on
set(gca, 'fontsize', 20)
xlabel("year");
ylabel("Concentration (nM)");
title("CYP1A2 mRNA");


%plotting free CYP1A2 concentration
V = 14;
figure(V);
plot(year, y(:,(N+1+M*L)), '-', 'color', color_scheme(:,1), 'linewidth', 2);
hold on
set(gca, 'fontsize', 20)
xlabel("year");
ylabel("Concentration (nM)");
title("free CYP1A2");


%plotting total CYP1A2 concentration
V1 = 12;  %Dioxin_CYP1A2 amount
V2 = 14;  %free CYP1A2 concentration
figure(V1+V2);
CYP1A2_total = sum(y(:,N+11 : L : N+(M-1)*L+11)./VL,2) + y(:,(N+1+M*L));
plot(year, CYP1A2_total, '-', 'color', color_scheme(:,1), 'linewidth', 2);
hold on
set(gca, 'fontsize', 20)
xlabel("year");
ylabel("Concentration (nM)");
title("total CYP1A2");

mixture_mode_data_output_combined = [year', CA, CF, CRB, CL_free, CL_total, CYP1A2_total];

mixture_mode_data_output.CA = CA;
mixture_mode_data_output.CRB = CRB;
mixture_mode_data_output.CF = CF;
mixture_mode_data_output.CL_total = CL_total;


%Calculate mixture-mode tissue concentration TEQ (ng/L)
mixture_mode_TEQ.CA = CA*(g.Systemic_TEF.*g.MW);
mixture_mode_TEQ.CF = CF*(g.Systemic_TEF.*g.MW);
mixture_mode_TEQ.CRB = CRB*(g.Systemic_TEF.*g.MW);
mixture_mode_TEQ.CL_total = CL_total*(g.Systemic_TEF.*g.MW);

figure(1001)
plot(year, mixture_mode_TEQ.CA, '-', 'color', color_scheme(:,1), 'linewidth', 2);
hold on
set(gca, 'fontsize', 20)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("Concentration (ng/L)");
title("Mixture mode CA TEQ");

figure(1006)
plot(year, mixture_mode_TEQ.CF, '-', 'color', color_scheme(:,1), 'linewidth', 2);
hold on
set(gca, 'fontsize', 20)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("Concentration (ng/L)");
title("Mixture mode CF TEQ");

figure(1008)
plot(year, mixture_mode_TEQ.CRB, '-', 'color', color_scheme(:,1), 'linewidth', 2);
hold on
set(gca, 'fontsize', 20)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("Concentration (ng/L)");
title("Mixture mode CRB TEQ");

figure(1033)
plot(year, mixture_mode_TEQ.CL_total, '-', 'color', color_scheme(:,1), 'linewidth', 2);
hold on
set(gca, 'fontsize', 20)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("Concentration (ng/L)");
title("Mixture mode CL_total TEQ");
toc
%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------Single-compound simultions----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.MSTOT = []; %reset
g.MSTOT = [0.02]; %Low dose: 7e-4, high dose: 0.02
num_compounds = 11;
M = 1; %number of compounds in the mixture, needed here for the ode file

for p = 1:num_compounds
    p
 %if p==1 || p==5 %For specific single compounds only 
    %assign parameter values to local parameter names respectively
    for k=1:1:num_param
      temp = param_values(p,k);
      eval(strcat('g.',param_names{k}, '=temp;'));  %param_names{k} returns string; param_names(k) returns cell
    end
    
     
    % ======================= INITIAL CONDITION ======================= %
    init = struct();
    init.MST = zeros(1, 1);
    init.LYMLUM = zeros(1, 1);
    init.LIMLUM = zeros(1, 1);
    init.AURI = zeros(1, 1);
    init.AFB = zeros(1, 1);
    init.AF = zeros(1, 1);
    init.ARBB = zeros(1, 1);
    init.ARB = zeros(1, 1);
    init.ALB = zeros(1, 1);
    init.AL = zeros(1, 1);
    init.A_Dioxin_AHR = zeros(1, 1);
    init.A_Dioxin_CYP1A2 = zeros(1, 1);
    init.CYP1A2_mRNA = g.CYP1A2_mRNA_1BASAL; %initial CYP1A2_mRNA concentration
    init.CYP1A2 = g.CYP1A2_1BASAL;     %initial CYP1A2 concentration

    y0 = [];
    y0(1) = init.MST(1);
    y0(2) = init.LYMLUM(1);
    y0(3) = init.LIMLUM(1);
    y0(4) = init.AURI(1);
    y0(5) = init.AFB(1);
    y0(6) = init.AF(1);
    y0(7) = init.ARBB(1);
    y0(8) = init.ARB(1);
    y0(9) = init.ALB(1);
    y0(10) = init.AL(1);
    y0(11) = init.A_Dioxin_AHR(1);
    y0(12) = init.A_Dioxin_CYP1A2(1);
    y0(13) = init.CYP1A2_mRNA;
    y0(14) = init.CYP1A2;
      
    % ======================= RUN SIMULATION ======================= %
    [y, istate, msg] = lsode("DLC_human_ode", y0, tspan); %Note: use @TCDD_Emond_ode3 does not work for passing argument

    single_compound_simultion_output.CA(:,p) = (QF.*(y(:,5)./VFB) + QRB.*(y(:,7)./VRBB) + QL.*(y(:,9)./VLB) + y(:,1).*g.KABS.*g.a) ./ (QC+g.CLURI);
    single_compound_simultion_output.CF(:,p) = y(:,6)./VF;
    single_compound_simultion_output.CRB(:,p) = y(:,8)./VRB;
    single_compound_simultion_output.CL_total(:,p) = (y(:,10)+y(:,11)+y(:,12))./VL;
    
    %TEQ calculation
    single_mode_TEQ.CA(:,p) = single_compound_simultion_output.CA(:,p)*(g.Systemic_TEF.*g.MW);
    single_mode_TEQ.CF(:,p) = single_compound_simultion_output.CF(:,p)*(g.Systemic_TEF.*g.MW);
    single_mode_TEQ.CRB(:,p) = single_compound_simultion_output.CRB(:,p)*(g.Systemic_TEF.*g.MW);
    single_mode_TEQ.CL_total(:,p) = single_compound_simultion_output.CL_total(:,p)*(g.Systemic_TEF.*g.MW);
    
    %plotting CA (nM)
    figure(1000)
    plot(year, single_compound_simultion_output.CA(:,p), '-', 'color', color_scheme(:,p), 'linewidth', 2, 'linestyle','--');
    set(gca, 'fontsize', 20)
    pbaspect([1.25,1,1])
    hold on
    xlabel("year");
    ylabel("Concentration (nM)");
    title("CA");
    
    %plotting compound concentration(nM) in Fat tissue
    V = 6;
    figure(V);
    plot(year, single_compound_simultion_output.CF(:,p), '-', 'color', color_scheme(:,p), 'linewidth', 2, 'linestyle','--');
    set(gca, 'fontsize', 20)
    pbaspect([1.25,1,1])
    hold on
    xlabel("year");
    ylabel("Concentration (nM)");
    title("Fat");
    
    %plotting compound concentration(nM) in RB tissue
    V = 8;
    figure(V);
    plot(year, single_compound_simultion_output.CRB(:,p), '-', 'color', color_scheme(:,p), 'linewidth', 2, 'linestyle','--');
    set(gca, 'fontsize', 20)
    pbaspect([1.25,1,1])
    hold on
    xlabel("year");
    ylabel("Concentration (nM)");
    title("RB");
    
    %plotting compound concentration(nM) in Liver tissue
    V = 33;
    figure(V);
    plot(year, single_compound_simultion_output.CL_total(:,p), '-', 'color', color_scheme(:,p), 'linewidth', 2, 'linestyle','--');
    set(gca, 'fontsize', 20)
    pbaspect([1.25,1,1])
    hold on
    xlabel("year");
    ylabel("Concentration (nM)");
    title("Liver total");
  %endif
end



%Calculate single-mode tissue concentration TEQ (ng/L)
single_mode_TEQ.CA   = sum(single_mode_TEQ.CA,2);
single_mode_TEQ.CF   = sum(single_mode_TEQ.CF,2);
single_mode_TEQ.CRB  = sum(single_mode_TEQ.CRB,2);
single_mode_TEQ.CL_total = sum(single_mode_TEQ.CL_total,2);

figure(1001)
plot(year, single_mode_TEQ.CA, '-', 'color', color_scheme(:,1), 'linewidth', 2, 'linestyle','--');
hold on
set(gca, 'fontsize', 20)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("Concentration (ng/L)");
title("Single mode CA TEQ");

figure(1006)
plot(year, single_mode_TEQ.CF, '-', 'color', color_scheme(:,1), 'linewidth', 2, 'linestyle','--');
hold on
set(gca, 'fontsize', 20)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("Concentration (ng/L)");
title("Single mode CF TEQ");

figure(1008)
plot(year, single_mode_TEQ.CRB, '-', 'color', color_scheme(:,1), 'linewidth', 2, 'linestyle','--');
hold on
set(gca, 'fontsize', 20)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("Concentration (ng/L)");
title("Single mode CRB TEQ");

figure(1033)
plot(year, single_mode_TEQ.CL_total, '-', 'color', color_scheme(:,1), 'linewidth', 2, 'linestyle','--');
hold on
set(gca, 'fontsize', 20)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("Concentration (ng/L)");
title("Single mode CL total TEQ");

%Calculate mixture:single congener concentration ratio

ratio_CA = mixture_mode_data_output.CA./single_compound_simultion_output.CA;
ratio_CRB = mixture_mode_data_output.CRB./single_compound_simultion_output.CRB;
ratio_CF = mixture_mode_data_output.CF./single_compound_simultion_output.CF;
ratio_CL_total = mixture_mode_data_output.CL_total./single_compound_simultion_output.CL_total;

min_CA = min(ratio_CA);
min_CRB = min(ratio_CRB);
min_CF = min(ratio_CF);
min_CL_total = min(ratio_CL_total);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------Pair-wise mixture simultions----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.MSTOT = []; %reset
g.MSTOT = [0.02;0.02]; %Low dose: 7e-4, high dose: 0.02
num_compounds = 11;
M = 2; %number of compounds in the mixture, needed here for the ode file

square_matrix_output = [];

for p = 1:(num_compounds-1)
    p
  for q = (p+1):num_compounds
    q
   %if p==1 && q==5 %For specific binary combination    
    %Select specific congener pairs from the congener listed in "DLC_human_parameters.xlsx" for mixture simulations
    congener_pair = [p,q];
    %assign parameter values to local parameter names respectively
    for k=1:1:num_param
      temp = param_values(congener_pair,k);
      eval(strcat('g.',param_names{k}, '=temp;'));  %param_name{k} returns string; param_name(k) returns cell
    end
    
    
    % ======================= INITIAL CONDITION ======================= %
    init = struct();
    init.MST = zeros(2, 1);
    init.LYMLUM = zeros(2, 1);
    init.LIMLUM = zeros(2, 1);
    init.AURI = zeros(2, 1);
    init.AFB = zeros(2, 1);
    init.AF = zeros(2, 1);
    init.ARBB = zeros(2, 1);
    init.ARB = zeros(2, 1);
    init.ALB = zeros(2, 1);
    init.AL = zeros(2, 1);
    init.A_Dioxin_AHR = zeros(2, 1);
    init.A_Dioxin_CYP1A2 = zeros(2, 1);
    init.CYP1A2_mRNA = g.CYP1A2_mRNA_1BASAL; %initial CYP1A2_mRNA concentration
    init.CYP1A2 = g.CYP1A2_1BASAL;     %initial CYP1A2 concentration

    y0 = [];
    for i = 1 : 1 : 2
      y0(N+(i-1)*L) = init.MST(i);
      y0(N+1+(i-1)*L) = init.LYMLUM(i);
      y0(N+2+(i-1)*L) = init.LIMLUM(i);
      y0(N+3+(i-1)*L) = init.AURI(i);
      y0(N+4+(i-1)*L) = init.AFB(i);
      y0(N+5+(i-1)*L) = init.AF(i);
      y0(N+6+(i-1)*L) = init.ARBB(i);
      y0(N+7+(i-1)*L) = init.ARB(i);
      y0(N+8+(i-1)*L) = init.ALB(i);
      y0(N+9+(i-1)*L) = init.AL(i);
      y0(N+10+(i-1)*L) = init.A_Dioxin_AHR(i);
      y0(N+11+(i-1)*L) = init.A_Dioxin_CYP1A2(i);
    end

    %non-repeating variables
      y0(N+2*L) = init.CYP1A2_mRNA;
      y0(N+1+2*L) = init.CYP1A2;
      
    % ======================= RUN SIMULATION ======================= %
    [y, istate, msg] = lsode("DLC_human_ode", y0, tspan); %Note: use @TCDD_Emond_ode3 does not work for passing argument


    two_compounds_CA = [];
    two_compounds_CF = [];
    two_compounds_CRB = [];
    two_compounds_CL_total = [];
    for i = 1 : 1 : 2 
        %plotting CA (nM)
        figure(1000)
        two_compounds_CA(:,i) = (QF.*(y(:, N+4+(i-1)*L)./VFB) + QRB.*(y(:, N+6+(i-1)*L)./VRBB) + QL.*(y(:, N+8+(i-1)*L)./VLB) + y(:,N+(i-1)*L).*g.KABS(i)'.*g.a(i)) ./ (QC+g.CLURI(i)); 
        plot(year, two_compounds_CA(:,i), '-', 'color', color_scheme(:,i), 'linewidth', 2);
        set(gca, 'fontsize', 20)
        hold on
        xlabel("year");
        ylabel("Concentration (nM)");
        title("CA");        
        
        %plotting compound concentration(nM) in Fat tissue
        V = 6;
        figure(V);
        two_compounds_CF(:,i) = y(:,(N-1)+V+(i-1)*L)./VF;
        plot(year, two_compounds_CF(:,i), '-', 'color', color_scheme(:,i), 'linewidth', 2);
        set(gca, 'fontsize', 20)
        hold on
        xlabel("year");
        ylabel("Concentration (nM)");
        title("Fat");
        
        %plotting compound concentration(nM) in RB tissue
        V = 8;
        figure(V);
        two_compounds_CRB(:,i) = y(:,(N-1)+V+(i-1)*L)./VRB;
        plot(year, two_compounds_CRB(:,i), '-', 'color', color_scheme(:,i), 'linewidth', 2);
        set(gca, 'fontsize', 20)
        hold on
        xlabel("year");
        ylabel("Concentration (nM)");
        title("Rest of Body");

        %plotting total compound concentration (nM) in liver tissue proper
        V1 = 10;  %AL (Amount of free PCDD and nonspecific bound PCDD in liver tissue proper)
        V2 = 11;  %A_Dioxin_AHR (Amount of Dioxin_AHR complex)
        V3 = 12;  %A_Dioxin_CYP1A2 (Amount of Dioxin_CYP1A2)
        two_compounds_CL_total(:,i) = (y(:,(N-1)+V1+(i-1)*L) + y(:,(N-1)+V2+(i-1)*L) + y(:,(N-1)+V3+(i-1)*L))./VL;  
        figure(33);
        plot(year, two_compounds_CL_total(:,i), '-', 'color', color_scheme(:,i), 'linewidth', 2);
        set(gca, 'fontsize', 20)
        hold on
        xlabel("year");
        ylabel("Concentration (nM)");
        title("Liver Total");   
    end
    
    square_matrix_output.CA{p,q} = two_compounds_CA; 
    square_matrix_output.CF{p,q} = two_compounds_CF;
    square_matrix_output.CRB{p,q} = two_compounds_CRB;
    square_matrix_output.CL_total{p,q} = two_compounds_CL_total;
   %endif
  end
end




%Merge single-compound and two-compound results
merged_matrix_output = square_matrix_output;

for i = 1:1:num_compounds
  merged_matrix_output.CA{i,i} = single_compound_simultion_output.CA(:,i);
  merged_matrix_output.CF{i,i} = single_compound_simultion_output.CF(:,i);
  merged_matrix_output.CRB{i,i} = single_compound_simultion_output.CRB(:,i);
  merged_matrix_output.CL_total{i,i} = single_compound_simultion_output.CL_total(:,i);
end



%Calculate ratio of tissue concnetrations in pair-mixture vs. single compounds models
num_time_points = length(tspan);
ratio_matrix_output = {};
for p = 1:(num_compounds-1)
    p
  for q = (p+1):num_compounds
    q
    
    %CA
    two_compounds_CA_p  = cell2mat(merged_matrix_output.CA(p,q))(:,1);
    two_compounds_CA_q  = cell2mat(merged_matrix_output.CA(p,q))(:,2);
    single_compound_CA_p   = cell2mat(merged_matrix_output.CA(p,p));
    single_compound_CA_q   = cell2mat(merged_matrix_output.CA(q,q));
    ratio_matrix_output.CA(p,q) = {[mean(two_compounds_CA_p(2:num_time_points)./single_compound_CA_p(2:num_time_points)), mean(two_compounds_CA_q(2:num_time_points)./single_compound_CA_q(2:num_time_points))]};
  
    %CF
    two_compounds_CF_p  = cell2mat(merged_matrix_output.CF(p,q))(:,1);
    two_compounds_CF_q  = cell2mat(merged_matrix_output.CF(p,q))(:,2);
    single_compound_CF_p   = cell2mat(merged_matrix_output.CF(p,p));
    single_compound_CF_q   = cell2mat(merged_matrix_output.CF(q,q));
    ratio_matrix_output.CF(p,q) = {[mean(two_compounds_CF_p(2:num_time_points)./single_compound_CF_p(2:num_time_points)), mean(two_compounds_CF_q(2:num_time_points)./single_compound_CF_q(2:num_time_points))]};
    
    %CRB
    two_compounds_CRB_p  = cell2mat(merged_matrix_output.CRB(p,q))(:,1);
    two_compounds_CRB_q  = cell2mat(merged_matrix_output.CRB(p,q))(:,2);
    single_compound_CRB_p   = cell2mat(merged_matrix_output.CRB(p,p));
    single_compound_CRB_q   = cell2mat(merged_matrix_output.CRB(q,q));
    ratio_matrix_output.CRB(p,q) = {[mean(two_compounds_CRB_p(2:num_time_points)./single_compound_CRB_p(2:num_time_points)), mean(two_compounds_CRB_q(2:num_time_points)./single_compound_CRB_q(2:num_time_points))]};
   
    %CL_total
    two_compounds_CL_total_p  = cell2mat(merged_matrix_output.CL_total(p,q))(:,1);
    two_compounds_CL_total_q  = cell2mat(merged_matrix_output.CL_total(p,q))(:,2);
    single_compound_CL_total_p   = cell2mat(merged_matrix_output.CL_total(p,p));
    single_compound_CL_total_q   = cell2mat(merged_matrix_output.CL_total(q,q));
    ratio_matrix_output.CL_total(p,q) = {[mean(two_compounds_CL_total_p(2:num_time_points)./single_compound_CL_total_p(2:num_time_points)), mean(two_compounds_CL_total_q(2:num_time_points)./single_compound_CL_total_q(2:num_time_points))]};
  end
end

%Adding 1 to diagonal
for i = 1:1:num_compounds
  ratio_matrix_output.CA{i,i} = [1,1];
  ratio_matrix_output.CF{i,i} = [1,1];
  ratio_matrix_output.CRB{i,i} = [1,1];
  ratio_matrix_output.CL_total{i,i} = [1,1];
end

%Combine CA and CRB matrix, and combine CF and CL_total matrix
ratio_matrix_output.CA_CRB = ratio_matrix_output.CA;
ratio_matrix_output.CF_CL_total = ratio_matrix_output.CF;
for p = 2:(num_compounds)
    p
  for q = 1:(p-1)
    q
    ratio_matrix_output.CA_CRB(p,q) = ratio_matrix_output.CRB'(p,q);
    ratio_matrix_output.CF_CL_total(p,q) = ratio_matrix_output.CL_total'(p,q);
  end
end

%Plot matrix result

matrix_CA_CRB       = cell2mat(ratio_matrix_output.CA_CRB);
matrix_CF_CL_total  = cell2mat(ratio_matrix_output.CF_CL_total);
num_row     = size(matrix_CA_CRB)(1);
num_column  = size(matrix_CA_CRB)(2);

lower_bound = min([min(matrix_CA_CRB), min(matrix_CF_CL_total)]);
upper_bound = max([max(matrix_CA_CRB), max(matrix_CF_CL_total)]);

%Plot CA and CRB
figure(11111)
%imagesc(matrix_CA_CRB, [lower_bound 1])
imagesc(matrix_CA_CRB, [0 1])
pbaspect([1,1,1])
axis("tic[]") %turn off tick
%colorbar('Location','eastoutside')
colorbar('Location','southoutside')
%Adding horizontal grid
x = [0.5:1:22.5];
for i=1.5:1:10.5
  y = i+0*x;
  line(x,y,'color','k')
end

%Adding vertical grid
y = [0.5:1:11.5];
for i=2.5:2:20.5
  x = i+0*y;
  line(x,y,'color','k')
end


%Plot CF and CL_total
figure(11113)
imagesc(matrix_CF_CL_total, [0 1.2])
pbaspect([1,1,1])
axis("tic[]") %turn off tick
colorbar('Location','southoutside')
%Adding horizontal grid
x = [0.5:1:22.5];
for i=1.5:1:10.5
  y = i+0*x;
  line(x,y,'color','k')
end

%Adding vertical grid
y = [0.5:1:11.5];
for i=2.5:2:20.5
  x = i+0*y;
  line(x,y,'color','k')
end

%
aaaa=[];
for p = 1:11
  for q=p:11
    aaaa=[aaaa,cell2mat(ratio_matrix_output.CL_total(p,q))];
  end
end
min(aaaa)
max(aaaa)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------Population Monte Carlo all-mixture simultions------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_individuals = 1000;

selected_congeners = [1:11]; %in the form of [a, b , c...] or [a:b, c...]
num_compounds = length(selected_congeners); %number of selected compounds
M = num_compounds; %number of compounds in the mixture, needed here for the ode file

%Dutch population data in 1994 (Liem et al., 2020)
[dose,txt1] = xlsread('DLC_human_dose.xlsx');
dose_mean = dose(:,1);
dose_CV = dose(:,2);

        
        
%assign parameter values to local parameter names respectively
for j=1:1:num_param
  %Select specific congeners for mixture simulations
  temp = param_values(selected_congeners,j);
  eval(strcat('g.',param_names{j}, '=temp;'));  %param_name{j} returns string; param_name(j) returns cell
end




BW_array = [];
AHRtot_array = [];
ktranscription_1A2_array = [];
CYP1A2_1EMAX_array = [];
dose_array = [];
population_mixture_CA_array = [];
population_mixture_CRB_array = [];
population_mixture_CF_array = [];
population_mixture_CL_total_array = [];
population_mixture_TEQ_CA_array = [];
population_mixture_TEQ_CRB_array = [];
population_mixture_TEQ_CF_array = [];
population_mixture_TEQ_CL_total_array = [];
population_mixture_TEQ_intake_per_kg_array = [];
population_mixture_TEQ_intake_array = [];

for n=1:num_individuals
    n
    % ======================= MC sampling ======================= %
    
      %-------DLC dose--------%
        CV = dose_CV; %coefficient of variance
        m = dose_mean; %mean
        v = (CV.*m).^2; %variance
        mu = log((m.^2)./sqrt(v+m.^2)); 
        sigma = sqrt(log(v./(m.^2)+1));
        g.MSTOT = []; %reset;
        g.MSTOT(1:num_compounds,1) = lognrnd(mu, sigma);
        
              
      %-------Growth curve--------%
        CV = 0.01; %coefficient of variance
        m = 1; %mean
        v = (CV*m)^2; %variance
        mu = log((m^2)/sqrt(v+m^2));
        sigma = sqrt(log(v/(m^2)+1));
        g.RN1 = lognrnd(mu, sigma);
        g.RN2 = lognrnd(mu, sigma);
        g.RN3 = lognrnd(mu, sigma);
        g.RN4 = lognrnd(mu, sigma);

        %Varying gender
        g.male = unidrnd(2)-1; %unidrnd(2) generates random number of either 1 or 2

        %Varying Body weight
        ageInYears = g.init_age + g.growswitch*year';    %time in term of year
        BWmale = g.male*(g.RN1*0.00058*ageInYears.^3 - g.RN2*0.0948*ageInYears.^2 + g.RN3*4.8434*ageInYears + g.RN4*2.2785);      %male body weight in kg as a function of ageInYears
        BWfemale = (1 - g.male)* (g.RN1*0.0006*ageInYears.^3 - g.RN2*0.0912*ageInYears.^2 + g.RN3*4.3200*ageInYears + g.RN4*3.6520);  %female body weight in kg as a function of ageInYears
        BW =  BWmale + BWfemale;
        BWgram = BW*1.0e3;    %body weight in grams

        %Liver volume fraction
        VLC = (3.59e-2 -(4.76e-7*BWgram)+(8.50e-12*BWgram.^2.0)-(5.45e-17*BWgram.^3.0));        %Liver fraction of BW (Wang et al (1997)), luecke paper (2007)

        %Fat volume fraction
        VFCmale = g.male*(-5.26e-20*BWgram.^4.0 +1.09e-14*BWgram.^3.0 -6.99e-10*BWgram.^2.0 +1.59e-5*BWgram+3.95e-2);   %male fat fraction of BW
        VFCfemale = (1-g.male)*(-6.36e-20*BWgram.^4.0 +1.12e-14*BWgram.^3.0 -5.8e-10*BWgram.^2.0 +1.2e-5*BWgram+5.91e-2);   %female fat fraction of BW
        VFC = VFCfemale + VFCmale;    %fat fraction of BW

        %Rest of Body (RB) fraction
        VRBC = (0.91 - (g.VLBC*VLC + g.VFBC*VFC + VLC + VFC))/(1+g.VRBBC);      %CHECK LATER SHOULD THERE BE BODY BLOOD VOLUME

        VL = VLC .* BW;       %Liver volume
        VF = VFC .* BW;       %Fat volume
        VRB = VRBC .* BW;     %RB volume

        VLB = g.VLBC * VL;    %Liver blood volume
        VFB = g.VFBC * VF;    %Fat blood volume
        VRBB = g.VRBBC * VRB; %RB blood volume

        QC = g.QCC*(BW).^0.75;   %Cardiac output (L/HR)
        QF = g.QFC*QC;    %Fat blood flow rate (L/HR)
        QL = g.QLC*QC;    %Liver blood flow rate (L/HR)
        QRB = g.QRBC*QC;    %RB blood flow rate (L/HR)
       
      %-------Varying biochemical parameters--------% 
        CV = 0.3; %coefficient of variance
        m = 1; %mean
        v = (CV*m)^2; %variance
        mu = log((m^2)/sqrt(v+m^2));
        sigma = sqrt(log(v/(m^2)+1));
        
        g.AHRtot = lognrnd(mu, sigma) * parameters_default.AHRtot;
        
        %The following 3 parameters are changed by the same scale
        RN = lognrnd(mu, sigma);
        g.CYP1A2_mRNA_1BASAL  = RN * parameters_default.CYP1A2_mRNA_1BASAL;
        g.CYP1A2_1BASAL       = RN * parameters_default.CYP1A2_1BASAL;
        g.ktranscription_1A2  = RN * parameters_default.ktranscription_1A2;
        
        
        g.CYP1A2_1EMAX = lognrnd(mu, sigma) * parameters_default.CYP1A2_1EMAX;
        
    
    % ======================= INITIAL CONDITION ======================= %
    init = struct();
    init.MST = zeros(M, 1);
    init.LYMLUM = zeros(M, 1);
    init.LIMLUM = zeros(M, 1);
    init.AURI = zeros(M, 1);
    init.AFB = zeros(M, 1);
    init.AF = zeros(M, 1);
    init.ARBB = zeros(M, 1);
    init.ARB = zeros(M, 1);
    init.ALB = zeros(M, 1);
    init.AL = zeros(M, 1);
    init.A_Dioxin_AHR = zeros(M, 1);
    init.A_Dioxin_CYP1A2 = zeros(M, 1);
    init.CYP1A2_mRNA = g.CYP1A2_mRNA_1BASAL; %initial CYP1A2_mRNA concentration
    init.CYP1A2 = g.CYP1A2_1BASAL;     %initial CYP1A2 concentration

    y0 = [];
    for i = 1 : 1 : M
      y0(N+(i-1)*L) = init.MST(i);
      
      y0(N+1+(i-1)*L) = init.LYMLUM(i);
      
      y0(N+2+(i-1)*L) = init.LIMLUM(i);
      
      y0(N+3+(i-1)*L) = init.AURI(i);
      
      y0(N+4+(i-1)*L) = init.AFB(i);
      
      y0(N+5+(i-1)*L) = init.AF(i);
      
      y0(N+6+(i-1)*L) = init.ARBB(i);
      
      y0(N+7+(i-1)*L) = init.ARB(i);
      
      y0(N+8+(i-1)*L) = init.ALB(i);
      
      y0(N+9+(i-1)*L) = init.AL(i);
      
      y0(N+10+(i-1)*L) = init.A_Dioxin_AHR(i);
      
      y0(N+11+(i-1)*L) = init.A_Dioxin_CYP1A2(i);
    end

    %non-repeating variables
      y0(N+M*L) = init.CYP1A2_mRNA;
      
      y0(N+1+M*L) = init.CYP1A2;
     
    % ======================= RUN SIMULATION ======================= %
      [y, istate, msg] = lsode("DLC_human_MC_ode", y0, tspan); 

    
    % ======================= Display simulation result  ======================= %  
      for i = 1 : 1 : M  
        
        %plotting CA
          CA(:,i) = (QF.*(y(:, N+4+(i-1)*L)./VFB) + QRB.*(y(:, N+6+(i-1)*L)./VRBB) + QL.*(y(:, N+8+(i-1)*L)./VLB) + y(:,N+(i-1)*L).*g.KABS(i)'.*g.a(i)) ./ (QC+g.CLURI(i)); 
##          figure(1000);
##          plot(year, CA(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 1);
##          set(gca, 'fontsize', 20)
##          hold on
##          xlabel("year");
##          ylabel("CA(nM)");
        
        %plotting CRB
          V = 8;
          CRB(:,i) = y(:,(N-1)+V+(i-1)*L)./VRB;
##          figure(V);
##          plot(year, CRB(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 1);
##          set(gca, 'fontsize', 20)
##          hold on
##          xlabel("year");
##          ylabel("CRB(nM)");
        
        
        %plotting CF
          V = 6;
          CF(:,i) = y(:,(N-1)+V+(i-1)*L)./VF;
##          figure(V);
##          plot(year, CF(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 1);
##          set(gca, 'fontsize', 20)
##          hold on
##          xlabel("year");
##          ylabel("CF (nM)");
        
        
        %plotting total compound concentration (nM) in liver tissue proper
          V1 = 10;  %AL (Amount of free PCDD and nonspecific bound PCDD in liver tissue proper)
          V2 = 11;  %A_Dioxin_AHR (Amount of Dioxin_AHR complex)
          V3 = 12;  %A_Dioxin_CYP1A2 (Amount of Dioxin_CYP1A2)
          CL_total(:,i) = (y(:,(N-1)+V1+(i-1)*L) + y(:,(N-1)+V2+(i-1)*L) + y(:,(N-1)+V3+(i-1)*L))./VL;  
##          figure(V1+V2+V3);
##          plot(year, CL_total(:,i), '-', 'color', color_scheme(:,selected_congeners(i)), 'linewidth', 1);
##          set(gca, 'fontsize', 20)
##          hold on
##          xlabel("year");
##          ylabel("CL total(nM)");
   
      end
      
      
      %Calculate mixture intake TEQ per kg bw(ng/kg/day TEQ)
        population_mixture_TEQ.intake_per_kg = sum(g.MSTOT .* g.TEF);
        
      %Calculate mixture intake TEQ (ng/day TEQ)
        population_mixture_TEQ.intake = BW * population_mixture_TEQ.intake_per_kg;
      
      
      
      %Calculate mixture tissue concentration TEQ (ng/L TEQ)
        population_mixture_TEQ.CA = CA * (g.MW .* g.Systemic_TEF);
        population_mixture_TEQ.CRB = CRB * (g.MW .* g.Systemic_TEF);
        population_mixture_TEQ.CF = CF * (g.MW .* g.Systemic_TEF);
        population_mixture_TEQ.CL_total = CL_total * (g.MW .* g.Systemic_TEF);

        figure(1001)
        plot(year, population_mixture_TEQ.CA, '-', 'color', color_scheme(:,1), 'linewidth', 1);
        hold on
        set(gca, 'fontsize', 20)
        xlabel("year");
        ylabel("CA TEQ (ng/L)");
        
        figure(1002)
        plot(year, population_mixture_TEQ.CRB, '-', 'color', color_scheme(:,1), 'linewidth', 1);
        hold on
        set(gca, 'fontsize', 20)
        xlabel("year");
        ylabel("CRB TEQ (ng/L)");
        
        figure(1003)
        plot(year, population_mixture_TEQ.CF, '-', 'color', color_scheme(:,1), 'linewidth', 1);
        hold on
        set(gca, 'fontsize', 20)
        xlabel("year");
        ylabel("CF TEQ (ng/L)");
        
        figure(1004)
        plot(year, population_mixture_TEQ.CL_total, '-', 'color', color_scheme(:,1), 'linewidth', 1);
        hold on
        set(gca, 'fontsize', 20)
        xlabel("year");
        ylabel("CL total TEQ (ng/L)");

              
      
      %Collect data
        BW_array = [BW_array, BW];
        AHRtot_array = [AHRtot_array, g.AHRtot];
        ktranscription_1A2_array = [ktranscription_1A2_array, g.ktranscription_1A2];
        CYP1A2_1EMAX_array = [CYP1A2_1EMAX_array, g.CYP1A2_1EMAX];
        dose_array = [dose_array, g.MSTOT];
        population_mixture_CA_array = [population_mixture_CA_array, CA];
        population_mixture_CRB_array = [population_mixture_CRB_array, CRB];
        population_mixture_CF_array = [population_mixture_CF_array, CF];
        population_mixture_CL_total_array = [population_mixture_CL_total_array, CL_total];
        
        population_mixture_TEQ_CA_array = [population_mixture_TEQ_CA_array, population_mixture_TEQ.CA];
        population_mixture_TEQ_CRB_array = [population_mixture_TEQ_CRB_array, population_mixture_TEQ.CRB];
        population_mixture_TEQ_CF_array = [population_mixture_TEQ_CF_array, population_mixture_TEQ.CF];
        population_mixture_TEQ_CL_total_array = [population_mixture_TEQ_CL_total_array, population_mixture_TEQ.CL_total];
        
        population_mixture_TEQ_intake_per_kg_array = [population_mixture_TEQ_intake_per_kg_array, population_mixture_TEQ.intake_per_kg];        
        population_mixture_TEQ_intake_array = [population_mixture_TEQ_intake_array, population_mixture_TEQ.intake];

end

figure(1)
plot(year, BW_array)
set(gca, 'fontsize', 24)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("BW (kg)");

figure(2)
hist(AHRtot_array, 'facecolor', [0.8500, 0.3250, 0.0980])
xlim([0,1.201])
set(gca, 'fontsize', 24)
pbaspect([1.25,1,1])
xlabel("AHRtot (nM)");
ylabel("Frequency");

figure(3)
hist(ktranscription_1A2_array, 'facecolor', [0.4940, 0.1840, 0.5560])
set(gca, 'fontsize', 24)
pbaspect([1.25,1,1])
xlabel("ktranscription_1A2 (nM/h)");
ylabel("Frequency");

figure(4)
hist(CYP1A2_1EMAX_array, 'facecolor', [0.4660, 0.6740, 0.1880])
set(gca, 'fontsize', 24)
pbaspect([1.25,1,1])
xlabel("CYP1A2_1EMAX");
ylabel("Frequency");


%Plot average dialy intake mixture dose distribution (ng/kg/day)
figure(11)
dose_array_mean = mean(dose_array, 2);
dose_array_mean_percentage=dose_array_mean/sum(dose_array_mean);
pie(dose_array_mean_percentage);
title("Mixture Daily Intake Dose Distribution");


%Plot average daily intake mixture TEQ distribution (ng/kg/day)
figure(12)
intake_TEQ = dose_array.*g.TEF;
intake_TEQ_mean = mean(intake_TEQ,2);
intake_TEQ_mean_percentage=intake_TEQ_mean/sum(intake_TEQ_mean);
pie(intake_TEQ_mean_percentage);
title("Mixture Daily Intake TEQ Distribution");


%Plot lifetime cumulative intake mixture dose distribution(ng)
figure(13)
dose_cumulative = [];
for n=1:num_individuals
  dose_inidivual_time_course = BW_array(:,n)*(dose_array(:,n))' * output_interval_days; %ng/half year
  dose_individual_cumulative = sum(dose_inidivual_time_course(1:(end-1),:), 1); %ng, 70-year time point not used
  dose_cumulative = [dose_cumulative, dose_individual_cumulative'];
end
dose_cumulative_mean = mean(dose_cumulative,2);
dose_cumulative_mean_percentage=dose_cumulative_mean/sum(dose_cumulative_mean);
pie(dose_cumulative_mean_percentage);
title("Mixture Cumulative Intake Dose Distribution");


%Plot lifetime cumulative intake mixture TEQ distribution (ng)
figure(14)
intake_TEQ_cumulative = dose_cumulative.*g.TEF;
intake_TEQ_cumulative_mean = mean(intake_TEQ_cumulative,2);
intake_TEQ_cumulative_mean_percentage=intake_TEQ_cumulative_mean/sum(intake_TEQ_cumulative_mean);
pie(intake_TEQ_cumulative_mean_percentage);
title("Mixture Cumulative Intake TEQ Distribution");
%legend(compound_names);


%Plot lifetime AUC CA TEQ distribution(ng/L*year)
figure(15)
CA_TEQ_AUC = [];
for n=1:num_individuals
  CA_inidivual_time_course = population_mixture_CA_array(:,(n-1)*num_compounds+1:n*num_compounds);
  CA_inidivual_AUC = sum(CA_inidivual_time_course(1:(end-1),:), 1) * output_interval_years; %nM*year
  CA_TEQ_inidivual_AUC = CA_inidivual_AUC' .* (g.MW .* g.Systemic_TEF); %ng/L*year
  CA_TEQ_AUC = [CA_TEQ_AUC, CA_TEQ_inidivual_AUC];
end
CA_TEQ_AUC_mean = mean(CA_TEQ_AUC,2);
CA_TEQ_AUC_mean_percentage = CA_TEQ_AUC_mean/sum(CA_TEQ_AUC_mean);
pie(CA_TEQ_AUC_mean_percentage);
title("Mixture CA TEQ AUC Distribution");
%legend(compound_names);

%Plot lifetime AUC CRB TEQ distribution(ng/L*year)
figure(16)
CRB_TEQ_AUC = [];
for n=1:num_individuals
  CRB_inidivual_time_course = population_mixture_CRB_array(:,(n-1)*num_compounds+1:n*num_compounds);
  CRB_inidivual_AUC = sum(CRB_inidivual_time_course(1:(end-1),:), 1) * output_interval_years; %nM*year
  CRB_TEQ_inidivual_AUC = CRB_inidivual_AUC' .* (g.MW .* g.Systemic_TEF); %ng/L*year
  CRB_TEQ_AUC = [CRB_TEQ_AUC, CRB_TEQ_inidivual_AUC];
end
CRB_TEQ_AUC_mean = mean(CRB_TEQ_AUC,2);
CRB_TEQ_AUC_mean_percentage = CRB_TEQ_AUC_mean/sum(CRB_TEQ_AUC_mean);
pie(CRB_TEQ_AUC_mean_percentage);
title("Mixture CRB TEQ AUC Distribution");
%legend(compound_names);


%Plot lifetime AUC CF TEQ distribution(ng/L*year)
figure(17)
CF_TEQ_AUC = [];
for n=1:num_individuals
  CF_inidivual_time_course = population_mixture_CF_array(:,(n-1)*num_compounds+1:n*num_compounds);
  CF_inidivual_AUC = sum(CF_inidivual_time_course(1:(end-1),:), 1) * output_interval_years; %nM*year
  CF_TEQ_inidivual_AUC = CF_inidivual_AUC' .* (g.MW .* g.Systemic_TEF); %ng/L*year
  CF_TEQ_AUC = [CF_TEQ_AUC, CF_TEQ_inidivual_AUC];
end
CF_TEQ_AUC_mean = mean(CF_TEQ_AUC,2);
CF_TEQ_AUC_mean_percentage = CF_TEQ_AUC_mean/sum(CF_TEQ_AUC_mean);
pie(CF_TEQ_AUC_mean_percentage);
title("Mixture CF TEQ AUC Distribution");
%legend(compound_names);


%Plot lifetime AUC CL_total TEQ distribution(ng/L*year)
figure(18)
CL_total_TEQ_AUC = [];
for n=1:num_individuals
  CL_total_inidivual_time_course = population_mixture_CL_total_array(:,(n-1)*num_compounds+1:n*num_compounds);
  CL_total_inidivual_AUC = sum(CL_total_inidivual_time_course(1:(end-1),:), 1) * output_interval_years; %nM*year
  CL_total_TEQ_inidivual_AUC = CL_total_inidivual_AUC' .* (g.MW .* g.Systemic_TEF); %ng/L*year
  CL_total_TEQ_AUC = [CL_total_TEQ_AUC, CL_total_TEQ_inidivual_AUC];
end
CL_total_TEQ_AUC_mean = mean(CL_total_TEQ_AUC,2);
CL_total_TEQ_AUC_mean_percentage = CL_total_TEQ_AUC_mean/sum(CL_total_TEQ_AUC_mean);
pie(CL_total_TEQ_AUC_mean_percentage);
title("Mixture CL_total TEQ AUC Distribution");
%legend(compound_names);


%Tissue TEQ timecourse
figure(2001)
plot(year, population_mixture_TEQ_CA_array(:,1:100), 'Color', [0.5 0.5 0.5])
hold on
population_mixture_TEQ_CA_mean = mean(population_mixture_TEQ_CA_array,2);
population_mixture_TEQ_CA_percentile = prctile(population_mixture_TEQ_CA_array,[2.5 97.5],2);
plot(year, population_mixture_TEQ_CA_mean, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
plot(year, population_mixture_TEQ_CA_percentile, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
set(gca, 'fontsize', 24)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("CA TEQ (ng/L)");
title("Population Mixture CA TEQ");


figure(2002)
plot(year, population_mixture_TEQ_CRB_array(:,1:100), 'Color', [0.5 0.5 0.5])
hold on
population_mixture_TEQ_CRB_mean = mean(population_mixture_TEQ_CRB_array,2);
population_mixture_TEQ_CRB_percentile = prctile(population_mixture_TEQ_CRB_array,[2.5 97.5],2);
plot(year, population_mixture_TEQ_CRB_mean, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
plot(year, population_mixture_TEQ_CRB_percentile, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
set(gca, 'fontsize', 24)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("CRB TEQ (ng/L)");
title("Population Mixture CRB TEQ");


figure(2003)
plot(year, population_mixture_TEQ_CF_array(:,1:100), 'Color', [0.5 0.5 0.5])
hold on
population_mixture_TEQ_CF_mean = mean(population_mixture_TEQ_CF_array,2);
population_mixture_TEQ_CF_percentile = prctile(population_mixture_TEQ_CF_array,[2.5 97.5],2);
plot(year, population_mixture_TEQ_CF_mean, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
plot(year, population_mixture_TEQ_CF_percentile, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
set(gca, 'fontsize', 24)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("CF TEQ (ng/L)");
title("Population Mixture CF TEQ");


figure(2004)
plot(year, population_mixture_TEQ_CL_total_array(:,1:100), 'Color', [0.5 0.5 0.5])
hold on
population_mixture_TEQ_CL_total_mean = mean(population_mixture_TEQ_CL_total_array,2);
population_mixture_TEQ_CL_total_percentile = prctile(population_mixture_TEQ_CL_total_array,[2.5 97.5],2);
plot(year, population_mixture_TEQ_CL_total_mean, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
plot(year, population_mixture_TEQ_CL_total_percentile, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3)
set(gca, 'fontsize', 24)
pbaspect([1.25,1,1])
xlabel("year");
ylabel("CL_total TEQ (ng/L)");
title("Population Mixture CL_total TEQ");

