%function [T_s, TSFC] = JetJoe(delta_eta)

clear all
data0 = [0 43000 79000 111500 161000;
 14.79 14.79 14.79 14.79 14.79;
 0 0.140 0.475 1.105 2.391;
 11.4 11.5 11.9 11.1 11.1;
 14.5 15.7 19.0 24.6 37.9;
 13.7 24.7 44.8 76.2 149.8;
 0 0.7 2.6 5.9 15.0;
 0 0.74 1.54 2.34 4.3;
 10 596 586 600 802]; 
%----------------------------------------------------------
% PAST YEAR'S DATA - 10 DATA SETS
% UNI2016 Set1
data1 = [0 42900 65500 104000 159500;
 14.79 14.79 14.79 14.79 14.79;
 0 0.116 0.334 0.91 2.44;
 12.1 12.0 11.3 11.5 11.6;
 14.5 15.7 17.6 22.7 38.0;
 15.5 24 35.85 66.4 148.7;
 0 0.8 2.0 5.2 15.5;
 0 0.8 1.66 2.64 4.18;
 10 654 600 602 834];
% UNI2015 Set 2
data2 = [0 42500 83200 114600 160000;
 14.64 14.64 14.64 14.64 14.64;
 0 0.11 0.471 1.082 2.28;
 23 22.7 22.7 22.8 22.8;
 14.4 15.6 18.9 24.2 36.5;
 25 36.4 58.5 89.9 160.7;
 0 0.6 2.7 6.0 14.5;
 0 0.80 1.48 2.27 4.3;
 22 641 586 604 822];
% UNI2015 Set 3
data3 = [0 42000 80300 110400 160300;
 14.64 14.64 14.64 14.64 14.64;
 0 0.111 0.456 1.000 2.300;
 23.0 22.8 22.7 22.7 22.6;
 14.4 15.6 18.9 23.6 36.8;
 27.3 36.9 57.4 86.6 160.8;
 0 0.7 2.7 5.7 14.6;
 0 0.80 1.48 2.21 4.30;
 24 642 586 600 810];
% UNI2015 Set 4
data4 = [0 43000 81000 113800 161000;
 14.65 14.65 14.65 14.65 14.65;
 0 0.112 0.470 1.080 2.280;
 23.2 22.9 23.0 22.8 23.1;
 14.4 15.6 19.0 24.3 36.5;
 26.2 37.4 58.7 90.7 161.9;
 0 0.7 2.7 6.0 14.6;
 0 0.80 1.29 1.91 3.81;
 22 630 575 597 808];
% UNI2015 Set 5
data5 = [0 45000 80000 120000 160000;
 14.71 14.71 14.71 14.71 14.71;
 0 0.120 0.477 1.020 2.307;
 23.1 22.8 22.9 23.1 23.11;
 14.5 15.8 19.1 24.0 36.4;
 25.1 36.2 57.6 88.0 160.0;
 0 0.8 2.8 5.9 14.5;
 0 0.86 1.54 2.21 4.30;
 23 632 590 608 812]; 
% UNI2015 Set 6
data6 = [0 44000 80400 112500 159400;
 14.71 14.71 14.71 14.71 14.71;
 0 0.120 0.468 1.050 2.270;
 23.4 23.1 22.9 23.1 23.3;
 14.5 15.7 19.0 24.2 36.5;
 27.2 36.8 58.7 89.7 160.8;
 0 0.8 2.8 6.1 14.4;
 0 0.86 1.54 2.21 4.30;
 24 634 586 606 814];
% UNI2015 Set 7
data7 = [0 43800 82000 110900 159500;
 14.71 14.71 14.71 14.71 14.71;
 0 0.110 0.482 1.014 2.270;
 23.2 22.9 22.9 23.0 23.0;
 14.5 15.7 19.1 24.2 37.0;
 27.4 37.2 58.7 89.4 160.6;
 0 0.8 3.0 5.9 14.7;
 0 0.80 1.54 2.27 4.43;
 24 638 590 632 854];
% UNI2015 Set 8
data8 = [0 42000 79800 112500 160500;
 14.71 14.71 14.71 14.71 14.71;
 0 0.113 0.453 1.056 2.285;
 23.4 23.2 23.3 23.5 23.6;
 14.5 15.7 18.9 24.2 36.9;
 27.4 37.4 57.4 89.6 161.5;
 0 0.8 2.7 6.0 14.6;
 0 0.86 1.54 2.27 4.43;
 24 684 636 646 848];
% UNI2015 Set 9
data9 = [0 42600 85700 112200 159200;
 14.79 14.79 14.79 14.79 14.79;
 0 0.115 0.544 1.053 2.258;
 23.7 23.3 23.6 23.4 23.6;
 14.6 15.8 19.7 24.3 37.1;
 25.6 37.5 63.1 89.8 158.8;
 0 0.7 3.1 5.9 14.3;
 0 0.86 1.60 2.27 4.30;
 20 688 634 650 860];
% UNI2015 Set 10
data10 = [0 43000 84400 111800 157700;
 14.79 14.79 14.79 14.79 14.79;
 0 0.114 0.520 1.048 2.250;
 23.7 23.5 23.5 23.7 23.8;
 14.6 15.8 19.8 24.2 36.9;
 26.6 37.7 61.0 89.7 157.8;
 0 0.7 3.0 5.9 14.1;
 0 0.86 1.60 2.27 4.12;
 22 696 636 648 840]; 

data = (data0 + data1 + data2 + data3 + data4 + data5 + data6 + data7 + data8 + data9 + data10)/11;

%Constants
gamma_c = 1.4;
gamma_h = 1.3;
R = 287;
rho_f = 6.843; %lb/gallon
A_1 = (3.5*0.0254/2)^2*pi; %cm-> m ->m^2
c_pc = R*gamma_c/(gamma_c-1);
c_ph = R*gamma_h/(gamma_h-1);

%Data
rpm = data(1,:);
p_0 = data(2,:)*6894.8; %psi -> Pa
T_0 = data(4,:)+273.15; %C -> K
p_1 = p_0 - data(3,:)*0.03613*6894.8; %inH2O -> psi -> Pa
p_3 = data(5,:)*6894.8; %psi -> Pa
T_3 = data(6,:)+273.15; %C -> K
T_5 = data(9,:)+273.15; %C -> K
T = data(7,:)*4.4482;   %lb -> N
mdot_f = data(8,:)*rho_f/3600*0.4536; %gal/hr * lb/gal * hr/s * kg/lb = kg/s
a_0 = sqrt(gamma_c*R*T_0);

%Part 4
M_1 = sqrt(((p_0./p_1).^(1-1/gamma_c)-1)/(gamma_c/2-1/2));
mdot_corrected = (1+(gamma_c-1)/2*M_1.^2).^(1/2*(1+gamma_c)/(1-gamma_c)).* M_1;
T_1 = T_0.*(1+(gamma_c-1)/2*M_1.^2).^-1;
mdot = p_1*A_1.*M_1.*sqrt(gamma_c./(R*T_1));
f = mdot_f./mdot;

%Part 5
pi_c = p_3./p_0;
T_3s = T_0.*pi_c.^(1-1/gamma_c);
eta_c = T_0.*(pi_c.^(1-1/gamma_c)-1)./(T_3-T_0);
T_4 = c_pc/c_ph.*(T_3 - T_0)./(1+f)+T_5;
%p_6 = p_0

%Part 6
delta_eta = [0 0.5 -.05 -.1 -.1];
eta_h = eta_c + delta_eta;
T_5s = (T_5 - T_4)./eta_h + T_4;
p_5 = p_3.*(T_5s./T_4).^(gamma_h/(gamma_h-1));
T_6_static = T_5.*(p_0./p_5).^(1-1/gamma_h);
M_6 = sqrt(((p_5./p_0).^(1-1/gamma_h)-1)/(gamma_h/2-1/2));
c_6 = M_6.*sqrt(gamma_h*R*T_6_static);
tau_s = pi_c.^(1-1/gamma_c);
eta_t = (1-1./tau_s).*(eta_c.*eta_h.*T_4./T_0-tau_s)./(1+eta_c.*(T_4./T_0-1)-tau_s);


%Specific thrust and TSFC

% subplot(2,1,1)
% specific_thrust = c_6./a_0;
% plot(rpm, specific_thrust, 'b--o', rpm, T./mdot./a_0, 'r-o')
% legend('Model Specific Thrust', 'Measured Specific Thrust')
% title('Model vs. Measured Specific Thrust')
% xlabel('Engine RPM')
% ylabel('Specific Thrust [unitless]')
% subplot(2,1,2)
% thrust = mdot.*a_0.*c_6./(sqrt(gamma_c*R*T_0));
% plot(rpm, mdot_f./thrust, 'b--o', rpm, mdot_f./T, 'r-o')
% title('Model vs. Measured Thrust Specific Fuel Consumption')
% legend('Model TSFC', 'Measured TSFC')
% xlabel('Engine RPM')
% ylabel('TSFC [kg/N-s]')

%Compressor and Turbine Adiabatic Efficiency, Thermal Efficiency

% subplot(3,1,1)
% plot(rpm,eta_c, 'r-o')
% title('Compressor Adiabatic Efficiency')
% xlabel('Engine RPM')
% ylabel('Compressor Adiabatic Efficiency [unitless]')
% subplot(3,1,2)
% plot(rpm,eta_h, 'r-o')
% title('Turbine Adiabatic Efficiency')
% xlabel('Engine RPM')
% ylabel('Turbine Adiabatic Efficiency [unitless]')
% subplot(3,1,3)
% plot(rpm,eta_t, 'r-o')
% title('Thermal Efficiency')
% xlabel('Engine RPM')
% ylabel('Thermal Efficiency [unitless]')

%Compressor and overall engine pressure ratio 

% subplot(2,1,1)
% plot(rpm,pi_c, 'r-o')
% title('Compressor Pressure Ratio')
% xlabel('Engine RPM')
% ylabel('Pressure Ratio [unitless]')
% subplot(2,1,2)
% plot(rpm,p_5./p_0, 'r-o')
% title('Overall Engine Pressure Ratio')
% xlabel('Engine RPM')
% ylabel('Pressure Ratio [unitless]')

%Compressor exit to inlet and turbine inlet to compressor inlet temperature
%ratio

% subplot(2,1,1)
% plot(rpm,T_3./T_0, 'r-o')
% title('Compressor Exit to Inlet Temperature Ratio')
% xlabel('Engine RPM')
% ylabel('Temperature Ratio [unitless]')
% subplot(2,1,2)
% plot(rpm,T_4./T_0, 'r-o')
% title('Turbine Inlet to Compressor Inlet Temperature Ratio')
% xlabel('Engine RPM')
% ylabel('Temperature Ratio [unitless]')

%Fuel mass flow and engine air mass flow

% subplot(2,1,1)
% plot(rpm,mdot_f, 'r-o')
% title('Fuel Mass Flow')
% xlabel('Engine RPM')
% ylabel('Mass Flow [kg/s]')
% subplot(2,1,2)
% plot(rpm,mdot,'r-o')
% title('Engine Air Mass Flow')
% xlabel('Engine RPM')
% ylabel('Mass Flow [kg/s]')

shg()
