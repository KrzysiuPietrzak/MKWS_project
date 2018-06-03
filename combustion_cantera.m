% kod spalania kerozyny z podtlenkiem azotu w silniku raketowym na ciekly
% material pedny z wykorzystaniem Cantery

%zmienne do podania
p_ox=50*oneatm; %cisnienie utleniacza w zbiorniku atm
d_ox=1/1000; %srednica wtrysku utleniacza w mm
p_cham=20*oneatm; %cisnienie w komorze spalania w atm
T0=298; %temperatura poczatkowowa w K
p_e=oneatm; %cisnienie na wylocie z dyszy

%zmienne wystepujace w kodzie
M_ful=170.34; % masa molowa paliwa w g/mol
M_ox=44.0129; %masa molowa utleniacza
ro_ox=1.799; %gestosc n2o w kg/m^3
R=8314.46; %uniwersalna stala gazowa w J/(kmol*K)
g=9.81;
%entalpie tworzenia zwiazkow
h_ful=-131060; % entalpia twozenia rp-1 w J/kg
h_ox=1860000; %entalpia tworzenia n2o w J/kg
h_co=-110.509*1e6/28.01;
h_n2=0;
h_h2=0;
h_h2o=-241.818*1e6/18.02;
h_co2=-393.509*1e6/44.01;

%obliczanie wydatku ulteniacza
dm_ox=pi*d_ox^2/4*sqrt(2*1.31/(1.31-1)*p_ox*ro_ox*((p_cham/p_ox)^(2/1.31)...
    -(p_cham/p_ox)^((1.31+1)/1.31)));

%przedial o do f
o_f=[1 10];
h=1;%krok petli
xx=1;
T_cham=zeros(1,(o_f(2)-o_f(1))/h+1);
Isp=zeros(1,(o_f(2)-o_f(1))/h+1);

for ii=o_f(1):h:o_f(2)
    dm_ful=dm_ox/ii;
    if ii<3.1
        % a1*C12H26 + a2*N2O -> a3*CO + a4*H2 + a5*N2 + a6*C12H26
        c=1;
        a=zeros(1,6);
        a(1)=dm_ful*1000/M_ful;
        a(2)=dm_ox*1000/M_ox;
        a(3)=a(2);
        a(5)=a(2);
        a(6)=a(1)-a(3)/12;
        a(4)=0.5*(26*a(1)-26*a(6));
        dm_co=a(3)*28.01/1000;
        dm_h2=a(4)*2/1000;
        dm_n2=a(5)*28/1000;
        dm_c12h26=a(6)*M_ful/1000;
        dm=dm_ox+dm_ful; % calkowity wydatek gazow wylotowych, kg/s
        H_subs=h_ful*dm_ful+h_ox*dm_ox;
        H_prod=h_co*dm_co+dm_c12h26*h_ful;
        H=-(H_prod-H_subs);
        %obliczanie parametrow gazow wylowych (wykorzystanie Cantery)
        %T_cham -temp w komorze spalania [K] ; k- wyk³adnik adiabaty
        % M_gas- masa molowa gazow wylotowych 
        [T_cham(xx), k, M_gas]=temp_calc(c,H,a,p_cham,T0,dm);
        M_gas=M_gas+M_ful*a(6)/sum(a(3:5));
        Isp(xx)=sqrt(k*R/M_gas*T_cham(xx)/(k-1)*(1-(p_e/p_cham)^((k-1)/k)))/g;
        clear a
    elseif ii==3.1
        % a1*C12H26 + a2*N2O -> a3*CO + a4*H2 + a5*N2 
        c=2;
        a=zeros(1,5);
        a(1)=dm_ful*1000/M_ful;
        a(2)=dm_ox*1000/M_ox;
        a(3)=a(2);
        a(4)=a(3)*26/24;
        a(5)=a(2);
        dm_co=a(3)*28.01/1000;
        dm_h2=a(4)*2/1000;
        dm_n2=a(5)*28/1000;
        dm=dm_ox+dm_ful; % calkowity wydatek gazow wylotowych, kg/s
        H_subs=h_ful*dm_ful+h_ox*dm_ox;
        H_prod=h_co*dm_co;
        H=-(H_prod-H_subs);
        %obliczanie parametrow gazow wylowych (wykorzystanie Cantery)
        %T_cham -temp w komorze spalania [K] ; k- wyk³adnik adiabaty
        % M_gas- masa molowa gazow wylotowych 
        [T_cham(xx), k, M_gas]=temp_calc(c,H,a,p_cham,T0,dm);
        Isp(xx)=sqrt(k*R/M_gas*T_cham(xx)/(k-1)*(1-(p_e/p_cham)^((k-1)/k)))/g;
        clear a
    elseif ii>3.1 && ii<6.2
        % a1*C12H26 + a2*N2O -> a3*CO + a4*H2 + a5*H2O + a6*N2 
        c=3;
        a=zeros(1,6);
        a(1)=dm_ful*1000/M_ful;
        a(2)=dm_ox*1000/M_ox;
        a(3)=12*a(1);
        a(5)=a(2)-a(3);
        a(4)=0.5*(26*a(1)-2*a(5));
        a(6)=a(2);
        dm_co=a(3)*28.01/1000;
        dm_h2=a(4)*2/1000;
        dm_n2=a(6)*28/1000;
        dm_h2o=a(5)*18/1000;
        dm=dm_ox+dm_ful; % calkowity wydatek gazow wylotowych, kg/s
        H_subs=h_ful*dm_ful+h_ox*dm_ox;
        H_prod=h_co*dm_co+dm_h2o*h_h2o;
        H=-(H_prod-H_subs);
        %obliczanie parametrow gazow wylowych (wykorzystanie Cantery)
        %T_cham -temp w komorze spalania [K] ; k- wyk³adnik adiabaty
        % M_gas- masa molowa gazow wylotowych 
        [T_cham(xx), k, M_gas]=temp_calc(c,H,a,p_cham,T0,dm);
        Isp(xx)=sqrt(k*R/M_gas*T_cham(xx)/(k-1)*(1-(p_e/p_cham)^((k-1)/k)))/g;
        clear a
    elseif ii>=6.2 && ii<9.56
        % a1*C12H26 + a2*N2O -> a3*CO2+ a4*H2O + a5*N2 + a6*H2
        c=4;
        a=zeros(1,5);
        a(1)=dm_ful*1000/M_ful;
        a(2)=dm_ox*1000/M_ox;
        a(3)=12*a(1);
        a(4)=a(2)-2*a(3);
        a(5)=a(2);
        a(6)= 0.5*(a(1)*26-2*a(4));
        dm_co2=a(3)*44.01/1000;
        dm_h2o=a(4)*18/1000;
        dm_n2=a(5)*28/1000;
        dm_h2=a(6)*2/1000;
        dm=dm_ox+dm_ful; % calkowity wydatek gazow wylotowych, kg/s
        H_subs=h_ful*dm_ful+h_ox*dm_ox;
        H_prod=h_co2*dm_co2+dm_h2o*h_h2o;
        H=-(H_prod-H_subs);
        %obliczanie parametrow gazow wylowych (wykorzystanie Cantery)
        %T_cham -temp w komorze spalania [K] ; k- wyk³adnik adiabaty
        % M_gas- masa molowa gazow wylotowych 
        [T_cham(xx), k, M_gas]=temp_calc(c,H,a,p_cham,T0,dm);
        Isp(xx)=sqrt(k*R/M_gas*T_cham(xx)/(k-1)*(1-(p_e/p_cham)^((k-1)/k)))/g;
        clear a
    else
        % a1*C12H26 + a2*N2O -> a3*CO2+ a4*H2O + a5*N2 + a6*N2O
        c=5;
        a=zeros(1,5);
        a(1)=dm_ful*1000/M_ful;
        a(2)=dm_ox*1000/M_ox;
        a(3)=12*a(1);
        a(4)=13*a(1);
        a(6)=a(2)-(2*a(3)+a(4));
        a(5)= a(2)-a(6);
        dm_co2=a(3)*44.01/1000;
        dm_h2o=a(4)*18/1000;
        dm_n2=a(5)*28/1000;
        dm_n2o=a(6)*44.0129/1000;
        dm=dm_ox+dm_ful; % calkowity wydatek gazow wylotowych, kg/s
        H_subs=h_ful*dm_ful+h_ox*dm_ox;
        H_prod=h_co2*dm_co2+dm_h2o*h_h2o+dm_n2o*h_ox;
        H=-(H_prod-H_subs);
        %obliczanie parametrow gazow wylowych (wykorzystanie Cantery)
        %T_cham -temp w komorze spalania [K] ; k- wyk³adnik adiabaty
        % M_gas- masa molowa gazow wylotowych 
        [T_cham(xx), k, M_gas]=temp_calc(c,H,a,p_cham,T0,dm);
        Isp(xx)=sqrt(k*R/M_gas*T_cham(xx)/(k-1)*(1-(p_e/p_cham)^((k-1)/k)))/g;
        clear a
    end
    
    xx=xx+1; 
end
%making plots and printing results
figure(1);
subplot(1,2,1);
plot([o_f(1):h:o_f(2)],Isp);
grid minor
xlabel('Oxidizer-Fuel Ratio');
ylabel('Specific impulse [s]');

subplot(1,2,2);
plot([o_f(1):h:o_f(2)],T_cham);
grid minor
xlabel('Oxidizer-Fuel Ratio');
ylabel('Temperature in combustion chamber [K]');
[Isp_max, ind]=max(Isp);

fprintf('Max specific impulse: Isp = %.2fs for o/f = %.2f\n',Isp_max,(ind-1)*h+o_f(1));
fprintf('Temperature in combustion chamber for max Isp: T_cham = %0.f K\n',T_cham(ind)); 
fprintf('Fuel consumption for max Isp: dm_ful=%.5f kg/s\n',dm_ox/((ind-1)*h+o_f(1)));
