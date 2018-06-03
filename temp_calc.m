function [T_cham, gamma,M] = temp_calc(c,H,a,p_cham,T0,dm)
% Kod wykorzystujacy Cantere w clu obliczen temeratury w komorze spalania z
% uzwglednieniem zmiennego ciepla wlasciwego gazow
% c- przyapadek do obliczen
% H- eneragia spalania
% a- wektor wspolczynnikow reakcji
% dm- calkowity wydatek gazow wylotowych
switch c
    case 1
        % a1*C12H26 + a2*N2O -> a3*CO + a4*H2 + a5*N2 + a6*C12H26
        gas=Solution('gri30.cti','gri30');
        ih2=speciesIndex(gas,'H2');
        ico=speciesIndex(gas,'CO');
        in2=speciesIndex(gas,'N2');
        x=zeros(nSpecies(gas),1);
        x(ih2,1) = a(4);
        x(ico,1) = a(3);
        x(in2,1) = a(5);
        set(gas,'P',p_cham,'T',T0,'X',x);
        cp=cp_mass(gas);
        T_cham=T0;
        dT=100;
        while H>0
            dH=cp*dm*dT;
            H=H-dH;
            T_cham=T_cham+dT;
            set(gas,'P',p_cham,'T',T_cham);
            cp=cp_mass(gas);
        end
        T_cham=T_cham-dT;
        dH=cp*dT*dm;
        H=H+dH;
        T_cham=T_cham+H/(cp*dm);
        set(gas,'P',p_cham,'T',T_cham);
        equilibrate(gas,'HP');
        gamma=cp_mass(gas)/cv_mass(gas);
        M=meanMolecularWeight(gas);
        T_cham=temperature(gas);
    case 2
        % a1*C12H26 + a2*N2O -> a3*CO + a4*H2 + a5*N2 
        gas=Solution('gri30.cti','gri30');
        ih2=speciesIndex(gas,'H2');
        ico=speciesIndex(gas,'CO');
        in2=speciesIndex(gas,'N2');
        x=zeros(nSpecies(gas),1);
        x(ih2,1) = a(4);
        x(ico,1) = a(3);
        x(in2,1) = a(5);
        set(gas,'P',p_cham,'T',T0,'X',x);
        cp=cp_mass(gas);
        T_cham=T0;
        dT=100;
        while H>0
            dH=cp*dm*dT;
            H=H-dH;
            T_cham=T_cham+dT;
            set(gas,'P',p_cham,'T',T_cham);
            cp=cp_mass(gas);
        end
        T_cham=T_cham-dT;
        dH=cp*dT*dm;
        H=H+dH;
        T_cham=T_cham+H/(cp*dm);
        set(gas,'P',p_cham,'T',T_cham);
        equilibrate(gas,'HP');
        gamma=cp_mass(gas)/cv_mass(gas);
        M=meanMolecularWeight(gas);
        T_cham=temperature(gas);
    case 3
        % a1*C12H26 + a2*N2O -> a3*CO + a4*H2 + a5*H2O + a6*N2
        gas=Solution('gri30.cti','gri30');
        ih2=speciesIndex(gas,'H2');
        ico=speciesIndex(gas,'CO');
        in2=speciesIndex(gas,'N2');
        ih2o=speciesIndex(gas,'H2O');
        x=zeros(nSpecies(gas),1);
        x(ih2,1) = a(4);
        x(ico,1) = a(3);
        x(in2,1) = a(6);
        x(ih2o,1)= a(5);
        set(gas,'P',p_cham,'T',T0,'X',x);
        cp=cp_mass(gas);
        T_cham=T0;
        dT=100;
        while H>0
            dH=cp*dm*dT;
            H=H-dH;
            T_cham=T_cham+dT;
            set(gas,'P',p_cham,'T',T_cham);
            cp=cp_mass(gas);
        end
        T_cham=T_cham-dT;
        dH=cp*dT*dm;
        H=H+dH;
        T_cham=T_cham+H/(cp*dm);
        set(gas,'P',p_cham,'T',T_cham);
        equilibrate(gas,'HP');
        gamma=cp_mass(gas)/cv_mass(gas);
        M=meanMolecularWeight(gas);
        T_cham=temperature(gas);
    case 4
        % a1*C12H26 + a2*N2O -> a3*CO2+ a4*H2O + a5*N2 
        gas=Solution('gri30.cti','gri30');
        ih2o=speciesIndex(gas,'H2O');
        ico2=speciesIndex(gas,'CO2');
        in2=speciesIndex(gas,'N2');
        ih2=speciesIndex(gas,'H2');
        x=zeros(nSpecies(gas),1);
        x(ico2,1) = a(3);
        x(in2,1) = a(5);
        x(ih2o,1)= a(4);
        x(ih2,1)= a(6);
        set(gas,'P',p_cham,'T',T0,'X',x);
        cp=cp_mass(gas);
        T_cham=T0;
        dT=100;
        while H>0
            dH=cp*dm*dT;
            H=H-dH;
            T_cham=T_cham+dT;
            set(gas,'P',p_cham,'T',T_cham);
            cp=cp_mass(gas);
        end
        T_cham=T_cham-dT;
        dH=cp*dT*dm;
        H=H+dH;
        T_cham=T_cham+H/(cp*dm);
        set(gas,'P',p_cham,'T',T_cham);
        equilibrate(gas,'HP');
        gamma=cp_mass(gas)/cv_mass(gas);
        M=meanMolecularWeight(gas);
        T_cham=temperature(gas);
    case 5
        % a1*C12H26 + a2*N2O -> a3*CO2+ a4*H2O + a5*N2 + a6*N2O
        gas=Solution('gri30.cti','gri30');
        ih2o=speciesIndex(gas,'H2O');
        ico2=speciesIndex(gas,'CO2');
        in2=speciesIndex(gas,'N2');
        in2o=speciesIndex(gas,'N2O');
        x=zeros(nSpecies(gas),1);
        x(ico2,1) = a(3);
        x(in2,1) = a(5);
        x(ih2o,1)= a(4);
        x(in2o,1)= a(6);
        set(gas,'P',p_cham,'T',T0,'X',x);
        cp=cp_mass(gas);
        T_cham=T0;
        dT=100;
        while H>0
            dH=cp*dm*dT;
            H=H-dH;
            T_cham=T_cham+dT;
            set(gas,'P',p_cham,'T',T_cham);
            cp=cp_mass(gas);
        end
        T_cham=T_cham-dT;
        dH=cp*dT*dm;
        H=H+dH;
        T_cham=T_cham+H/(cp*dm);
        set(gas,'P',p_cham,'T',T_cham);
        equilibrate(gas,'HP');
        gamma=cp_mass(gas)/cv_mass(gas);
        M=meanMolecularWeight(gas);
        T_cham=temperature(gas);
end
