%PROPEP32 obliczenia parametrow spalania w silniku rakietowym na ciekly
%material pedny

%czyszczenie okna 
clc, clear;


%dane do podania:
p_cham=input('podaj cisnienie w komorze spalania w atm ');
p_e=input('podaj cisnienie na wyjsciu z dyszy w atm ');
p_ox=input('podaj cisnienie utleniacza w zbiorniku w atm ');
T0=input('podaj temperature substratow w K: ');

%zmienne wystepujace w kodzie
ro_ox=1.799; %gestosc n2o w kg/m^
d_ox=1/1000; %srednica wtrysku w m
n_ox=16; %liczba wtryskiwaczy

% oblicznaie wydatku utleniacza na podstawie cisnienie w instalacji
dm_ox=n_ox*pi*d_ox^2/4*sqrt(2*1.31*p_ox*oneatm/((1.31-1)*ro_ox)*((p_cham/p_ox)^(2/1.31)...
    -(p_cham/p_ox)^((1.31+1)/1.31)));

%zakres petli o_f
n_p=5;% poczatek
n_k=7;% koniec
h=0.1;%krok
xx=1;

%petla obliczniowa
for ii=n_p:h:n_k 
    
    %wprowadzanie danych do pliku input
    f_in=fopen('input.dat','w');
    fprintf(f_in,'print.txt\n');%zdefiniowanie pliku out
    fprintf(f_in,'\n');
    string2=sprintf('    %d    %d    %d ',2,0,1);
    fprintf(f_in,'%s\n',string2);
    fprintf(f_in,'1.\n');
    fprintf(f_in,'298\n');
    fprintf(f_in,'\n');
    fprintf(f_in,'0001000000\n');
    fprintf(f_in,' 1095  708\n');
    string3=sprintf('%d.0        %d.0         100.        %d',p_cham,p_e,ii*100);
    string3=[string3 '.'];
    fprintf(f_in,'%s\n',string3);
    fclose(f_in);   %zamkniecie pliku input
    system('skrypt.vbs');   %odpalenie skryptu uruchamiajacego program
    pause(0.1);
    status=dos('taskkill /F /IM propep32.exe');
    while status==0
        system('skrypt.vbs');   %odpalenie skryptu uruchamiajacego program
        pause(0.1);
        status=dos('taskkill /F /IM propep32.exe');
    end
    %otwarcie pliku output jako str
    f_out=fileread('print.txt');
    %sczytanie wynikow z output
    index=strfind(f_out,'RT/V');

    T_cham(xx)=sscanf(f_out(index+7:end),'%g',1);%temperatura w komorze
    index2=strfind(f_out,'EX-T');
    Isp(xx)=sscanf(f_out(index2+9:end),'%g',1);%impuls wlasciwy
    xx=xx+1;
    ii+ii+h;
    fclose('all');
    end
 
%rysowanie wykresow
figure(1);
subplot(1,2,1);
plot([n_p:h:n_k],Isp);
grid minor
axis([n_p n_k 226 228]);
xlabel('Oxidizer-Fuel Ratio');
ylabel('Specific impulse [s]');

subplot(1,2,2);
plot([n_p:h:n_k],T_cham);
grid minor
xlabel('Oxidizer-Fuel Ratio');
ylabel('Temperature in combustion chamber [K]');
[Isp_max, ind]=max(Isp);

fprintf('Max specific impulse: Isp = %.2fs for o/f = %.2f\n',Isp_max,(ind-1)*h+n_p);
fprintf('Temperature in combustion chamber for max Isp: T_cham = %0.f K\n',T_cham(ind)); 
fprintf('Fuel consumption for max Isp: dm_ful=%.5f kg/s\n',dm_ox/((ind-1)*h+n_p));
