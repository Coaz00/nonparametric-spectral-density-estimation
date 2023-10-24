clear
close all
clc

% Q = 1 -> x1.csv

% ucitavanje podataka
data = load('x1.csv');

% noramlizovane ucestanosti
Nf = 100;
f = linspace(-1/2, 1/2, Nf);

% dimenzije ucitane matrice
s = size(data);
R = s(1); 
N = s(2);

%% Periodogram

Pper = periodogram(data(1,:), f);

figure(1)
plot(f, Pper)
title('Periodogram')
xlabel('normalizovane ucestanosti')
ylabel('$\hat{P}_{xx}(f)$', 'Interpreter', 'Latex')

%% Average periodogram

K = 8;
Pavper = average_periodogram(data(1,:), f, K);

figure(2)
plot(f, Pavper)
title('Usrednjeni periodogram: K = 8')
xlabel('normalizovane ucestanosti')
ylabel('$\hat{P}_{xx}(f)$', 'Interpreter', 'Latex')

%% Welch method

K = 20;
p = 0.5;                      

% koristi se hamming-ov prozor
Pwelch = welch(data(1,:),f,K,p);

figure(3)
plot(f,Pwelch)
title('Welch-ova metoda: K = 20, p = 0.5')
xlabel('normalizovane ucestanosti')
ylabel('$\hat{P}_{xx}(f)$', 'Interpreter', 'Latex')

%% Blackman-Tukey method

window = "bartlet";

Pbt = blackman_tukey(data(1,:),f,window,N/2-1);

figure(4)
plot(f,Pbt)
title('BT metoda: Barlett-ov prozor, M = 127')
xlabel('normalizovane ucestanosti')
ylabel('$\hat{P}_{xx}(f)$', 'Interpreter', 'Latex')

%% 2. zad

% Izbor 5 razlicitih nasumicnih realizacija
while 1
    idx = randi([1 R],1,5);
    
    if length(idx) == length(unique(idx))
        break
    end
end

for i = 1:5
    figure(5)
    hold all
    plot(f,periodogram(data(idx(i),:),f))
    title('Periodogram za razlicite realizacije')
    xlabel('normalizovane ucestanosti') 
    ylabel('$\hat{P}_{xx}(f)$', 'Interpreter', 'Latex')
end

%% 3. zad

% Izbor optimalne duzine prozora za BT metodu

window = "bartlet";

for M = 1:128 
    figure(6)
    plot(f,blackman_tukey(data(1,:),f,window,2*M-1))
    2*M-1
    pause()
end 

%%     

% M = 25; 51 - vece, 13 manje
M_opt = 25;
M_big = 51;
M_small = 13;

figure(21)
plot(f,blackman_tukey(data(1,:),f,window,M_small))
hold all
plot(f,blackman_tukey(data(1,:),f,window,M_opt))
plot(f,blackman_tukey(data(1,:),f,window,M_big))
hold off
legend('M = 13','M = 25', 'M = 51')
title('BT metoda za razlicite duzine prozora')
xlabel('normalizovane ucestanosti')
ylabel('$\hat{P}_{xx}(f)$', 'Interpreter', 'Latex')

%% 4. zad

% Izbor optimalnih parametara za Welch-ov metod
p = [0.1 0.25 0.33 0.5 0.67 0.75];

for i = 1:length(p)
    p(i)
    for K = 1:70
        figure(7)
        plot(f,welch(data(1,:),f,K,p(i)))
        K
        pause()
    end
end

%%

% Optimalno: K=17, p=0.25

figure(22)
subplot(3,3,1)
plot(f,welch(data(1,:),f,11,0.1))
title('K = 11')
ylabel('p = 0.1','fontweight','bold')
subplot(3,3,4)
plot(f,welch(data(1,:),f,11,0.25))
ylabel('p = 0.25','fontweight','bold')
subplot(3,3,7)
plot(f,welch(data(1,:),f,11,0.67))
ylabel('p = 0.67','fontweight','bold')

subplot(3,3,2)
plot(f,welch(data(1,:),f,17,0.1))
title('K = 17')
subplot(3,3,5)
plot(f,welch(data(1,:),f,17,0.25))
subplot(3,3,8)
plot(f,welch(data(1,:),f,17,0.67))

subplot(3,3,3)
plot(f,welch(data(1,:),f,44,0.1))
title('K = 44')
subplot(3,3,6)
plot(f,welch(data(1,:),f,44,0.25))
subplot(3,3,9)
plot(f,welch(data(1,:),f,44,0.67))

%% 5. zad

Ppers = zeros(R,Nf);
Pwelchs = zeros(R,Nf);
Pbts = zeros(R,Nf);

% Estimacija SGS za sve realizacije
for i = 1:R
    Ppers(i,:) = periodogram(data(i,:),f);
    Pwelchs(i, :) = welch(data(i,:),f,17,0.25);
    Pbts(i,:) = blackman_tukey(data(i,:),f,'bartlet',25);
end

% Varijansa dobijena usrednjavanjem po realizacijama
var_per = var(Ppers);
var_welch = var(Pwelchs);
var_bt = var(Pbts);

% Varijanse po ucestanostima sve 3 metode
figure(8)
plot(f,var_per)
hold all
plot(f,var_welch)
plot(f,var_bt)
hold off
title('Varijanse za sve 3 metode')
xlabel('normalizovane ucestanosti')
legend('Per','Welch','BT')

% Uvelicane varijanse za BT i Welch-ov metod
figure(20)
plot(f,var_welch)
hold all
plot(f,var_bt)
hold off
title('Uvelican grafik varijansi za Welch-ov i BT metod')
xlabel('normalizovane ucestanosti')
legend('Welch','BT')

% Medijane varijansi
var_med_per = median(var_per);
var_med_welch = median(var_welch);
var_med_bt = median(var_bt);

%% 6. zad

per_aprox = mean(Ppers);

% Teorijski (aproksimativno) izracunata varijansa periodograma
var_per_theory = per_aprox.^2*(1 + (sin(2*pi*N.*f)/N/sin(2*pi.*f)).^2);

figure(9)
plot(f, var_per)
hold all
plot(f, var_per_theory)
title('Varijansa periodograma')
xlabel('normalizovane ucestanosti')
legend('Experiment', 'Theory')

%% 7. zad

% interval poverenja 1 - a = 95%
a = 0.05;
K = 12;
M = 25;

Pbt = blackman_tukey(data(1,:),f,'bartlet',M);
%v = 2*N/sum(bartlett(floor(M)).^2);
v = 3*N/(floor(M));

% Donja i gornja granica intervala poverenja
pr_down = Pbt./(chi2inv(1 - a/2,v))*v;
pr_up = Pbt./(chi2inv(a/2,v))*v;

% Iscrtavanje intervala poverenja
figure(10)
plot(f,10*log10(pr_up),'Color','r','LineStyle','--','LineWidth',1.5)
hold all
plot(f, 10*log10(pr_down),'Color','r','LineStyle','--','LineWidth',1.5)
hold off

% Usrednjeni periodogrami za sve realizacije
for i = 1:R
    Pavper = average_periodogram(data(i,:), f, K);
    figure(10)
    hold all
    plot(f,10*log10(Pavper),'Color',[0 0 1 0.2])
end

figure(10)
title('Usrednjeni periodogrami za sve realizacije')
xlabel('normalizovane ucestanosti')
legend('Interval poverenja')

%% 8. zad

Ppers_short = zeros(R,Nf);

% Periodogrami za skracenu sekvencu od 64 merenja
for i = 1:R
    Ppers_short(i,:) = periodogram(data(i,1:64),f);
end

% medijane varijansi periodograma pune i skracene sekvence
var_per_short = var(Ppers_short);
var_med_per_short = median(var_per_short);

figure(11)
plot(f,var_per,f,var_per_short)
title('Varijanse periodograma za razlicite duzine sekvence')
legend('N = 256','N = 64')
xlabel('normalizovane ucestanosti')
