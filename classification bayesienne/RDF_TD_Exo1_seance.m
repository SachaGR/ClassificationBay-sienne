% TD Reconnaissance de Formes - Exercice 1

clc;
close all;
clear;

% Compléter les informations manquantes (notées ...)
Max_x = 4;
Min_x = -4;
Pas = 0.01;
x= [Min_x:Pas:Max_x];
L = length(x);

% Considèrons que 2 classes sont modélisées par les densités de
% probabilité gaussiennes suivantes :
pxw1 = exp(-(x.*x))./sqrt(pi);
mu2 = 1;
sigma2 = sqrt(0.5); % ecart-type
%sigma2 = 2;
pxw2 = exp(-(x-mu2).*(x-mu2)/(2*sigma2^2))./sqrt(2*pi*sigma2^2);

%% 1. 1er cas:  Pw1=Pw2=0,5
Pw1 = 0.5;
Pw2 = 1 - Pw1;

px = pxw1 * Pw1 + pxw2 * Pw2; % calcul de la probabilité totale

Pw1x = pxw1 * Pw1 ./ px;
Pw2x = pxw2 * Pw2 ./ px;


% Seuil
i = 1;
g1 = log(pxw1(1) * Pw1);
g2 = log(pxw2(1) * Pw2);
for k = 2 : length(x)
    g1m = log(pxw1(k) * Pw1);
    g2m = log(pxw2(k) * Pw2);
    
    if (g1 - g2)*(g1m - g2m) < 0
        xB(i) = x(k);
        i = i + 1;
    end
    
    g1 = log(pxw1(k) * Pw1);
    g2 = log(pxw2(k) * Pw2);
    
end

% Calcul de l'erreur
%indices = find(x <= xB(1) | xB(2) <= x);
indices = find(x <= xB(1));
M = length(indices);
Perreur = Pas * sum(pxw2(indices) * Pw2);
Perreur = Perreur + Pas * sum(pxw1(indices(M) + 1:L) * Pw1);

%  Tracés des densités de probabilités des 2 classes
figure(1)
hold on;
plot(x,pxw1,'color','red');
plot(x,pxw2,'color','green');
plot(x,Pw1x,'color','blue');
plot(x,Pw2x,'color','magenta');
plot([xB(1) xB(1)],[0 1],'--','color','black');
%plot([xB(2) xB(2)],[0 1],'--','color','black');
hold off;
legend('Classe 1 : densité de probabilité','Classe 2 : densité de probabilité','Classe 1 : propabilité à postériori','Classe 2 : propabilité à postériori','x_B');








%% 2. Modification des valeurs de probabilité à priori
Pw1 = 9/10;
Pw2 = 1/10;

px = pxw1 * Pw1 + pxw2 * Pw2; % calcul de la probabilité totale

Pw1x = pxw1 * Pw1 ./ px;
Pw2x = pxw2 * Pw2 ./ px;


% Seuil
i = 1;
g1 = log(pxw1(1) * Pw1);
g2 = log(pxw2(1) * Pw2);
for k = 2 : length(x)
    g1m = log(pxw1(k) * Pw1);
    g2m = log(pxw2(k) * Pw2);
    
    if (g1 - g2)*(g1m - g2m) < 0
        xB(i) = x(k);
        i = i + 1;
    end
    
    g1 = log(pxw1(k) * Pw1);
    g2 = log(pxw2(k) * Pw2);
    
end

% Calcul de l'erreur
indices = find(x <= xB);
M = length(indices);
Perreur = Pas * sum(pxw2(indices) * Pw2);
Perreur = Perreur + Pas * sum(pxw1(indices(M) + 1:L) * Pw1);

%  Tracés des densités de probabilités des 2 classes
figure(2) 
hold on
plot(x,pxw1,'color','red');
plot(x,pxw2,'color','green');
plot(x,Pw1x,'color','blue');
plot(x,Pw2x,'color','magenta');
plot([xB xB],[0 1],'--','color','black');
hold off;
legend('Classe 1 : densité de probabilité','Classe 2 : densité de probabilité','Classe 1 : propabilité à postériori','Classe 2 : propabilité à postériori','x_B');

    
% Prise en compte du numérateur uniquement
Pw1x_modif = pxw1 * Pw1;
Pw2x_modif = pxw2 * Pw2;


figure(3) 
hold on
plot(x,pxw1,'color','red');
plot(x,pxw2,'color','green');
plot(x,Pw1x_modif,'color','blue');
plot(x,Pw2x_modif,'color','magenta');
plot([xB xB],[0 0.6],'--','color','black');
hold off;
legend('Classe 1 : densité de probabilité','Classe 2 : densité de probabilité','Classe 1 : propabilité à postériori','Classe 2 : propabilité à postériori','x_B');

