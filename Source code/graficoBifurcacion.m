% diagrama de bifurcacion hacia adelante (forward bifuraction)

mu=1/(82*365); %mortalidad natural
alfa = 1/7; %periodo de latencia -> 7 dias
gamma = 1/17; %periodo de recuperacion -> 17 dias

beta1 = linspace(0.001,3,2000); % valores de beta1
beta2 = linspace(0.001,3,2000); % valores de beta2
beta3 = linspace(0.001,3,2000); % valores de beta3

cita1 = 0.1; %tasa de deteccion superpropagadores -> supuestos
cita2 = 0.1; %tasa de deteccion asintomaticos -> supuestos
omega1 = 0.6; %probabilidad de convertirse en asintomatico
omega2 = 0.01; %probabilidad de convertirser en infeccioso detectado 
omega3 = 1-0.6-0.01; %probabildiad de convertirse en superpropagador (10 % gente que no hace caso)

N=10000;

%variables
% S - susceptibles 
% E - infectados no infecciosos  
% I1 - infectados asintomaticos 
% I2 - infectados diagnosticados 
% I3 - infectados superpropagadores
% R - recuperados


%numero basico reproductivo
R01 = (alfa.*beta1.*omega1)./((cita2+gamma+mu).*(mu+alfa));
R02 = (beta2).*alfa*(  omega2+  (cita1.*omega3)./(cita1+mu+gamma) +   (cita2.*omega1)/(cita2+mu+gamma)  )./((mu+gamma).*(mu+alfa));
R03 = (alfa.*beta3.*omega3)./((cita1+gamma+mu).*(mu+alfa));
R0 = R01 + R02 + R03;

%R0=beta./(mu+gamma); % número básico reproductivo

EE1 = R01.*(R0-1).*mu.*N./(R0.*beta1); % punto de equilibrio endémico
EE2 = R02.*(R0-1).*mu.*N./(R0.*beta2);
EE3 = R03.*(R0-1).*mu.*N./(R0.*beta3);
% gráfica de bifurcación
figure
plot(R0,EE1);
hold on
plot(R0,EE2);
hold on 
plot(R0,EE3);

xlabel('R_0')
ylabel('I1^*, I2^*, I3^*')
legend('Asintomatico','Diagnosticado','Superpropagador')
axis([0 4 0 4]);