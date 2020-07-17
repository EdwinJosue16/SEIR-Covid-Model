% modelo superpropagadores

%function seriesDeTiempo(tiempoFinal)
tiempoFinal=365;
%parametros
mu=1/(85*365); %mortalidad natural
alfa = 1/7; %periodo de latencia -> 7 dias
gamma = 1/17; %periodo de recuperacion -> 17 dias
beta1 = 0.26; %tasa de transmision asintomatico
beta2 = 0.21; %tasa de transmision infectado diagnosticado (normal)
beta3 = 0.31; %tasa de transmision super propagador 
cita1 = 0.2; %tasa de deteccion superpropagadores -> supuestos
cita2 = 0.05; %tasa de deteccion asintomaticos -> supuestos
omega1 = 0.2; %probabilidad de convertirse en asintomatico
omega2 = 0.7; %probabilidad de convertirser en infeccioso detectado 
omega3 = 0.1; %probabildiad de convertirse en superpropagador (10 % gente que no hace caso)


%variables
% S - susceptibles 
% E - infectados no infecciosos  
% I1 - infectados asintomaticos 
% I2 - infectados diagnosticados 
% I3 - infectados superpropagadores
% R - recuperados

%tiempo
tspan=[0 tiempoFinal];

%numero basico reproductivo
R01 = (alfa*beta1*omega1)/((cita2+gamma+mu)*(mu+alfa));
R02 = (beta2)*alfa*(  omega2+  (cita1*omega3)/(cita1+mu+gamma) +   (cita2*omega1)/(cita2+mu+gamma)  )/((mu+gamma)*(mu+alfa));
R03 = (alfa*beta3*omega3)/((cita1+gamma+mu)*(mu+alfa));
R0 = R01 + R02 + R03;

%condiciones iniciales
condInit=[9997,0,1,1,1,0];
N=sum(condInit); %poblacion total

%sistema de ecuaciones diferenciales 
f = @(t,x) [
                 mu*N - x(1)*(beta1*x(3)  + beta2*x(4) + beta3*x(5))/N - mu*x(1); 
                 x(1)*(beta1*x(3)  + beta2*x(4) + beta3*x(5))/N - (alfa+mu)*x(2); 
                 alfa*omega1*x(2) - (cita2 + mu + gamma)*x(3); 
                 alfa*omega2*x(2) + cita1 * x(5) + cita2*x(3) - (mu+gamma)*x(4);
                 alfa*omega3*x(2) - (cita1 + mu + gamma)*x(5); 
                 (x(3) + x(4) + x(5))*gamma - mu*x(6)
           ];

[t,xa] = ode45(f,tspan,condInit); % paquete de matlab que hace la integracion

%Tres tipos de propagadores gráficados a la vez

plot(t,xa(:,3),'g') %grafica 
hold on
plot(t,xa(:,4),'b') %grafica 
hold on
plot(t,xa(:,5),'r')
xlabel('tiempo (dias)')
ylabel('Infectados')
legend('Asintomatico','Diagnosticado','Superpropagador')
title(['R0 = ' num2str(round(R0,2))])

