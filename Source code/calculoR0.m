%Esto calcula el jacobiano (matriz V)
%X=I1
%Y=I2
%Z=I3

syms S E X Y Z alfa omegaU omegaD omegaT citaU citaD gamma mu N bU bD bT

F=[0 bU bD bT; 0 0 0 0; 0 0 0 0; 0 0 0 0]; %Jacobiana evaluada en el DFE
V=jacobian([(alfa*omegaU*E+alfa*omegaD*E+alfa*omegaT*E + mu*E),(-alfa*omegaU*E+X*(citaD+mu+gamma)),(-alfa*omegaD*E-citaU*Z-citaD*X+(mu+gamma)*Y),(-alfa*omegaT*E+Z*(citaU+mu+gamma))],[E,X,Y,Z]);
Vinversa=inv(V); %invertir matriz
F*Vinversa %el producto de matrices para encontrar el R0

