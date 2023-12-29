function [w]=Biot_Savart(p_1,p_2,p_c,g,toll)


% calcola la velocità indotta nel punto p_c da un vortice esteso da p_1 a
% p_2

% INPUT

% p_1 : 3 coord punto iniziale filamento vorticoso 
% p_2 : 3 coord punto finale filamento vorticoso
% p_c : 3 coord punto controllo
% g circolazione (1 se omessa)
% toll : distanza sotto la quale il vortice non induce velocità 

% OUTPUT
% w : velocità indotta nel punto p_c dal segmento vorticoso di estremi p_1
% e p_2 



r1=p_c-p_1;
r2=p_c-p_2;
r0=p_2-p_1;

dist= norm(cross(r1,r2))/norm(r0);

if dist < toll
  w= 0*r0;
else
  w = ((g/(4*pi))  *  cross(r1,r2)/(norm(cross(r1,r2)))^2)  *  dot(r0,((r1/norm(r1))-(r2/norm(r2))));
end 