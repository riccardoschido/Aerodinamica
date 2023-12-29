clc
close all
clear

%% dati in ingresso al problema

[x1,y1]=importvector("./naca_7512_vector"); %profilo 1
[x2,y2]=importvector("./naca_6512_vector"); %profilo 2
U=1;  %velocità indisturbata
a1=2; %angolo di incidenza in gradi profilo 1
a2=-4; %angolo di incidenza in gradi profilo 2

x12=1.3; %distanza lungo x profilo 1-profilo 2
% x12=2; %distanza lungo x profilo 1-profilo 2

y12=0.1; %distanza lungo y profilo 1 - profilo 2

%% geometrizzazione e rotazione pannelli profilo 1

x1=flipud(x1); 
y1=flipud(y1);
N1=length(x1)-1; % numero di pannelli selezionato
alpha1=a1*pi/180;
chord1 = x1(1)-x1(N1/2); % lunghezza della corda del secondo profilo
ca1= chord1/4;
vec_ca1 = [ca1 0]; % vettore origine del centro aerodinamico del secondo 
                        % profilo 
R1=[cos(alpha1), -sin(alpha1);
    sin(alpha1), cos(alpha1)]; 
%matrice di rotazione per considerare profilo...
... inclinato rispetto orizzontale. (risultati...
...confrontati con xfoil e corretti)    

V1=[]; 
%rotazione dei nodi

for i=1:N1+1
    vet=[x1(i)-vec_ca1(1); y1(i)-vec_ca1(2)];
  vett1=R1'*vet;
  V1=[V1,vett1];
end

%vettori punti ruotati di angolo incidenza
x1=V1(1,:)'; 
y1=V1(2,:)';

%% geometrizzazione e rotazione pannelli profilo 2

x2=flipud(x2); 
y2=flipud(y2);
points2 = [x2 y2];

N2=length(x2)-1; % numero di pannelli selezionato
alpha2=a2*pi/180;
chord2 = x2(1)-x2(N2/2); % lunghezza della corda del secondo profilo
ca2= chord2/4;
vec_ca2 = [ca2 0]; % vettore origine del centro aerodinamico del secondo 
                        % profilo 

R=[cos(alpha2), -sin(alpha2);
    sin(alpha2), cos(alpha2)]; 
%matrice di rotazione per considerare profilo...
... inclinato rispetto orizzontale. (risultati...
...confrontati con xfoil e corretti)    

V=[];
%rotazione dei nodi

for i=1:N2+1
  vet=[x2(i)-vec_ca2(1); y2(i)-vec_ca2(2)];
  vett1=R'*vet;
  V=[V,vett1];
end
%vettori punti ruotati di angolo incidenza
x2=x12+V(1,:)'; 
y2=y12+V(2,:)';

%% definizione vettore inclinazione pannelli profilo 1
theta1=zeros(N1,1);
for i=1:N1
    Y1=y1(i+1)-y1(i);
    X1=x1(i+1)-x1(i);
    theta1(i)=pi+atan2(Y1,X1);   
end

%% definizione vettore inclinazione pannelli profilo 2
theta2=zeros(N2,1);
for i=1:N2
    Y2=y2(i+1)-y2(i);
    X2=x2(i+1)-x2(i);
    theta2(i)=pi+atan2(Y2,X2);   
end

%% definizione normali e tangenti profilo 1

nor1=zeros(N1,2);
tau1=zeros(N1,2);

for i=1:N1
nor1(i,1)=-sin(theta1(i));
nor1(i,2)=cos(theta1(i));
tau1(i,1)=cos(theta1(i));
tau1(i,2)=sin(theta1(i));
end

%% definizione normali e tangenti profilo 2

nor2=zeros(N2,2);
tau2=zeros(N2,2);

for i=1:N2
nor2(i,1)=-sin(theta2(i));
nor2(i,2)=cos(theta2(i));
tau2(i,1)=cos(theta2(i));
tau2(i,2)=sin(theta2(i));
end

%% definizione centro pannello profilo 1

centro1=zeros(2,N1);

for i=1:N1
    centro1(1,i)=(x1(i+1)+x1(i))/2;
    centro1(2,i)=(y1(i+1)+y1(i))/2;
end

%% definizione centro pannello profilo 2

centro2=zeros(2,N2);

for i=1:N2
    centro2(1,i)=(x2(i+1)+x2(i))/2;
    centro2(2,i)=(y2(i+1)+y2(i))/2;
end

%% definizione estremi pannelli profilo 1
 %estremo_i_j  i=estremo j=pannello
estremo_1_1=zeros(2,N1);
estremo_2_1=zeros(2,N1);

for i=1:N1
    estremo_1_1(1,i)=x1(i);
    estremo_1_1(2,i)=y1(i);
    estremo_2_1(1,i)=x1(i+1);
    estremo_2_1(2,i)=y1(i+1);
end

%% definizione estremi pannelli profilo 2

estremo_1_2=zeros(2,N2);
estremo_2_2=zeros(2,N2);

for i=1:N2
    estremo_1_2(1,i)=x2(i);
    estremo_1_2(2,i)=y2(i);
    estremo_2_2(1,i)=x2(i+1);
    estremo_2_2(2,i)=y2(i+1);
end

%% componenti velocità asintotica

U_inf(1)=U;  %velocità in direzione x
U_inf(2)=0;  %non ho componente y perche l'incidenza la...
             ...definisco sul profilo

%% condizioni al contorno

b=zeros(N1+N2+2,1);

%condizione di non penetrabilità

%profilo 1

for j=1:N1
    b(j)=-U_inf(1)*nor1(j,1)-U_inf(2)*nor1(j,2);
end

%profilo 2
b(N1+1:N2+N1)=-U_inf(1)*nor2(:,1)-U_inf(2)*nor2(:,2);

%condizione kutta 

%profilo 1
b(end-1,1)=-U_inf(1)*tau1(1,1)-U_inf(1)*tau1(N1,1)...
    -U_inf(2)*tau1(1,2)-U_inf(2)*tau1(N1,2);
%profilo 2
b(end,1)=-U_inf(1)*tau2(1,1)-U_inf(1)*tau2(N2,1)...
    -U_inf(2)*tau2(1,2)-U_inf(2)*tau2(N2,2);
 
%% scrittura matrice A 
A=zeros(N1+N2+2,N1+N2+2);

%sottomatrice A_11

A_11=zeros(N1,N1);

for i=1:N1
    for j=1:N1
        M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
    M_g_l=M_l_g';
   u_s = ViSorgente(centro1(:,i), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
    A_11(i,j)=u_s(1)*nor1(i,1)+u_s(2)*nor1(i,2);
    
    end
end

%sotto matrice A_12

A_12=zeros(N1,N2); %sottomatrice di induzione su pannello j...
                  ... del profilo 2 da centro su profilo 1
for i=1:N1
    for j=1:N2
         M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
         M_g_l=M_l_g';
         u_s = ViSorgente(centro1(:,i), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
          A_12(i,j)=u_s(1)*nor1(i,1)+u_s(2)*nor1(i,2);
    end
end

%sottomatrice A_21

A_21=zeros(N2,N1);
for i=1:N2
    for j=1:N1
        M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
         M_g_l=M_l_g';
         u_s = ViSorgente(centro2(:,i), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
          A_21(i,j)=u_s(1)*nor2(i,1)+u_s(2)*nor2(i,2);
    end
end
 
%sottomatrice A_22

A_22=zeros(N2,N2);

for i=1:N2
    for j=1:N2
        M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
         M_g_l=M_l_g';
         u_s = ViSorgente(centro2(:,i), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
          A_22(i,j)=u_s(1)*nor2(i,1)+u_s(2)*nor2(i,2);
    end
end

%sottomatrice a_11

a_11=zeros(N1,1);

for i=1:N1
    for j=1:N1
     M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
         M_g_l=M_l_g';
u_v= ViVortice(centro1(:,i), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
a_11(i,1)=a_11(i,1)+u_v(1)*nor1(i,1)+u_v(2)*nor1(i,2);
    end
end

%sottomatrice a_12

a_12=zeros(N1,1);

for i=1:N1
    for j=1:N2
     M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
         M_g_l=M_l_g';
         u_v= ViVortice(centro1(:,i), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
a_12(i,1)=a_12(i,1)+u_v(1)*nor1(i,1)+u_v(2)*nor1(i,2);
    end
end

%sottomatrice a_21

a_21=zeros(N2,1);

for i=1:N2
    for j=1:N1
        M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
         M_g_l=M_l_g';
         u_v= ViVortice(centro2(:,i), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
a_21(i,1)=a_21(i,1)+u_v(1)*nor2(i,1)+u_v(2)*nor2(i,2);
    end
end

%sottomatrice a_22

a_22=zeros(N2,1);

for i=1:N2
    for j=1:N2
        M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
         M_g_l=M_l_g';
         u_v= ViVortice(centro2(:,i), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
a_22(i,1)=a_22(i,1)+u_v(1)*nor2(i,1)+u_v(2)*nor2(i,2);
    end
end
  
%sottomatrice c_11_s
 
 c_11_s=zeros(1,N1);
 
 for j=1:N1
      M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
         M_g_l=M_l_g';
         u_s_1 = ViSorgente(centro1(:,1), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
          c_11_s(j)=u_s_1(1)*tau1(1,1)+u_s_1(2)*tau1(1,2);
           u_s_n = ViSorgente(centro1(:,N1), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
          c_11_s(j)=c_11_s(j)+u_s_n(1)*tau1(N1,1)+u_s_n(2)*tau1(N1,2);
 end
 
 %sottomatrice c_12_s
 
 c_12_s=zeros(1,N2);
 
 for j=1:N2
      M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
         M_g_l=M_l_g';
         u_s_1 = ViSorgente(centro1(:,1), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
          c_12_s(j)=u_s_1(1)*tau1(1,1)+u_s_1(2)*tau1(1,2);
           u_s_n = ViSorgente(centro1(:,N1), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
          c_12_s(j)=c_12_s(j)+u_s_n(1)*tau1(N1,1)+u_s_n(2)*tau1(N1,2);
 end
 
 %sottomatrice c_21_s
 
 c_21_s=zeros(1,N1);
 
 for j=1:N1
     M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
         M_g_l=M_l_g';
         u_s_1 = ViSorgente(centro2(:,1), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
          c_21_s(j)=u_s_1(1)*tau2(1,1)+u_s_1(2)*tau2(1,2);
           u_s_n = ViSorgente(centro2(:,N2), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
          c_21_s(j)=c_21_s(j)+u_s_n(1)*tau2(N2,1)+u_s_n(2)*tau2(N2,2);
 end
 
 %sottomatrice c_22_s
 
 c_22_s=zeros(1,N1);
 
 for j=1:N2
     M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
         M_g_l=M_l_g';
         u_s_1 = ViSorgente(centro2(:,1), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
          c_22_s(j)=u_s_1(1)*tau2(1,1)+u_s_1(2)*tau2(1,2);
           u_s_n = ViSorgente(centro2(:,N2), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
          c_22_s(j)=c_22_s(j)+u_s_n(1)*tau2(N2,1)+u_s_n(2)*tau2(N2,2);
 end

 %sottomatrice c_11_v
 for j=1:N1  
     M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
    M_g_l=M_l_g';
   u_v_1= ViVortice(centro1(:,1), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
   A(end-1,end-1)=A(end-1,end-1)+u_v_1(1)*tau1(1,1)+u_v_1(2)*tau1(1,2);
   u_v_n= ViVortice(centro1(:,N1), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
     A(end-1,end-1)=A(end-1,end-1)+u_v_n(1)*tau1(N1,1)+u_v_n(2)*tau1(N1,2);  
 end
    
 %sottomatrice c_12_v
 
 for j=1:N2 
      M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
    M_g_l=M_l_g';
   u_v_1= ViVortice(centro1(:,1), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
   A(end-1,end)=A(end-1,end)+u_v_1(1)*tau1(1,1)+u_v_1(2)*tau1(1,2);
   u_v_n= ViVortice(centro1(:,N1), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
     A(end-1,end)=A(end-1,end)+u_v_n(1)*tau1(N1,1)+u_v_n(2)*tau1(N1,2);  
 end
 
 %sottomatrice c_21_v
 
 for j=1:N1 
        M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
    M_g_l=M_l_g';
   u_v_1= ViVortice(centro2(:,1), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
   A(end,end-1)=A(end,end-1)+u_v_1(1)*tau2(1,1)+u_v_1(2)*tau2(1,2);
   u_v_n= ViVortice(centro2(:,N2), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
     A(end,end-1)=A(end,end-1)+u_v_n(1)*tau2(N2,1)+u_v_n(2)*tau2(N2,2);  
 end
 
 %sottomatrice c_22_v
 
 for j=1:N2 
     M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
    M_g_l=M_l_g';
   u_v_1= ViVortice(centro2(:,1), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
   A(end,end)=A(end,end)+u_v_1(1)*tau2(1,1)+u_v_1(2)*tau2(1,2);
   u_v_n= ViVortice(centro2(:,N2), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
     A(end,end)=A(end,end)+u_v_n(1)*tau2(N2,1)+u_v_n(2)*tau2(N2,2);  
 end
  
 %costruzione matrice 
 
 A(1:N1,1:N1)=A_11;
 A(1:N1,N1+1:N2+N1)=A_12;
 A(N1+1:N2+N1,1:N1)=A_21;
 A(N1+1:N2+N1,N1+1:N2+N1)=A_22;
 A(1:N1,end-1)=a_11;
 A(1:N1,end)=a_12;
 A(N1+1:N2+N1,end-1)=a_21;
 A(N1+1:N2+N1,end)=a_22;
 A(end-1,1:N1)=c_11_s;
 A(end-1,N1+1:N1+N2)=c_12_s;
 A(end,1:N1)=c_21_s;
 A(end,N1+1:N1+N2)=c_22_s;
 
 %risoluzione sistema Ax=b
 solution=linsolve(A,b);
 q1=solution(1:N1); %sorgenti sugli N pannelli profilo 1
 q2=solution(N1+1:N1+N2);
 gamma1=solution(end-1); %intensità del vortice
 gamma2=solution(end);
 
 %% calcolo coefficienti
 
 %lunghezze pannelli
 
 length_1=zeros(N1,1);
for i=1:N1
    length_1(i)=sqrt((y1(i+1)-y1(i))^2+(x1(i+1)-x1(i))^2);
    
    
end
length_2=zeros(N2,1);
for i=1:N2
    length_2(i)=sqrt((y2(i+1)-y2(i))^2+(x2(i+1)-x2(i))^2);
    
end


%coeff portanza Kutta

GAMMA_1=sum(length_1,'all')*gamma1;
GAMMA_2=sum(length_2,'all')*gamma2;
C_l_1=-2*GAMMA_1/U;
C_l_2=-2*GAMMA_2/U;


%calcolo della portanza con andamento cp...
%calcolo portanza con cp profilo 1


 u_y_1=zeros(N1,1);
u_x_1=zeros(N1,1);

for i=1:N1
   u_x_1(i)=U_inf(1);
   u_y_1(i)=U_inf(2);
   for j=1:N1
          M_l_g=[cos(theta1(j)) -sin(theta1(j));
               sin(theta1(j))   cos(theta1(j))];
    M_g_l=M_l_g';
    u_s_1 = ViSorgente(centro1(:,i), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
   u_v_1= ViVortice(centro1(:,i), estremo_1_1(:,j), estremo_2_1(:,j), M_l_g, M_g_l);
   u_x_1(i)=u_x_1(i)+u_s_1(1)*q1(j)+u_v_1(1)*gamma1;
   u_y_1(i)=u_y_1(i)+u_s_1(2)*q1(j)+u_v_1(2)*gamma1;
end
u_x_1(i)=u_x_1(i);
u_y_1(i)=u_y_1(i);
for k=1:N2
    
       M_l_g=[cos(theta2(k)) -sin(theta2(k));
               sin(theta2(k))   cos(theta2(k))];
    M_g_l=M_l_g';
   u_s_2 = ViSorgente(centro1(:,i), estremo_1_2(:,k), estremo_2_2(:,k), M_l_g, M_g_l);
   u_v_2= ViVortice(centro1(:,i), estremo_1_2(:,k), estremo_2_2(:,k), M_l_g, M_g_l);
   u_x_1(i)=u_x_1(i)+u_s_2(1)*q2(k)+u_v_2(1)*gamma2;
   u_y_1(i)=u_y_1(i)+u_s_2(2)*q2(k)+u_v_2(2)*gamma2;
  
   end
end

u_tang_1=u_x_1.*tau1(:,1)+u_y_1.*tau1(:,2);
modu=sqrt(U_inf(1)^2+U_inf(2)^2);
Cp_1=1-(((u_tang_1).^2)/(modu)^2);
cl_1 =sum(Cp_1.*length_1.*cos(theta1))

%calcolo portanza con cp profilo 2

u_y_2=zeros(N2,1);
u_x_2=zeros(N2,1);
       
for i=1:N2
   u_x_2(i)=U_inf(1);
   u_y_2(i)=U_inf(2);
   for j=1:N2
          M_l_g=[cos(theta2(j)) -sin(theta2(j));
               sin(theta2(j))   cos(theta2(j))];
    M_g_l=M_l_g';
    u_s_2 = ViSorgente(centro2(:,i), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
   u_v_2= ViVortice(centro2(:,i), estremo_1_2(:,j), estremo_2_2(:,j), M_l_g, M_g_l);
   u_x_2(i)=u_x_2(i)+u_s_2(1)*q2(j)+u_v_2(1)*gamma2;
   u_y_2(i)=u_y_2(i)+u_s_2(2)*q2(j)+u_v_2(2)*gamma2;
   end
u_x_2(i)=u_x_2(i);
u_y_2(i)=u_y_2(i);
   for k=1:N1
     
       M_l_g=[cos(theta1(k)) -sin(theta1(k));
               sin(theta1(k))   cos(theta1(k))];
    M_g_l=M_l_g';
    u_s_1 = ViSorgente(centro2(:,i), estremo_1_1(:,k), estremo_2_1(:,k), M_l_g, M_g_l);
   u_v_1= ViVortice(centro2(:,i), estremo_1_1(:,k), estremo_2_1(:,k), M_l_g, M_g_l);
   u_x_2(i)=u_x_2(i)+u_s_1(1)*q1(k)+u_v_1(1)*gamma1;
   u_y_2(i)=u_y_2(i)+u_s_1(2)*q1(k)+u_v_1(2)*gamma1;
   
   end
end

u_tang_2=u_x_2.*tau2(:,1)+u_y_2.*tau2(:,2);
modu=sqrt(U_inf(1)^2+U_inf(2)^2);
Cp_2=1-(((u_tang_2).^2)/(modu)^2);
cl_2 =sum(Cp_2.*length_2.*cos(theta2))

%% grafici
% configurazione profilo pannellizati
figure(1)
axis equal
grid on
grid minor
hold on
plot(estremo_1_1(1,:),estremo_1_1(2,:),'|b', 'LineWidth', 0.6, 'MarkerSize', 8)
plot(estremo_2_1(1,:),estremo_2_1(2,:),'|b', 'LineWidth', 0.6, 'MarkerSize', 8)
plot(centro1(1,:),centro1(2,:),'+r', 'LineWidth', 0.6, 'MarkerSize', 8)
plot(x1,y1,'k', 'LineWidth', 3)

plot(estremo_1_2(1,:),estremo_1_2(2,:),'|b', 'LineWidth', 0.6, 'MarkerSize', 8)
plot(estremo_2_2(1,:),estremo_2_2(2,:),'|b', 'LineWidth', 0.6, 'MarkerSize', 8)
plot(centro2(1,:),centro2(2,:),'+r', 'LineWidth', 0.6, 'MarkerSize', 8)
plot(x2,y2,'k', 'LineWidth', 3)
xlabel('x')
ylabel('y')
axis([-0.3 2.1 -0.2 0.3])
set(gca,'FontWeight','bold', 'FontSize', 26)

% grafico cp
figure(2)
plot(centro1(1,1:N1/2-3)+vec_ca2(1),-Cp_1(1:N1/2-3),'LineWidth',3,'Color','r')
hold on
plot(centro1(1,N1/2-3:end)+vec_ca2(1),-Cp_1(N1/2-3:end),'LineWidth',3,'Color','b')
hold on
centro_2=zeros(2,N2);
centro_2(1,:)=centro2(1,:)-x12 + vec_ca1(1);
centro_2(2,:)=centro2(2,:)-y12;
plot(centro_2(1,1:N2/2+1),-Cp_2(1:N2/2+1), '-.', 'LineWidth',3, 'color', [0.7 0.4470 0.7410])
hold on
plot(centro_2(1,N2/2+1:end),-Cp_2(N2/2+1:end), '-.', 'LineWidth', 3, 'color', [0.4 0.3 0.3])
legend(' Ventre ala principale ',' Dorso  ala principale ',' Ventre ala di coda ',' Dorso ala di coda')
grid on 
grid minor
xlabel('x')
ylabel('-C_p')
set(gca,'FontWeight','bold', 'FontSize', 26)
