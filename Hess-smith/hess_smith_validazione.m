clc
close all
clear

%% dati in ingresso problema

[x1,y1]=importvector("./naca_0012_101"); %import da xfoil profilo naca
U=1;
a=2; %angolo di incidenza in gradi

%% geometrizzazione del profilo

x1=flipud(x1);
y1=flipud(y1);
N=length(x1)-1; % numero di pannelli selezionato
alpha=a*pi/180;

R=[cos(alpha), -sin(alpha);
    sin(alpha), cos(alpha)]; 
%matrice di rotazione per considerare profilo...
... inclinato rispetto orizzontale. (risultati...
...confrontati con xfoil e corretti)    
V=[];
for i=1:N+1
   
    vet=[x1(i); y1(i)];
    
  vett1=R'*vet;
  V=[V,vett1];
end
%vettori punti ruotati di angolo incidenza
x=V(1,:)'; 
y=V(2,:)';


%% definizione vettore inclinazione pannelli
theta=zeros(N,1);
for i=1:N
    Y=y(i+1)-y(i);
    X=x(i+1)-x(i);
    theta(i)=pi+atan2(Y,X);   
end
%% definizione normali e tangenti

nor=zeros(N,2);
tau=zeros(N,2);

for i=1:N
nor(i,1)=-sin(theta(i));
nor(i,2)=cos(theta(i));
tau(i,1)=cos(theta(i));
tau(i,2)=sin(theta(i));
end

%% definizione centro pannello

centro=zeros(2,N);

for i=1:N
    centro(1,i)=(x(i+1)+x(i))/2;
    centro(2,i)=(y(i+1)+y(i))/2;
end

%% definizione estremi pannelli

estremo_1=zeros(2,N);
estremo_2=zeros(2,N);

for i=1:N
    estremo_1(1,i)=x(i);
    estremo_1(2,i)=y(i);
    estremo_2(1,i)=x(i+1);
    estremo_2(2,i)=y(i+1);
end

%% componenti velocità asintotica
U_inf(1)=U;
U_inf(2)=0;

%% condizioni al contorno
b=zeros(N+1,1);

%cond non penetrabilità
for j=1:N
    b(j)=-U_inf(1)*nor(j,1)-U_inf(2)*nor(j,2);
end
%cond di Kutta
b(N+1)=-U_inf(1)*tau(1,1)-U_inf(1)*tau(N,1)...
    -U_inf(2)*tau(1,2)-U_inf(2)*tau(N,2);

%% matrice A

A=zeros(N+1,N+1);

for i=1:N
    for j=1:N
    
    M_l_g=[cos(theta(j)) -sin(theta(j));
               sin(theta(j))   cos(theta(j))];
    M_g_l=M_l_g';
   u_s = ViSorgente(centro(:,i), estremo_1(:,j), estremo_2(:,j), M_l_g, M_g_l);
   u_v= ViVortice(centro(:,i), estremo_1(:,j), estremo_2(:,j), M_l_g, M_g_l);
    A(i,j)=u_s(1)*nor(i,1)+u_s(2)*nor(i,2);
    A(i,N+1)=A(i,N+1)+u_v(1)*nor(i,1)+u_v(2)*nor(i,2);
    end
   
end

    for j=1:N
        M_l_g=[cos(theta(j)) -sin(theta(j));
               sin(theta(j))   cos(theta(j))];
    M_g_l=M_l_g';
    u_s_1 = ViSorgente(centro(:,1), estremo_1(:,j), estremo_2(:,j), M_l_g, M_g_l);
   u_v_1= ViVortice(centro(:,1), estremo_1(:,j), estremo_2(:,j), M_l_g, M_g_l);
    A(N+1,j)=u_s_1(1)*tau(1,1)+u_s_1(2)*tau(1,2);
    A(N+1,N+1)=A(N+1,N+1)+u_v_1(1)*tau(1,1)+u_v_1(2)*tau(1,2);
    u_s_n = ViSorgente(centro(:,N), estremo_1(:,j), estremo_2(:,j), M_l_g, M_g_l);
   u_v_n= ViVortice(centro(:,N), estremo_1(:,j), estremo_2(:,j), M_l_g, M_g_l);
   A(N+1,j)=A(N+1,j)+u_s_n(1)*tau(N,1)+u_s_n(2)*tau(N,2);
   A(N+1,N+1)=A(N+1,N+1)+u_v_n(1)*tau(N,1)+u_v_n(2)*tau(N,2);  
    end
 solution=linsolve(A,b);
 q=solution(1:N); %sorgenti sugli N pannelli
 gamma=solution(N+1); %intensità del vortice
 
 %% calcolo dei coefficienti
u_y=zeros(N,1);
u_x=zeros(N,1);

for i=1:N
   u_x(i)=U_inf(1);
   u_y(i)=U_inf(2);
   for j=1:N
          M_l_g=[cos(theta(j)) -sin(theta(j));
               sin(theta(j))   cos(theta(j))];
    M_g_l=M_l_g';
    u_s = ViSorgente(centro(:,i), estremo_1(:,j), estremo_2(:,j), M_l_g, M_g_l);
   u_v= ViVortice(centro(:,i), estremo_1(:,j), estremo_2(:,j), M_l_g, M_g_l);
   u_x(i)=u_x(i)+u_s(1)*q(j)+u_v(1)*gamma;
   u_y(i)=u_y(i)+u_s(2)*q(j)+u_v(2)*gamma;
   end
end

if (max(u_x.*nor(:,1) + u_y.*nor(:,2))>10^(-14))
    disp('There is a bug in the program!')
    
end

u_tang=u_x.*tau(:,1)+u_y.*tau(:,2);
modu=sqrt(U_inf(1)^2+U_inf(2)^2);
Cp=1-(((u_tang).^2)/(modu)^2);

length=zeros(N,1);
for i=1:N
    length(i)=sqrt((y(i+1)-y(i))^2+(x(i+1)-x(i))^2);
    
end
GAMMA=sum(length,'all')*gamma;
C_l_k=-2*GAMMA/U; %Kutta zukowsky
cl = sum(Cp.*length.*cos(theta)) %pressure integration


%% plot grafico profilo pannellizzato

figure(1)
axis equal
grid on
grid minor
hold on
plot(estremo_1(1,:),estremo_1(2,:),'|r', 'LineWidth', 1.5, 'MarkerSize', 10)
plot(estremo_2(1,:),estremo_2(2,:),'|r', 'LineWidth', 1.5, 'MarkerSize', 10)
plot(centro(1,:),centro(2,:),'+b', 'LineWidth', 1.5, 'MarkerSize', 10)
plot(x,y,'k', 'LineWidth', 3)
axis ([0 1 -0.2 0.2])
xlabel('x')
ylabel('y')
set(gca,'FontWeight','bold', 'FontSize', 26)

%% plot grafico andamento CP
figure(2)
plot(centro(1,1:(N+1)/2),-Cp(1:(N+1)/2), 'LineWidth', 3, 'Color','r')
hold on
plot(centro(1,(N+1)/2-2:end),-Cp((N+1)/2-2:end), 'LineWidth', 3,'Color','b')
legend(' Ventre ',' Dorso ')
grid on
grid minor
xlabel('x')
ylabel('-C_p')
set(gca,'FontWeight','bold', 'FontSize', 26)


   
