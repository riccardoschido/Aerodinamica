clear
close all 
clc

%% dati wing

rootchord = 1.50;
tipchord = 1.20;
span = 8.6; % in pianta
sweep = 5/180*pi; % angolo con congiuntente 25% corda
dihedral = 0/180*pi ;

%% dati tail

rootchord_t = 0.70;
tipchord_t = 0.60;
span_t = 2; % in pianta
sweep_t = 0/180*pi; % angolo con congiuntente 25% corda
dihedral_t = 0/180*pi ;

%% discretizzazione wing

Nchord = 6;
Nspan = 22; % solo di metà ala 


%% discretizzazione tail

Nchord_t = 3;
Nspan_t = 22; % solo di metà ala 

%% altri dati

alfa = 2*pi/180;
beta = 0*pi/180;
alfa_t = -4*pi/180; 
mac = 0.5*(rootchord+tipchord);
dx = 3; % cambiare per ottimizzazione, adimensionalizzato
dz_0 =0.5;
dz = dz_0 - mac*sin(alfa) +0.1;

%% adimensionalizzazione secondo la corda media aerodianamica wing

dz_0 = dz_0/mac;
dz = dz/mac;
span = span/mac;
rootchord = rootchord/mac;
tipchord = tipchord/mac;
s = 0.5*(rootchord+tipchord)*span; % superficie in pianta adimensionalizzata
coord_polo = [0.25*rootchord, 0, 0]';

%% adimensionalizzazione secondo la corda media aerodianamica tail

mac_t = 0.5*(rootchord_t+tipchord_t);
span_t = span_t/mac;
rootchord_t = rootchord_t/mac;
tipchord_t = tipchord_t/mac;
s_t = 0.5*(rootchord_t+tipchord_t)*span_t; % superficie in pianta adimensionalizzata
coord_polo_t = [0.25*rootchord_t, 0, 0]';

%% dati velocità

v0 = 80000; % mm/s
v = v0/mac;
v_inf = [v 0 0];
rho = 1;

%% generazione geometria wing

[X,Y,Z,X_V,Y_V,Z_V,coor_C,coor_V_sn,coor_V_dx,N] = wing_geo_rotaz(rootchord,tipchord,span,sweep,dihedral ,Nchord, Nspan, alfa, beta);

Z = Z + dz_0; 
Z_V = Z_V + dz_0;

coor_C(3,:) = coor_C(3,:) + dz_0;

coor_V_sn(3,:) = coor_V_sn(3,:) + dz_0;
coor_V_dx(3,:) = coor_V_dx(3,:) + dz_0;

% generazione geometria wing specchio

X_s = X;
Z_s = -Z;
Y_s = Y;
X_V_s = X_V;
Z_V_s = -Z_V;
Y_V_s = Y_V;

coor_C_s(1,:) = coor_C(1,:);
coor_C_s(3,:) = -coor_C(3,:);
coor_C_s(2,:) = coor_C(2,:);

coor_V_sn_s(1,:) = coor_V_sn(1,:);
coor_V_sn_s(3,:) = -coor_V_sn(3,:);
coor_V_sn_s(2,:) = coor_V_sn(2,:);
coor_V_dx_s(1,:) = coor_V_dx(1,:);
coor_V_dx_s(3,:) = -coor_V_dx(3,:);
coor_V_dx_s(2,:) = coor_V_dx(2,:);

%% generazione geometria tail

[X_t,Y_t,Z_t,X_V_t,Y_V_t,Z_V_t,coor_C_t,coor_V_sn_t,coor_V_dx_t,N_t] = wing_geo_rotaz(rootchord_t,tipchord_t,span_t,sweep_t,dihedral_t ,Nchord_t, Nspan_t, alfa_t, beta);

X_t = X_t + dx;
Z_t = Z_t + dz;
X_V_t = X_V_t + dx;
Z_V_t = Z_V_t + dz;

coor_C_t(1,:) = coor_C_t(1,:) + dx;
coor_C_t(3,:) = coor_C_t(3,:) + dz ;

coor_V_sn_t(1,:) = coor_V_sn_t(1,:) + dx;
coor_V_sn_t(3,:) = coor_V_sn_t(3,:) + dz;
coor_V_dx_t(1,:) = coor_V_dx_t(1,:) + dx;
coor_V_dx_t(3,:) = coor_V_dx_t(3,:) + dz;

% generazione geometria tail specchio

X_t_s = X_t;
Z_t_s = -Z_t;
Y_t_s = Y_t;
X_V_t_s = X_V_t;
Z_V_t_s = -Z_V_t;
Y_V_t_s = Y_V_t;

coor_C_t_s(1,:) = coor_C_t(1,:);
coor_C_t_s(3,:) = -coor_C_t(3,:);
coor_C_t_s(2,:) = coor_C_t(2,:);

coor_V_sn_t_s(1,:) = coor_V_sn_t(1,:);
coor_V_sn_t_s(3,:) = -coor_V_sn_t(3,:);
coor_V_sn_t_s(2,:) = coor_V_sn_t(2,:);
coor_V_dx_t_s(1,:) = coor_V_dx_t(1,:);
coor_V_dx_t_s(3,:) = -coor_V_dx_t(3,:);
coor_V_dx_t_s(2,:) = coor_V_dx_t(2,:);


%% tolleranza Biot-Savart

toll = 10^-4;

%% matrice A delle incidenze, b termine noto e soluzione

l_w = length(coor_C(1,:));
l_t = length(coor_C_t(1,:));

A = zeros(l_w + l_t);

%Aww

for ii = 1 : l_w
    for jj = 1 : l_w
      A(ii,jj) = (ind_vort(coor_V_sn(:,jj),coor_V_dx(:,jj),coor_C(:,ii),1,toll))'*N(:,ii)-...
          (ind_vort(coor_V_sn_s(:,jj),coor_V_dx_s(:,jj),coor_C(:,ii),1,toll))'*N(:,ii);
    end
end

%Atw

for ii = 1 : l_w
    for jj = (l_w + 1) : (l_w + l_t)
      A(ii,jj) = (ind_vort(coor_V_sn_t(:,jj - l_w),coor_V_dx_t(:,jj - l_w),coor_C(:,ii),1,toll))'*N(:,ii)-...
          (ind_vort(coor_V_sn_t_s(:,jj - l_w),coor_V_dx_t_s(:,jj - l_w),coor_C(:,ii),1,toll))'*N(:,ii);
    end
end

%Awt

for ii = (l_w + 1) : (l_w + l_t)
    for jj = 1 : l_w
      A(ii,jj) = (ind_vort(coor_V_sn(:,jj),coor_V_dx(:,jj),coor_C_t(:,ii - l_w),1,toll))'*N_t(:,ii - l_w)-...
          (ind_vort(coor_V_sn_s(:,jj),coor_V_dx_s(:,jj),coor_C_t(:,ii - l_w),1,toll))'*N_t(:,ii - l_w);
    end
end

%Att

for ii = (l_w + 1) : (l_w + l_t)
    for jj = (l_w + 1) : (l_w + l_t)
      A(ii,jj) = (ind_vort(coor_V_sn_t(:,jj - l_w),coor_V_dx_t(:,jj - l_w),coor_C_t(:,ii - l_w),1,toll))'*N_t(:,ii - l_w)-...
          (ind_vort(coor_V_sn_t_s(:,jj - l_w),coor_V_dx_t_s(:,jj - l_w),coor_C_t(:,ii - l_w),1,toll))'*N_t(:,ii - l_w);
    end
end

  
b = -(v_inf*N)';
b_t = -(v_inf*N_t)';

b_tot =[b;b_t];
gamma = A\b_tot;


gamma_w = gamma(1 : l_w);
gamma_t = gamma(l_w + 1 : l_w + l_t);
GAMMA_w = (reshape(gamma_w, 2*Nspan, Nchord))';
GAMMA_t = (reshape(gamma_t, 2*Nspan_t, Nchord_t))';
gamma_span_w = sum(GAMMA_w);
gamma_span_t = sum(GAMMA_t);
    
%% Calcolo della portanza

dy_w = abs(Y(1,1)-Y(1,2));
dy_t = abs(Y_t(1,1)-Y_t(1,2));

portanza1d_w = rho*v*gamma_w*cos(dihedral);
portanza1d_t = rho*v*gamma_t*cos(dihedral);
l_1d_w = rho*v*GAMMA_w*cos(dihedral);
l_1d_t = rho*v*GAMMA_t*cos(dihedral);
l_2d_w = rho*v*gamma_span_w*cos(dihedral);
l_2d_t = rho*v*gamma_span_t*cos(dihedral_t);
L_w = dy_w*sum(l_2d_w);
L_t = dy_t*sum(l_2d_t);

cl_w = 2*L_w/(rho*s*v^2);
cl_t = 2*L_t/(rho*s*v^2);
cl_2d_w = 2*l_2d_w/(rho*s*v^2);

%% Calcolo dell'incidenza indotta e resistenza wing

gamma_scia_w = [gamma_span_w(1),-gamma_span_w(1:(end-1))+gamma_span_w(2:end),-gamma_span_w(end)];

y_w = Y(1,:);
z_w = Z(1,:);
y_c_w = y_w(1:end-1)+0.5*dy_w;
z_c_w = 0.5*(z_w(1:end-1)+z_w(2:end));

v_ind_w = zeros(size(y_c_w));

for ii = 1:length(y_c_w)
   for jj = 1:length(gamma_scia_w)
   dist1 = sqrt((y_c_w(ii)-y_w(jj))^2+(z_c_w(ii)-z_w(jj))^2);        
   v_ind_w(ii) = v_ind_w(ii) + gamma_scia_w(jj)*(-y_c_w(ii)+y_w(jj))/(4*pi*dist1^2);
   end
end

alfa_ind_w = -atan(v_ind_w/v); % vettore incidenze indotte
d_2d_w = l_2d_w .* alfa_ind_w; 
D_w = sum (d_2d_w)*dy_w;

%% calcolo dell'incidenza indotta e resistenza tail

gamma_scia_t = [gamma_span_t(1),-gamma_span_t(1:(end-1))+gamma_span_t(2:end),-gamma_span_t(end)];

y_t = Y_t(1,:);
z_t = Z_t(1,:);
y_c_t = y_t(1:end-1)+0.5*dy_t;
z_c_t = 0.5*(z_t(1:end-1)+z_t(2:end));

v_ind_t = zeros(size(y_c_t));

for ii = 1:length(y_c_t)
   for jj = 1:length(gamma_scia_t)
   dist2 = sqrt((y_c_t(ii)-y_t(jj))^2+(z_c_t(ii)-z_t(jj))^2);        
   v_ind_t(ii) = v_ind_t(ii) + gamma_scia_t(jj)*(-y_c_t(ii)+y_t(jj))/(4*pi*dist2^2);
   end
end

alfa_ind_t = -atan(v_ind_t/v); % vettore incidenze indotte
d_2d_t = l_2d_t .* alfa_ind_t; 
D_t = sum (d_2d_t)*dy_t;

%% plot
figure(1)
M=mesh(X,Y,Z);
set(M,'facealpha',0)
set(M,'edgecolor',[.2 .4 .9])
hold on
plot3(coor_C(1,:),coor_C(2,:),coor_C(3,:),'og')
plot3(coor_V_sn(1,:),coor_V_sn(2,:),coor_V_sn(3,:),'^r')
axis equal
hold on
M=mesh(X_t,Y_t,Z_t);
set(M,'facealpha',0)
set(M,'edgecolor',[.2 .4 .9])
hold on
plot3(coor_C_t(1,:),coor_C_t(2,:),coor_C_t(3,:),'og')
plot3(coor_V_sn_t(1,:),coor_V_sn_t(2,:),coor_V_sn_t(3,:),'^r')
axis equal
hold on
M=mesh(X_s,Y_s,Z_s);
set(M,'facealpha',0)
set(M,'edgecolor',[.2 .4 .9])
hold on
plot3(coor_C_s(1,:),coor_C_s(2,:),coor_C_s(3,:),'og')
plot3(coor_V_sn_s(1,:),coor_V_sn_s(2,:),coor_V_sn_s(3,:),'^r')
axis equal
hold on
M=mesh(X_t_s,Y_t_s,Z_t_s);
set(M,'facealpha',0)
set(M,'edgecolor',[.2 .4 .9])
hold on
plot3(coor_C_t_s(1,:),coor_C_t_s(2,:),coor_C_t_s(3,:),'og')
plot3(coor_V_sn_t_s(1,:),coor_V_sn_t_s(2,:),coor_V_sn_t_s(3,:),'^r')
axis equal
grid minor
set(gca,'FontWeight','bold', 'FontSize', 26)
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
surf(X,Y,Z,l_1d_w),colorbar, axis equal
title (['Portanza con ', '\alpha',' = ', num2str(alfa*180/pi),'° ',...
    '\beta',' = ', num2str(beta*180/pi),'°']),xlabel 'x / x_{M_{AC}}', ylabel 'y / y_{M_{AC}}', zlabel 'z / z_{M_{AC}}'
hold on
surf(X_t,Y_t,Z_t,l_1d_t),colorbar, axis equal
title (['Portanza con ', '\alpha',' = ', num2str(alfa*180/pi),'° ',...
    '\beta',' = ', num2str(beta*180/pi),'°']),xlabel 'x / x_{M_{AC}}', ylabel 'y / y_{M_{AC}}', zlabel 'z / z_{M_{AC}}'
grid minor
set(gca,'FontWeight','bold', 'FontSize', 26)

figure(3)
plot(coor_C(2,1:2*Nspan),GAMMA_w,'o', 'LineWidth', 3, 'MarkerSize', 12), grid on
title (['\Gamma wing', ' sulle diverse corde con ', '\alpha',' = ',...
    num2str(alfa*180/pi),'°', '\beta',' = ', num2str(beta*180/pi),'°'])
xlabel 'y / y_{M_{AC}}', ylabel '\Gamma'
grid minor
set(gca,'FontWeight','bold', 'FontSize', 26)

figure(4)
plot(coor_C_t(2,1:2*Nspan_t),GAMMA_t,'o', 'LineWidth', 3, 'MarkerSize', 12), grid on
title (['\Gamma tail', ' sulle diverse corde con ', '\alpha',' = ',...
    num2str(alfa_t*180/pi),'°', '\beta',' = ', num2str(beta*180/pi),'°'])
xlabel 'y / y_{M_{AC}}', ylabel '\Gamma'
grid minor
set(gca,'FontWeight','bold', 'FontSize', 26)

Pmed_vort=0.5*(coor_V_sn+coor_V_dx);
Pmed_vort_t=0.5*(coor_V_sn_t+coor_V_dx_t);



Q_w = zeros(1,length(d_2d_w));
Q_t = zeros(1,length(d_2d_t));
dx_w = abs(X(1,1)-X(2,1));
dx_t = abs(X_t(1,1)-X_t(2,1));

figure(5)
surf(X,Y,Z,'FaceColor','black', 'EdgeColor','None')
colormap hsv
alpha(.4)
hold on 
quiver3(X(2,1:end-1)+ 0.5*dx_w,(Y(2,1:end-1) + 0.5*dy_w),(Z(2,1:end-1) - sin(alfa)*0.5*dx_w),d_2d_w,Q_w,Q_w,'r', 'LineWidth', 2)
axis equal
hold on
surf(X_t,Y_t,Z_t,'FaceColor','black', 'EdgeColor','None')
colormap hsv
alpha(.4)
hold on 
quiver3(X_t(1,1:end-1)+0.75*dx_t,(Y_t(2,1:end-1) + 0.5*dy_t),(Z_t(1,1:end-1) - sin(alfa_t)*0.75*dx_t),d_2d_t,Q_t,Q_t,'r', 'LineWidth', 2)
axis equal
hold on 
quiver3((X(2,1:end-1)+ 0.5*dx_w),(Y(2,1:end-1) + 0.5*dy_w),(Z(2,1:end-1) - sin(alfa)*0.5*dx_w),Q_w,Q_w,l_2d_w,'b', 'LineWidth', 2)
axis equal
hold on 
quiver3(X_t(1,1:end-1)+0.75*dx_t,(Y_t(2,1:end-1) + 0.5*dy_t),(Z_t(1,1:end-1) - sin(alfa_t)*0.75*dx_t),Q_t,Q_t,l_2d_t,'b', 'LineWidth', 2)
axis equal
grid minor
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontWeight','bold', 'FontSize', 26)