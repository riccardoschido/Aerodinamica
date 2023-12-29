function [X,Y,Z,X_V,Y_V,Z_V,coordC,coordV_sn,coordV_dx,N]=wing_geo(rootchord,tipchord,span,sweep,dihedral ,Nchord, Nspan, alfa, beta) 


% dati i valori geometrici e il numero di suddivisioni desiderate
% ricostruisce l'ala

% INPUT
% rootchord : lunghezza della corda alla radice
% tipchord : lunghezza della corda all'estremità
% span : apertura alare
% sweep : angolo di freccia 
% dihedral : angolo di diero
% Nchord : numero di suddivisioni lungo la corda
% Nspan : numero di suddivisioni lungo la SEMI-apertura

% OUTPUT
% X,Y,Z : matrici contenenti tutte le coordinate x,y,z dei punti nodali (la
% posizione riga e colonna rappresenta la posizione riga-colonna sull'ala
% vista dall'alto
% X_V,Y_V,Z_V : matrici contenenti le coordinate dei punti di estremità
% dei vortici a staffa. (la posizione riga e colonna rappresenta 
% la posizione riga-colonna sull'ala vista dall'alto)
% coordC : matrice contenente le coordinate di tutti i punti di controllo ordinati
% dal'estremo in alto a sinistra dell'ala fino in basso a dstra (per righe). Le 3
% righe della matrice rappresentano le 3 coordinate
% coordV_sn,coordV_dx : matrici contenenti le coordinate dei punti si estremità dei
% vortici a staffa. coordV_sn contiene tutti i punti "a sinistra"
% dei vortici, la matrice coordV_dx solo quelli "a destra"
% N : matrice contenente tutte le normali ai punti di controllo. ordinati
% dal'estremo in alto a sinistra dell'ala fino in basso a dstra (per righe). Le 3
% righe della matrice rappresentano le 3 coordinate

%% definizioni
dy=span/2/Nspan;
dx_root=rootchord/Nchord;
dx_tip=tipchord/Nchord;
Yhalf=ones(Nchord+1,1)*[0:dy:span/2];
d=span/2*tan(sweep)+0.25*rootchord-0.25*tipchord;

x_root=0:dx_root:rootchord;
x_tip=d+[0:dx_tip:tipchord];

%% coordinate punti nodali
%coefficienti angolari delle rette congiungenti i nodi delle corde 
%di radice ed estremità
m=2*(x_tip-x_root)/span;
coeff=[m' ,x_root'];
Xhalf=zeros(Nchord+1,Nspan+1);
y_half=0:dy:(span/2);
for k=1: (Nchord+1)
    Xhalf(k,:)=polyval(coeff(k,:),y_half);
end 

% ruotiamo l'ala del diedro 
%Xhalf rimane la stessa 
Yhalf=Yhalf*cos(dihedral);
Zhalf=Yhalf*sin(dihedral);

X=[fliplr(Xhalf(:,2:end)),Xhalf];
Y=[-fliplr(Yhalf(:,2:end)),Yhalf];
Z=[fliplr(Zhalf(:,2:end)),Zhalf];


%% coordinate estremi vortici a staffa

X_V = X(1:end-1,:)-0.25*(X(1:end-1,:)-X(2:end,:));
Y_V = Y(1:end-1,:)-0.25*(Y(1:end-1,:)-Y(2:end,:));
Z_V = Z(1:end-1,:)-0.25*(Z(1:end-1,:)-Z(2:end,:));

X_V_sn = X_V(:,1:end-1);
Y_V_sn = Y_V(:,1:end-1);
Z_V_sn = Z_V(:,1:end-1);

X_V_dx = X_V(:,2:end);
Y_V_dx = Y_V(:,2:end);
Z_V_dx = Z_V(:,2:end);

%% coordinate punti controllo

X_C = X(1:end-1,:)-0.75*(X(1:end-1,:)-X(2:end,:));
X_C = 0.5*(X_C(:,1:end-1)+X_C(:,2:end));

Y_C = 0.5*(Y(:,1:end-1)+Y(:,2:end));
Y_C = Y_C(1:end-1,:);

Z_C = 0.5*(Z(:,1:end-1)+Z(:,2:end));
Z_C = Z_C(1:end-1,:);

%%  matrice per i punti di controllo:
% le tre righe sono le x,le y e le z dei punti di controllo
% ordine sn-dx alto-basso

[m_c,n_c] = size(X_C);
coordC = [reshape(X_C',1,m_c*n_c);reshape(Y_C',1,m_c*n_c);reshape(Z_C',1,m_c*n_c)];


coordV_sn = [reshape(X_V_sn',1,m_c*n_c);reshape(Y_V_sn',1,m_c*n_c);reshape(Z_V_sn',1,m_c*n_c)];
coordV_dx = [reshape(X_V_dx',1,m_c*n_c);reshape(Y_V_dx',1,m_c*n_c);reshape(Z_V_dx',1,m_c*n_c)];

%% Rotazioni 

R=[cos(alfa)*cos(beta) -sin(beta) sin(alfa)*cos(beta) ;...
   cos(alfa)*sin(beta)  cos(beta) sin(alfa)*sin(beta) ;...
   -sin(alfa)               0         cos(alfa)];

for j = 1 : length (X(1,:))
    for i = 1 : length(X(:,1))

        v = [X(i,j) Y(i,j) Z(i,j)];
        v_r = R * v';
        X(i,j) = v_r(1);
        Y(i,j) = v_r(2);
        Z(i,j) = v_r(3);
    end
end


for j = 1 : length (X_V(1,:))
    for i = 1 : length(X_V(:,1))

        v = [X_V(i,j) Y_V(i,j) Z_V(i,j)];
        v_r = R * v';
        X_V(i,j) = v_r(1);
        Y_V(i,j) = v_r(2);
        Z_V(i,j) = v_r(3);
    end
end

for i = 1 : length(coordC)
    v = coordC(:,i);
    v_r = R * v;
    coordC(:,i) = v_r;

end

for i = 1 : length(coordV_sn)
    v = coordV_sn(:,i);
    v_r = R * v;
    coordV_sn(:,i) = v_r;

end

for i = 1 : length(coordV_dx)
    v = coordV_dx(:,i);
    v_r = R * v;
    coordV_dx(:,i) = v_r;

end

%% matrice delle normali
% ogni colonna per un pannello

n = [ones(Nchord,Nspan),-ones(Nchord,Nspan)];
n = [reshape(n',1,2*Nchord*Nspan)];
N = [zeros(1,2*Nspan*Nchord) ; n*sin(dihedral) ; ones(1,2*Nspan*Nchord)*cos(dihedral)];

for i = 1 : length(N)
    v = N(:,i);
    v_r = R * v;
    N(:,i) = v_r;

end


