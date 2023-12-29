
function [v_ind] = ind_vort(x_v_sn,x_v_dx,x_c,g,toll)


% calcola la velocità indotta nel punto di controllo del pannello i dal vortice a staffa del pannello j 
% applicando 3 volte biot - savart per la testa del vortice e le due code
% INPUT
% x_v_sn : vettore posizione estremo sinistro vortice staffa
% x_v_dx : vettore posizione estremo destro vortice staffa
% x_c : vettore posizione punto di controllo
% g : valore della circolazione 
% OUTPUT
% v_ind : vettore velocità



x_inf_sn = x_v_sn+[1000 0 0]';
x_inf_dx = x_v_dx+[1000 0 0]';


v_ind = Biot_Savart(x_inf_sn,x_v_sn,x_c,g,toll)+Biot_Savart(x_v_sn,x_v_dx,x_c,g,toll)+Biot_Savart(x_v_dx,x_inf_dx,x_c,g,toll);