function [F,G,Tcheck] = dioph_mtx(A,B,T);
% Alternative version for polynomial Diophantine equation solving
% Solves AF + BG = T for F & G.
da = length(A)-1; db=length(B)-1; % da = deg(A) etc.
T = [T, zeros(1,da+db-length(T)+1)]; % pad with zeros
dt=length(T)-1; % Convert to full length
dg =da-1; df=dt-da ; % assuming for RAB to be square
B = [ zeros(1, df-db+1), B]; % pad with leading zeros
Rac = [convmtx(A,df+1);convmtx(B, dg+1)]; % Construct RAB = [RA;RB]
FG = T/Rac; % [F,G] = T R^-1
F = FG(1:df+1); G = FG(df+2:df+dg+2);% Note F & G are row vectors.
Tcheck = polyadd(conv(A,F), conv(B,G)); % Verify solution
return