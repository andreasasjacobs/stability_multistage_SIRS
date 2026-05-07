%Created by Jacobs Somnic on 20260227
%Latest revision: 
%Purpose: Revision of multistage SIRS model
%This version currently only works for n=1,2,3

clc;clear;

%==================NUMBER of STAGES=====================
prompt="How many stages? n=";
n= input(prompt);

prompt="DFE or EE? (DFE-->1,EE-->2)";
type= input(prompt);

%==================COMPARTMENTS=========================
syms S R
assume([S R],'real');assumeAlso([S R],'positive')
I=sym('I', [1 n],'real'); assumeAlso(I,'positive')


%==================PARAMETERS===========================
%Defining the general parameters
syms alpha beta lambda mu R0
assume([alpha beta lambda mu R0],'real');assumeAlso([alpha beta lambda mu R0],'positive');

%Defining the stages parameters
delta=sym('delta', [1 n],'real'); assumeAlso(delta,'positive')
gamma= sym('gamma', [1 n],'real'); assumeAlso(gamma,'positive')
theta= sym('theta', [1 n],'real'); assumeAlso(theta,'positive')

%Defining Positive coefficients
m=(n^2+3*n)/2;
C= sym('C', [1 m],'real'); assumeAlso(gamma,'positive')
%===========================================================


%==================SYSTEM of ODE=============================
dS =lambda+alpha*R-mu*S-beta*S*I(1) ;

if n == 1
    % Single-stage model (no gamma transitions)
    dI(1) = beta*S*I(1) - (mu + delta(1) + theta(1))*I(1);
else
    % Multistage model
    for k = 1:n
        if k == 1
            dI(k) = beta*S*I(k) - (mu + gamma(k) + delta(k) + theta(k))*I(k);
        elseif k == n
            dI(k) = gamma(k-1)*I(k-1) - (mu + delta(k) + theta(k))*I(k);
        else
            dI(k) = gamma(k-1)*I(k-1) - (mu + gamma(k) + delta(k) + theta(k))*I(k);
        end
    end
end
dR = sum(theta.*I) - (alpha + mu)*R;


%==================DELTA VARIABLES======================
syms DS DR real
DI = sym('DI', [1 n], 'real');

Delta = [DS, DI, DR];
Delta_dot = [dS, dI, dR];

%==================EQUILIBRIUM=========================
syms Sstar Rstar
assume([Sstar Rstar],'real');assumeAlso([Sstar Rstar],'positive')
Istar = sym('Istar', [1 n], 'real'); assumeAlso(Istar,'positive')

dS_star = subs(dS, [S, I, R], [Sstar, Istar, Rstar]);
dI_star = subs(dI, [S, I, R], [Sstar, Istar, Rstar]);
dR_star = subs(dR, [S, I, R], [Sstar, Istar, Rstar]);

%==================DELTA VERSION OF THE SYSTEM==========
% S = Sstar + DS, I = Istar + DI, R = Rstar + DR
dS_delta = simplify(subs(dS-dS_star, [S, I, R], [Sstar + DS, Istar + DI, Rstar + DR]));
dI_delta = simplify(subs(dI-dI_star, [S, I, R], [Sstar + DS, Istar + DI, Rstar + DR]));
dR_delta = simplify(subs(dR-dR_star, [S, I, R], [Sstar + DS, Istar + DI, Rstar + DR]));

Delta_dot = [dS_delta, dI_delta, dR_delta];

%==================DEFINE LYAPUNOV FUNCTION=============
if type==1
    if n==1
        V = sym(1/2)*(DS + DI(1) + DR)^2 ...
         +C(1)*(DI(1)) ...
         +C(2)/2*DR^2;
    elseif n==2
        V = sym(1/2)*(DS + DI(1)+ DI(2) + DR)^2 ...
         +C(1)*(DI(1)) ...
         +C(2)/2*(DS + DI(1)+ DR)^2 ...
         +C(3)/2*(DI(2)+ DR)^2 ...
         +C(4)/2*DI(2)^2 ...
         +C(5)/2*DR^2;
    
    elseif n==3
        V = sym(1/2)*(DS + DI(1)+ DI(2)+ DI(3) + DR)^2 ...
         +C(1)*DI(1) ...
         +C(2)/2*(DS + DI(1)+ DR)^2 ...
         +C(3)/2*(DS + DI(1)+ DI(2) + DR)^2 ...
         +C(4)/2*(DI(2)+ DR)^2 ...
         +C(5)/2*(DI(2)+DI(3)+DR)^2 ...
         +C(6)/2*(DI(2)+DI(3))^2 ...
         +C(7)/2*DI(2)^2 ...
         +C(8)/2*DI(3)^2 ...
         +C(9)/2*DR^2;
    end
elseif type==2
    if n==1
        V = sym(1/2)*(DS + DI(1) + DR)^2 ...
         +C(1)*(DI(1) - Istar(1)*log((DI(1) + Istar(1))/Istar(1))) ...
         +C(2)/2*DR^2;
    elseif n==2
        V = sym(1/2)*(DS + DI(1)+ DI(2) + DR)^2 ...
         +C(1)*(DI(1) - Istar(1)*log((DI(1) + Istar(1))/Istar(1))) ...
         +C(2)/2*(DS + DI(1)+ DR)^2 ...
         +C(3)/2*(DI(2)+ DR)^2 ...
         +C(4)/2*DI(2)^2 ...
         +C(5)/2*DR^2;
    elseif n==3
        V = sym(1/2)*(DS + DI(1)+ DI(2)+ DI(3) + DR)^2 ...
         +C(1)*(DI(1) - Istar(1)*log((DI(1) + Istar(1))/Istar(1))) ...
         +C(2)/2*(DS + DI(1)+ DR)^2 ...
         +C(3)/2*(DS + DI(1)+ DI(2) + DR)^2 ...
         +C(4)/2*(DI(2)+ DR)^2 ...
         +C(5)/2*(DI(2)+DI(3)+DR)^2 ...
         +C(6)/2*(DI(2)+DI(3))^2 ...
         +C(7)/2*DI(2)^2 ...
         +C(8)/2*DI(3)^2 ...
         +C(9)/2*DR^2;
    end
end
%==================DERIVATIVE of LYAPUNOV FUNCTION======
dV = jacobian(V, Delta) * Delta_dot.';

%Some necessarry subtitutions
if type==1
    dV=subs(dV, Sstar,S-DS);dV=subs(dV, Istar,I-DI);dV=subs(dV,Rstar,R-DR);
    dV=subs(dV,I,DI);dV=subs(dV,S,DS+Sstar);
    if n==1
        dV=subs(dV,Sstar,R0*(mu+delta(1)+theta(1))/beta);
    else
        dV=subs(dV,Sstar,R0*(mu+delta(1)+gamma(1)+theta(1))/beta);
    end
elseif type==2
    if n==1
        dV=subs(dV, Sstar,(mu+delta(1)+theta(1))/beta);
    else
        dV=subs(dV, Sstar,(mu+delta(1)+theta(1)+gamma(1))/beta);
    end
    dV=subs(dV, Istar(1),I(1)-DI(1));
end



dV = simplify(expand(dV));
dV = collect(dV, [DS, DI, DR])


%================EXTRACT THE PT COEFFICIENTS with second derivative=============
vars_cross = [DS, DI, DR];
H = jacobian(jacobian(dV, vars_cross), vars_cross);

CrossVec = H(triu(true(size(H)),1));


% remove DS*DR term
%all DS*DI(k) terms: n terms---all DI(i)*DI(j) terms among the infected classes: n(n-1)/2---so n+n(n-1)/2 = n(n+1)/2
pos_DSDR=n*(n+1)/2 +1; 
CrossVec(pos_DSDR) = [];
%================SOLVE LINEAR SYSTEMS to OBTAIN THE POSITIVE COEFFICIENTS=========
[A_coeff,b_coeff]=equationsToMatrix(CrossVec,C);
X=linsolve(A_coeff,b_coeff)
%=======================================================================================


