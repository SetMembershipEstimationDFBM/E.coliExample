% Shen, X., & Budman, H. (2021). Set Membership Estimation with Dynamic 
% Flux Balance Models. Processes, 9(10), 1762.
% https://www.mdpi.com/2227-9717/9/10/1762

% Copyright 2021 Xin Shen
% Copyright 2021 Hector Budman

% Multi-Parametric Toolbox 3.0 is required
% Algorithm has been tested on MATLAB R2018a.

%% Parameter calculations
clear;
clc;
startTime = datetime('now');

% stoichiometry coefficient
A = [0     -9.46   -9.84   -19.23;...   %glc
    -35    -12.92  -12.73  0;...       %o2
    -39.43   0       1.24    12.12;...    %ace
    1        1      1        1];       %X
% objective vector
c = ones(4,1);

% true plant initial concentration
z0 = [0.4 0.21 0.2 0.001]'; %(glc [mM], o2[mM], ace[mM], X[mM])
vol = 0.3;    % initla volume V[L]

% lower bound and upper bound
lb = [0.38; 0.1995; 0.19; 0.00095];
ub = [0.42; 0.2205; 0.21; 0.00105];

z = lb + (ub-lb).*rand(length(z0),1); 

% define initla set
Box = Polyhedron('lb',lb,'ub',ub);

Km = 0.015;        % [mM]
GUR_max = 6.5;     % [mM/g-dw/hr]
OUR_max = 12;      % [mM/g-dw/hr]
kla = 4;           % [hr^-1]
z_fd = 5;

nstep = 160;
tend = 8; % h, total time of cell culture growth, guess
dt = tend/nstep; % h, time, guess
Timeseries = (1:nstep)*dt;

Constr = [-A(2,:); A(3,:); -A(1:3,:)*dt;-A(1,:);-A(4,:)*dt];
MptOptions.fastbreak = true;

b0 = [OUR_max; 100];
flux = sdpvar(4,1);     % decision variables
theta = sdpvar(5,1);  % parameters variables
obj = -c.'*flux;        % LP objective
G = [zeros(2,5);eye(5)];  
B = [b0;zeros(5,1)];
Con1 = [Constr*flux <= B+G*theta, flux>=0];    
plp = Opt(Con1, obj, theta, flux);
% solve the mpLP to obtain different critical regions
Solution = plp.solve();
% Only critical region 12 and 14 are related to this research.
% when concentration is very low, MPT toolbox may locate the state coexist
% in several critical regions because of tolerance issue. 
% 12 and 14 are the true critical regions.

%% Observer Initialization
% rng('default');
% sigma = [0.025,0.0033,0.0200,0.004].';
nstate = length(z0);
mindex = [4];    %measurements
nmeas = length(mindex);

oindex = cell(14,1);    %observable state(subsystem), include bio and volume
oindex{12} = [1,4];
oindex{14} = [1,3,4];

modmsindex = cell(24,1);    %model measurements index
modmsindex{12} = [2];
modmsindex{14} = [3];
CRindex = [12];
C = [0 0 0 1];

% noise
Q0 = diag(zeros(4,1));                     % NO process noise covariance (5 states)!!!
P0 = diag((z0*0.05/3).^2);                  %initial state convariance (5 states)
R = diag((sigma(mindex)).^2);         %measurement noise covariance

Q = Q0(oindex{12},oindex{12});             %process noise covariance (3 states for CR12)
P = P0(oindex{12},oindex{12});             %initial state convariance (3 states for CR12)

%range
% This research assume no process noise.
Qlb = sigma*0;                             %No process noise !!!
Qub = sigma*0;                             %No process noise !!!
Mlb = -sigma*0.1;
Mub = sigma*0.1;
RBox = Polyhedron('lb',Mlb(mindex),'ub',Mub(mindex));

plantstate = z;
Obserstate = z0(oindex{CRindex});
xPlant = zeros(nstate,nstep);
xPlant(:,1) = plantstate;
switchpoint = [];

xlb = zeros(nstate,nstep);
xub = zeros(nstate,nstep);
xlb(:,1) = lb;
xub(:,1) = ub;
y = zeros(nmeas,nstep);

xPre = cell(nstep,1);
PPre = cell(nstep,1);
xCrr = cell(nstep,1);
PCrr = cell(nstep,1);
Amatrix = cell(nstep,1);
PhdSet = cell(nstep,1);
xCrr{12} = Obserstate;
xCrr{14} = [];
PCrr{1} = P0(oindex{CRindex},oindex{CRindex}); 

BoxSet = cell(nstep,2);
BoxSet{1,1} = Box;
BoxSet{1,2} = [];

PlantStateWrap = @(z,rf,rp)PlantState(z,rf,rp,OUR_max,z_fd,dt,GUR_max,Km,kla,A,c,Constr,Q0,Qlb,Qub,vol);
SensorWrap = @(plantstate)Sensor(plantstate,mindex,R,Mlb,Mub);

tol = 0.08;

IsSwitched = false;
method = 'vrep';

fig = figure;
fig.Color = [1 1 1];
subplot(2,2,1)
ylabel('Glucose');
hold on;
h1 = animatedline;
h1l = animatedline('Color','r');
h1u = animatedline('Color','r');
h1e = animatedline('Color','g');
xlim([0 8]);

subplot(2,2,2)
ylabel('Oxygen');
hold on;
h2 = animatedline;
h2l = animatedline('Color','r');
h2u = animatedline('Color','r');
xlim([0 8]);

subplot(2,2,3)
ylabel('Acetate');
hold on;
h3 = animatedline;
h3l = animatedline('Color','r');
h3u = animatedline('Color','r');
h3e = animatedline('Color','g');
xlim([0 8]);

subplot(2,2,4)
ylabel('Biomass');
hold on;
h4 = animatedline;
h4l = animatedline('Color','r');
h4u = animatedline('Color','r');
h4e = animatedline('Color','g');
xlim([0 8]);

flag = false;

%%
for k=2:nstep
 if IsSwitched
     % check observablility
     if rank(ObFun14(Obserstate))~=3 || rank(AFun14(Obserstate))~=3
         disp(Obserstate);
     end
 end
    
%batch process
    rf = 0;
    rp = 0;
    
%--------------------PLANT----------------------
    plantstate = PlantStateWrap(plantstate,rf,rp);
    xPlant(:,k) = plantstate;
    ym = SensorWrap(plantstate);
    y(:,k) = ym;
        
 %------------------Montinor--------------------
    if IsSwitched
    else
        [IsSwitched,CRindex,est] = Monitor(Obserstate,GUR_max,Km,Solution,tol,CRindex,xlb(:,k-1),xub(:,k-1),rf,dt,vol);

        if IsSwitched && (flag == false) 
            disp('Montior switch critical region!');
            switchpoint = [switchpoint k];
            temp = zeros(3,1);
            temp(1) = Obserstate(1);
            temp(3) = Obserstate(2);
            temp(2) = est;
            Obserstate = temp;
            Ptemp = zeros(3);
            Ptemp([1 3],[1 3]) = P;
            Ptemp(2,2) = ((xub(2,k-1)-xlb(2,k-1))/2/3)^2;
            P = Ptemp;
            P = P*1.1;    %If necessary, the covariance can be initilized larger to deal with the early or late swtich issue
            Q = Q0(oindex{CRindex},oindex{CRindex});     %no process noise has been set as 0
            flag = true;
            xCrr{14} = [xCrr{14} Obserstate];
        end
    end   
    
    
%----------Set-Membership-Estimation-----------------    
    [Box] = BoxPropogate(Box,ym,Obserstate,rf,rp,z_fd,dt,GUR_max,Km,kla,A,...
        Solution,Qlb,Qub,RBox,C,CRindex,method,sigma,P,vol);
    BoxSet{k,1} = Box;
    BoxSet{k,2} = [];  %scaling
    
    xlb(:,k) = Box.Internal.lb;
    xub(:,k) = Box.Internal.ub;
               
%------------------Kalman Filter--------------------    


    ModelStateWrap = @(z,rf,rp)ModelState(z,rf,rp,z_fd,dt,GUR_max,Km,A,Solution,CRindex,vol);
    ModelMeasureWrap = @(state)ModelMeasure(state,modmsindex{CRindex});

    [xPre{k},PPre{k},Obserstate,P,Amatrix{k}] = ekf(@(z)ModelStateWrap(z,rf,rp),...
            Obserstate,P,ModelMeasureWrap,ym,Q,R);
    if ~IsSwitched
        xCrr{12} = [xCrr{12} Obserstate];
    else
        xCrr{14} = [xCrr{14} Obserstate];
    end
    PCrr{k} = P;     

    
%---------------------Plot--------------------        
    addpoints(h1,k*dt,plantstate(1));
    addpoints(h2,k*dt,plantstate(2));
    addpoints(h3,k*dt,plantstate(3));
    addpoints(h4,k*dt,plantstate(4));

    addpoints(h1l,k*dt,xlb(1,k));
    addpoints(h2l,k*dt,xlb(2,k));
    addpoints(h3l,k*dt,xlb(3,k));
    addpoints(h4l,k*dt,xlb(4,k));

    addpoints(h1u,k*dt,xub(1,k));
    addpoints(h2u,k*dt,xub(2,k));
    addpoints(h3u,k*dt,xub(3,k));
    addpoints(h4u,k*dt,xub(4,k));
   
    drawnow;

end
endTime = datetime('now');

function stateobserved = Sensor(plantstate,mindex,R,Mlb,Mub)
    stateobserved = plantstate(mindex);
    % add bounded Gaussian noise
    if length(mindex)>1
        noise = mvrandn(Mlb(mindex),Mub(mindex),R,1);
    else
        while true
            noise = normrnd(0,R);
            if noise>=Mlb(mindex) && noise<=Mub(mindex)
                break ;
            end
        end
    end
    stateobserved = plantstate(mindex) + noise;
end
function stateobserved = ModelMeasure(state,index)
    stateobserved = state(index);
end

function z = ModelState(z,rf,rp,z_fd,dt,GUR_max,Km,A,Solution,CRindex,vol)
switch CRindex    
    case 12
        b = zeros(5,1);
        b(4) = GUR_max*z(1)/(Km+z(1));
        flux = round(Solution.xopt.Set(12).Functions('primal').F,6)*b + round(Solution.xopt.Set(12).Functions('primal').g,6);
        z(1) = z(1) - rf*dt*z(1)/vol + A(1,:)*flux*z(2)*dt + rf*dt/vol*z_fd;
        z(2) = z(2) - (rf - rp)*dt/vol*z(2) + A(4,:)*flux*z(2)*dt;
    case 14
        b = zeros(5,1);
        b(3) = z(2)/z(3) - rf/vol*z(2)/z(3)*dt;
        b(4) = GUR_max*z(1)/(Km+z(1));
        flux = round(Solution.xopt.Set(14).Functions('primal').F,6)*b + round(Solution.xopt.Set(14).Functions('primal').g,6);
        z(1) = z(1) - rf*dt*z(1)/vol + A(1,:)*flux*z(3)*dt + rf*dt/vol*z_fd;
        z(2) = z(2) - rf*dt*z(2)/vol + A(3,:)*flux*z(3)*dt;
        z(3) = z(3) - (rf - rp)*dt/vol*z(3) + A(4,:)*flux*z(3)*dt;
end
        
end



function z = PlantState(z,rf,rp,OUR_max,z_fd,dt,GUR_max,Km,kla,A,c,Constr,Q,Qlb,Qub,vol)
    global b1
    b0 = [OUR_max; 100];
    b1 = z(1:3)/z(4) - rf/vol*z(1:3)/z(4)*dt;
    b1(1) = b1(1) + rf/vol*z_fd/z(4)*dt;
    b1(2) = b1(2) + kla*(0.21 - z(2))/z(4)*dt;
    b2 = GUR_max*(z(1)/(Km + z(1)));
    b3 = 1 - (rf-rp)*dt/vol;
    
    global nu 
    problem = Opt('f',-c,'A',Constr,'b',[b0;b1;b2;b3],'lb',zeros(4,1));
    solution = problem.solve;
    nu = round(solution.xopt,6);
    
    z = z + ODE(0,z,nu,rf,z_fd,rp,A,kla,vol)*dt;
    z(z<=0)=0;
    
    function dzdt = ODE(t,z,nu,rf,z_fd,rp,A,kla,vol)
        dzdt = zeros(4,1);
        dzdt(1:4) = A*nu*z(4) - rf/vol*z(1:4);
        dzdt(1) = dzdt(1) + rf/vol*z_fd;
        dzdt(2) = dzdt(2) + kla*(0.21 - z(2));
        dzdt(4) = dzdt(4) + z(4)*rp/vol;
    end
end

function [BoxNew] = BoxPropogate(Box,ym,Obserstate,rf,rp,z_fd,dt,GUR_max,...
                       Km,kla,A,Solution,Qlb,Qub,RBox,C,CRindex,method,sigma,P,vol)
switch CRindex 
    case 12

        b = zeros(5,1);
        bhat = zeros(5,2);
        [b(4),bhat(4,1)]=jaccsd(@(x)(GUR_max*x/(Km+x)),Obserstate(1));
        flux = round(Solution.xopt.Set(12).Functions('primal').F,6)*b + round(Solution.xopt.Set(12).Functions('primal').g,6);
        bffHat = zeros(4,2);        
        translate = [z_fd*rf*dt/vol; 0.21*kla*dt; 0; 0];
        
        bff = A*flux*Obserstate(2)*dt + translate;

        Aff = [1 - rf*dt/vol; 1 - rf*dt/vol-kla*dt;...
               1 - rf*dt/vol; 1 - (rf-rp)*dt/vol];
        T = diag(Aff);
        
        termB = (3*P(2,2)^0.5 + Obserstate(2))*dt*A*round(Solution.xopt.Set(12).Functions('primal').F,6)*bhat + bffHat;
        termD = dt*A*flux;
        
        var = 3*diag(P).^0.5;
        Err = Polyhedron('lb',-var,'ub',var);
        Errbio = Polyhedron('lb',-var(2),'ub',var(2));
        term1 = PolyUnion(termB*Err).outerApprox;
        term2 = PolyUnion(termD*Errbio).outerApprox;
        safetyl = term1.Internal.lb + term2.Internal.lb;
        safetyu = term1.Internal.ub + term2.Internal.ub;

        Prior0 = Box.affineMap(T,method) + bff; 
        Prior1 = PolyUnion(Prior0);  
        Prior1Box = Prior1.outerApprox;
        
        Prior1Boxlb = Prior1Box.Internal.lb + Qlb + safetyl;
        Prior1Boxub = Prior1Box.Internal.ub + Qub + safetyu;
        BoxPrior = Polyhedron('lb',Prior1Boxlb,'ub',Prior1Boxub);
        
    case 14
        b = zeros(5,1);
        bhat = zeros(5,3);

        [b(3),bhat(3,:)]=jaccsd(@(x)x(2)/x(3) - rf/vol*x(2)/x(3)*dt,Obserstate);
        [b(4),bhat(4,1)]=jaccsd(@(x)(GUR_max*x/(Km+x)),Obserstate(1));

        flux = round(Solution.xopt.Set(14).Functions('primal').F,6)*b + round(Solution.xopt.Set(14).Functions('primal').g,6);


        bffHat = zeros(4,3);
                
        translate = [z_fd*rf*dt/vol; 0.21*kla*dt; 0; 0];
        
        bff = A*flux*Obserstate(3)*dt + translate;
                
        Aff = [1 - rf*dt/vol; 1 - rf*dt/vol-kla*dt;...
               1 - rf*dt/vol; 1 - (rf-rp)*dt/vol];
        T = diag(Aff);
        
        termB = (3*P(3,3)^0.5 + Obserstate(3))*dt*A*round(Solution.xopt.Set(14).Functions('primal').F,6)*bhat + bffHat;
        termD = dt*A*flux;
        
        var = 3*diag(P).^0.5;
        Err = Polyhedron('lb',-var,'ub',var);
        Errbio = Polyhedron('lb',-var(3),'ub',var(3));
        term1 = PolyUnion(termB*Err).outerApprox;
        term2 = PolyUnion(termD*Errbio).outerApprox;
        safetyl = term1.Internal.lb + term2.Internal.lb;
        safetyu = term1.Internal.ub + term2.Internal.ub;
        
        Prior0 = Box.affineMap(T,method) + bff; 
        Prior1 = PolyUnion(Prior0);  
        Prior1Box = Prior1.outerApprox;

        Prior1Boxlb = Prior1Box.Internal.lb + Qlb + safetyl;
        Prior1Boxub = Prior1Box.Internal.ub + Qub + safetyu;
        BoxPrior = Polyhedron('lb',Prior1Boxlb,'ub',Prior1Boxub);
        
end

    R = Polyhedron('A',-RBox.A*C,'b',RBox.b-RBox.A*ym);
    BoxPost = intersect(R,BoxPrior);

    U = PolyUnion(BoxPost);    
    OutBox = U.outerApprox;
    
    % When the state interval is very small. Numerical issue can happen. So
    % a very tiny interval is given artifically to let the algorithm
    % continue.
    lb = OutBox.Internal.lb - sigma.*0;
    lb(lb<0) = 0;
   
    ub = OutBox.Internal.ub + sigma.*0 ;
     ub(ub<0) = 0;
    lb1 = lb;
    ub1 = ub;
    distance = 0.000001;
    
    conver = abs(ub-lb)<=distance & lb>0;
    lb1(conver) = (lb(conver) + ub(conver))*0.5 - distance*0.01;
    ub1(conver) = (lb(conver) + ub(conver))*0.5 + distance*0.01;
    
    conver1 = abs(ub-lb)<=distance & lb<=0;
    lb1(conver1) = 0;
    ub1(conver1) = distance*0.01;

    BoxNew = Polyhedron('lb',lb1,'ub',ub1);
    BoxNew.computeVRep();       
end

function [IsSwitched,CRindex,est] = Monitor(z,GUR_max,Km,Solution,tol,CRindex,xlb,xub,rf,dt,vol)
    theta3lb = xlb(3)/z(2) - rf/z(2)*xlb(3)/vol*dt;
    theta3ub = xub(3)/z(2) - rf/z(2)*xub(3)/vol*dt;
    theta3p = (theta3lb+theta3ub)/2;
    
    b = zeros(5,1);
    b(3) = theta3p;
    b(4) = GUR_max*z(1)/(Km+z(1));
    vLeft = round(Solution.xopt.Set(12).Functions('primal').F,6)*b + round(Solution.xopt.Set(12).Functions('primal').g,6);
    vRight = round(Solution.xopt.Set(14).Functions('primal').F,6)*b + round(Solution.xopt.Set(14).Functions('primal').g,6);
    diff = norm(vLeft-vRight);
    
    % detect switch
    if norm(diff)<=tol
        CRindex = 14;
        IsSwitched = true;

        Theta = sym('Theta',[5,1]);
        F = (round(Solution.xopt.Set(12).Functions('primal').F,6)-round(Solution.xopt.Set(14).Functions('primal').F,6))*Theta...
        +(round(Solution.xopt.Set(12).Functions('primal').g,6)-round(Solution.xopt.Set(14).Functions('primal').g,6));
        theta4 = b(4);
        F0 = subs(F,Theta(4),theta4);
        theta3 = eval(solve(F0==0));
        % estimate the unobserved concetration
        est = theta3./(1/z(2)- rf/z(2)/vol*dt);
    else
        IsSwitched = false;
        est = [];
    end
end









