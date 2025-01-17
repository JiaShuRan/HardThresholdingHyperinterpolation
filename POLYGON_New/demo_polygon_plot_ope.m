function demo_polygon_plot_ope
% Codes based on Alvise Sommariva (University of Padova)
% Date: 18 Sep, 2023
% demo_interval_plot_new is is used to verify Theorem 4.5 and
% lower degree approximation on the polygon

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%--------------------------------------------------------------------------
LV=15;     % Hyperinterpolant tot degree.
NV=30;   % Degree of the rule used in hyperinterpolation.
NR=40;     % Degree of precision of the reference rule (estimate L2 error).

[xv,yv]=example_polygon(2);

% * Function to approximate:
% 1. degree L poly., 2. degree floor(L/2)-1 poly. 3. test functions
% (see line 265 approx).
funct_example=3;


% ........ Numerical approximation, varying the degree in "nV" ............

AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics
JzMV=[]; HzMV=[];

% Define quadrature rule for hyperinterpolation at L, with  ADE=NV.
XYW=cub_polygon(NV,xv,yv); X=XYW(:,1); Y=XYW(:,2); W=XYW(:,3);

% Test points
XYWR=cub_polygon(NR,xv,yv); XR=XYWR(:,1); YR=XYWR(:,2); WR=XYWR(:,3);

% Compute bounding box
xmin=min(xv); xmax=max(xv);
ymin=min(yv); ymax=max(yv);
dbox=[xmin xmax ymin ymax];

% Compute orthonormal basis matrix at nodes.
jvec=1:(LV+1)*(LV+2)/2;
[U,~,Q,R,~,degs] = dORTHVAND(LV,[X Y],W,jvec,[],dbox);

% .. testing AE_L2err hyperinterpolation error for each "f" at "deg" ..
poly_coeffs=[]; lambdaV=[];

% ... define function to approximate ...
g=define_function(funct_example,LV);

% ... evaluate function to approximate ...
gXY=feval(g,X,Y);

% ... add noise (if present) ...

a = 0; sigma = 0.1;
pert_impulse=a*(1-2*rand(length(gXY),1))*binornd(1,0.5);

var=sigma^2;pert_gauss=sqrt(var)*randn(size(gXY));

% add gaussian + impulse noise
pert=pert_impulse+pert_gauss;

% perturbed values
gXY_pert=gXY+pert;

% ... determine polynomial hyperinterpolant ...
% compute hyperinterpolant coefficients
coeff0=Q'*(sqrt(W).*gXY_pert);

% ... test hyperinterpolant with or withour filters ...

lambdas=sort(abs(coeff0),'descend'); coeff0=coeff0';
lambdaL=lambdas(8);

zeta0 = Q'*(sqrt(W).*pert); zeta0=zeta0';

Val_f = sum((gXY.^2).*W);
Val_f = ones(1,length(coeff0))*Val_f;

for k = 1:length(zeta0)
    coeff0_new = coeff0.*(abs(coeff0) > lambdas(k));
    Jz(k) = 2*coeff0_new*zeta0' - coeff0_new*coeff0_new';
    sgn_coeff = sgnfun(coeff0_new);
    zeta0_new = zeta0.*(abs(coeff0) > lambdas(k));
    Hz(k) = sum(abs(sgn_coeff))*lambdas(k)^2 - 2*lambdas(k)*sgn_coeff*zeta0_new';
end

for p=1:length(Hz)
    if Hz(p) >= 0
        Hz_great_0(p) = Hz(p);
    else
        Hz_great_0(p) = NaN;
    end
end


AErr_lasso = sqrt(Jz+Val_f+Hz);
AErr_hard = sqrt(Jz+Val_f);

coeff0 = coeff0';


Candiate = [8,16];

LV_initial = 7;
LV_final = 15;
LV_step = 1;

for j = 1:length(Candiate)

    for LV = LV_initial:LV_step:LV_final

        % Compute orthonormal basis matrix at nodes.
        jvec=1:(LV+1)*(LV+2)/2;
        [U,~,Q,R,~,degs] = dORTHVAND(LV,[X Y],W,jvec,[],dbox);

        % compute hyperinterpolant coefficients
        coeff0=Q'*(sqrt(W).*gXY_pert);

        % ... test hyperinterpolant with or withour filters ...

        lambdaL=lambdas(Candiate(j)); %coeff0=coeff0';
        jtest = LV-6;

        for ktest=1:6

            switch ktest

                case 1
                    hypermode='tikhonov';
                    parms.lambda=lambdaL;
                    parms.mu=[];
                    parms.b=ones(size(coeff0));
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 2
                    hypermode='filtered';
                    parms.lambda=[];
                    parms.mu=[];
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 3
                    hypermode='lasso';
                    parms.lambda=lambdaL;
                    parms.mu=ones(size(coeff));
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 4
                    hypermode='hybrid';
                    parms.lambda=lambdaL;
                    parms.mu=ones(size(coeff0));
                    parms.b=ones(size(coeff0));
                    parms.w=W;
                    parms.pert=pert;
                    parms.hybrid=0; % establishes it is a pre-choosen parameter.
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 5
                    hypermode='hard';
                    parms.lambda=lambdaL;
                    parms.mu=[];
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 6
                    hypermode='hyperinterpolation';
                    parms.lambda=[];
                    parms.mu=[];
                    parms.b=[];
                    coeff=coeff0;
            end

            % evaluate polynomial at reference points.
            gXYR=feval(g,XR,YR);
            VR=chebvand(LV,[XR YR],dbox);
            pXYR(:,ktest) = (VR(:,jvec)/R)*coeff;

            % errors
            AEinfV(jtest,ktest)=norm(gXYR-pXYR(:,ktest),inf); % absolute error (inf norm)
            AE2V(jtest,ktest)=sqrt(WR'*((gXYR-pXYR(:,ktest)).^2)); % absolute error (2 norm)
            beta0V(jtest,ktest)=sum(abs(coeff) > 0);
        end

    end

    InfV(j).matrix = AEinfV; L2Err(j).matrix = AE2V; 
end

%% Plot 

GG = LV_initial:LV_step:LV_final;

A1 = L2Err(1).matrix; A2 = L2Err(2).matrix; %A3 = L2Err(3).matrix; A4 = L2Err(4).matrix; 

figure(2)
subplot(2,2,1)
semilogy(GG,A1(:,5),'d','linewidth',1,'MarkerSize',22,'color','blue'), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
hold on

semilogy(GG,A1(:,3),'pentagram','linewidth',1,'MarkerSize',20,'color','red','MarkerFaceColor','red'), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
hold on

xlabel({'\textbf{Degree} $(a)$'},'interpreter','latex','fontsize',30);ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',30);
legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',25);
%title({'1'},'interpreter','latex','fontsize',35);

subplot(2,2,2)
semilogy(GG,A2(:,5),'d','linewidth',1,'MarkerSize',22,'color','blue'), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
hold on

semilogy(GG,A2(:,3),'pentagram','linewidth',1,'MarkerSize',20,'color','red','MarkerFaceColor','red'), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
hold on

xlabel({'\textbf{Degree} $(b)$'},'interpreter','latex','fontsize',30);ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',30);
%legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',25);
%title({'2'},'interpreter','latex','fontsize',35);

% subplot(2,2,3)
% semilogy(GG,A3(:,5),'d','linewidth',1,'MarkerSize',22,'color','blue'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
% hold on
% 
% semilogy(GG,A3(:,3),'pentagram','linewidth',1,'MarkerSize',20,'color','red','MarkerFaceColor','red'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
% hold on
% 
% xlabel({'\textbf{Degree} $(c)$'},'interpreter','latex','fontsize',30);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',30);
% %legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',25);
% %title({'3'},'interpreter','latex','fontsize',35);
% 
% subplot(2,2,4)
% semilogy(GG,A4(:,5),'d','linewidth',1,'MarkerSize',22,'color','blue'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
% hold on
% 
% semilogy(GG,A4(:,3),'pentagram','linewidth',1,'MarkerSize',20,'color','red','MarkerFaceColor','red'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
% hold on
% 
% xlabel({'\textbf{Degree} $(d)$'},'interpreter','latex','fontsize',30);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',30);
% %legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',25);
% %title({'4'},'interpreter','latex','fontsize',35);
subplot(2,2,3)

plot(1:length(zeta0_new),Hz,'-.o','linewidth',1,'MarkerSize',10,'color','k'), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

plot(1:length(zeta0),Hz_great_0,'o','linewidth',1,'MarkerSize',10,'color','red',"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

xlabel({'\textbf{Sparsity} $(c)$'},'interpreter','latex','fontsize',35);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',35);
legend({'\textbf{Negative}','\textbf{Nonnegative}'},'interpreter','latex','fontsize',30);
%title({'$J(z)$ \textbf{for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',35);

subplot(2,2,4)
semilogy(1:length(zeta0),sqrt(Jz+Val_f),'bd','linewidth',1,'markersize',18), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

semilogy(1:length(zeta0),sqrt(Hz+Jz+Val_f),'rpentagram','linewidth',1,'markersize',15,"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

xlabel({'\textbf{Sparsity} (d)'},'interpreter','latex','fontsize',35);ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',35);
%legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',30);
%title({'$L_2$ \textbf{Error estimate for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',35);


end


%% Function used in this programm

function sgn_coeff = sgnfun(coeff0)

for k=1:length(coeff0)
    if coeff0(k) > 0
        sgn_coeff(k) = 1;
    elseif coeff0(k) < 0
        sgn_coeff(k) = -1;
    else
        sgn_coeff(k) = 0;
    end
end

end

function [xv,yv,iv]=example_polygon(example)

switch example
    case 1
        polygon_sides=[0.1 0; 0.7 0.2; 1 0.5; 0.75 0.85; 0.5 1; 0 0.25; 0.1 0];
        xv=polygon_sides(:,1); yv=polygon_sides(:,2);
        iv=length(xv); % This variable depends on the holes or not connected domain.
        % In these simple cases the domains are without holes and
        % domains are connected.
        
    case 2
        polygon_sides=(1/4)*[1 2; 1 0; 3 2; 3 0; 4 2; 3 3; 3 0.85*4; 2 4;
            0 3; 1 2];
        xv=polygon_sides(:,1); yv=polygon_sides(:,2);
        iv=length(xv); % This variable depends on the holes or not connected domain.
        % In these simple cases the domains are without holes and
        % domains are connected.
    case 3
        % fprintf('\n \t [POLYGON]: QUATERFOIL LIKE');
        warning off;
        M=129;
        th=linspace(0,2*pi,M);
        %th=(th(1:end-1))';
        polygon_sides=[cos(th').*(sin(2*th')) sin(th').*(sin(2*th'))];
        polygon_sides=polygon_sides(1:end-1,:);
        xv=polygon_sides(:,1); yv=polygon_sides(:,2);
        iv=length(xv); % This variable depends on the holes or not connected domain.
        % In these simple cases the domains are without holes and
        % domains are connected.
        
    case 4 % domain not simply connected (optics)
        Nsides=100;
        y=[0         0   -0.1184   -0.1184   -0.3761];
        r=[1.0000    0.6120    0.5663    1.0761    1.2810];
        th=linspace(0,2*pi,Nsides); th=(th(1:end-1))';
        C1=[0 y(1)]; P1v=C1+r(1)*[cos(th) sin(th)]; P1=polyshape(P1v);
        C2=[0 y(2)]; P2v=C2+r(2)*[cos(th) sin(th)]; P2=polyshape(P2v);
        C3=[0 y(3)]; P3v=C3+r(3)*[cos(th) sin(th)]; P3=polyshape(P3v);
        C4=[0 y(4)]; P4v=C4+r(4)*[cos(th) sin(th)]; P4=polyshape(P4v);
        C5=[0 y(5)]; P5v=C5+r(5)*[cos(th) sin(th)]; P5=polyshape(P5v);
        Pout=intersect(P1,P4);
        Pout=intersect(Pout,P5);
        Pin=union(P2,P3);
        xv=subtract(Pout,Pin);
        yv=[]; iv=[];
        
        
end

end

function g=define_function(funct_example,L)

% function to test

switch funct_example

    case 1 % test exactness hyperinterpolation
        nexp=L;
        c0=rand(1); c1=rand(1); c2=rand(1);
        g=@(x,y) (c0+c1*x+c2*y).^nexp;

    case 2 % test exactness filt. hyperinterpolation
        nexp=max(floor(L/2)-1,0);
        c0=rand(1); c1=rand(1); c2=rand(1);
        g=@(x,y) (c0+c1*x+c2*y).^nexp;

    case 3 % function of that type

        funct_example_sub=1;

        fstring='Not available';

        switch funct_example_sub
            case 1
                g=@(x,y) (1-x.^2-y.^2).*exp(x.*cos(y));
                fstring='(1-x.^2-y.^2).*exp(x.*cos(y))';
            case 2
                g=@(x,y) exp(-(x.^2+y.^2));
                fstring='exp(-(x.^2+y.^2))';
            case 3
                g=@(x,y) sin(-(x.^2+y.^2));
                fstring='sin(-(x.^2+y.^2))';
            case 4
                g=@(x,y) 1+0*x+0*y;
                fstring='1+0*x+0*y';
            case 5
                g=@(x,y) sqrt((x-0.5).^2+(y-0.5).^2);
            case 6
                g=@(x,y) (0.2*x+0.5*y).^19;
            case 7
                g=@(x,y) exp((x-0.5).^2+(y-0.5).^2);
            case 8
                g=@(x,y) exp(-100*((x-0.5).^2+(y-0.5).^2));
            case 9
                g=@(x,y) cos(30*(x+y));
            case 10
                g=@(x,y) cos(5*(x+y));
            case 11
                g=@(x,y) exp((x-0.5).^1+(y-0.5).^1);
            case 12
                g=@(x,y) exp((x-0.5).^3+(y-0.5).^3);
            case 13
                g=@(x,y) (0.2*x+0.5*y).^15;
            case 14
                g=@(x,y) 1./(x.^2+y.^2);
            case 15
                % g=@(x,y) (x+y).^ade;
                g=@(x,y) (1+x+0.5*y).^ade;
            case 16
                x0=0.5; y0=0.5;
                g=@(x,y) exp(-((x-x0).^2+(y-y0).^2));
            case 17
                x0=0.5; y0=0.5;
                g=@(x,y) ((x-x0).^2 + (y-y0).^2).^(3/2);
            case 18 % franke
                g=@(x,y) .75*exp(-((9*x-2).^2 + (9*y-2).^2)/4) + ...
                    .75*exp(-((9*x+1).^2)/49 - (9*y+1)/10) + ...
                    .5*exp(-((9*x-7).^2 + (9*y-3).^2)/4) - ...
                    .2*exp(-(9*x-4).^2 - (9*y-7).^2);
        end


end


end

function [pXYR_re, abs_err] = evaluate_disk(f,coeff,n)

t = -1:0.025:1;
[XR, YR] = meshgrid(t);

fXYR = feval(f, XR, YR);

% Compute orthonormal basis matrix at nodes.
jvec=1:(n+1)*(n+2)/2;
[U,~,Q,R,~,degs] = dORTHVAND(n,[XR YR],W,jvec,[],dbox);

VR=chebvand(LV,[XR YR],dbox);
pXYR0 = (VR(:,jvec)/R)*coeff;


pXYR = reshape(pXYR0,size(XR,1),size(XR,2));

abs_err = abs(fXYR-pXYR);
pXYR_re = pXYR;
end

function [noisy_function,real_function ]= noise_fun(f,sigma)

t = -1:0.025:1;
[XR, YR] = meshgrid(t);

fXYR = feval(f, XR, YR);

%val = (ones(size(XR))- XR.^2 - YR.^2) >=0;

%val = ()

%val = val./val;
var=sigma^2; pert_gauss=sqrt(var)*randn(size(fXYR));

noisy_function = (fXYR + pert_gauss);
real_function = fXYR;
end

function [Mx,My]= indomain_polygon(XR,YR,NR)

g1 = @(x) -x + 0.75;
%g2
g3 = @(x) x - 0.25;
%g4
g5 = @(x) 2*(x-0.75);
g6 = @(x) -(x-1) + 0.5; 
%g7 = @(x) 
g8 = @(x) -0.6*(x-0.75) + 0.85;
g9 = @(x) 0.5*x + 0.75;
    for j = 1:length(XR) 
        if XR(j) <= 0.25
            if YR(j) >= feval(g1,XR(j)) && YR(j) <= feval(g9,XR(j))
                YR(j) = YR(j);
            else
                YR(j) =NaN;
            end
        elseif XR(j) <= 0.5
            if YR(j) >= feval(g3,XR(j)) && YR(j) <= feval(g9,XR(j))
                YR(j) = YR(j);
            else
                YR(j) = NaN;
            end
        elseif XR(j) <= 0.75
            if YR(j) >= feval(g3,XR(j)) && YR(j) <= feval(g8,XR(j))
                YR(j) = YR(j);
            else
                YR(j) = NaN;
            end
        else
            if YR(j) >= feval(g5,XR(j)) && YR(j) <= feval(g6,XR(j))
                YR(j) = YR(j);
            else
                YR(j) = NaN;
            end
        end

    end
Mx = reshape(XR,NR,NR);
My = reshape(YR,NR,NR);
end
