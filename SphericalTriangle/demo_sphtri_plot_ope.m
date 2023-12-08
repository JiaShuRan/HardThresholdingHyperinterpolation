function demo_sphtri_plot_ope

% This function is used to verify Theorem 4.5 and
% lower degree approximation over the spherical triangle 
% with vertices A=[1,0,0], B=[0,1,0] and C=[0,0,1].
% Date: 26 Sep, 2023
% Codes based on Alvise Sommariva (University of Padova)


clear; clf;

domain_example=0;

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%--------------------------------------------------------------------------
%nV=1:10;
LLV=8;
NV=11;
NR=30;
%--------------------------------------------------------------------------
% Noise and lasso parameter.
%--------------------------------------------------------------------------

noise=1;          % 0: no noise, 1: noise
a=0.1;              % defining impulse noise (in experiment 2)
sigma=0.2; %0.02;   % defining gaussian noise (in experiment 2)
%lambda=10^(-2);   % defining lasso parameter
pos=0;            % extraction type.
domain_structure.domain='spherical-triangle';

% ....... Special settings .......

% Approximation type parameter "pts_type".
%     case 1, pts_type='Hyperinterpolation full set';
%     case 2, pts_type='Hyperinterpolation compressed set';
%
% Note: the case "2" should be used for mild values of "n", say at most 15.
pts_type=1;

% Plot domain and nodes: do_plot=1 (yes), do_plot=0 (no).
do_plot=1;

% testing functions:
% 1. polynomial of degree n,
% 2. polynomial of degree floor(n/2)-1
% 3. gaussian like exponential
%funct_example=7;

funct_example=7;

% ....... Apply settings to define domain, pointsets, and functions .......

% Domain
vertices=define_domain(domain_example);

% Test points
%[XYZWR,dbox]=define_cub_rule(domain_structure,30);
P1=vertices(1,:); P2=vertices(2,:); P3=vertices(3,:);
XYZWR = cub_sphtri(NR,P1',P2',P3',pos);

XR=XYZWR(:,1:end-1); WR=XYZWR(:,end);

% ........ Numerical approximation, varying the degree in "nV" ............
fprintf('\n \t ');
AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics

n = NV;
dimpoly=(n+1)^2;

% ... extract hyperinterpolation set (notice that ade=2*n) ...
if pts_type == 2 % compressed set
    % fprintf('\n \t * Compressed set')
    XYZW = cub_sphtri(2*n,P1',P2',P3',pos);
    [pts,weights,momerr,dbox] =...
        dCATCH(2*n,XYZW(:,1:3),XYZW(:,4));
else % full set
    % fprintf('\n \t * Full set')
    XYZW = cub_sphtri(2*n,P1',P2',P3',pos);
    pts=XYZW(:,1:end-1); weights=XYZW(:,end);
end

% .. testing AE_L2err hyperinterpolation error for each "f" at "deg" ..
g = define_function(funct_example);

% ... evaluate function to approximate ...
g_pts=feval(g,pts(:,1),pts(:,2),pts(:,3));

% ... Add noise (if present) ...
[g_pts_pert,pert] = add_noise(noise,a,sigma,g_pts);

% ... determine polynomial hyperinterpolant ...
[coeff0,R,jvec,dbox,degs] = dHYPERFIT2(LLV,pts,weights,g_pts_pert,...
    [],[],domain_structure,dimpoly);

[zeta0,R,jvec,dbox,degs] = dHYPERFIT2(LLV,pts,weights,pert,...
    [],[],domain_structure,dimpoly);

Val_f = sum((g_pts.^2).*weights);
Val_f = ones(1,length(coeff0))*Val_f;


if iscell(jvec), degs=degs(jvec{1}); else, degs=degs(jvec); end

lambdas=sort(abs(coeff0),'descend');


for k = 1:length(zeta0)
    coeff0_new = coeff0.*(abs(coeff0) > lambdas(k));
    Jz(k) = 2*coeff0_new'*zeta0 - coeff0_new'*coeff0_new;
    sgn_coeff = sgnfun(coeff0_new)';
    zeta0_new = zeta0.*(abs(coeff0) > lambdas(k));
    Hz(k) = sum(abs(sgn_coeff))*lambdas(k)^2 - 2*lambdas(k)*sgn_coeff'*zeta0_new;
end


AErr_lasso = sqrt(Jz+Val_f+Hz);
AErr_hard = sqrt(Jz+Val_f);

for p=1:length(Hz)
    if Hz(p) >= 0
        Hz_great_0(p) = Hz(p);
    else
        Hz_great_0(p) = NaN;
    end
end



LV_initial = 5;
LV_final = 11;
LV_step = 1;


for LV = LV_initial:LV_step:LV_final
    jtest = LV-4;
    % ... determine polynomial hyperinterpolant ...
    [coeff0,R,jvec,dbox,degs] = dHYPERFIT2(LV,pts,weights,g_pts_pert,...
        [],[],domain_structure,dimpoly);

    % ... compute the number of plus and multiplication for obtaining hyper coeff...
    complexity_pm(jtest) = (3*length(g_pts_pert) - 1 )*length(coeff0);

    lambdas=sort(abs(coeff0),'descend');
    lambdaL=lambdas(6);
   
    % test hyperinterpolant with or withour filters.
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
                parms.w=weights;
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

        % ... evaluate hyperinterpolant at initial pointset ...
        p_XR(:,ktest)=dPOLYVAL2(LV,coeff,XR,R,jvec,dbox,domain_structure,dimpoly);

        % ... estimating hyperinterpolant error ...
        g_XR=feval(g,XR(:,1),XR(:,2),XR(:,3));

        AEinfV(jtest,ktest)=norm(g_XR-p_XR(:,ktest),inf); % absolute error (inf norm)
        AE2V(jtest,ktest)=sqrt(WR'*((g_XR-p_XR(:,ktest)).^2)); % absolute error (2 norm)
        beta0V(jtest,ktest)=sum(abs(coeff) > 0);

    end

end


%% Plot nodes
% domain_parms=vertices;
% figure(1)
% plot_s2(domain_example,domain_parms,pts,[],'',[]);


%% Step 3

GG = LV_initial:LV_step:LV_final;

figure(1)

% subplot(2,2,1)
% semilogy(GG,AE2V(:,5),'d','linewidth',1,'MarkerSize',22,'color','blue'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
% hold on
% 
% semilogy(GG,AE2V(:,3),'pentagram','linewidth',1,'MarkerSize',20,'color','red','MarkerFaceColor','red'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
% hold on
% 
% xlabel({'\textbf{Degree} $(a)$'},'interpreter','latex','fontsize',30);ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',30);
% legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',25);
% %title({'1'},'interpreter','latex','fontsize',35);
% 
% subplot(2,2,2)
% semilogy(GG,complexity_pm,'o','linewidth',1,'MarkerSize',18,'color','blue'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
% hold on
% 
% % semilogy(GG,A2(:,3),'pentagram','linewidth',1,'MarkerSize',20,'color','red','MarkerFaceColor','red'), box on, set(gca,'fontsize',16),
% % set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),%axis([0,TT-1,-1.5,0.5]),
% % hold on
% 
% xlabel({'\textbf{Degree} $(b)$'},'interpreter','latex','fontsize',30);ylabel({'\textbf{Total Operations}'},'interpreter','latex','fontsize',30);
% %legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',25);
% %title({'2'},'interpreter','latex','fontsize',35);
% subplot(2,2,3)
% 
% plot(1:length(zeta0),Hz,'-.o','linewidth',1,'MarkerSize',10,'color','k'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
% hold on
% 
% plot(1:length(zeta0),Hz_great_0,'o','linewidth',1,'MarkerSize',10,'color','red',"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
% hold on
% 
% xlabel({'\textbf{Sparsity} $(c)$'},'interpreter','latex','fontsize',35);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',35);
% legend({'\textbf{Negative}','\textbf{Nonnegative}'},'interpreter','latex','fontsize',30);
% 
% xlabel({'\textbf{Sparsity} $(c)$'},'interpreter','latex','fontsize',24);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',24);
% %legend({'J function','H function'},'fontsize',16.5);
% %title({'$J(z)$ \textbf{for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',30);
% 
% subplot(2,2,4)
% semilogy(1:length(zeta0),sqrt(Jz+Val_f),'bd','linewidth',1,'markersize',18), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
% hold on
% 
% semilogy(1:length(zeta0),sqrt(Hz+Jz+Val_f),'rpentagram','linewidth',1,'markersize',15,"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
% hold on
% 
% xlabel({'\textbf{Sparsity} $(d)$'},'interpreter','latex','fontsize',24);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',24);
% %legend({'Hard','Lasso'},'fontsize',16.5);
% %title({'\textbf{Error estimate for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',30);


subplot(1,2,1)

plot(1:length(zeta0),Hz,'-.o','linewidth',1,'MarkerSize',10,'color','k'), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

plot(1:length(zeta0),Hz_great_0,'o','linewidth',1,'MarkerSize',10,'color','red',"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

xlabel({'\textbf{Sparsity} $(a)$'},'interpreter','latex','fontsize',35);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',35);
legend({'\textbf{Negative}','\textbf{Nonnegative}'},'interpreter','latex','fontsize',30);

%legend({'J function','H function'},'fontsize',16.5);
%title({'$J(z)$ \textbf{for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',30);

subplot(1,2,2)
semilogy(1:length(zeta0),sqrt(Jz+Val_f),'bd','linewidth',1,'markersize',18), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

semilogy(1:length(zeta0),sqrt(Hz+Jz+Val_f),'rpentagram','linewidth',1,'markersize',15,"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

xlabel({'\textbf{Sparsity} $(b)$'},'interpreter','latex','fontsize',35);ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',35);
legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',30);
%legend({'Hard','Lasso'},'fontsize',16.5);
%title({'\textbf{Error estimate for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',30);


end



%% Function used in this program

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

function vertices=define_domain(example)


switch example

    case 0
        % large domain
        vertices=[ 1 0 0;
            0 1 0;
            0 0 1];

    case 1
        a=0.5; b=0.1; % MEDIUM-LARGE SIZE (Australia)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];

end

end


function g = define_function(funct_example)

switch funct_example

    case 1 % test exactness hyperinterpolation
        nexp=n;
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

    case 2 % test exactness filt. hyperinterpolation
        nexp=floor(n/2);
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

    case 3 % exponential type
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) exp(-(c0*x.^2+c1*y.^2+c2*z.^2+c3));

    case 4
        g=@(x,y,z) 1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2;

    case 5
        g=@(x,y,z) cos(10*(x+y+z));

    case 6
        g=@(x,y,z) exp(-(x.^2+y.^2+(z-1).^2));

    case 7
        g=@(x,y,z) exp(-((x-1/sqrt(3)).^2 + (y-1/sqrt(3)).^2+ (z-1/sqrt(3)).^2 ) );

    case 8 % Test function used in the accompanying paper.
        g=@(x,y,z) 0.75*exp(-((9*x-2).^2)/4- ((9*y-2).^2)/4- ((9*z-2).^2)/4)+...
                   0.75*exp(-((9*x+1).^2)/49- ((9*y+1).^2)/10- ((9*z+1).^2)/10)+...
                   0.5*exp(-((9*x-7).^2)/4-((9*y-3).^2)/4-((9*z-5).^2)/4)-...
                   0.2*exp(-((9*x-4).^2) - ((9*y-7).^2)-((9*z-5).^2));

end
end

function [g_pts_pert,pert] = add_noise(noise,a,sigma,g_pts);
if noise
    % add impulse noise
    pert_impulse=0;
    if a > 0
        pert_impulse=a*(1-2*rand(length(g_pts),1))*binornd(1,0.5);
        while norm(pert_impulse) == 0
            pert_impulse=a*(1-2*rand(length(g_pts),1))*binornd(1,0.5);
        end
    end

    % add gaussian noise
    pert_gauss=0;
    if sigma > 0
        var=sigma^2;
        pert_gauss=sqrt(var)*randn(size(g_pts));
        while norm(pert_gauss) == 0
            pert_gauss=sqrt(var)*randn(size(g_pts));
        end
    end

    % add gaussian + impulse noise
    pert=pert_impulse+pert_gauss;

    % perturbed values
    g_pts_pert=g_pts+pert;
else
    g_pts_pert=g_pts;
end
end
