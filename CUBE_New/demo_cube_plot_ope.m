function demo_cube_plot_ope
% Codes based on Alvise Sommariva (University of Padova)
% Date: 26 Sep, 2023
% demo_cube_plot_ope is used to verify Theorem 4.5 and
% lower degree approximation

LV=5;      % Hyperinterpolant tot degree.
%NV=2*LV;    % Degree of the rule.
NV=40;      % Quadrature exactness
NR=50;      % Reference rule for computing L2 errors.

% * Function to approximate:
% 1. degree L poly., 2. degree floor(L/2)-1 poly. 3. test functions
% (see line 274 approx).
funct_example=3;

% The degree d of the polynomial space is d = (LV+3)(LV+2)(LV+1)/6.
% The quadrature points N = (NV+2)^3/4.



% ........ Numerical approximation, varying the degree in "nV" ............

AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics

% Define quadrature rule for hyperinterpolation at degree N.
XYZW=cub_cube(NV); X=XYZW(:,1); Y=XYZW(:,2); Z=XYZW(:,3); W=XYZW(:,4);

% Test points
XYZWR=cub_cube(NR); XR=XYZWR(:,1); YR=XYZWR(:,2); ZR=XYZWR(:,3); WR=XYZWR(:,4);

% Vandermonde matrix at nodes and polynomial degrees.
[V1,dbox,duples]=dCHEBVAND0(LV,[X Y Z]);
degs=sum(duples,2);

% define function (see attached file at the bottom)
g=choose_function(funct_example,LV);

% ... evaluate function to approximate ...
gXYZ=feval(g,X,Y,Z);


sigma = 0.2; var=sigma^2; pert=sqrt(var)*randn(size(gXYZ));

% add Gaussian noise

gXYZ_pert=gXYZ+pert;

% ... determine polynomial hyperinterpolant ...
coeff1=(gXYZ_pert.*W)'*V1; coeff1=coeff1';

% ... determine polynomial hyperinterpolant (noise) ...
noise1=(pert.*W)'*V1; noise1=noise1';

lambdak=sort(abs(coeff1),'descend');
for k=1:length(coeff1)
    coeff1_new = coeff1.*(abs(coeff1) > lambdak(k));
    Jz(k) = 2*coeff1_new'*noise1 - coeff1_new'*coeff1_new;
    sgn_coeff = sgnfun(coeff1_new)';
    noise1_new = noise1.*(abs(coeff1) > lambdak(k));
    Hz(k) = sum(abs(sgn_coeff))*lambdak(k)^2 - 2*lambdak(k)*sgn_coeff'*noise1_new;
end

%Hz_great_0 = Hz.*(Hz >=0);

for p=1:length(Hz)
    if Hz(p) >= 0
        Hz_great_0(p) = Hz(p);
    else
        Hz_great_0(p) = NaN;
    end
end


Val_f = sum((gXYZ.^2).*W);
Val_f = ones(1,length(coeff1))*Val_f;

AErr_lasso = sqrt(Jz+Val_f+Hz);
AErr_hard = sqrt(Jz+Val_f);

% test hyperinterpolant with or withour filters.


LV_initial = 5;
LV_final = 20;
LV_step = 1;

for LV = LV_initial:LV_step:LV_final

    jtest = LV-4;

    % Vandermonde matrix at nodes and polynomial degrees.
    [V,dbox,duples]=dCHEBVAND0(LV,[X Y Z]);
    degs=sum(duples,2);


    % ... determine polynomial hyperinterpolant ...
    coeff0=(gXYZ_pert.*W)'*V; coeff0=coeff0';

    % ... compute the number of plus and multiplication for obtaining hyper coeff...
    complexity_pm(jtest) = (3*length(gXYZ_pert) - 1 )*length(coeff0);

    lambdas=sort(abs(coeff0),'descend');
    lambdaL=lambdas(8);

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

        gXYZR=feval(g,XR,YR,ZR);
        [VR]=dCHEBVAND0(LV,[XR YR ZR]);
        pXYZR(:,ktest)=VR*coeff;

        % errors
        AEinfV(jtest,ktest)=norm(gXYZR-pXYZR(:,ktest),inf); % absolute error (inf norm)
        AE2V(jtest,ktest)=sqrt(WR'*((gXYZR-pXYZR(:,ktest)).^2)); % absolute error (2 norm)
        beta0V(jtest,ktest)=sum(abs(coeff) > 0);

    end

end



%% Plot
% Step1
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
% xlabel({'\textbf{Degree} $(b)$'},'interpreter','latex','fontsize',30);ylabel({'\textbf{Total Operations}'},'interpreter','latex','fontsize',30);
% %legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',25);
% %title({'2'},'interpreter','latex','fontsize',35);
% 
% subplot(2,2,3)
% 
% plot(1:length(noise1),Hz,'-.o','linewidth',1,'MarkerSize',10,'color','k'), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
% hold on
% 
% plot(1:length(noise1),Hz_great_0,'o','linewidth',1,'MarkerSize',10,'color','red',"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
% hold on
% 
% xlabel({'\textbf{Sparsity} $(c)$'},'interpreter','latex','fontsize',35);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',35);
% legend({'\textbf{Negative}','\textbf{Nonnegative}'},'interpreter','latex','fontsize',30);
% %title({'$J(z)$ \textbf{for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',35);
% 
% subplot(2,2,4)
% semilogy(1:length(noise1),sqrt(Jz+Val_f),'bd','linewidth',1,'markersize',18), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
% hold on
% 
% semilogy(1:length(noise1),sqrt(Hz+Jz+Val_f),'rpentagram','linewidth',1,'markersize',15,"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
% set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
% hold on
% 
% xlabel({'\textbf{Sparsity} (d)'},'interpreter','latex','fontsize',35);ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',35);
% %legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',30);
% %title({'$L_2$ \textbf{Error estimate for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',35);

subplot(1,2,1)

plot(1:length(noise1),Hz,'-.o','linewidth',1,'MarkerSize',10,'color','k'), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

plot(1:length(noise1),Hz_great_0,'o','linewidth',1,'MarkerSize',10,'color','red',"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

xlabel({'\textbf{Sparsity} $(a)$'},'interpreter','latex','fontsize',35);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',35);
legend({'\textbf{Negative}','\textbf{Nonnegative}'},'interpreter','latex','fontsize',30);
%title({'$J(z)$ \textbf{for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',35);

subplot(1,2,2)
semilogy(1:length(noise1),sqrt(Jz+Val_f),'bd','linewidth',1,'markersize',18), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

semilogy(1:length(noise1),sqrt(Hz+Jz+Val_f),'rpentagram','linewidth',1,'markersize',15,"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

xlabel({'\textbf{Sparsity} (b)'},'interpreter','latex','fontsize',35);ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',35);
legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',30);
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

function g=choose_function(funct_example,L)

switch funct_example

    case 1 % test exactness hyperinterpolation
        nexp=L;
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

    case 2 % test exactness filt. hyperinterpolation
        nexp=max(floor(L/2)-1,0);
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

    case 3 % function of that type

        funct_example_sub=0;

        switch funct_example_sub
            case 0 % Test function used in the accompanying paper.
                g=@(x,y,z) exp(-1./(x.^2+y.^2+z.^2));
                % fstring='exp(-1./(x.^2+y.^2+z.^2))';
            case 1
                g=@(x,y,z) (1-x.^2-y.^2-z.^2).*exp(x.*cos(y));
                % fstring='(1-x.^2-y.^2-z.^2).*exp(x.*cos(y))';
            case 2
                g=@(x,y,z) exp((x.^6).*cos(y+2*z));
                % fstring='exp((x.^6).*cos(y+2*z))';
        end

end        
end

function [pXYZR_3D, abs_err_3D] = evaluate_cube(f,coeff,n)

    TT = -1:0.1:1;
     
    [XR_3, YR_3, ZR_3] = meshgrid(TT);

    fXYZR_3D = feval(f, XR_3, YR_3, ZR_3);

    [V,dbox,duples]=dCHEBVAND0(n,[XR_3(:) YR_3(:) ZR_3(:)]);

    pXYZR0 = V*coeff;

    pXYZR_3D = reshape(pXYZR0,size(XR_3,1),size(YR_3,1),size(ZR_3,1));

    abs_err_3D = abs(fXYZR_3D-pXYZR_3D);
end