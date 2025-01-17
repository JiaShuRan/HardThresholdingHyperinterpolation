function demo_cube_hard_0818
% Codes based on Alvise Sommariva (University of Padova)
% Date: 18 August, 2023
% demo_cube_plot is used to plot original function
% and its recovery function via different hyperinterpolants
% in the cube

clear all; clc;

LV=20;      % Hyperinterpolant tot degree.
NV=2*LV;    % Degree of the rule.
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
[V,dbox,duples]=dCHEBVAND0(LV,[X Y Z]);
degs=sum(duples,2);

% define function (see attached file at the bottom)
g=choose_function(funct_example,LV);

% ... evaluate function to approximate ...
gXYZ=feval(g,X,Y,Z);


sigma = 0.2; pert=sigma*randn(size(gXYZ));

% add Gaussian noise

gXYZ_pert=gXYZ+pert;

% ... determine polynomial hyperinterpolant ...
coeff0=(gXYZ_pert.*W)'*V; coeff0=coeff0';

% ... determine polynomial hyperinterpolant (noise) ...
noise0=(pert.*W)'*V; noise0=noise0';

% we reference to the method in "RBF approximation of noisy scattered data on the sphere"
lambdas = -15:0.1:7;

lambdak=2.^lambdas;

for k = 1:length(lambdas)
    lambdaL = lambdak(k);

    criteria(k) = lambdaL*ones(1,length(coeff0))*( abs(coeff0) > lambdaL) - 2*(sign(coeff0))'*(noise0.*( abs(coeff0) > lambdaL));


for ktest=1:3
            switch ktest
                case 1
                    hypermode='hyperinterpolation';
                    parms.lambda=[];
                    parms.mu=[];
                    parms.b=[];
                    coeff=coeff0;
                case 2
                    hypermode='lasso';
                    parms.lambda=lambdaL;
                    parms.mu=ones(size(coeff));
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 3
                    hypermode='hard';
                    parms.lambda=lambdaL;
                    parms.mu=[];
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
            end

            gXYZR=feval(g,XR,YR,ZR);
            [VR]=dCHEBVAND0(LV,[XR YR ZR]);
            pXYZR(:,ktest)=VR*coeff;

            % errors
            AEinfV(k,ktest)=norm(gXYZR-pXYZR(:,ktest),inf); % absolute error (inf norm)
            AE2V(k,ktest)=sqrt(WR'*((gXYZR-pXYZR(:,ktest)).^2)); % absolute error (2 norm)
            beta0V(k,ktest)=sum(abs(coeff) > 0);

            [pXYZR_3D, abs_err_3D] = evaluate_cube(g,coeff,LV);

            pXYZR_3D_plot(k,ktest).matrix = pXYZR_3D;
            abs_err_3D_plot(k,ktest).matrix = abs_err_3D;
end

end


%[xmin_1, ] = find(AE2V(:,1) == min(AE2V(:,1))); % Hyperinterpolation
[xmin_2, ] = find(AE2V(:,2) == min(AE2V(:,2))); % Lasso
[xmin_3, ] = find(AE2V(:,3) == min(AE2V(:,3))); % Hard

%[x_1, ] = find(AEinfV(:,1) == min(AEinfV(:,1))); % Hyperinterpolation
[x_2, ] = find(AEinfV(:,2) == min(AEinfV(:,2))); % Lasso
[x_3, ] = find(AEinfV(:,3) == min(AEinfV(:,3))); % Hard






AA =  pXYZR_3D_plot(1,1).matrix; 
AA_1 = abs_err_3D_plot(1,1).matrix;
BB =  pXYZR_3D_plot(xmin_2,2).matrix; lambda_2 = lambdak(xmin_2); Lambda_2 = lambdak(x_2);
BB_1 = abs_err_3D_plot(xmin_2,2).matrix;
CC =  pXYZR_3D_plot(min(xmin_3),3).matrix; lambda_3 = lambdak(min(xmin_3)); Lambda_3 = lambdak(min(x_3));
CC_1 = abs_err_3D_plot(min(xmin_3),3).matrix;


fprintf("....Variants.....lambda.....smallest $L_2$ error.....sparsity.....lambda.....smallest maximum error.....sparsity.....\n")
fprintf("{\\em{Hyperint.}} & $ - $ & $ %.6f $ & $ %d $ & $ - $ & $%.6f$ & $%d$ \n",AE2V(1,1),beta0V(1,1),AEinfV(1,1),beta0V(1,1))
fprintf("\n")
fprintf("{\\em{Lasso}} & $ %.6f $ & $ %.6f $ & $ %d $ & $ %.6f $ & $%.6f$ & $%d$ \n",lambda_2,AE2V(xmin_2,2),beta0V(xmin_2,2),Lambda_2,AEinfV(x_2,2),beta0V(x_2,2))
fprintf("\n")
fprintf("{\\em{Hard}} & $ %.6f $ & $ %.6f $ & $ %d $ & $ %.6f $ & $%.6f$ & $%d$ \n",lambda_3,AE2V(min(xmin_3),3),beta0V(min(xmin_3),3),Lambda_3,AEinfV(min(x_3),3),beta0V(min(x_3),3))
fprintf("\n")

It_matrix = [0,AE2V(1,1),beta0V(1,1),0,AEinfV(1,1),beta0V(1,1);lambda_2,AE2V(xmin_2,2),beta0V(xmin_2,2),Lambda_2,AEinfV(x_2,2),beta0V(x_2,2);lambda_3,AE2V(min(xmin_3),3),beta0V(min(xmin_3),3),Lambda_3,AEinfV(min(x_3),3),beta0V(min(x_3),3)];

L2_matrix = zeros(length(lambdak),5);
L2_matrix(:,1:3) = AE2V;  L2_matrix(:,4) = lambdak; L2_matrix(:,5) = criteria;

zz = find(criteria >0);

x_zero = zz(1);

save("cube_error_lambda_criteria.mat","L2_matrix");
save("cube_data.mat","It_matrix");

[XR_3, YR_3, ZR_3] = meshgrid(-1:0.1:1);


gXYZR_3D = feval(g, XR_3, YR_3, ZR_3);
gXYZR_pert_3D = gXYZR_3D + sigma*randn(size(gXYZR_3D));

%% Plot


figure(1)
loglog(lambdak,AE2V(:,2),'linewidth',3,'color','blue'), box on, 
set(gca, 'FontSize', 35, 'XMinorGrid', 'on'), set(gca, 'FontSize', 35, 'YMinorGrid', 'on'),
hold on,


loglog(lambdak,AE2V(:,3),'-.','linewidth',3,'color','black'), box on, 
set(gca, 'FontSize', 35, 'XMinorGrid', 'on'), set(gca, 'FontSize', 35, 'YMinorGrid', 'on'),
hold on,

loglog(lambda_2,min(AE2V(:,2)),'pentagram','linewidth',2,'MarkerSize',35,'color','red','MarkerFaceColor','red'), hold on,
loglog(lambdak(xmin_3),min(AE2V(:,3)),'diamond','linewidth',2,'MarkerSize',25,'color','red','MarkerFaceColor','red'), hold on,


xlabel({'\textbf{Regularization parameter} $\lambda$'},'interpreter','latex','fontsize',35);
ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',35);
legend({'$\|\mathcal{L}_n^{\lambda} f^{\epsilon} - f\|_2$','$\|\mathcal{H}_n^{\lambda} f^{\epsilon} - f\|_2$',['{\textbf{Minimum at}} $\lambda=$',num2str(lambda_2)],['{\textbf{Minimum at}} $\lambda\in [$',num2str(lambda_3), '$,$', num2str(lambdak(max(xmin_3))),'$]$']},'interpreter','latex','fontsize',40,'Location','northwest');


figure(2)

fontsize_baselinea = 10;
fontsize_baseline = 10;
fontsize_baselinet = 25;
marksize = 80;

colormap(jet)


xslice = [-.25,.5,1]; 
yslice = [0,1]; 
zslice = [-1,0];

%Primal and noisy function
axes('position',[0.03,0.55,0.2,0.4]) 
slice(XR_3,YR_3,ZR_3,gXYZR_3D,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('\textbf{Original function} $f$','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside')

axes('position',[0.03,0.05,0.2,0.4])
slice(XR_3,YR_3,ZR_3,gXYZR_pert_3D,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('\textbf{Noisy function} $f^{\epsilon}$','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'),caxis([0,0.7])

% Hyper. and its error
axes('position',[0.27,0.55,0.2,0.4])
slice(XR_3,YR_3,ZR_3,AA,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('$\mathcal{L}_nf^{\epsilon}$','interpreter','latex','fontsize',fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'), view(-36,15), colorbar('eastoutside'), caxis([0,0.7])

axes('position',[0.27,0.05,0.2,0.4])
slice(XR_3,YR_3,ZR_3,AA_1,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('\textbf{Absolute error}','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'), caxis([0,0.25])

% Lasso hyper. and its error
axes('position',[0.51,0.55,0.2,0.4])
slice(XR_3,YR_3,ZR_3,BB,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title(['$\mathcal{L}_n^{\lambda} f^{\epsilon}$ \textbf{at} $\lambda=$',num2str(lambda_2)], 'interpreter','latex','fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'), view(-36,15),colorbar('eastoutside'),caxis([0,.7])
 
axes('position',[0.51,0.05,0.2,0.4])
slice(XR_3,YR_3,ZR_3,BB_1,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),... 
title('\textbf{Absolute error}','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15),colorbar('eastoutside'), caxis([0,0.25])

% Hard hyper. and its error
axes('position',[0.75,0.55,0.2,0.4])
slice(XR_3,YR_3,ZR_3,CC,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title(['$\mathcal{H}_n^{\lambda} f^{\epsilon}$ \textbf{at} $\lambda=$',num2str(lambda_3)], 'interpreter','latex','fontsize', fontsize_baselinet),...
grid on,set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'),caxis([0,.7])

axes('position',[0.75,0.05,0.2,0.4])
slice(XR_3,YR_3,ZR_3,CC_1,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),... 
title('\textbf{Absolute error}','interpreter','latex', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15),colorbar('eastoutside'),caxis([0,0.25])


 end

%% Function used in this programm

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