function demo_sphtri_hard_0816

% This function is used to verify Theorem 4.5 and
% lower degree approximation over the spherical triangle 
% with vertices A=[1,0,0], B=[0,1,0] and C=[0,0,1].
% Date: 26 Sep, 2023
% Codes based on Alvise Sommariva (University of Padova)

clear all; clc;


% This function is used to plot the denoising effect of hyper. and its 
% variants over the spherical triangle with vertices A=[1,0,0], B=[0,1,0]
% and C=[0,0,1].
% Date: 20 Sep, 2023
% Codes based on Alvise Sommariva (University of Padova)


domain_example=0;

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%--------------------------------------------------------------------------
%nV=1:10;
LV=10;
NV=2*LV;
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

n = LV;
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

% ............................ Plot nodes .................................
domain_parms=vertices;
plot_s2(domain_example,domain_parms,pts,[],'',[]);

% ... Add noise (if present) ...
[g_pts_pert,pert,Noise_matrix] = add_noise(noise,a,sigma,g_pts);

% ... determine polynomial hyperinterpolant ...
[coeff0,R,jvec,dbox,degs] = dHYPERFIT2(LV,pts,weights,g_pts_pert,...
    [],[],domain_structure,dimpoly);

% ... determine polynomial hyperinterpolant (noise) ...
[noise0,R,jvec,dbox,degs] = dHYPERFIT2(LV,pts,weights,pert,...
    [],[],domain_structure,dimpoly);

if iscell(jvec), degs=degs(jvec{1}); else, degs=degs(jvec); end



% we reference to the method in "RBF approximation of noisy scattered data on the sphere"
lambdas = -15:0.1:7;

lambdak=2.^lambdas;

for k = 1:length(lambdas)
    lambdaL = lambdak(k);

    criteria(k) = lambdaL*ones(1,length(coeff0))*( abs(coeff0) > lambdaL) - 2*(sign(coeff0))'*(noise0.*( abs(coeff0) > lambdaL));

    %LV=NV;

    % test hyperinterpolant with or withour filters.
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

        % ... evaluate hyperinterpolant at initial pointset ...
        p_XR(:,ktest)=dPOLYVAL2(LV,coeff,XR,R,jvec,dbox,domain_structure,dimpoly);

        % ... estimating hyperinterpolant error ...
        g_XR=feval(g,XR(:,1),XR(:,2),XR(:,3));
        
        pt_eqphyper(k,ktest).matrix  = dPOLYVAL2(LV,coeff,XR,R,jvec,dbox,domain_structure,dimpoly);
        AEinfV(k,ktest)=norm(g_XR-p_XR(:,ktest),inf); % absolute error (inf norm)
        AE2V(k,ktest)=sqrt(WR'*((g_XR-p_XR(:,ktest)).^2)); % absolute error (2 norm)
        beta0V(k,ktest)=sum(abs(coeff) > 0);

    end

end


%% Plot denoising effect

%[xmin_1, ] = find(AE2V(:,1) == min(AE2V(:,1))); % Hyperinterpolation
[xmin_2, ] = find(AE2V(:,2) == min(AE2V(:,2))); % Lasso
[xmin_3, ] = find(AE2V(:,3) == min(AE2V(:,3))); % Hard

%[x_1, ] = find(AEinfV(:,1) == min(AEinfV(:,1))); % Hyperinterpolation
[x_2, ] = find(AEinfV(:,2) == min(AEinfV(:,2))); % Lasso
[x_3, ] = find(AEinfV(:,3) == min(AEinfV(:,3))); % Hard

AA =  pt_eqphyper(1,1).matrix; 
BB =  pt_eqphyper(xmin_2,2).matrix; lambda_2 = lambdak(xmin_2); Lambda_2 = lambdak(x_2);
CC_1 =  pt_eqphyper( min(xmin_3), 3).matrix; lambda_3 = lambdak(min(xmin_3)); Lambda_3 = lambdak(min(x_3));


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

save("sphere_impulse_gasuu_pert.mat","Noise_matrix");
save("sphtri_error_lambda_criteria.mat","L2_matrix");
save("sphtri_data.mat","It_matrix");


figure(2)
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


figure(3)

fontsize_baselinet = 35;

% C. compute triangulation (valid only for sph.poly. lying in some hemisphere)

% 1. rotate to North Pole

CC=mean(vertices); CC=CC/norm(CC);

% ................ rotation matrix centroid to north pole .................

[az,el,r] = cart2sph(CC(1),CC(2),CC(3));
phi=az; theta=pi/2-el;
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta);
R1=[ct 0 -st; 0 1 0; st 0 ct]; R2=[cp sp 0; -sp cp 0; 0 0 1];
rotmat=R1*R2; inv_rotmat=rotmat';

% ........................ rotate vertices to north pole...................

vertices_NP=(rotmat*vertices')';

% ................... stereographic map from south pole ...................

% ....... vertices .......

XX_SP=vertices_NP(:,1); YY_SP=vertices_NP(:,2); ZZ_SP=vertices_NP(:,3);
rat=1./(1+ZZ_SP);

XX_SPm=rat.*XX_SP; YY_SPm=rat.*YY_SP; % vertices on the plane

% ....... points .......

rat2=1./(1+XR(:,3));
X_SP=rat2.*XR(:,1); % points on the plane
Y_SP=rat2.*XR(:,2);

% ...... triangulation on the plane ....

tri = delaunay(X_SP,Y_SP);

% exact function and its noisy function
axes('position',[0.025,0.55,0.2,0.47])
fg = trisurf(tri,XR(:,1), XR(:,2), XR(:,3), g_XR,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('\textbf{Original function} $f$','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'), %caxis([min(g_XR)-0.02,max(g_XR)+0.02])
axis off

% ....... points for noisy function .......

rat21=1./(1+pts(:,3));
X_SP1=rat21.*pts(:,1); % points on the plane
Y_SP1=rat21.*pts(:,2);

% ...... triangulation on the plane ....

tri1 = delaunay(X_SP1,Y_SP1);

axes('position',[0.025 0.05 0.2 0.47])
fg = trisurf(tri1, pts(:,1), pts(:,2), pts(:,3), g_pts_pert,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('\textbf{Noisy function}  $f^{\epsilon}$','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),%caxis([min(g_pts_pert)-0.02,max(g_pts_pert)+0.02]) %caxis([-3.0,3.0]),
axis off

% Hyper. and its error

axes('position',[0.275 0.55 0.2 0.47])

fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), AA,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('$\mathcal{L}_nf^{\epsilon}$','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1]),
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(g_XR),max(g_XR)]),
axis off

axes('position',[0.275,0.05,0.2,0.47])
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), (AA-g_XR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('\textbf{Absolute error}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(BB-g_XR)),max(abs(BB-g_XR))])
axis off

% Lasso hyper. and its error
axes('position',[0.525 0.55 0.2 0.47])
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), BB,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title(['$\mathcal{L}_n^{\lambda} f^{\epsilon}$ \textbf{at} $\lambda=$',num2str(lambda_2)], 'interpreter','latex','fontsize', fontsize_baselinet, 'position',[0,0,1]),
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(g_XR),max(g_XR)]),
axis off

axes('position',[0.525,0.05,0.2,0.47])
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), abs(BB-g_XR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('\textbf{Absolute error}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(BB-g_XR)),max(abs(BB-g_XR))])
axis off

% Hard thresholding hyper. and its error
axes('position',[0.775 0.55 0.2 0.47])    
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), CC_1,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title(['$\mathcal{H}_n^{\lambda} f^{\epsilon}$ \textbf{at} $\lambda=$',num2str(lambda_3)], 'interpreter','latex','fontsize', fontsize_baselinet, 'position',[0,0,1]),
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(g_XR),max(g_XR)]),
axis off

axes('position',[0.775,0.05,0.2,0.47])    
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), abs(CC_1-g_XR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('\textbf{Absolute error}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(BB-g_XR)),max(abs(BB-g_XR))])
axis off


end



%% Function used in this program


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

function [g_pts_pert,pert,Noise_matrix] = add_noise(noise,a,sigma,g_pts)
if noise
    % add impulse noise
    pert_impulse=0;
    if a > 0
        pert_impulse=a*(1-2*rand(length(g_pts),1)).*binornd(1,0.5,length(g_pts),1);
        while norm(pert_impulse) == 0
            pert_impulse=a*(1-2*rand(length(g_pts),1)).*binornd(1,0.5,length(g_pts),1);
        end
    end

    % add gaussian noise
    pert_gauss=0;
    if sigma > 0
        pert_gauss=sigma*randn(size(g_pts));
        while norm(pert_gauss) == 0
            pert_gauss=sigma*randn(size(g_pts));
        end
    end

    % add gaussian + impulse noise
    pert=pert_impulse+pert_gauss;

    % perturbed values
    g_pts_pert=g_pts+pert;

    Noise_matrix = zeros(length(pert_gauss),3);
    Noise_matrix(:,1) = pert_impulse; Noise_matrix(:,2) = pert_gauss;   Noise_matrix(:,3) = pert;

else
    g_pts_pert=g_pts;
end
end
