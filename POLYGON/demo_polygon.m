
function [AEinfMV,AE2MV,beta0MV,lambdaV,XYW,XYWR,JzMV,HzMV,Wg_norm2]=...
    demo_polygon(lambda_index,a,sigma,XYW,XYWR)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% Numerical experiment in "Hybrid hyperinterpolation over general regions".
% Region: polygon.
%--------------------------------------------------------------------------
% Usage:
% >> demo_polygon
%--------------------------------------------------------------------------
% Note:
% The routine uses 'binornd' that requires Statistics and Machine Learning
% Toolbox.
%--------------------------------------------------------------------------
% Dates:
% Written on January 1, 2023: A. Sommariva.
% Modified on April 22, 2023: A. Sommariva.
%--------------------------------------------------------------------------
% COPYRIGHT
%--------------------------------------------------------------------------
% Copyright (C) 2023- 
%
% Authors:
% Alvise Sommariva
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%--------------------------------------------------------------------------
LV=15;     % Hyperinterpolant tot degree.
NV=2*LV;   % Degree of the rule used in hyperinterpolation.
NR=40;     % Degree of precision of the reference rule (estimate L2 error).


%--------------------------------------------------------------------------
% Noise and choice of lasso, hybrid, hard threshold parameter.
%--------------------------------------------------------------------------

noise=1;

if noise
    if nargin <2,a=0.02;end     % defining impulse noise (in experiment 2)
    if nargin <3,sigma=0.02;end % defining gaussian noise (in experiment 2)
else
    a=0; sigma=0; % no noise.
end

% * Function to approximate:
% 1. degree L poly., 2. degree floor(L/2)-1 poly. 3. test functions
% (see line 265 approx).
funct_example=3;

% No table or stats.
display_stats=0;

%--------------------------------------------------------------------------
% Special settings.
%--------------------------------------------------------------------------

% Plot domain and nodes: do_plot=1 (yes), do_plot=0 (no).
do_plot=0;

% Number of tests for reconstructing functions of a type on this domain.
ntests=1;

if nargin < 1, lambda_index=10; end



% ........................ Main code below ................................


% ......................... Define domain .................................

[xv,yv]=example_polygon(2);


% ........ Numerical approximation, varying the degree in "nV" ............

AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics
JzMV=[]; HzMV=[];

for k=1:length(NV)
    N=NV(k); % Quadrature points.
    L=LV(k); % Hyperinterpolant degree.


    % Define quadrature rule for hyperinterpolation at L, with  ADE=N.
    if nargin < 4, XYW=cub_polygon(N,xv,yv); end
    if isempty(XYW), XYW=cub_polygon(N,xv,yv); end
    X=XYW(:,1); Y=XYW(:,2); W=XYW(:,3);

    % Test points
    if nargin < 5, XYWR=cub_polygon(NR,xv,yv); end
    if isempty(XYWR), XYWR=cub_polygon(NR,xv,yv); end
    XR=XYWR(:,1); YR=XYWR(:,2); WR=XYWR(:,3);


    % Compute bounding box
    xmin=min(xv); xmax=max(xv);
    ymin=min(yv); ymax=max(yv);
    dbox=[xmin xmax ymin ymax];

    % Compute orthonormal basis matrix at nodes.
    jvec=1:(L+1)*(L+2)/2;
    [U,~,Q,R,~,degs] = dORTHVAND(L,[X Y],W,jvec,[],dbox);

    % .. testing AE_L2err hyperinterpolation error for each "f" at "deg" ..
    poly_coeffs=[];
    lambdaV=[];

    for j=1:ntests

        % ... define function to approximate ...
        g=define_function(funct_example,L);

        % ... evaluate function to approximate ...
        gXY=feval(g,X,Y);

        % ... add noise (if present) ...

        % a) add impulse noise
        pert_impulse=0;
        if a > 0
            pert_impulse=a*(1-2*rand(length(gXY),1))*binornd(1,0.5);
            while norm(pert_impulse) == 0
                pert_impulse=a*(1-2*rand(length(gXY),1))*binornd(1,0.5);
            end
        end

        % b) add gaussian noise
        pert_gauss=0;
        if sigma > 0
            var=sigma^2;
            pert_gauss=sqrt(var)*randn(size(gXY));
            while norm(pert_gauss) == 0
                pert_gauss=sqrt(var)*randn(size(gXY));
            end
        end

        % add gaussian + impulse noise
        pert=pert_impulse+pert_gauss;

        % perturbed values
        gXY_pert=gXY+pert;

        % ... determine polynomial hyperinterpolant ...
        % compute hyperinterpolant coefficients
        coeff0=Q'*(sqrt(W).*gXY_pert);

        % ... test hyperinterpolant with or withour filters ...

        lambdas=sort(abs(coeff0),'descend');
        lambdaL=lambdas(lambda_index);

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
            VR=chebvand(L,[XR YR],dbox);
            pXYR = (VR(:,jvec)/R)*coeff;

            % errors
            AEinfV(ktest,j)=norm(gXYR-pXYR,inf); % absolute error (inf norm)
            AE2V(ktest,j)=sqrt(WR'*((gXYR-pXYR).^2)); % absolute error (2 norm)
            beta0V(ktest,j)=sum(abs(coeff) > 0);
                        % evaluate J(coeff) and H(coeff), that are error relevant 
            % parameters, as observed in Thm 5.1.

            JzV(ktest,j)=evaluate_J(coeff,coeff0);
            HzV(ktest,j)=evaluate_H(coeff,W,pert,U);

        end

        lambdaV=[lambdaV lambdas];

    end

    % averages of the errors (vectors 5 x 1)
    AEinfM=mean(AEinfV,2);
    AE2M=mean(AE2V,2);
    beta0M=mean(beta0V,2);
    JzM=mean(JzV,2);
    HzM=mean(HzV,2);

    if display_stats
    fprintf('\n       ........ table at degree: %2.0f ........ \n \n ',N);
    HypType=categorical({'tikhonov'; 'filtered'; 'lasso'; 'hybrid'; ...
        'hard'; 'hyperint.'});
    T = table(HypType,AEinfM,AE2M,beta0M,JzM,HzM); disp(T)
    end
    
    AEinfMV=[AEinfMV AEinfM]; AE2MV=[AE2MV AE2M]; beta0MV=[beta0MV beta0M];
    JzMV=[JzMV JzM]; HzMV=[HzMV HzM];

end

Wg_norm2=(norm(sqrt(W).*gXY,2))^2;







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





function Jz=evaluate_J(z,alpha)

%--------------------------------------------------------------------------
% Object:
% Evaluate function J(z)=sum( (z(l))^2-2*z(l)*alpha(l) )
%--------------------------------------------------------------------------
% Input:
% z    : vector of dimension d x 1
% alpha: vector of dimension d x 1
%--------------------------------------------------------------------------
% Output:
% Jz: value of J(z)=sum( (z(l))^2-2*z(l)*alpha(l) )
%--------------------------------------------------------------------------
% Reference:
% Quantity relevant in Thm. 5.1 of the paper
% "Hybrid hyperinterpolation over general regions"
% Congpei An · Alvise Sommariva · Jia-Shu Ran
%--------------------------------------------------------------------------

% Jz=sum(z.^2-2*z.*alpha);
Jz=z'*z -2*z'*alpha;




function Hz=evaluate_H(z,w,err,V)

%--------------------------------------------------------------------------
% Object:
% Evaluate function H(z)=2*sum_l( z(l) * sum_j( w(j)*err(j)*V(l,j) ) )
%--------------------------------------------------------------------------
% Input:
% z    : vector of dimension d x 1
% w    : vector of dimension N x 1
% err  : vector of dimension N x 1
% V    : matrix of dimension d x N
%--------------------------------------------------------------------------
% Output:
% Hz: value of H(z)=2*sum_l( z(l) * sum_j( w(j)*err(j)*V(l,j) ) )
%--------------------------------------------------------------------------
% Reference:
% Quantity relevant in Thm. 5.1 of the paper
% "Hybrid hyperinterpolation over general regions"
% Congpei An · Alvise Sommariva · Jia-Shu Ran
%--------------------------------------------------------------------------

inner_term=V'*(w.*err);
outer_term=z'*inner_term;
Hz=2*outer_term;





