
function [AEinfMV,AE2MV,beta0MV,lambdaV,XYZW,XYZWR,JzMV,HzMV,Wg_norm2]=...
    demo_cube(lambda_index,a,sigma,XYZW,XYZWR)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% Numerical experiment in "Hybrid hyperinterpolation over general regions".
% Region: unit-cube: [-1,1]^3.
%--------------------------------------------------------------------------
% Usage:
% >> demo_disk
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
LV=20;      % Hyperinterpolant tot degree.
NV=2*LV;    % Degree of the rule.
NR=50;      % Reference rule for computing L2 errors.

%--------------------------------------------------------------------------
% Noise and choice of lasso, hybrid, hard thresh. parameter.
%--------------------------------------------------------------------------

% Defining elastic parameter (knowning hyp. coeffs, it is the k-th in
% descending absolute value!)
if nargin<1,lambda_index=1;end

noise=1;            % 0: no noise, 1: noise (see parameters below)

% In the reference paper of Lasso (Table 1): a=0, sigma=0.2.
if noise
    if nargin <2,a=0.5;end     % defining impulse noise (in experiment 2)
    if nargin <3,sigma=0;end   % defining gaussian noise (in experiment 2)
else
    a=0; sigma=0; % no noise.
end

% * Function to approximate:
% 1. degree L poly., 2. degree floor(L/2)-1 poly. 3. test functions
% (see line 274 approx).
funct_example=3;

%--------------------------------------------------------------------------
% Special settings.
%--------------------------------------------------------------------------

% Number of tests for reconstructing functions of a type on this domain.
ntests=30;

% Show statistics table
display_stats=0;



% ........................ Main code below ................................



% ........ Numerical approximation, varying the degree in "nV" ............

AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics
JzMV=[]; HzMV=[];

for k=1:length(NV)
    N=NV(k); % Degree of the rule.
    L=LV(k); % Hyperinterpolant degree.

    % Define quadrature rule for hyperinterpolation at degree N.
    if nargin < 4, XYZW=cub_cube(N); end
    if isempty(XYZW), XYZW=cub_cube(N); end
    X=XYZW(:,1); Y=XYZW(:,2); Z=XYZW(:,3); W=XYZW(:,4);

    % Test points
    if nargin < 5, XYZWR=cub_cube(NR); end
    if isempty(XYZWR), XYZWR=cub_cube(NR); end
    XR=XYZWR(:,1); YR=XYZWR(:,2); ZR=XYZWR(:,3); WR=XYZWR(:,4);

    % Vandermonde matrix at nodes and polynomial degrees.
    [V,dbox,duples]=dCHEBVAND0(L,[X Y Z]);
    degs=sum(duples,2);

    % Testing hyperinterpolation errors for each "f" at "deg".

    poly_coeffs=[];
    lambdaV=[];

    for j=1:ntests

        % define function (see attached file at the bottom)
        g=choose_function(funct_example,L);

        % ... evaluate function to approximate ...
        gXYZ=feval(g,X,Y,Z);

        % ... Add noise (if present) ...

        % add impulse noise
        pert_impulse=0;
        if a > 0
            pert_impulse=a*(1-2*rand(length(gXYZ),1))*binornd(1,0.5);
            while norm(pert_impulse) == 0
                pert_impulse=a*(1-2*rand(length(gXYZ),1))*binornd(1,0.5);
            end
        end

        % add gaussian noise
        pert_gauss=0;
        if sigma > 0
            var=sigma^2;
            pert_gauss=sqrt(var)*randn(size(gXYZ));
            while norm(pert_gauss) == 0
                pert_gauss=sqrt(var)*randn(size(gXYZ));
            end
        end

        % add gaussian + impulse noise
        pert=pert_impulse+pert_gauss;

        % perturbed values
        gXYZ_pert=gXYZ+pert;

        % ... determine polynomial hyperinterpolant ...
        coeff0=(gXYZ_pert.*W)'*V; coeff0=coeff0';

        % test hyperinterpolant with or withour filters.

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

            gXYZR=feval(g,XR,YR,ZR);
            [VR]=dCHEBVAND0(L,[XR YR ZR]);
            pXYZR=VR*coeff;

            % errors
            AEinfV(ktest,j)=norm(gXYZR-pXYZR,inf); % absolute error (inf norm)
            AE2V(ktest,j)=sqrt(WR'*((gXYZR-pXYZR).^2)); % absolute error (2 norm)
            beta0V(ktest,j)=sum(abs(coeff) > 0);

            % evaluate J(coeff) and H(coeff), that are error relevant
            % parameters, as observed in Thm 5.1.

            JzV(ktest,j)=evaluate_J(coeff,coeff0);
            HzV(ktest,j)=evaluate_H(coeff,W,pert,V);

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
        fprintf('\n       ........ table at degree: %2.0f ........ \n \n ',...
            N);
        HypType=categorical({'tikhonov'; 'filtered'; 'lasso'; 'hybrid'; ...
            'hard'; 'hyperint.'});
        T = table(HypType,AEinfM,AE2M,beta0M,JzM,HzM); disp(T)
    end

    AEinfMV=[AEinfMV AEinfM]; AE2MV=[AE2MV AE2M]; beta0MV=[beta0MV beta0M];
    JzMV=[JzMV JzM]; HzMV=[HzMV HzM];

end

Wg_norm2=(norm(sqrt(W).*gXYZ,2))^2;


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






