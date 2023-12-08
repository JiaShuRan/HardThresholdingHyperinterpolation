
function [xyw,xywc,xvc,yvc,pgon,tri]=cub_polygon(ade,xv,yv,iv)

%--------------------------------------------------------------------------
% Important: this "cub_polygon" version requires at least Matlab 9.3.0.
%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
% ade: algebraic degree of precision of the rule
% xvc: abscissae of the vertices.
% yvc: ordinates of the vertices.
% iv: index of the polygon components. Example: if part of the boundary is
%     made by the first 5 coordinates and two other non connected
%     components by other 10 and 12 coordinates, then iv=[5,10,12].
%
% Important: 
% 1. differently from classical codes on polygons, the first and
% last vertex must NOT be equal.
% 2. if the polygon is of polyshape class, set "xv" such polygon.
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% xyw: N x 3 matrix, where (x,y) with x=xyw(:,1), y=xyw(:,2) are the nodes
%      and w=xyw(:,3) are the weights.
% xvc, yvc: if "xv" is a polyshape object then "xvc","yvc" are column vectors,
%    otherwise they are cell arrays.
% pgon: polygon coded in polyshape form.
% tri: polygon triangulation in polyshape form.
%--------------------------------------------------------------------------
% Examples:
%--------------------------------------------------------------------------
% >> P=[0 0; 1 0; 1 1; 0 1];
% >> pgon=polyshape(P);
% >> ade=2;
% >> [xyw,xvc,yvc,pgon,tri]=polygauss_2018(ade,pgon)
% xyw =
%     0.1667    0.6667    0.1667
%     0.1667    0.1667    0.1667
%     0.6667    0.1667    0.1667
%     0.3333    0.8333    0.1667
%     0.8333    0.3333    0.1667
%     0.8333    0.8333    0.1667
%
% xvc =
%      0
%      0
%      1
%      1
%      0
% 
% yvc =
%      0
%      1
%      1
%      0
%      0
%
% pgon = 
%   polyshape with properties:
%       Vertices: [4×2 double]
%     NumRegions: 1
%       NumHoles: 0
% tri = 
%   triangulation with properties:
%               Points: [4×2 double]
%     ConnectivityList: [2×3 double]
%--------------------------------------------------------------------------
% >> xv=[0 1 1 0]; % square as vertices
% >> yv=[0 0 1 1];
% >> ade=2;
% >> [xyw,xvc,yvc,pgon,tri]=polygauss_2018(ade,xv,yv)
% xyw =
%     0.1667    0.6667    0.1667
%     0.1667    0.1667    0.1667
%     0.6667    0.1667    0.1667
%     0.3333    0.8333    0.1667
%     0.8333    0.3333    0.1667
%     0.8333    0.8333    0.1667
% 
% xvc =
% 
%   1×1 cell array
%     {4×1 double}
% 
% yvc =
%   1×1 cell array
%     {4×1 double}
% 
% pgon = 
% 
%   polyshape with properties:
%       Vertices: [4×2 double]
%     NumRegions: 1
%       NumHoles: 0
% 
% tri = 
% 
%   triangulation with properties:
%               Points: [4×2 double]
%     ConnectivityList: [2×3 double]
%--------------------------------------------------------------------------
% Routines required..
%--------------------------------------------------------------------------
% 1. cub_triangle (external)
% 2. rule_conversion (attached here)
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2007-2018 Alvise Sommariva, Marco Vianello.
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Author:  Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: December 07, 2018.
%--------------------------------------------------------------------------


S=class(xv);

S_double = strcmp(S,'double');

% TROUBLESHOOTING.
if S_double % distinguishing polyshape class from double.
    if nargin < 3
        pts=xv;
        if size(pts,1) <= size(pts,2) % working with column vectors.
            pts=xv';
        end
        if pts(1,:) == pts(end,:)
            pts=pts(1:end-1,:);
        end
        xv=pts(:,1); yv=pts(:,2); iv=length(xv);
    else
        if size(xv,1) <= size(xv,2) % working with column vectors.
            xv=xv';
        end
        if size(yv,1) <= size(yv,2) % working with column vectors.
            yv=yv';
        end
        %     if [xv(1) yv(1)] == [xv(end) yv(end)]
        %         xv=xv(1:end-1); yv=yv(1:end-1);
        %     end
        if nargin < 4
            iv=length(xv);
        else
            if isempty(iv)
                iv=length(xv);
            end
        end
    end
    
    % Here we convert the given information into cells of coordinates, where
    % each component represent a geometrical component of the polygon, not
    % connected or at most tangent to other geometrical components.
    
    xvc={xv(1:iv(1))}; yvc={yv(1:iv(1))};
    ii1=1; ii2=iv(1);
    for ii=2:length(iv)
        ii1=ii2+1;
        ii2=ii2+iv(ii);
        xvl=xv(ii1:ii2);
        yvl=yv(ii1:ii2);
        
        xvc{end+1}=xvl; yvc{end+1}=yvl;
    end
    
    pgon = polyshape(xvc,yvc);
else
    pgon=xv;
    [xvc,yvc]=boundary(pgon);
end

tri = triangulation(pgon);


TCL=tri.ConnectivityList;
X = tri.Points(:,1);
Y = tri.Points(:,2);

N=size(TCL,1);

% fprintf('\n \t SIDES: %5.0f TRIANGLES: %5.0f CPU: %2.2e',length(xv),N,t1);

% clf;
% hold on;
% triplot(tri);
% for ii=1:N
%     TTloc=TCL(ii,:);
%     Xloc=X(TTloc'); Xloc=[Xloc;Xloc(1)];
%     Yloc=Y(TTloc'); Yloc=[Yloc;Yloc(1)];
%
% end


% rule in barycentric coordinates
[xyw,xyw_bar]=cub_triangle(ade);


xyw=[];
for ii=1:N
    TTloc=TCL(ii,:);
    Xloc=X(TTloc'); Yloc=Y(TTloc');
    xywloc=rule_conversion(xyw_bar,[Xloc Yloc]);
    xyw=[xyw; xywloc];
end

if nargout > 1
[pts,w,momerr] = comprexcub(ade,xyw(:,1:2),xyw(:,3),1); xywc=[pts w];
end


function xyw=rule_conversion(xyw_bar,vertices)

% INPUT:
% ade: ALGEBRAIC DEGREE OF EXACTNESS.
% vertices: 3 x 2 MATRIX OF VERTICES OF THE SIMPLEX.

% OUTPUT:
% xw: NODES AND WEIGHTS OF STROUD CONICAL RULE TYPE OF ADE ade ON THE SIMPLEX
%     WITH VERTICES vertices.

bar_coord=xyw_bar(:,1:3);
xx=bar_coord*vertices;

A=polyarea(vertices(:,1),vertices(:,2));
ww=A*xyw_bar(:,4);

xyw=[xx ww];




