% This is an ordinary state-based peridynamics(OSBPD) programming
% By Jin-Hu Pan at WHU
clear
clc
close all
format long
addpath('Common Files')
% Input information
[Geome,Mater] = AAInput();
dt = Geome.dt;
dx = Geome.dx;
dof = Geome.dof;
vol = Geome.vol;
delta = Geome.delta;
thick = Geome.thick;
radij = Geome.radij;
nt = Geome.nt;
a1 = Mater.a1;
a2 = Mater.a2;
a3 = Mater.a3;
b = Mater.b;
d = Mater.d;
alpha = Mater.alpha;
dtemp = Mater.dtemp;
pratio = Mater.pratio;
sc = Mater.sc;

% Generate material points
[coord,totint] = Gen_nodes(Geome);

% Find displacement and force boundary
[coord,totnode,fnodes,inr,dire,inb] = AAFindBoundary(Geome,coord,totint);

% Visualize nodes and boundaries
plot_node(coord,totint,Geome,fnodes)

% Find neighboor nodes
tic
[nodefam,numfam,pointfam,fail] = neighboornodes(coord,Geome);
toc
% Surface correction
%fncst = Surc(coord,totnode,Geome,Mater,numfam,nodefam,pointfam);
svcr = Surc_vm(Geome,coord,totnode,numfam,pointfam,nodefam);
% Stable mass vector computation
if Geome.ADR == 1
    if dof == 2
        mass = 0.25 * dt * dt * (pi * (delta)^2 * thick) * 4*delta*b / dx * 1.1;
        massvec = ones(totnode,dof).*mass;
%     elseif dof == 3
%         mass = 0.25 * dt * dt * ((4.0 / 3.0) * pi * (delta)^3) * bc / dx;
%         massvec = ones(totnode,dof).*mass;
    end
end
D = reshape((massvec(1:totint,:))',[],1);

% Applied force
bforce = zeros(totnode,dof);
if isfield(Geome,'P')
    bforce(fnodes(1:inr),dire) = Geome.P/dx;
    bforce(fnodes(inr+1:end),dire) = -Geome.P/dx;
end

% Initialization of displacements and velocities
velhalf = zeros(totnode,dof);
velhalfold = zeros(totnode,dof);
vel = zeros(totnode,dof);
% vel: velocity of a material point
disp = zeros(totnode,dof);
pforce = zeros(totnode,dof);
% pforce: total peridynamic force acting on a material point
pforceold = zeros(totnode,dof);
% pforceold: total peridynamic force acting on a material point in the previous time step 1:x-coord,2:y-coord
acc = zeros(totnode,dof);
% acc: acceleration of a material point
dmg = zeros(totint, 1);

coord_disp_pd_nt_675 = zeros(totint, 5);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
coord_disp_pd_nt_750 = zeros(totint, 5);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
coord_disp_pd_nt_825 = zeros(totint, 5);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
coord_disp_pd_nt_1000 = zeros(totint, 5);
% Peridynamic displacement and Analytical displacement of all points at time step of nt

cn = 0.0;
cn1 = 0.0;
cn2 = 0.0;  % These should be put in loop.
V = pi*delta^2*thick;
%V = max(svcr);
%% Time integration
for tt = 1:nt
    fprintf("%d/%d\n",tt,nt);
    tic
    % Application of boundary conditions
    if isfield(Geome,'uinc') && (sum(Geome.uinc) ~= 0)
        disp(totint+1:totint+bot,:) = -ones(bot,1)*Geome.uinc(1:dof).*tt*dt;
        disp(totint+bot+1:totnode,:) = ones(totnode-bot-totint,1)*Geome.uinc(1:dof).*tt*dt;
    end
    acoord = coord+disp;
    for i = 1:totint
        dmgpar1 = 0.0;
        dmgpar2 = 0.0;
        pforce(i, :) = 0.0;
        dilai = dilation(numfam,i,nodefam,pointfam,coord,acoord,Geome,Mater);
        for j = 1:numfam(i,1)
            cnode = nodefam(pointfam(i, 1)+j, 1);
            dilaj = dilation(numfam,cnode,nodefam,pointfam,coord,acoord,Geome,Mater);
            if dof == 2
                idist = sqrt((coord(cnode,1)-coord(i,1))^2+(coord(cnode,2)-coord(i, 2))^2);
                nlength = sqrt((acoord(cnode,1)-acoord(i,1))^2+(acoord(cnode,2)-acoord(i,2))^2);
                dot_xy = (coord(cnode,1)-coord(i,1))*(acoord(cnode,1)-acoord(i,1))+...
                         (coord(cnode,2)-coord(i,2))*(acoord(cnode,2)-acoord(i,2));
            elseif dof == 3
                idist = sqrt((coord(cnode,1)-coord(i,1))^2+(coord(cnode,2)-coord(i, 2))^2+(coord(cnode,3)-coord(i, 3))^2);
                nlength = sqrt((acoord(cnode,1)-acoord(i,1))^2+(acoord(cnode,2)-acoord(i,2))^2+(acoord(cnode,3)-acoord(i,3))^2);
                dot_xy = (coord(cnode,1)-coord(i,1))*(acoord(cnode,1)-acoord(i,1))+...
                         (coord(cnode,2)-coord(i,2))*(acoord(cnode,2)-acoord(i,2))+...
                         (coord(cnode,3)-coord(i,3))*(acoord(cnode,3)-acoord(i,3));
            end
            % Volume correction
            if (idist <= delta - radij)
                fac = 1.0;
            elseif (idist <= delta + radij)
                fac = (delta + radij - idist) / (2.0 * radij);
            else
                fac = 0.0;
            end
            Vi = svcr(i);
            Vj = svcr(cnode);
            scr = 2*V/(Vi+Vj);
            scr = max(1,scr);
            
            fn1 = (coord(cnode,1)+disp(cnode,1)-coord(i,1) - disp(i,1))/nlength;
            fn2 = (coord(cnode,2)+disp(cnode,2)-coord(i,2) - disp(i,2))/nlength;
            weight = delta/idist;
            A = 4*weight*(d*dot_xy/idist/nlength*(a1*dilai-1/2*a2*dtemp)+b*(nlength-idist)-alpha*dtemp*idist);
            B = 4*weight*(d*dot_xy/idist/nlength*(a1*dilaj-1/2*a2*dtemp)+b*(nlength-idist)-alpha*dtemp*idist);
            dforce1A = 1/2*A*fn1;
            dforce2A = 1/2*A*fn2;
            dforce1B = -1/2*B*fn1;
            dforce2B = -1/2*B*fn2;
            dforce1 = dforce1A-dforce1B;
            dforce2 = dforce2A-dforce2B;
            pforce(i,1) = pforce(i,1)+dforce1*vol*fac*scr;
            pforce(i,2) = pforce(i,2)+dforce2*vol*fac*scr;
            if dof == 3
                fn3 = (coord(cnode,3)+disp(cnode,3)-coord(i,3) - disp(i,3))/nlength;
                dforce3 = bc*((nlength-idist)/idist-(alpha*dtemp))*vol*scr*fac*fn3*fail(i,j);
                pforce(i,3) = pforce(i,3)+dforce3;
            end
            % Definition of a no-fail zone 
            if (abs((nlength - idist)/idist) > sc)
                fail(i, j) = 0;
            end
            dmgpar1 = dmgpar1 + fail(i, j)*vol*fac;
            dmgpar2 = dmgpar2 + vol*fac;
            if dforce1 ~= 0
                j;
            end
        end
        % Calculation of the damage parameter 
        dmg(i, 1) = 1.0 - dmgpar1/dmgpar2;
    end
    toc
    
    % Adaptive dynamic relaxation ⬇⬇⬇
    tic
    for i = 1:totint
        if (velhalfold(i, 1) ~= 0.0)
            cn1 = cn1 - disp(i, 1) * disp(i, 1) * (pforce(i, 1) / massvec(i, 1) - pforceold(i, 1) / massvec(i, 1)) / (dt * velhalfold(i, 1));
        end
        if (velhalfold(i, 2) ~= 0.0)
            cn1 = cn1 - disp(i, 2) * disp(i, 2) * (pforce(i, 2) / massvec(i, 2) - pforceold(i, 2) / massvec(i, 2)) / (dt * velhalfold(i, 2));
        end
        cn2 = cn2 + disp(i, 1) * disp(i, 1);
        cn2 = cn2 + disp(i, 2) * disp(i, 2);
        if dof == 3 && (velhalfold(i, 3) ~= 0.0)
            cn1 = cn1 - disp(i, 3) * disp(i, 3) * (pforce(i, 3) / massvec(i, 3) - pforceold(i, 3) / massvec(i, 3)) / (dt * velhalfold(i, 3));
            cn2 = cn2 + disp(i, 3) * disp(i, 3);
        end
    end
    
    if (cn1/cn2 > 0) && (cn2 ~= 0.0)
        cn = 2.0 * sqrt(cn1/cn2);
    else
        cn = 0.0;
    end
    
    if (cn > 2.0) % cn需要<2
        cn = 1.9;
    end
    
    Fb_int = reshape((bforce(1:totint,:)+pforce(1:totint,:))',[],1);
    velhalfold_int = reshape((velhalfold(1:totint,:))',[],1);
    if tt == 1
        velhalf_int = dt./D.*Fb_int/2;
    else
        velhalf_int = ((2-cn*dt).*velhalfold_int+2*dt./D.*Fb_int)/(2+cn*dt);
    end
    velhalf(1:totint,:) = (reshape(velhalf_int,dof,totint))';
    vel = 0.5*(velhalf+velhalfold);
    disp = disp+dt.*velhalf;
    % Adaptive dynamic relaxation ⬆⬆⬆
    % Update
    pforceold = pforce;
    velhalfold = velhalf;
    if mod(tt,50) == 0
        Processing(Geome,Mater,totint,coord,disp,dmg,tt)
        drawnow
    end
end

Processing(Geome,Mater,totint,coord,disp,dmg,tt)
