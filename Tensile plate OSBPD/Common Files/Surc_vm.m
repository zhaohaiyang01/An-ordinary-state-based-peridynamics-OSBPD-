function svcr = Surc_vm(Geome,coord,totnode,numfam,pointfam,nodefam)
% original point locates in center
delta = Geome.delta;
thick = Geome.thick;
area = Geome.area;
radij = Geome.radij;
dv = area*thick;
dx = Geome.dx;
V = pi*delta^2*thick;
svcr = ones(totnode,1).*V;
for i =1:totnode
    Vs = 0; % The volume of support domain
    for k = 1:numfam(i,1)
        cnode = nodefam(pointfam(i, 1)+k, 1);
        idist = sqrt((coord(cnode,1)-coord(i,1))^2+(coord(cnode,2)-coord(i, 2))^2);
        if (idist <= delta - radij)
            fac = 1.0;
        elseif (idist <= delta + radij)
            fac = (delta + radij - idist) / (2.0 * radij);
        else
            fac = 0.0;
        end
        Vs = Vs+fac*dv;
    end
        svcr(i) = Vs+dv;
end





