function dila = dilation(numfam,node,nodefam,pointfam,coord,acoord,Geome,Mater)
delta = Geome.delta;
radij = Geome.radij;
alpha = Mater.alpha;
dtemp = Mater.dtemp;
dof = Geome.dof;
dila = 0;
d = Mater.d;
vol = Geome.vol;
for j = 1:numfam(node,1)
    cnodej = nodefam(pointfam(node, 1)+j, 1);
    if dof == 2
        idist = sqrt((coord(cnodej,1)-coord(node,1))^2+(coord(cnodej,2)-coord(node, 2))^2);
        nlength = sqrt((acoord(cnodej,1)-acoord(node,1))^2+(acoord(cnodej,2)-acoord(node,2))^2);
        dot_xy = (coord(cnodej,1)-coord(node,1))*(acoord(cnodej,1)-acoord(node,1))+...
            (coord(cnodej,2)-coord(node,2))*(acoord(cnodej,2)-acoord(node,2));
    elseif dof == 3
        idist = sqrt((coord(cnodej,1)-coord(node,1))^2+(coord(cnodej,2)-coord(node, 2))^2+(coord(cnodej,3)-coord(node, 3))^2);
        nlength = sqrt((acoord(cnodej,1)-acoord(node,1))^2+(acoord(cnodej,2)-acoord(node,2))^2+(acoord(cnodej,3)-acoord(node,3))^2);
        dot_xy = (coord(cnodej,1)-coord(node,1))*(acoord(cnodej,1)-acoord(node,1))+...
            (coord(cnodej,2)-coord(node,2))*(acoord(cnodej,2)-acoord(node,2))+...
            (coord(cnodej,3)-coord(node,3))*(acoord(cnodej,3)-acoord(node,3));
    end
    sij = (nlength-idist)/idist;
    if (idist <= delta - radij)
        fac = 1.0;
    elseif (idist <= delta + radij)
        fac = (delta + radij - idist) / (2.0 * radij);
    else
        fac = 0.0;
    end
    weight = delta/idist;
    dila = dila+d*weight*(sij-alpha*dtemp)*dot_xy/nlength*vol*fac+3*alpha*dtemp;
end

% Plot nodes in support domain
% figure(2)
% plot(coord(:,1),coord(:,2),'r.')
% hold on
% plot(coord(node,1),coord(node,2),'ko',coord(nodefam(pointfam(node, 1)+1:pointfam(node, 1)+numfam(node,1), 1),1),...
%    coord(nodefam(pointfam(node, 1)+1:pointfam(node, 1)+numfam(node,1), 1),2),'k*')
% para = [coord(node,1)-delta,coord(node,2)-delta, 2*delta, 2*delta];
% rectangle('Position', para, 'Curvature', [1 1])%,'FaceColor','b');
% axis equal
