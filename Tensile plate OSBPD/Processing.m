function Processing(Geome,Mater,totint,coord,disp,dmg,tt)
% Analytical Solution

% % plane stress
x = (-0.5*Geome.length:0.005:0.5*Geome.length);
E = Mater.emod;
p1 = Geome.P;
nu = Mater.pratio;
u = p1/E*x;
y = (-0.5*Geome.width:0.005:0.5*Geome.width);
v = -p1*nu/E*y;

% plane strain
% x = (-0.025:0.005:0.025);
% E = Mater.emod;
% p1 = Geome.P;
% nu = Mater.pratio;
% E = E/(1-nu^2);
% nu = nu/(1-nu);
% u = p1/E*x;
% y = (-0.025:0.005:0.025);
% v = -p1*nu/E*y;


coord_disp_pd_nt = zeros(totint, 4);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
horiCnt = 0;
horizontal_disps = zeros(Geome.ndivx, 4);
% Peridynamic displacement and Analytical displacement of points at y = 0;
vertiCnt = 0;
vertical_disps = zeros(Geome.ndivy, 4);
% Peridynamic displacement and Analytical displacement of points at x = 0;
dx = Geome.dx;
for i = 1:totint
    coord_disp_pd_nt(i, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
    if (abs(coord(i, 2)-(dx / 2.0)) <= 1e-8)
        horiCnt = horiCnt + 1;
        horizontal_disps(horiCnt, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
    end
    if (abs(coord(i, 1)-(dx / 2.0)) <= 1e-8)
        vertiCnt = vertiCnt + 1;
        vertical_disps(vertiCnt, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
    end
end

colormap jet;
subplot(221)
scale = 1;
scatter(coord_disp_pd_nt(:, 1)+scale*coord_disp_pd_nt(:, 3), coord_disp_pd_nt(:, 2)+scale*coord_disp_pd_nt(:, 4), [], coord_disp_pd_nt(:, 4), "filled")
% subplot(222)
% plot(steady_check(:, 1), steady_check(:, 2), steady_check(:, 1), steady_check(:, 3))
subplot(223)
plot(horizontal_disps(1:horiCnt, 1), horizontal_disps(1:horiCnt, 3),'-ro',x,u,'-b*')
tit = ['step = ',num2str(tt)];
title(tit)
subplot(224)
plot(vertical_disps(1:vertiCnt, 2), vertical_disps(1:vertiCnt, 4),'-go',y,v,'-k*')
if tt == Geome.nt
    ink = 0;
    for i = 1:totint
        if (abs(coord(i, 2)-max(coord(:, 2))) <= 1.0e-8)
            ink = ink + 1;
            kkk(ink, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
        end
    end
    figure(2)
    plot(horizontal_disps(1:horiCnt, 1), horizontal_disps(1:horiCnt, 3),'-ro',x,u,'-b*',kkk(1:ink,1),kkk(1:ink,3),'-go')
    
    
    ink = 0;
    for i = 1:totint
        if (abs(coord(i, 1)-max(coord(:, 1))) <= 1.0e-8)
            ink = ink + 1;
            kkk(ink, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
        end
    end
    figure(3)
    plot(vertical_disps(1:vertiCnt, 2), vertical_disps(1:vertiCnt, 4),'-ro',y,v,'-k*',kkk(1:ink,2),kkk(1:ink,4),'-go')
end