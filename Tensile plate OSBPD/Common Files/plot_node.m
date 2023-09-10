function plot_node(coord,totint,Geome,fnodes)
scrsz = get(groot, 'ScreenSize');
figure('Position',[scrsz(4)/15 scrsz(4)/12 9*scrsz(3)/10 4.*scrsz(4)/5]);
subplot(2,2,1)
if Geome.dof == 2
    scatter(coord(1:totint,1),coord(1:totint,2),1,[1,0,0])
    if size(coord,1) > totint
        hold on
        scatter(coord(totint+1:end,1),coord(totint+1:end,2),1,[0,0,1])
    end
    if fnodes ~= 0 
        hold on
        scatter(coord(fnodes,1),coord(fnodes,2),1,[0,1,0])
    end
    axis equal ; axis on;
    xlim(1.1*[min(coord(:,1)) max(coord(:,1))]);ylim(1.1*[min(coord(:,2)) max(coord(:,2))])
elseif Geome.dof == 3
    scatter3(coord(1:totint,1),coord(1:totint,2),coord(1:totint,3),1,[1 0 0])
    if size(coord,1) > totint
        hold on
        scatter3(coord(totint+1:end,1),coord(totint+1:end,2),coord(totint+1:end,3),1,[0 0 1])
    end
    if fnodes ~= 0 
        hold on
        scatter3(coord(fnodes,1),coord(fnodes,2),coord(fnodes,3),1,[0,1,0])
    end
    axis equal ; axis on;
    xlim(1.1*[min(coord(:,1)) max(coord(:,1))]);ylim(1.1*[min(coord(:,2)) max(coord(:,2))])
    zlim(1.1*[min(coord(:,3)) max(coord(:,3))])
end
drawnow