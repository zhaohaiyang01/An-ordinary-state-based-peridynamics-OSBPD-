function [coord,totint] = Gen_nodes(Geome)
% Generate interior material points
% coord - the coordinate (only contain the interior material points)
% totint - the total account of interior material points
dx = Geome.dx;
dof = Geome.dof;
if dof == 2
    totint = Geome.ndivx*Geome.ndivy;
elseif dof == 3
    totint = Geome.ndivx*Geome.ndivy*Geome.ndivz;
else
    error('Please enter the correct number of dof')
end
coord = zeros(totint,dof);
int = 0;
if dof == 2
    for i = 1:Geome.ndivy
        for j = 1:Geome.ndivx
            int = int+1;
            coord(int,1) = -0.5*Geome.length+0.5*dx+(j-1)*dx;
            coord(int,2) = -0.5*Geome.width+0.5*dx+(i-1)*dx;
        end
    end
elseif dof == 3
    for i = 1:Geome.ndivz
        for j = 1:Geome.ndivy
            for k =1:Geome.ndivx
                int = int+1;
                coord(int,1) = -0.5*Geome.length+0.5*dx+(k-1)*dx;
                coord(int,2) = -0.5*Geome.width+0.5*dx+(j-1)*dx;
                coord(int,3) = -0.5*Geome.thick+0.5*dx+(i-1)*dx;
            end
        end
    end
end
    % for i = 1:Geome.ndivy
    %     for j = 1:Geome.ndivx
    %         if dof == 2
    %             int = int+1;
    %             coord(int,1) = 0.5*dx+(j-1)*dx;
    %             coord(int,2) = -0.5*Geome.width+0.5*dx+(i-1)*dx;
    %         elseif dof == 3
    %             for k = 1:Geome.ndivz
    %                 int = int+1;
    %                 coord(int,1) = 0.5*dx+(i-1)*dx;
    %                 coord(int,2) = -0.5*Geome.width+0.5*dx+(j-1)*dx;
    %                 coord(int,3) = -0.5*Geome.thick+0.5*dx+(k-1)*dx;
    %             end
    %         end
    %     end
    % end
    
