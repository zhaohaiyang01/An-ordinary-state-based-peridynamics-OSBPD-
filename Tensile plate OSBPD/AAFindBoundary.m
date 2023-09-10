function [coord,totnode,fnodes,inr,dire,inb] = AAFindBoundary(Geome,coord,totint)
int = 0;
inb = 0;
dx = Geome.dx;
dof = Geome.dof;
fnodes = 0;
dire = 0;
%% Find displacement boundary 
if isfield(Geome,'uinc')
    if dof == 2
        for i = 1:Geome.nbnd
            for j = 1:Geome.ndivx
                inb = inb+1;
                coord(totint+inb,1) = -0.5*Geome.length+0.5*dx+(j-1)*dx; 
                coord(totint+inb,2) = -0.5*Geome.width-0.5*dx-(i-1)*dx;   % bottom boundary
            end
        end
        
        for i = 1:Geome.nbnd
            for j = 1:Geome.ndivx
                int = int+1;
                coord(totint+int+inb,1) = -0.5*Geome.length+0.5*dx+(j-1)*dx; 
                coord(totint+int+inb,2) = 0.5*Geome.width+0.5*dx+(i-1)*dx;   % top boundary
            end
        end
    elseif dof == 3
        for i = 1:Geome.ndivz
            for j = 1:Geome.ndivy
                for k = 1:Geome.nbnd
                    int = int+1;
                coord(totint+int,1) = -0.5*dx-(k-1)*dx; % left boundary
                coord(totint+int,2) = -0.5*Geome.width+0.5*dx+(j-1)*dx;
                coord(totint+int,3) = -0.5*Geome.thick+0.5*dx+(i-1)*dx;
                end
            end
        end
    end
end
totnode = totint+int+inb;

%% Find force boundary
if isfield(Geome,'P')
   fnodes1 = find(coord(:,1)==max(coord(:,1)));  % The nodes applied force 
   inr = size(fnodes1,1);
   fnodes2 = find(coord(:,1)==min(coord(:,1)));  % The nodes applied force 
   fnodes = [fnodes1,fnodes2];
   dire = 1;
end
        
                
            
