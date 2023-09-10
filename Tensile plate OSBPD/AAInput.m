function [Geome,Mater] = AAInput()
%% Geometry information
Geome.length = 0.05;
Geome.width = 0.05;
Geome.thick = 0.0005;
Geome.ndivx = 150;
Geome.ndivy = 150;
Geome.ndivz = 1;
Geome.dof = 2;
Geome.nbnd = 3; 
Geome.nt = 1000;                    
% nt: Total number of time step
Geome.dt = 1;
Geome.ADR = 1;
%Geome.uinc = [0,2.7541e-7,0];                  
% Displacement BC,uinc:velcoty
Geome.P = 10e6;                       
% P:Force
Geome.prec = 0;                        
% Prec:Whether there is a pre-crack 1 - exist,0 - non-exist
% *** The precrack process applies only to 2D problems ***
Geome.dx = Geome.length/Geome.ndivx;
Geome.radij = Geome.dx/2;
Geome.area = Geome.dx^2;
%Geome.delta = 3.015 * Geome.dx;
Geome.delta = 0.0015075;
%% Material information
Mater.alpha = 23.0e-6;                 
% alpha: Coefficient of thermal expansion
Mater.dens = 8000.0;                   
% dens: Density
Mater.emod = 10.0e9;                  
% emod: Elastic modulus
Mater.pratio = 1.0/8.0;              
% pratio: Poisson's ratio

if Geome.dof ==  2
    Geome.vol = Geome.area*Geome.thick;
    vmod = Mater.emod/(2*(1-Mater.pratio));% Volume modulus
    smod = Mater.emod/(2*(1+Mater.pratio));% Shear modulus
    Mater.a1 = 1/2*(vmod-2*smod);
    Mater.a2 = 4*Mater.alpha*Mater.a1;
    Mater.a3 = Mater.alpha*Mater.a2;
    Mater.b = 6*smod/pi/Geome.thick/Geome.delta^4;
    Mater.d = 2/pi/Geome.thick/Geome.delta^3;
% elseif Geome.dof == 3
%     Geome.vol = Geome.area*Geome.dx;
%     Mater.bc = 12.0 * Mater.emod / (pi * (Geome.delta^4));
end

Mater.dtemp = 0.0;                     
% dtemp: Temperature change
Mater.sc = 1;
