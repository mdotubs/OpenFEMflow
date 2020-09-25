function [Cl, Cd, Cm, dCl_dAlpha, dCd_dAlpha, dCl_dMach, dCd_dMach, dCm_dMach, dCl_dU, dCd_dU, dCm_dU, dCl_dXnodes, dCd_dXnodes, dCm_dXnodes, ...
    dCl_dYnodes, dCd_dYnodes, dCm_dYnodes] = Forces(U,Alpha, Mach)

global nodes P Surf XM YM
global Adjoint


Xu = nodes(Surf(Surf(:,2)==1),1);
Yu = nodes(Surf(Surf(:,2)==1),2);

Xl = nodes(Surf(Surf(:,2)==2),1);
Yl = nodes(Surf(Surf(:,2)==2),2);

Pu = U((Surf(Surf(:,2)==1)-1)*4+4);
Pl = U((Surf(Surf(:,2)==2)-1)*4+4);

Cpu = 2*(Pu-P)/Mach^2;
Cpl = 2*(Pl-P)/Mach^2;

Cn = sum((Xl(2:end)-Xl(1:end-1)).*0.5.*(Cpl(2:end)+Cpl(1:end-1))) - sum((Xu(2:end)-Xu(1:end-1)).*0.5.*(Cpu(2:end)+Cpu(1:end-1)));
Ca = sum((Yu(2:end)-Yu(1:end-1)).*0.5.*(Cpu(2:end)+Cpu(1:end-1))) - sum((Yl(2:end)-Yl(1:end-1)).*0.5.*(Cpl(2:end)+Cpl(1:end-1)));

Cmx = sum((Xu(2:end)-Xu(1:end-1)).*0.5.*(Cpu(2:end)+Cpu(1:end-1)).*(0.5.*(Xu(2:end)+Xu(1:end-1))-XM)) - sum((Xl(2:end)-Xl(1:end-1)).*0.5.*(Cpl(2:end)+Cpl(1:end-1)).*(0.5.*(Xl(2:end)+Xl(1:end-1))-XM)); 
Cmy = sum((Yu(2:end)-Yu(1:end-1)).*0.5.*(Cpu(2:end)+Cpu(1:end-1)).*(0.5.*(Yu(2:end)+Yu(1:end-1))-YM)) - sum((Yl(2:end)-Yl(1:end-1)).*0.5.*(Cpl(2:end)+Cpl(1:end-1)).*(0.5.*(Yl(2:end)+Yl(1:end-1))-YM));

Cl = Cn*cosd(Alpha) - Ca*sind(Alpha);
Cd = Cn*sind(Alpha) + Ca*cosd(Alpha);
Cm = Cmx+Cmy;

if Adjoint ==1
    
   dCl_dAlpha = -Cn*sind(Alpha) - Ca*cosd(Alpha);
   dCd_dAlpha =  Cn*cosd(Alpha) - Ca*sind(Alpha);
    
   dCl_dCn =  cosd(Alpha);  dCl_dCa =  -sind(Alpha);
   dCd_dCn =  sind(Alpha);  dCd_dCa =  cosd(Alpha);
   
   %  ************** d_dU *********************
   dCn_dCpu = sparse(size(Cpu));
   dCn_dCpl = sparse(size(Cpl));
   dCmx_dCpu = sparse(size(Cpu));
   dCmx_dCpl = sparse(size(Cpl));
   dCmy_dCpu = sparse(size(Cpu));
   dCmy_dCpl = sparse(size(Cpl));
 
      
   dCn_dCpu(1) = -0.5*(Xu(2)-Xu(1));  
   for i = 2:length(Cpu)-1
        dCn_dCpu(i) = -0.5*(Xu(i)-Xu(i-1)) - 0.5*(Xu(i+1)-Xu(i));
   end
   dCn_dCpu(i+1) = -0.5*(Xu(i+1)-Xu(i));

   
   dCn_dCpl(1) = 0.5*(Xl(2)-Xl(1));
   for i = 2:length(Cpl)-1
        dCn_dCpl(i) = 0.5*(Xl(i)-Xl(i-1)) + 0.5*(Xl(i+1)-Xl(i));
   end
   dCn_dCpl(i+1) = 0.5*(Xl(i+1)-Xl(i));
   
      
   dCa_dCpu = sparse(size(Cpu));
   dCa_dCpl = sparse(size(Cpl));
 
   dCa_dCpu(1) = 0.5*(Yu(2)-Yu(1));
   for i = 2:length(Cpu)-1
        dCa_dCpu(i) = 0.5*(Yu(i)-Yu(i-1)) + 0.5*(Yu(i+1)-Yu(i));
   end
   dCa_dCpu(i+1) = 0.5*(Yu(i+1)-Yu(i));
    
   
   dCa_dCpl(1) = -0.5*(Yl(2)-Yl(1));
   for i = 2:length(Cpl)-1
        dCa_dCpl(i) = -0.5*(Yl(i)-Yl(i-1)) - 0.5*(Yl(i+1)-Yl(i));
   end
   dCa_dCpl(i+1) = -0.5*(Yl(i+1)-Yl(i));
 
Cmx = sum((Xu(2:end)-Xu(1:end-1)).*0.5.*(Cpu(2:end)+Cpu(1:end-1)).*(0.5.*(Xu(2:end)+Xu(1:end-1))-XM)) - sum((Xl(2:end)-Xl(1:end-1)).*0.5.*(Cpl(2:end)+Cpl(1:end-1)).*(0.5.*(Xl(2:end)+Xl(1:end-1))-XM)); 
% Cmy = sum((Yu(2:end)-Yu(1:end-1)).*0.5.*(Cpu(2:end)+Cpu(1:end-1)).*(0.5.*(Yu(2:end)+Yu(1:end-1))-YM)) - sum((Yl(2:end)-Yl(1:end-1)).*0.5.*(Cpl(2:end)+Cpl(1:end-1)).*(0.5.*(Yl(2:end)+Yl(1:end-1))-YM));


   dCmx_dCpu(1) = 0.5*(Xu(2)-Xu(1))*(0.5.*(Xu(2)+Xu(1))-XM); 
   for i = 2:length(Cpu)-1
        dCmx_dCpu(i) = 0.5*(Xu(i)-Xu(i-1))*(0.5.*(Xu(i)+Xu(i-1))-XM)+0.5*(Xu(i+1)-Xu(i))*(0.5.*(Xu(i+1)+Xu(i))-XM);
   end
   dCmx_dCpu(i+1) = 0.5*(Xu(i+1)-Xu(i))*(0.5.*(Xu(i+1)+Xu(i))-XM); 
   
    dCmx_dCpl(1) = -0.5*(Xl(2)-Xl(1))*(0.5.*(Xl(2)+Xl(1))-XM); 
   for i = 2:length(Cpl)-1
        dCmx_dCpl(i) = -0.5*(Xl(i)-Xl(i-1))*(0.5.*(Xl(i)+Xl(i-1))-XM) - 0.5*(Xl(i+1)-Xl(i))*(0.5.*(Xl(i+1)+Xl(i))-XM);
   end
   dCmx_dCpl(i+1) = -0.5*(Xl(i+1)-Xl(i))*(0.5.*(Xl(i+1)+Xl(i))-XM); 

   dCmy_dCpu(1) = 0.5*(Yu(2)-Yu(1))*(0.5.*(Yu(2)+Yu(1))-YM); 
   for i = 2:length(Cpu)-1
        dCmy_dCpu(i) = 0.5*(Yu(i)-Yu(i-1))*(0.5.*(Yu(i)+Yu(i-1))-YM)+0.5*(Yu(i+1)-Yu(i))*(0.5.*(Yu(i+1)+Yu(i))-YM);
   end
   dCmy_dCpu(i+1) = 0.5*(Yu(i+1)-Yu(i))*(0.5.*(Yu(i+1)+Yu(i))-YM); 
   
   dCmy_dCpl(1) = -0.5*(Yl(2)-Yl(1))*(0.5.*(Yl(2)+Yl(1))-YM); 
   for i = 2:length(Cpl)-1
        dCmy_dCpl(i) = -0.5*(Yl(i)-Yl(i-1))*(0.5.*(Yl(i)+Yl(i-1))-YM) - 0.5*(Yl(i+1)-Yl(i))*(0.5.*(Yl(i+1)+Yl(i))-YM);
   end
   dCmy_dCpl(i+1) = -0.5*(Yl(i+1)-Yl(i))*(0.5.*(Yl(i+1)+Yl(i))-YM); 
 

   dCpu_dU = sparse(length(Cpu),length(U));
   dCpl_dU = sparse(length(Cpl),length(U));
   du = (Surf(Surf(:,2)==1)-1)*4+4;
   dl = (Surf(Surf(:,2)==2)-1)*4+4;
   
   for i=1:length(Cpu)
       dCpu_dU(i,du(i)) = 2/Mach^2;
   end
   for i=1:length(Cpl)
       dCpl_dU(i,dl(i)) = 2/Mach^2;
   end
   
   dCn_dU = dCn_dCpu*dCpu_dU + dCn_dCpl*dCpl_dU;
   dCa_dU = dCa_dCpu*dCpu_dU + dCa_dCpl*dCpl_dU;

   dCl_dU = dCl_dCn*dCn_dU + dCl_dCa*dCa_dU;
   dCd_dU = dCd_dCn*dCn_dU + dCd_dCa*dCa_dU;
   dCm_dU = dCmx_dCpu*dCpu_dU + dCmx_dCpl*dCpl_dU + dCmy_dCpu*dCpu_dU + dCmy_dCpl*dCpl_dU;

   % *************** d_dMach ********************

   dCpu_dMach = -4*(Pu-P)/Mach^3;
   dCpl_dMach = -4*(Pl-P)/Mach^3;

   dCn_dMach = dCn_dCpu*dCpu_dMach + dCn_dCpl*dCpl_dMach;
   dCa_dMach = dCa_dCpu*dCpu_dMach + dCa_dCpl*dCpl_dMach;

   dCl_dMach = dCl_dCn*dCn_dMach + dCl_dCa*dCa_dMach;
   dCd_dMach = dCd_dCn*dCn_dMach + dCd_dCa*dCa_dMach;
   dCm_dMach = dCmx_dCpu*dCpu_dMach + dCmx_dCpl*dCpl_dMach + dCmy_dCpu*dCpu_dMach + dCmy_dCpl*dCpl_dMach; 
  
   %  ************** d_dX *********************
   
   dCn_dXu = sparse(size(Xu));
   dCn_dXl = sparse(size(Xl));

   dCmx_dXu = sparse(size(Xu));
   dCmx_dXl = sparse(size(Xl));

   dCn_dXu(1) = 0.5*(Cpu(2)+Cpu(1));
   for i = 2:length(Xu)-1
        dCn_dXu(i) = 0.5*(Cpu(i+1)+Cpu(i)) - 0.5*(Cpu(i)+Cpu(i-1));
   end
   dCn_dXu(i+1) = -0.5*(Cpu(end)+Cpu(end-1));
 
 
   dCmx_dXu(1) = -0.5*Xu(1)*(Cpu(2)+Cpu(1))+0.25*(Cpu(2)+Cpu(1))*XM;
   for i = 2:length(Xu)-1
        dCmx_dXu(i) = 0.5*Xu(i)*(Cpu(i)+Cpu(i-1)) - 0.25*(Cpu(i)+Cpu(i-1))*XM - 0.5*Xu(i)*(Cpu(i+1)+Cpu(i))+0.25*(Cpu(i+1)+Cpu(i))*XM;
   end
   dCmx_dXu(i+1) = 0.5*Xu(i+1)*(Cpu(end)+Cpu(end-1)) - 0.25*(Cpu(end)+Cpu(end-1))*XM;
   
   
   dCn_dXl(1) = -0.5*(Cpl(2)+Cpl(1));
   for i = 2:length(Xl)-1
        dCn_dXl(i) = -0.5*(Cpl(i+1)+Cpl(i)) + 0.5*(Cpl(i)+Cpl(i-1));
   end
   dCn_dXl(i+1) = 0.5*(Cpl(end)+Cpl(end-1));

   
   dCmx_dXl(1) = 0.5*Xl(1)*(Cpl(2)+Cpl(1)) - 0.25*(Cpl(2)+Cpl(1))*XM;
   for i = 2:length(Xl)-1
        dCmx_dXl(i) = -0.5*Xl(i)*(Cpl(i)+Cpl(i-1)) + 0.25*(Cpl(i)+Cpl(i-1))*XM + 0.5*Xl(i)*(Cpl(i+1)+Cpl(i)) - 0.25*(Cpl(i+1)+Cpl(i))*XM;
   end
   dCmx_dXl(i+1) = -0.5*Xl(i+1)*(Cpl(end)+Cpl(end-1)) + 0.25*(Cpl(end)+Cpl(end-1))*XM;

   
   dCa_dYu = sparse(size(Yu));
   dCa_dYl = sparse(size(Yl));
   
   dCmy_dYu = sparse(size(Yu));
   dCmy_dYl = sparse(size(Yl));
   
   dCa_dYu(1) = -0.5*(Cpu(2)+Cpu(1));
   for i = 2:length(Yu)-1
        dCa_dYu(i) = 0.5*(Cpu(i)+Cpu(i-1)) - 0.5*(Cpu(i+1)+Cpu(i));
   end
   dCa_dYu(i+1) = 0.5*(Cpu(end)+Cpu(end-1));
    
   dCa_dYl(1) = 0.5*(Cpl(2)+Cpl(1));
   for i = 2:length(Yl)-1
        dCa_dYl(i) = -0.5*(Cpl(i)+Cpl(i-1)) + 0.5*(Cpl(i+1)+Cpl(i));
   end
   dCa_dYl(i+1) = -0.5*(Cpl(end)+Cpl(end-1));

   dCmy_dYu(1) = -0.5*Yu(1)*(Cpu(2)+Cpu(1))+0.25*(Cpu(2)+Cpu(1))*YM;
   for i = 2:length(Yu)-1
        dCmy_dYu(i) = 0.5*Yu(i)*(Cpu(i)+Cpu(i-1)) - 0.25*(Cpu(i)+Cpu(i-1))*YM - 0.5*Yu(i)*(Cpu(i+1)+Cpu(i))+0.25*(Cpu(i+1)+Cpu(i))*YM;
   end
   dCmy_dYu(i+1) = 0.5*Yu(i+1)*(Cpu(end)+Cpu(end-1)) - 0.25*(Cpu(end)+Cpu(end-1))*YM;

   dCmy_dYl(1) = 0.5*Yl(1)*(Cpl(2)+Cpl(1)) - 0.25*(Cpl(2)+Cpl(1))*YM;
   for i = 2:length(Yl)-1
        dCmy_dYl(i) = -0.5*Yl(i)*(Cpl(i)+Cpl(i-1)) + 0.25*(Cpl(i)+Cpl(i-1))*YM + 0.5*Yl(i)*(Cpl(i+1)+Cpl(i)) - 0.25*(Cpl(i+1)+Cpl(i))*YM;
   end
   dCmy_dYl(i+1) = -0.5*Yl(i+1)*(Cpl(end)+Cpl(end-1)) + 0.25*(Cpl(end)+Cpl(end-1))*YM;

   
   dXu_dXnodes = sparse(length(Xu),size(nodes,1));
   dXl_dXnodes = sparse(length(Xl),size(nodes,1));
   
   dYu_dYnodes = sparse(length(Yu),size(nodes,1));
   dYl_dYnodes = sparse(length(Yl),size(nodes,1));
   
   du = Surf(Surf(:,2)==1);
   dl = Surf(Surf(:,2)==2);
   
   for i=1:length(Xu)
       dXu_dXnodes(i,du(i)) = 1;
   end
   for i=1:length(Xl)
       dXl_dXnodes(i,dl(i)) = 1;
   end
   
   for i=1:length(Yu)
       dYu_dYnodes(i,du(i)) = 1;
   end
   for i=1:length(Yl)
       dYl_dYnodes(i,dl(i)) = 1;
   end
   
   dCn_dXnodes = dCn_dXu*dXu_dXnodes + dCn_dXl*dXl_dXnodes;
   dCa_dYnodes = dCa_dYu*dYu_dYnodes + dCa_dYl*dYl_dYnodes;
   
   dCl_dXnodes = dCl_dCn*dCn_dXnodes;
   dCl_dYnodes = dCl_dCa*dCa_dYnodes;
   
   dCd_dXnodes = dCd_dCn*dCn_dXnodes;
   dCd_dYnodes = dCd_dCa*dCa_dYnodes;
   
   dCm_dXnodes = dCmx_dXu*dXu_dXnodes + dCmx_dXl*dXl_dXnodes;
   dCm_dYnodes = dCmy_dYu*dYu_dYnodes + dCmy_dYl*dYl_dYnodes;

   
   dCl_dU = dCl_dU';
   dCd_dU = dCd_dU';
   dCm_dU = dCm_dU';
end