function CP = PostProc(U,Mach,Options)

global nodes P Surf


%% u,v,p

u = U(2:4:end);
v = U(3:4:end);
p = U(4:4:end);

x = nodes(:,1);
y = nodes(:,2);


% xq = meshgrid(linspace(min(x),max(x),400));
% yq = meshgrid(linspace(min(y),max(y),400))';

xq = meshgrid(linspace(-2,3,500));
yq = meshgrid(linspace(-2,3,500))';

uq = griddata(x,y,full(u),xq,yq);
vq = griddata(x,y,full(v),xq,yq);
pq = griddata(x,y,full(p),xq,yq);

velq = sqrt(uq.^2+vq.^2);


%% wall data

Xu = nodes(Surf(Surf(:,2)==1),1);
Yu = nodes(Surf(Surf(:,2)==1),2);

Xl = nodes(Surf(Surf(:,2)==2),1);
Yl = nodes(Surf(Surf(:,2)==2),2);

Pu = U((Surf(Surf(:,2)==1)-1)*4+4);
Pl = U((Surf(Surf(:,2)==2)-1)*4+4);

Cpu = 2*(Pu-P)/Mach^2;
Cpl = 2*(Pl-P)/Mach^2;

CP.Xu = Xu;
CP.Xl = Xl;
CP.Cpu = Cpu;
CP.Cpl = Cpl;

%% plotting

if Options.Plot ==1
figure
subplot(2,2,1);
hold on
contour(xq,yq,velq,20);
% area(xg,yg);
area([Xu; Xl],[Yu;Yl]);
axis equal
colorbar; 
xlabel('x');
ylabel('y');
title('Countours of Mach');
hold off


subplot(2,2,2);
hold on
contour(xq,yq,pq,20);
area([Xu; Xl],[Yu;Yl]);
% area([Xu; Xl]-0.5,[Yu;Yl]);
axis equal
colorbar; 
xlabel('x');
ylabel('y');
title('Countours of Pressure');
hold off

subplot(2,2,3);
hold on
quiver(x,y,u,v)
axis equal
colorbar; 
xlabel('x');
ylabel('y');
title('Velocity vectors');
hold off

subplot(2,2,4);
hold on
plot(Xu,-Cpu,'-r');
plot(Xl,-Cpl,'-r');
xlabel('x');
ylabel('-C_p');
title('Wall pressure coefficient');

end
% save('CP.mat','CP');






