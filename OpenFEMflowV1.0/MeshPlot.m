function MeshPlot(nodes,element)


for e=1:size(element,1)
    for i=1:3
        xm(i,e) = nodes(element(e,i),1);
        ym(i,e) = nodes(element(e,i),2);
    end
end


figure
patch(xm,ym,ones(size(xm)))
axis equal
xlabel('x');
ylabel('y');
title('Mesh');