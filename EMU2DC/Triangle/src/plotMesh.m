function plotMesh

 load nodes.dat
 load elems.dat

 x = nodes(:,1);
 y = nodes(:,2);

 plot(x, y, 'ro'); hold on;
 for ii=1:length(elems)
  node1 = elems(ii,2)+1;
  node2 = elems(ii,3)+1;
  node3 = elems(ii,4)+1;
  plot([x(node1) x(node2) x(node3) x(node1)],... 
       [y(node1) y(node2) y(node3) y(node1)], 'b-'); 
 end
