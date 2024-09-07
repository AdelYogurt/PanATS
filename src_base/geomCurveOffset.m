function node_list=geomCurveOffset(node_list,offset)
% shifting line with offset
% direction is outside
%

% pgon=polyshape(node_list(:,1),node_list(:,2));
% plot(pgon);hold on;axis equal;

% calculate normal vector of edge
E_nmvctr_list=[diff(node_list,1,1);node_list(1,:)-node_list(end,:)];
E_nmvctr_list=E_nmvctr_list*[0,-1;1,0];
E_nmvctr_list=E_nmvctr_list./sqrt(sum(E_nmvctr_list.^2,2));

% calculate sin of half of angle between external normal vectors
cos_ang_list=sum(E_nmvctr_list.*circshift(E_nmvctr_list,1,1),2);
coshang_list=sqrt((1+cos_ang_list)/2);

% calculate vertex offset vector
P_nmvctr_list=(E_nmvctr_list+circshift(E_nmvctr_list,1,1))/2;
P_nmvctr_list=P_nmvctr_list./sqrt(sum(P_nmvctr_list.^2,2));

% offset vertex
P_offset_list=offset./coshang_list;
node_list=node_list+P_offset_list.*P_nmvctr_list;

% quiver(node_list(:,1),node_list(:,2),E_nmvctr_list(:,1),E_nmvctr_list(:,2))
% quiver(node_list(:,1),node_list(:,2),P_nmvctr_list(:,1),P_nmvctr_list(:,2))
% pgon=polyshape(node_list(:,1),node_list(:,2));

end
