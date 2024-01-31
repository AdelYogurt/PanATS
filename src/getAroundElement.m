function [Elem_arou,Vertex_arou]=getAroundElement(element_list,elem_idx,node_idx)
% base on element_list and node_idx
% find out vertex and element around node_idx
% base on half edge to search all around element
%
global user_model

vertex_list=user_model.vertex_list;
pnt_base_idx=element_list(elem_idx).Node_idx(node_idx);
Vertex_arou=vertex_list{pnt_base_idx,1};
Elem_arou=vertex_list{pnt_base_idx,2};
Vertex_arou=Vertex_arou(Vertex_arou ~= 0);
Elem_arou=Elem_arou(Elem_arou ~= 0);

% elem_idx_origin=elem_idx;
% elem=element_list(elem_idx);
% vertex_idx_orign=elem.Node_idx(node_idx);
% if node_idx == elem.node_num
%     vertex_idx=elem.Node_idx(1);
% else
%     vertex_idx=elem.Node_idx(node_idx+1);
% end
% Elem_arou=elem_idx;
% Vertex_arou=vertex_idx;
% 
% % anti-clock order search loop
% elem_idx=elem.Vertex_next(node_idx);
% while elem_idx_origin ~= elem_idx && elem_idx ~= 0 && length(Elem_arou) < 16
%     elem=element_list(elem_idx);
%     Node_idx=elem.Node_idx;
%     node_idx=find(Node_idx == vertex_idx_orign);
%     if node_idx == elem.node_num
%         vertex_idx=Node_idx(1);
%     else
%         vertex_idx=Node_idx(node_idx+1);
%     end
%     Elem_arou=[Elem_arou,elem_idx];
%     Vertex_arou=[Vertex_arou,vertex_idx];
% 
%     % search nex element
%     elem_idx=elem.Vertex_next(node_idx);
% 
%     if elem_idx == 0
%         break;
%     end
% end
% 
% % clock order search loop
% if elem_idx == 0
%     % begin clock order search loop
%     elem=element_list(elem_idx_origin);
%     Node_idx=elem.Node_idx;
%     node_idx=find(Node_idx == vertex_idx_orign);
%     if node_idx == 1
%         elem_idx=elem.Vertex_next(elem.node_num);
%     else
%         elem_idx=elem.Vertex_next(node_idx-1);
%     end
% end
% while elem_idx_origin ~= elem_idx && elem_idx ~= 0 && length(Elem_arou) < 16
%     elem=element_list(elem_idx);
%     Node_idx=elem.Node_idx;
%     node_idx=find(Node_idx == vertex_idx_orign);
%     if node_idx == elem.node_num
%         vertex_idx=Node_idx(1);
%     else
%         vertex_idx=Node_idx(node_idx+1);
%     end
%     Elem_arou=[elem_idx,Elem_arou];
%     Vertex_arou=[vertex_idx,Vertex_arou];
% 
%     % search nex element
%     if node_idx == 1
%         elem_idx=elem.Vertex_next(elem.node_num);
%     else
%         elem_idx=elem.Vertex_next(node_idx-1);
%     end
% 
%     if elem_idx == 0
%         break;
%     end
% end

end
