function preModelPanel()
% function to prepare model for panel method(only basical function)
% 
% calculate monitor element geomertry (element normal vector)
% project element center point to flow vector
% sort g_Element by element_center_point_project_flow_vector
% nearby element only search in nearby after sort g_Element
% calculate all elemenet geometry
%
% element tree is an ADT tree
% root is stagnation point element of all marker element(struct)
% nodes have properity parent nodes and children nodes(struct)
% parent nodes and children nodes can been more
% number of parent node and number of children nodes add is number of element edges
% a node can be children of same height of tree
%
% copyright Adel 2022.11
%
global g_model

% calculate element geometry parameter
volume_center_point=zeros(1,3);
center_point_list=zeros(element_number,3);
HATS_element_list=repmat(HATSElement(0,0),element_number,1); % all element list
for element_index=1:element_number 
    element_type=g_Element(element_index,1);
    
    switch element_type
        case 5 % triangle
            point_number=3;
        case 9 % quadrilateral
            point_number=4;
    end
    point_index_list=g_Element(element_index,2:1+point_number);
    center_point=sum(g_Point(point_index_list,1:3),1)/point_number;
    
    center_point_list(element_index,:)=center_point;
    HATS_element=HATSElement(element_type,point_index_list);
    
    % calculate geomertry properity
    dr12=g_Point(point_index_list(2),1:3)-g_Point(point_index_list(1),1:3);
    dr23=g_Point(point_index_list(3),1:3)-g_Point(point_index_list(2),1:3);
    
    % calculate norm_vector of element
    cross_vector=cross(dr12,dr23);
    length_cross_vector=norm(cross_vector,2);
    HATS_element.area=length_cross_vector/2;
    HATS_element.normal_vector=cross_vector/length_cross_vector;
    
    switch element_type
        case 5 % triangle
        case 9 % quadrilateral
            % quadrilateral include another triangle
            dr34=g_Point(point_index_list(4),:)-g_Point(point_index_list(3),:);
            dr41=g_Point(point_index_list(1),:)-g_Point(point_index_list(4),:);
            HATS_element.area=HATS_element.area+norm(cross(dr34,dr41)/2,2);
    end
    
    if length_cross_vector < eps
        disp('length_cross_vector less than 0');
    end
    
    HATS_element_list(element_index)=HATS_element;
    volume_center_point=volume_center_point+center_point;
end
volume_center_point=volume_center_point/element_number;

% project center_point to vector_flow
center_point_project_list=center_point_list*vector_flow;
[~,index_list]=sort(center_point_project_list);
center_point_list=center_point_list(index_list,:);
HATS_element_list=HATS_element_list(index_list);

% construct ADT by center point
ADtree_marker_element=ADTreeElement...
    (center_point_list,g_geometry.min_bou,g_geometry.max_bou);

% allocation element parent node and children node
addContactElement(1)
for element_index=2:element_number
    geoMakeConsistent(element_index,element_index-1)
    addContactElement(element_index)
end

% due to normal vector is continue, only need to judge one element normal
% if one element normal vector is incorrect, reverse all element
% volume_center_point_stagnation_point=ADtree_marker_element.center_point_list(1,:)-volume_center_point;
% normal_vector_project=volume_center_point_stagnation_point*HATS_element_list(1).normal_vector';
% if normal_vector_project < 0
%     % reverse all element
%     for element_index=1:element_number
%         HATS_element_list(element_index).normal_vector=-HATS_element_list(element_index).normal_vector;
%         HATS_element_list(element_index).point_index_list=fliplr(HATS_element_list(element_index).point_index_list);
%     end
% end

    function geoMakeConsistent(current_index,parent_index)
        % make normal vector consistent
        % check if normal_vector is consistent with parent normal_vector
        % if not, reverse normal_vector
        % base on public node index order to judge
        % normal vector of element and node order are right-handed rule
        point_index_list_current=HATS_element_list(current_index).point_index_list;
        point_index_list_parent=HATS_element_list(parent_index).point_index_list;
        
        % public point 1
        public_index=1;
        public_flag=0;
        while public_index <= length(point_index_list_current)
            public_parent_index=1;
            while public_parent_index <= length(point_index_list_parent)
                if point_index_list_current(public_index) == point_index_list_parent(public_parent_index)
                    public_flag=1;
                    break;
                end
                public_parent_index=public_parent_index+1;
            end
            if public_flag
                break;
            end
            public_index=public_index+1;
        end
        
        % public point 2
        public_index_sub=public_index-1;
        if public_index_sub == 0
            public_index_sub=length(point_index_list_current);
        end
        public_parent_index_sub=public_parent_index-1;
        if public_parent_index_sub == 0
            public_parent_index_sub=length(point_index_list_parent);
        end
        
        public_index_plus=public_index+1;
        if public_index_plus > length(point_index_list_current)
            public_index_plus=1;
        end
        public_parent_index_plus=public_parent_index+1;
        if public_parent_index_plus > length(point_index_list_parent)
            public_parent_index_plus=1;
        end
        
        % if edge direction do not reverse, normal should reverse
        if point_index_list_current(public_index_sub) == point_index_list_parent(public_parent_index_sub) || ...
                point_index_list_current(public_index_plus) == point_index_list_parent(public_parent_index_plus)
            HATS_element_list(current_index).normal_vector=-HATS_element_list(current_index).normal_vector;
            HATS_element_list(current_index).point_index_list=fliplr(HATS_element_list(current_index).point_index_list);
        end
    end
    function addContactElement(element_index)
        % add children to element
        %
        point_index_list__=HATS_element_list(element_index).point_index_list;
        % process symmetry
        on_symmetry_flag=0;
        if ~isempty(symmetry)
            on_symmetry_number=0;
            for point_index=1:point_number
                point_unit=point_index_list__(point_index);
                switch symmetry
                    case 'XOY'
                        dimension_index=3;
                    case 'YOZ'
                        dimension_index=1;
                    case 'ZOX'
                        dimension_index=2;
                    otherwise
                        error('preModel: unknown symmetry');
                end
                % check point if in symmetry face
                if abs(g_Point(point_unit,dimension_index)) < geometry_torlance
                    on_symmetry_number=on_symmetry_number+1;
                end
            end
            if on_symmetry_number == 2
                on_symmetry_flag=1;
            end
        end
        
        % find contact element index
        contact_index_list=findContactElement(element_index,on_symmetry_flag);
        
        % base on inflow vector judge which is inflow or outflow
        % if vertical, contact node list is out flow list
        if abs(HATS_element_list(element_index).normal_vector*vector_flow-1) < geometry_torlance
            out_index_list=contact_index_list;
            in_index_list=[];
        else
            out_index_list=[];
            in_index_list=[];
            for contact_index__=1:length(contact_index_list)
                contact_index=contact_index_list(contact_index__);
                
                if ((ADtree_marker_element.center_point_list(contact_index,:)-...
                        ADtree_marker_element.center_point_list(element_index,:))*vector_flow) > 0
                    % outflow
                    out_index_list=[out_index_list,contact_index];
                else
                    % inflow
                    in_index_list=[in_index_list,contact_index];
                end
            end
        end
        
        HATS_element_list(element_index).in_index_list=in_index_list;
        HATS_element_list(element_index).out_index_list=out_index_list;
    end
    function contact_index_list=findContactElement(element_index,on_symmetry_flag)
        % search nearby chilren element_index(node_index)
        % if on_symmetry_flag == 1, total contact element should subtract 1
        %
        point_index_list__=HATS_element_list(element_index).point_index_list;
        
        contact_index_list=[];
        search_start_index_old=element_index;
        search_end_index_old=element_index;
        search_range_unit=search_range;
        
        % loop to search contact element
        % if number of out_index_list plus parent_node_number do not equal element edge number will contiune
        % search range will increase if one search do no found enough element
        while (length(contact_index_list)+on_symmetry_flag) ~= length(point_index_list__) &&...
                search_range_unit < element_number
            search_start_index=search_start_index_old-search_range_unit;
            search_start_index=max(search_start_index,1);
            
            search_end_index=search_end_index_old+search_range_unit;
            search_end_index=min(search_end_index,element_number);
            
            % search nearby element
            element_search_index_list=[search_start_index:search_start_index_old-1,search_end_index_old+1:search_end_index];
            for search_index=1:length(element_search_index_list)
                element_search_index=element_search_index_list(search_index);
                element_point_index_list=HATS_element_list(element_search_index).point_index_list;
                
                % if have two point same, mean is nearby element
                same_point_number=0;
                for point_index__=1:length(point_index_list__)
                    for element_point_index__=1:length(element_point_index_list)
                        if element_point_index_list(element_point_index__) ==...
                                point_index_list__(point_index__)
                            same_point_number=same_point_number+1;
                        end
                    end
                end
                
                if same_point_number == 2
                    contact_index_list=[contact_index_list;element_search_index];
                end
            end
            
            % increase range
            search_range_unit=search_range_unit*2;
            search_start_index_old=search_start_index;
            search_end_index_old=search_end_index;
        end
    end
end

function Marker_Element=getMarkerElement(marker_name)
% return specified marker
% marker is include all element unit of marker
% element include element type and point index
%
global g_Marker
Marker_Element=[];
for marker_index=1:size(g_Marker,1)
    if strcmp(g_Marker{marker_index,1},marker_name)
        Marker_Element=g_Marker{marker_index,3};
    end
end
if isempty(Marker_Element)
    error('getMarkerElement: no marker found')
end
end