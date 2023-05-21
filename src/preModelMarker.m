function [marker_list,element_empty]=preModelMarker(part_list)
% base on mesh format part_list generate marker list
% each part will generater each marker
%
part_number=length(part_list);

marker_list=repmat(struct('name',[],'element_number',[],'element_list',[]),part_number,1);
element_empty=HATSElement([],[]);

for part_index=1:part_number
    part=part_list{part_index};

    mesh_list=part.mesh_list;
    name=part.name;
    element_number=part.element_number;
    
    % allocate memory
    element_list=repmat(element_empty,element_number,1);
    element_offset=0;

    % convert each mesh into HATSElement
    for mesh_index=1:length(mesh_list)
        mesh=mesh_list{mesh_index};
        mesh_element_list=mesh.element_list;
        
        for element_index=1:mesh.element_number
            % new element
            element=HATSElement(mesh.element_ID,mesh_element_list(element_index,:));
            element.marker_index=part_index;
            element.element_index=element_index+element_offset;

            % give element
            element_list(element_index+element_offset)=element;
        end

        element_offset=+element_offset+mesh.element_number;
    end

    marker_list(part_index).name=name;
    marker_list(part_index).element_number=element_number;
    marker_list(part_index).element_list=element_list;
end

end
