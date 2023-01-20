function [marker_Element,marker_index]=getMarkerElement(marker_name,marker_list)
% function to find specified marker element
%
% input marker_name(char list), marker_list(struct)
% return marker_element, marker_index(specified marker index in marker_list)
%
marker_Element=[];
for marker_index=1:size(marker_list,1)
    if strcmp(marker_list{marker_index,1},marker_name)
        marker_Element=marker_list{marker_index,3};
        return;
    end
end
if isempty(marker_Element)
    error('getMarkerElement: no marker found')
end
end