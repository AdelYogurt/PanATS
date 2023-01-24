function [marker_element,marker_index]=getMarkerElement(marker_name,marker_list)
% function to find specified marker element
%
% input marker_name(char list), marker_list(struct)
% return marker_element, marker_index(specified marker index in marker_list)
%
marker_element=[];
for marker_index=1:length(marker_list)
    if strcmp(marker_list(marker_index).name,marker_name)
        marker_element=marker_list(marker_index).element_list;
        return;
    end
end
if isempty(marker_element)
    error('getMarkerElement: no marker found')
end
end