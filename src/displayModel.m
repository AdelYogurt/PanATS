function displayModel(type)
% display solve result
%
global g_geometry g_Point g_Element g_Marker...
    ADtree_marker_element HATS_element_list...
    streamline_output inviscid_output heat_output viscid_output FEM_output

if length(HATS_element_list) < 5000
    DRAW_FACE_FLAG=1;
else
    DRAW_FACE_FLAG=0;
end

DRAW_SYM_FLAG=1;
DRAW_INITIAL_GEO=0;

element_number=length(HATS_element_list);
center_point_list=ADtree_marker_element.center_point_list;

% inflow_direction_list=streamline_output.inflow_direction_list;

figure_result=figure();
% colormap('jet');

switch type
    case 'P'
        P_list=inviscid_output.P_list;
        drawData(P_list);
        title('Pressure');
    case 'log_P'
        P_1=g_geometry.P_1;
        log_P_list=log(inviscid_output.P_list/P_1)/log(10)+1;
        drawData(log_P_list);
        title('log Pressure');
    case 'Cp'
        Cp_list=inviscid_output.Cp_list;
        drawData(Cp_list);
        title('Surface Pressure Coefficient');
    case 'Q'
        Q_list=heat_output.Q_list;
        drawData(Q_list);
        title('Heat Flux');
    case 'Re_x'
        Re_x_list=heat_output.Re_x_list;
        drawData(Re_x_list);
        title('Surface Reynolds number');
    case 'SL'
        SL_list=streamline_output.streamline_length_list;
        drawData(SL_list);
        title('Streamline Length');
    case 'FEM'
        DOF_node=6;
        
        U_list=FEM_output.U_list;
        R_list=FEM_output.R_list;
        surf_stress_list=FEM_output.surf_stress_list;
        surf_stress_point_list=FEM_output.surf_stress_point_list;
        
        abs_delta_max=max(abs(U_list([1:DOF_node:end,2:DOF_node:end,3:DOF_node:end])));
        abs_delta_min=min(abs(U_list([1:DOF_node:end,2:DOF_node:end,3:DOF_node:end])));
        [von_surf_stress_max,node_max]=max(surf_stress_list(:,7,:));
        [von_surf_stress_max,element_max]=max(von_surf_stress_max);
        [von_surf_stress_min,node_min]=min(surf_stress_list(:,7,:));
        [von_surf_stress_min,element_min]=min(von_surf_stress_min);
        fprintf('Max Surface Von Stress: %12.4e,Element index: %4d,node index: %4d\n',...
            von_surf_stress_max,element_max,g_Element(element_max,node_max(:,:,element_max)));
        fprintf('Min Surface Von Stress: %12.4e,Element index: %4d,node index: %4d\n',...
            von_surf_stress_min,element_min,g_Element(element_min,node_min(:,:,element_min)));
        
        % drawing displacements
        % calculate zoom in/out figure scale
        min_node=min(g_Point);
        max_node=max(g_Point);
        max_distance=sqrt(sum((max_node-min_node).^2));
        delta=sqrt(sum((abs_delta_max).^2));
        if delta == 0
            scale=1;
        else
            scale=0.1*max_distance/delta;
        end
        
        element_number=size(g_Element,1);
        draw_order=[1,2,3,1];
        for element_index=1:element_number
            xpt=g_Point(g_Element(element_index,draw_order+1),1);
            ypt=g_Point(g_Element(element_index,draw_order+1),2);
            zpt=g_Point(g_Element(element_index,draw_order+1),3);
            disx=U_list(DOF_node*g_Element(element_index,draw_order+1)-5);
            disy=U_list(DOF_node*g_Element(element_index,draw_order+1)-4);
            disz=U_list(DOF_node*g_Element(element_index,draw_order+1)-3);
            surf_stress_vm=sum(surf_stress_list(draw_order,7,element_index))/4;
            patch(xpt+scale*disx,ypt+scale*disy,zpt+scale*disz,surf_stress_vm,'LineStyle','none');
        end
        title('Deformed Shape and Surface Stress-von mises Distribution');
end

colorbar;
view(3);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
drawnow;

    function drawData(data_list)
        if DRAW_FACE_FLAG
            for element_index=1:element_number
                point_index_list=HATS_element_list(element_index).point_index_list;
                patch(g_Point(point_index_list,1),g_Point(point_index_list,2),g_Point(point_index_list,3),...
                    data_list(element_index),'LineStyle','none')
                
                if DRAW_INITIAL_GEO
                    line(g_geometry.initial_point(point_index_list,1),g_geometry.initial_point(point_index_list,2),g_geometry.initial_point(point_index_list,3),...
                        'LineStyle',':','Color','r','LineWidth',0.5);
                end
                
                if DRAW_SYM_FLAG
                    switch g_geometry.symmetry
                        case 'XOY'
                            patch(g_Point(point_index_list,1),g_Point(point_index_list,2),-g_Point(point_index_list,3),...
                                data_list(element_index),'LineStyle','none')
                            if DRAW_INITIAL_GEO
                                line(g_geometry.initial_point(point_index_list,1),g_geometry.initial_point(point_index_list,2),-g_geometry.initial_point(point_index_list,3),...
                                    'LineStyle',':','Color','r','LineWidth',0.5);
                            end
                        case 'YOZ'
                            patch(-g_Point(point_index_list,1),g_Point(point_index_list,2),g_Point(point_index_list,3),...
                                data_list(element_index),'LineStyle','none')
                            if DRAW_INITIAL_GEO
                                line(-g_geometry.initial_point(point_index_list,1),g_geometry.initial_point(point_index_list,2),g_geometry.initial_point(point_index_list,3),...
                                    'LineStyle',':','Color','r','LineWidth',0.5);
                            end
                        case 'ZOX'
                            patch(g_Point(point_index_list,1),-g_Point(point_index_list,2),g_Point(point_index_list,3),...
                                data_list(element_index),'LineStyle','none')
                            if DRAW_INITIAL_GEO
                                line(g_geometry.initial_point(point_index_list,1),-g_geometry.initial_point(point_index_list,2),g_geometry.initial_point(point_index_list,3),...
                                    'LineStyle',':','Color','r','LineWidth',0.5);
                            end
                    end
                end
            end
        else
            scatter3(center_point_list(:,1),center_point_list(:,2),center_point_list(:,3),...
                [],data_list,'Marker','.');
        end
    end
end