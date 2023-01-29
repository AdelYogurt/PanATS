function displayModel(type)
% display solve result
%
global user_model

DRAW_FACE_FLAG=1;
DRAW_SYM_FLAG=1;
DRAW_INITIAL_GEO=0;

point_list=user_model.point_list;
marker_list=user_model.marker_list;

% inflow_direction_list=streamline_output.inflow_direction_list;

figure_result=figure();
% colormap('jet');

switch type
    case 'P'
        P_list={user_model.inviscid_output.P_list};
        drawData(P_list);
        title('Pressure');
        figure_result.Children.set('ColorScale','log')
    case 'log_P'
        P_1=g_geometry.P_1;
        log_P_list={log(user_model.inviscid_output.P_list/P_1)/log(10)+1};
        drawData(log_P_list);
        title('log Pressure');
        figure_result.Children.set('ColorScale','linear')
    case 'Cp'
        Cp_list={user_model.inviscid_output.Cp_list};
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
        min_node=min(point_list);
        max_node=max(point_list);
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
            xpt=point_list(g_Element(element_index,draw_order+1),1);
            ypt=point_list(g_Element(element_index,draw_order+1),2);
            zpt=point_list(g_Element(element_index,draw_order+1),3);
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
            for marker_index=1:length(data_list)
                if isempty(data_list{marker_index})
                    continue;
                end

                % only have data will be draw
                data=data_list{marker_index};
                marker_element=marker_list(marker_index).element_list;
                element_number=marker_list(marker_index).element_number;
                for element_index=1:element_number
                    point_index_list=marker_element(element_index).point_index_list;
                    patch(point_list(point_index_list,1),point_list(point_index_list,2),point_list(point_index_list,3),...
                        data(element_index),'LineStyle','none')

                    if DRAW_INITIAL_GEO
                        line(g_geometry.initial_point(point_index_list,1),g_geometry.initial_point(point_index_list,2),g_geometry.initial_point(point_index_list,3),...
                            'LineStyle',':','Color','r','LineWidth',0.5);
                    end

                    if DRAW_SYM_FLAG
                        switch user_model.SYMMETRY
                            case 'XOY'
                                patch(point_list(point_index_list,1),point_list(point_index_list,2),-point_list(point_index_list,3),...
                                    data(element_index),'LineStyle','none')
                                if DRAW_INITIAL_GEO
                                    line(g_geometry.initial_point(point_index_list,1),g_geometry.initial_point(point_index_list,2),-g_geometry.initial_point(point_index_list,3),...
                                        'LineStyle',':','Color','r','LineWidth',0.5);
                                end
                            case 'YOZ'
                                patch(-point_list(point_index_list,1),point_list(point_index_list,2),point_list(point_index_list,3),...
                                    data(element_index),'LineStyle','none')
                                if DRAW_INITIAL_GEO
                                    line(-g_geometry.initial_point(point_index_list,1),g_geometry.initial_point(point_index_list,2),g_geometry.initial_point(point_index_list,3),...
                                        'LineStyle',':','Color','r','LineWidth',0.5);
                                end
                            case 'ZOX'
                                patch(point_list(point_index_list,1),-point_list(point_index_list,2),point_list(point_index_list,3),...
                                    data(element_index),'LineStyle','none')
                                if DRAW_INITIAL_GEO
                                    line(g_geometry.initial_point(point_index_list,1),-g_geometry.initial_point(point_index_list,2),g_geometry.initial_point(point_index_list,3),...
                                        'LineStyle',':','Color','r','LineWidth',0.5);
                                end
                        end
                    end
                end
            end
        else
        end
    end
end