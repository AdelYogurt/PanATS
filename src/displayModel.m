function displayModel(type)
% display solve result
%
global user_model

DRAW_FACE_FLAG=1;
DRAW_SYM_FLAG=1;

geometry=user_model.geometry;
element_list=user_model.element_list;
SYMMETRY=user_model.SYMMETRY;

dimension=geometry.dimension;
point_list=geometry.point_list;

% load geometry
center_point_list=geometry.center_point_list;
area_list=geometry.area_list;
normal_vector_list=geometry.normal_vector_list;

output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_heat=user_model.output_heat;
output_viscid=user_model.output_viscid;
output_FEM=user_model.output_FEM;
output_post=user_model.output_post;

figure_result=figure();
% colormap('jet');

switch type
    case 'P'
        P_list=output_inviscid.P_list;
        drawData(P_list);
        title('Pressure');
        figure_result.Children.set('ColorScale','log')
    case 'log_P'
        P_1=g_geometry.P_1;
        log_P_list={log(output_inviscid.P_list/P_1)/log(10)+1};
        drawData(log_P_list);
        title('log Pressure');
        figure_result.Children.set('ColorScale','linear')
    case 'Cp'
        Cp_list=output_inviscid.Cp_list;
        drawData(Cp_list);
        title('Surface Pressure Coefficient');
        figure_result.Children.set('ColorScale','log')
    case 'SL'
        SLL_list=output_streamline.streamline_len_list;
        drawData(SLL_list);
%         HF_list=output_heat.HF_list;
%         drawData(HF_list);
%         figure_result.Children.set('ColorScale','log')
        SL_list=output_streamline.streamline_list;
        for SL_index=[1,6:length(SL_list)]
            streamline=SL_list{SL_index};
            line(streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');
            if DRAW_SYM_FLAG
                switch user_model.SYMMETRY
                    case 'XOY'
                        line(streamline(:,1),streamline(:,2),-streamline(:,3),'Color','r');
                    case 'YOZ'
                        line(-streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');
                    case 'ZOX'
                        line(streamline(:,1),-streamline(:,2),streamline(:,3),'Color','r');
                end
            end
        end
        element_attach_list=output_streamline.element_attach_list;
        for attach_index=1:length(element_attach_list)
            elem=element_list(element_attach_list(attach_index));

            draw_index_list=elem.Node_idx;
            draw_index_list=[draw_index_list,draw_index_list(1)];

            line(point_list(draw_index_list,1),point_list(draw_index_list,2),point_list(draw_index_list,3),'Color','r');
            if DRAW_SYM_FLAG
                switch user_model.SYMMETRY
                    case 'XOY'
                        line(point_list(draw_index_list,1),point_list(draw_index_list,2),-point_list(draw_index_list,3),'Color','r');
                    case 'YOZ'
                        line(-point_list(draw_index_list,1),point_list(draw_index_list,2),point_list(draw_index_list,3),'Color','r');
                    case 'ZOX'
                        line(point_list(draw_index_list,1),-point_list(draw_index_list,2),point_list(draw_index_list,3),'Color','r');
                end
            end
        end
        title('Streamline Length');
    case 'rou_e'
        rou_e_list=output_heat.rou_e_list;
        drawData(rou_e_list);
        title('Edge of Boundary Layer Density');
    case 'V_e'
        V_e_list=output_heat.V_e_list;
        drawData(V_e_list);
        title('Edge of Boundary Layer Velocity');
    case 'miu_e'
        miu_e_list=output_heat.miu_e_list;
        drawData(miu_e_list);
        title('Edge of Boundary Layer Viscosity');
    case 'H_e'
        H_e_list=output_heat.H_e_list;
        drawData(H_e_list);
        title('Edge of Boundary Layer Enthalpy ');
    case 'Re_x'
        Re_x_list=output_heat.Re_x_list;
        drawData(Re_x_list);
        title('Surface Reynolds Number');
    case 'HF'
        HF_list=output_heat.HF_list;
        drawData(HF_list);
        title('Heat Flux');
        figure_result.Children.set('ColorScale','log')
    case 'Cf'
        Cf_list=output_viscid.Cf_list;
        drawData(Cf_list);
        title('Surface Friction Coefficient');
        figure_result.Children.set('ColorScale','log')
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

        elem_num=size(g_Element,1);
        draw_order=[1,2,3,1];
        for elem_idx=1:elem_num
            xpt=point_list(g_Element(elem_idx,draw_order+1),1);
            ypt=point_list(g_Element(elem_idx,draw_order+1),2);
            zpt=point_list(g_Element(elem_idx,draw_order+1),3);
            disx=U_list(DOF_node*g_Element(elem_idx,draw_order+1)-5);
            disy=U_list(DOF_node*g_Element(elem_idx,draw_order+1)-4);
            disz=U_list(DOF_node*g_Element(elem_idx,draw_order+1)-3);
            surf_stress_vm=sum(surf_stress_list(draw_order,7,elem_idx))/4;
            patch(xpt+scale*disx,ypt+scale*disy,zpt+scale*disz,surf_stress_vm,'LineStyle','none');
        end
        title('Deformed Shape and Surface Stress-von mises Distribution');
end

% colorbar;
x_range=xlim();
y_range=ylim();
z_range=zlim();

center=[mean(x_range),mean(y_range),mean(z_range)];
range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;

axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);

xlim([center(1)-range,center(1)+range]);
ylim([center(2)-range,center(2)+range]);
zlim([center(3)-range,center(3)+range]);
colorbar;
drawnow;

    function drawData(data_list)
        if DRAW_FACE_FLAG
            elem_num=length(element_list);
            for elem_idx=1:elem_num
                Node_idx=element_list(elem_idx).Node_idx;
                patch(point_list(Node_idx,1),point_list(Node_idx,2),point_list(Node_idx,3),...
                    data_list(elem_idx),'LineStyle','none')

                if DRAW_SYM_FLAG
                    switch user_model.SYMMETRY
                        case 'XOY'
                            patch(point_list(Node_idx,1),point_list(Node_idx,2),-point_list(Node_idx,3),...
                                data_list(elem_idx),'LineStyle','none')
                        case 'YOZ'
                            patch(-point_list(Node_idx,1),point_list(Node_idx,2),point_list(Node_idx,3),...
                                data_list(elem_idx),'LineStyle','none')
                        case 'ZOX'
                            patch(point_list(Node_idx,1),-point_list(Node_idx,2),point_list(Node_idx,3),...
                                data_list(elem_idx),'LineStyle','none')
                    end
                end
            end
        else
            
        end
    end
end

function drawFlowElement(element)
surface_flow=element.surface_flow;
point=element.center_point;
hold on;
quiver3(point(1),point(2),point(3),surface_flow(1),surface_flow(2),surface_flow(3),0.1);
hold off;
end