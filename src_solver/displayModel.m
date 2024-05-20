function fig_res=displayModel(type,full_salce)
% display model
%
global user_model
if nargin < 2
    full_salce=[];
    if nargin < 1
        type='';
    end
end

if isempty(full_salce), full_salce=false();end

config=user_model.config;
geometry=user_model.geometry;
element_list=user_model.element_list;
SYMMETRY=config.SYMMETRY;

% load geometry
dimension=geometry.dimension;
point_list=geometry.point_list;
center_point_list=geometry.center_point_list;
area_list=geometry.area_list;
normal_vector_list=geometry.normal_vector_list;

% load calculate result
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_boulay=user_model.output_boulay;
output_viscid=user_model.output_viscid;
output_heat=user_model.output_heat;

fig_res=figure();

% colormap('jet');

switch type
    case 'P'
        data_list=output_inviscid.P_list;
        title('Pressure');
        fig_res.Children.set('ColorScale','log')
    case 'CMz'
        data_list=output_inviscid.dMn_list(:,3);
        title('CMz');
        fig_res.Children.set('ColorScale','linear')
    case 'log_P'
        P_1=g_geometry.P_1;
        data_list={log(output_inviscid.P_list/P_1)/log(10)+1};
        title('log Pressure');
        fig_res.Children.set('ColorScale','linear')
    case 'Cp'
        data_list=output_inviscid.Cp_list;
        title('Surface Pressure Coefficient');
        fig_res.Children.set('ColorScale','log')
    case 'SL'
        data_list=output_streamline.streamline_len_list;
        SL_list=output_streamline.streamline_list;
        for SL_index=[1,6:length(SL_list)]
            streamline=SL_list{SL_index};
            line(streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');
            if full_salce
                switch config.SYMMETRY
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
            if full_salce
                switch config.SYMMETRY
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
    case 'rho_e'
        data_list=output_boulay.rho_e_list;
        title('Edge of Boundary Layer Density');
    case 'V_e'
        data_list=output_boulay.V_e_list;
        title('Edge of Boundary Layer Velocity');
    case 'miu_e'
        data_list=output_boulay.miu_e_list;
        title('Edge of Boundary Layer Viscosity');
    case 'H_e'
        data_list=output_boulay.H_e_list;
        title('Edge of Boundary Layer Enthalpy ');
    case 'Re_x'
        data_list=output_boulay.Re_x_list;
        title('Surface Reynolds Number');
    case 'HF'
        data_list=output_heat.HF_list;
        title('Heat Flux');
        fig_res.Children.set('ColorScale','log')
    case 'Cf'
        data_list=output_viscid.Cf_list;
        title('Surface Friction Coefficient');
        fig_res.Children.set('ColorScale','log')
    case 'FEM'
        title('Deformed Shape and Surface Stress-von mises Distribution');
    case ''
        data_list=[];
    otherwise
        data_list=[];
end

% convert CGNS mesh data structure to matlab mesh data struct
elem_num=length(element_list);
Elem_idx=zeros(elem_num,3);
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;
    if size(Elem_idx,2) < elem.node_num
        Elem_idx=[Elem_idx,nan(elem_num,elem.node_num-size(Elem_idx,2))];
    end
    Elem_idx(elem_idx,:)=Node_idx;
end

if isempty(data_list)
    patch('Faces',Elem_idx,'Vertices',point_list,'FaceColor','none')
else
    patch('CData',data_list,'FaceColor','flat','Faces',Elem_idx,'Vertices',point_list,'EdgeColor','none')
    colorbar;
end

if full_salce
    if ~isempty(SYMMETRY)
        point_list_sym=point_list;
        switch SYMMETRY
            case 'XOY'
                point_list_sym(:,3)=-point_list_sym(:,3);
            case 'YOZ'
                point_list_sym(:,1)=-point_list_sym(:,1);
            case 'ZOX'
                point_list_sym(:,2)=-point_list_sym(:,2);
        end
        if isempty(data_list)
            patch('Faces',Elem_idx,'Vertices',point_list_sym,'FaceColor','none')
        else
            patch('CData',data_list,'FaceColor','flat','Faces',Elem_idx,'Vertices',point_list_sym,'EdgeColor','none')
        end
    end
end

axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);

% x_range=xlim();
% y_range=ylim();
% z_range=zlim();
% center=[mean(x_range),mean(y_range),mean(z_range)];
% range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
% xlim([center(1)-range,center(1)+range]);
% ylim([center(2)-range,center(2)+range]);
% zlim([center(3)-range,center(3)+range]);

end
