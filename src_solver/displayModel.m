function varargout=displayModel(axe_hdl,model,type,full_salce)
% display model
%
if nargin < 4
    full_salce=[];
    if nargin < 3
        type='';
        if nargin < 2
            axe_hdl=[];
        end
    end
end

if isempty(axe_hdl),axe_hdl=gca();end
if isempty(full_salce), full_salce=false();end

config=model.config; % load config
geometry=model.geometry; % load geometry
numerics=model.numerics; % load numerics
output=model.output; % load output

% load config
SYMMETRY=config.SYMMETRY;

% load geometry
dim=geometry.dimension;
pnt_list=geometry.point_list;
elem_list=geometry.element_list;

%% draw data on axes

axe_hdl.set('FontName','times new roman','FontSize',12);
switch type
    case 'P'
        data_list=numerics.P_list;
        title('Pressure');
        axe_hdl.set('ColorScale','log')
    case 'CMz'
        data_list=numerics.dMn_list(:,3);
        title('CMz');
        axe_hdl.set('ColorScale','linear')
    case 'log_P'
        P_1=g_geometry.P_1;
        data_list={log(numerics.P_list/P_1)/log(10)+1};
        title('log Pressure');
        axe_hdl.set('ColorScale','linear')
    case 'Cp'
        data_list=numerics.Cp_list;
        title('Surface Pressure Coefficient');
        axe_hdl.set('ColorScale','log')
    case 'SL'
        data_list=numerics.streamline_len_list;
        sle_list=numerics.streamline_list;
        for sle_idx=[1,6:length(sle_list)]
            sle=sle_list{sle_idx};
            line(sle(:,1),sle(:,2),sle(:,3),'Color','r');
            if full_salce
                switch config.SYMMETRY
                    case 'XOY'
                        line(sle(:,1),sle(:,2),-sle(:,3),'Color','r');
                    case 'YOZ'
                        line(-sle(:,1),sle(:,2),sle(:,3),'Color','r');
                    case 'ZOX'
                        line(sle(:,1),-sle(:,2),sle(:,3),'Color','r');
                end
            end
        end
        elem_att_list=numerics.element_attachment_list;
        for att_idx=1:length(elem_att_list)
            elem_idx=elem_att_list(att_idx);

            if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);
            else,node_idx_list=elem_list(elem_idx,1:4);end
            node_idx_list=[node_idx_list,node_idx_list(1)];

            line(pnt_list(node_idx_list,1),pnt_list(node_idx_list,2),pnt_list(node_idx_list,3),'Color','g');
            if full_salce
                switch config.SYMMETRY
                    case 'XOY'
                        line(pnt_list(node_idx_list,1),pnt_list(node_idx_list,2),-pnt_list(node_idx_list,3),'Color','g');
                    case 'YOZ'
                        line(-pnt_list(node_idx_list,1),pnt_list(node_idx_list,2),pnt_list(node_idx_list,3),'Color','g');
                    case 'ZOX'
                        line(pnt_list(node_idx_list,1),-pnt_list(node_idx_list,2),pnt_list(node_idx_list,3),'Color','g');
                end
            end
        end
        title('Streamline Length');
    case 'rho_e'
        data_list=numerics.rho_e_list;
        title('Edge of Boundary Layer Density');
    case 'V_e'
        data_list=numerics.V_e_list;
        title('Edge of Boundary Layer Velocity');
    case 'miu_e'
        data_list=numerics.miu_e_list;
        title('Edge of Boundary Layer Viscosity');
    case 'H_e'
        data_list=numerics.H_e_list;
        title('Edge of Boundary Layer Enthalpy ');
    case 'Re_x'
        data_list=numerics.Re_x_list;
        title('Surface Reynolds Number');
    case 'HF'
        data_list=numerics.HF_list;
        title('Heat Flux');
        axe_hdl.set('ColorScale','log')
    case 'Cf'
        data_list=numerics.Cf_list;
        title('Surface Friction Coefficient');
        axe_hdl.set('ColorScale','log')
    case 'FEM'
        title('Deformed Shape and Surface Stress-von mises Distribution');
    case ''
        data_list=[];
    otherwise
        data_list=[];
end

% convert CGNS mesh data structure to matlab mesh data struct

if dim == 2

else
    if isempty(data_list)
        patch('Faces',elem_list,'Vertices',pnt_list,'FaceColor','w')
    else
        patch('CData',data_list,'FaceColor','flat','Faces',elem_list,'Vertices',pnt_list,'EdgeColor','none')
        bar_hdl=colorbar;
        bar_hdl.Location="northoutside";
        bar_hdl.FontName='times new roman';
        bar_hdl.FontSize=12;
    end

    if full_salce
        if ~isempty(SYMMETRY)
            pnt_sym_list=pnt_list;
            switch SYMMETRY
                case 'XOY'
                    pnt_sym_list(:,3)=-pnt_sym_list(:,3);
                case 'YOZ'
                    pnt_sym_list(:,1)=-pnt_sym_list(:,1);
                case 'ZOX'
                    pnt_sym_list(:,2)=-pnt_sym_list(:,2);
            end
            if isempty(data_list)
                patch('Faces',elem_list,'Vertices',pnt_sym_list,'FaceColor','none')
            else
                patch('CData',data_list,'FaceColor','flat','Faces',elem_list,'Vertices',pnt_sym_list,'EdgeColor','none')
            end
        end
    end
end

x_range=xlim();
y_range=ylim();
z_range=zlim();
center=[mean(x_range),mean(y_range),mean(z_range)];
range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
xlim([center(1)-range,center(1)+range]);
ylim([center(2)-range,center(2)+range]);
zlim([center(3)-range,center(3)+range]);

axis equal;view(3);
xlabel('\itX','FontName','times new roman','FontSize',12);
ylabel('\itY','FontName','times new roman','FontSize',12);
zlabel('\itZ','FontName','times new roman','FontSize',12);

varargout={};
if nargout > 0
    varargout={axe_hdl};
end

end
