function preModelFEM...
    (P_1,symmetry,thickness,Elastic_modulus,Poisson_ratio,density,y_cut)
% preprocess,define total model
% using global variables
% g_Point is coordinate of all node
% g_Element is element' node
% g_Material is meterial list
% g_Border is border displace
% g_Focus is concentrated force
% g_Focus is face force
% g_Volume is volume force
%
global g_geometry g_Point g_Element g_Marker g_Material g_Border g_Focus g_Face g_Volume...
    ADtree_marker_element HATS_element_list...
    streamline_output inviscid_output heat_output viscid_output...

% format:
% x, y, z, thickness
g_Point=[g_Point(:,1:3),ones(size(g_Point,1),1)*thickness];

% format:
% element_type, point_index, element_material
g_Element=[g_Element(:,1:4),ones(size(g_Element,1),1)];

% format:
% Elastic modulus, Poisson's ratio, density
g_Material=[Elastic_modulus,Poisson_ratio,density];

% format:
% point_index, direction, displacement
g_Border=[];
point_index_list=preSYM(symmetry);
g_Border=[g_Border;
    point_index_list,ones(length(point_index_list),1),zeros(length(point_index_list),1);
    point_index_list,ones(length(point_index_list),1)*3,zeros(length(point_index_list),1);
    point_index_list,ones(length(point_index_list),1)*5,zeros(length(point_index_list),1)];
point_index_list=preInnerFace(y_cut);

% format:
g_Focus=[];

% format: (face index: down(1) up(2) anti_clockwise_form_1 (3 4 5 6))
% (pressure)element_index, face_index, pressure
% (surface stress)element_index, face_index, stress vector(1x3)
g_Face=[[1:size(g_Element,1)]',1*ones(size(g_Element,1),1),ones(size(g_Element,1),1)*1e3];
preFace()
% load pressure from aerodynamic output
% g_Face=[
%     [1:size(g_Element,1)]',2*ones(size(g_Element,1),1),inviscid_output.P_list;
%     [1:size(g_Element,1)]',1*ones(size(g_Element,1),1),ones(size(g_Element,1),1)*P_1];

% format:
% element_index, acceleration vector
% g_Volume=[[1:size(g_Element,1)]',repmat([0 0 -9.8],size(g_Element,1),1)];
g_Volume=[];

    function point_index_list=preSYM(symmetry)
        % find point on SYM face
        %
        point_index_list=[];
        
        switch symmetry
            case 'XOY'
                dimension_index=3;
                direction_list=[3;4;5];
            case 'YOZ'
                dimension_index=1;
                direction_list=[1;5;6];
            case 'ZOX'
                dimension_index=2;
                direction_list=[2;6;4];
            otherwise
                error('preModel: unknown symmetry');
        end
        
        for point_index=1:size(g_Point,1)
            if g_Point(point_index,dimension_index) < 1e-6
                point_index_list=[point_index_list;point_index];
                g_Border=[g_Border;
                    [repmat(point_index,3,1),direction_list,zeros(3,1)]];
            end
        end
    end
    function point_index_list=preInnerFace(y_cut)
        for point_index=1:size(g_Point,1)
            point_index_list=[];
            direction_list=[1;2;3;4;5;6];
            if g_Point(point_index,2) < (y_cut+1e-6)
                point_index_list=[point_index_list;point_index];
                g_Border=[g_Border;
                    [repmat(point_index,6,1),direction_list,zeros(6,1)]];
            end
        end
    end
    function preFace()
       for element_index=1:size(ADtree_marker_element.center_point_list,1)
           if ADtree_marker_element.center_point_list(element_index,3) < 0
               g_Face=[g_Face;element_index,2,100e3];
           end
       end
    end
end