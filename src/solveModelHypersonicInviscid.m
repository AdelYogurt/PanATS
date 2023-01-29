function [Cl,Cd,LDratio,Cmy]=solveModelHypersonicInviscid()
% Newton method to calculate hypersonic aircraft
%
% point_list is coordinate of all node
% element_list contain element(element_type, node_index1, node_index2, ...)
% marker_list contain make{marker_name,marker_element_number,marker_element}
% marker_element contain contain element(element_type, node_index1, node_index2, ...)
%
% copyright Adel 2022.11
%
global user_model

dimension=user_model.dimension;
point_list=user_model.point_list;
marker_list=user_model.marker_list;
MARKER_MONITORING=user_model.MARKER_MONITORING;

% calculate inflow vector
free_flow_vector=[1;0;0];

AOA=user_model.AOA/180*pi;
cos_AOA=cos(AOA);
sin_AOA=sin(AOA);
rotation_AOA=[
    cos_AOA 0 -sin_AOA;
    0 1 0;
    sin_AOA 0 cos_AOA];

AOS=user_model.SIDESLIP_ANGLE/180*pi;
cos_AOS=cos(AOS);
sin_AOS=sin(AOS);
rotation_SIDESLIP_ANGLE=[
    cos_AOS -sin_AOS 0;
    sin_AOS cos_AOS 0;
    0 0 1];

free_flow_vector=rotation_AOA*rotation_SIDESLIP_ANGLE*free_flow_vector;
user_model.free_flow_vector=free_flow_vector;

% reference value
ref_point=[user_model.REF_ORIGIN_MOMENT_X,user_model.REF_ORIGIN_MOMENT_Y,user_model.REF_ORIGIN_MOMENT_Z];
ref_area=user_model.REF_AREA;
ref_length=user_model.REF_LENGTH;

T_1=user_model.FREESTREAM_TEMPERATURE;
P_1=user_model.FREESTREAM_PRESSURE;
Ma_1=user_model.MACH_NUMBER;
gama=user_model.GAMA_VALUE;
Re=user_model.REYNOLDS_NUMBER;

R=287.0955;
rou_1=P_1/R/T_1;
a_1=sqrt(gama*R*T_1);
V_1=a_1*Ma_1;
q_1=rou_1*V_1*V_1/2;

% solve prepare
Ma_1_sq=Ma_1*Ma_1;
gama_plus=gama+1;
gama_sub=gama-1;

% modified Dahlem Buck parameter
a=(6.0-0.3*Ma_1)+sin((log(Ma_1)-0.588)/1.2*pi);
n=-1.15-0.5*sin((log(Ma_1)-0.916)/3.29*pi);

% Dejarnetle parameter
G=2.6054749-0.465998*Ma_1+0.09309305*Ma_1*Ma_1+...
    -0.00817329*Ma_1*Ma_1*Ma_1+0.00026447*Ma_1*Ma_1*Ma_1*Ma_1;
D=1.0081057-0.0132323*Ma_1+0.00164956*Ma_1*Ma_1-0.00006797*Ma_1*Ma_1*Ma_1;

% Prandtl Mayer parameter
delta_max=(sqrt(gama_plus/gama_sub)-1)*pi/2-...
    sqrt(gama_plus/gama_sub)*atan(sqrt(gama_sub/gama_plus*(Ma_1_sq-1)))+...
    atan(sqrt(Ma_1_sq-1));

% modified Newton Cp_max
Cp_max=2/gama/Ma_1_sq*(...
    (gama_plus^2*Ma_1_sq/(4*gama*Ma_1_sq-2*gama_sub))^(gama/gama_sub)*((1-gama+2*gama*Ma_1_sq)/gama_plus)...
    -1);

if Ma_1 < 5
    warning('solveModelHypersonicAerodynamic: Mach number less than 5');
elseif Ma_1 < 8.5
    HIGH_HYPERSONIC_FLAG=0;
else
    HIGH_HYPERSONIC_FLAG=1;
end

% initialize data sort array
inviscid_output=repmat(...
    struct('delta_list',[],'Cp_list',[],'P_list',[],'dFn_list',[],'dMn_list',[]),...
    length(marker_list),1); % delta, Cp, P, dFn, dMn
force=zeros(1,3);
moment=zeros(1,3);

for monitor_index=1:length(MARKER_MONITORING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITORING(monitor_index),marker_list);

    delta_list=zeros(marker_list(marker_index).element_number,1);
    Cp_list=zeros(marker_list(marker_index).element_number,1);
    P_list=zeros(marker_list(marker_index).element_number,1);
    dFn_list=zeros(marker_list(marker_index).element_number,3);
    dMn_list=zeros(marker_list(marker_index).element_number,3);

    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);

        normal_vector=element.normal_vector;
        area=element.area;

        % theta is angle of negative normal vector and nomlz_vec_air
        cos_theta=-normal_vector*free_flow_vector;
        if cos_theta < -1
            theta=pi;
        elseif cos_theta > 1
            theta=0;
        else
            theta=acos(cos_theta);
        end
        % delta is attack angle
        delta=abs(pi/2-theta);

        % Cp
        if theta < pi/2
            % modified Newton
            Cp=Cp_max*cos_theta^2;
        else
            % ACM experantial
            Cp=max(-(delta*3.8197)*(1/Ma_1_sq),-(1/Ma_1_sq));
        end

        % p
        P=q_1*Cp+P_1;

        if P < 0
            disp('P < 0');
        end

        delta_list(element_index,:)=delta;
        Cp_list(element_index,:)=Cp;
        P_list(element_index,:)=P;
        dFn_list(element_index,:)=-normal_vector*area*P;
        dMn_list(element_index,:)=cross(element.center_point-ref_point,dFn_list(element_index,:));

        force=force+dFn_list(element_index,:);
        moment=moment+dFn_list(element_index,:);
    end

    inviscid_output(marker_index).delta_list=delta_list;
    inviscid_output(marker_index).Cp_list=Cp_list;
    inviscid_output(marker_index).P_list=P_list;
    inviscid_output(marker_index).dFn_list=dFn_list;
    inviscid_output(marker_index).dMn_list=dMn_list;
end

% calculate lift and drag coefficient
drag=force*free_flow_vector;
rotation_matrix=[0,0,1;
    0,1,0;
    -1,0,0]';
lift=force*(rotation_matrix*free_flow_vector);
Cl=lift/ref_area/q_1;
Cd=drag/ref_area/q_1;
LDratio=Cl/Cd;

% calculate moment
moment=sum(dMn_list,1);
moment_y=moment*[0;1;0];
Cmy=moment_y/ref_area/ref_length/q_1;

user_model.inviscid_output=inviscid_output;
end