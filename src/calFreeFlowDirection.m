function free_flow_vector=calFreeFlowDirection(AOA,SIDESLIP_ANGLE)
% input:
% AOA(deg), SIDESLIP_ANGLE(deg)
%

free_flow_vector=[1;0;0];

AOA=AOA/180*pi;
cos_AOA=cos(AOA);
sin_AOA=sin(AOA);
rotation_AOA=[
    cos_AOA 0 -sin_AOA;
    0 1 0;
    sin_AOA 0 cos_AOA];

SIDESLIP_ANGLE=SIDESLIP_ANGLE/180*pi;
cos_AOS=cos(SIDESLIP_ANGLE);
sin_AOS=sin(SIDESLIP_ANGLE);
rotation_SIDESLIP_ANGLE=[
    cos_AOS -sin_AOS 0;
    sin_AOS cos_AOS 0;
    0 0 1];

free_flow_vector=rotation_AOA*rotation_SIDESLIP_ANGLE*free_flow_vector;
end