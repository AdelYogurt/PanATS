function rotate_matrix=calFreeFlowDirection(AOA,SIDESLIP_ANGLE,dim)
% input:
% AOA(deg), SIDESLIP_ANGLE(deg)
%

if dim == 2
    AOA=AOA/180*pi;
    cos_AOA=cos(AOA);
    sin_AOA=sin(AOA);
    rotation_AOA=[
        cos_AOA -sin_AOA;
        sin_AOA cos_AOA];

    rotate_matrix=rotation_AOA*eye(dim);

elseif dim == 3
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

    rotate_matrix=rotation_AOA*rotation_SIDESLIP_ANGLE*eye(dim);
else
    error('calFreeFlowDirection: dimension only can be 2 and 2');
end

end