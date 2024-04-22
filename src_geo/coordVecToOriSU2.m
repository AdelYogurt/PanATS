function coord_vec=coordVecToOriSU2(AOA,SIDESLIP_ANGLE,dim)
% input:
% AOA(deg), SIDESLIP_ANGLE(deg)
%

if dim == 2
    AOA=AOA/180*pi;
    cos_AOA=cos(AOA);
    sin_AOA=sin(AOA);
    rot_AOA=[
        cos_AOA -sin_AOA;
        sin_AOA cos_AOA];

    coord_vec=rot_AOA*eye(dim);

elseif dim == 3
    AOA=AOA/180*pi;
    cos_AOA=cos(AOA);
    sin_AOA=sin(AOA);
    rot_AOA=[
        cos_AOA 0 -sin_AOA;
        0 1 0;
        sin_AOA 0 cos_AOA];

    SIDESLIP_ANGLE=SIDESLIP_ANGLE/180*pi;
    cos_AOS=cos(SIDESLIP_ANGLE);
    sin_AOS=sin(SIDESLIP_ANGLE);
    rot_SIDESLIP_ANGLE=[
        cos_AOS -sin_AOS 0;
        sin_AOS cos_AOS 0;
        0 0 1];

    coord_vec=rot_AOA*rot_SIDESLIP_ANGLE*eye(dim);
else
    error('coordVecToOriSU2: dimension only can be 2 and 3');
end

end