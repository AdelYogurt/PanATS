function [Cl,Cd,LDratio,Cmz,max_heat_flow]=solveModelCoupling(rou_1,V_1,FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
    ref_point,REF_LENGTH,REF_AREA,REYNOLDS_NUMBER,...
    MARKER_ISOTHERMAL,SYMMETRY,thickness,Elastic_modulus,Poisson_ratio,density,y_cut,D_x,D_z)
% solve model coupling structure deform and aerodynamic
%
global g_geometry g_Point g_Element g_Marker...
    ADtree_marker_element HATS_element_list...
    streamline_output inviscid_output heat_output viscid_output FEM_output

g_geometry.initial_point=g_Point;

LDratio_old=inf;

done=0;
iteration=1;
iteration_max=20;
torlance=1e-3;

while ~done
    if iteration ~= 1
        LDratio_old=LDratio;
    end
    
    preModelPanel(AOA,SIDESLIP_ANGLE,SYMMETRY)
    solveModelStreamline...
        ([],[],[],[],[],[],AOA,SIDESLIP_ANGLE);
    [~,~,~,~]=solveModelHypersonicInviscid...
        ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
        ref_point,REF_LENGTH,REF_AREA,REYNOLDS_NUMBER);
    max_heat_flow=solveModelHypersonicHeat...
        ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
        MARKER_ISOTHERMAL,D_x,D_z);
    [Cl,Cd,LDratio,Cmz]=solveModelHypersonicViscid...
        ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
        ref_point,REF_LENGTH,REF_AREA,REYNOLDS_NUMBER);
    preModelFEM...
        (FREESTREAM_PRESSURE,SYMMETRY,thickness,Elastic_modulus,Poisson_ratio,density,y_cut)
    solveModelFEM(0)
    
    U_list=reshape(FEM_output.U_list,6,size(g_Point,1))';
    g_Point(:,1:3)=g_Point(:,1:3)+U_list(:,1:3);
    
    if iteration >= iteration_max || abs((LDratio-LDratio_old)/LDratio) <= torlance
        done=1;
    end
    
    iteration=iteration+1;
end

end