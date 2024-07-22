function CFD=PanATS(fig_filestr)
% driver of PanATS
%
global user_model

%% pre model
config=PanATSConfig(fig_filestr);
preModelPanel(config);

%% solve model

config=user_model.config;

[area,area_x,area_y,area_z,volume]=solveGeometry();
if config.REF_AREA == 0, config.REF_AREA=area_z;end

if strcmp(config.SOLVER,'PANEL_INVISCID') || strcmp(config.SOLVER,'PANEL_VISCID')
    % inviscid
    [CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff]=solveModelHypersonicInviscid();

    if strcmp(config.SOLVER,'PANEL_VISCID')
        % viscid
        CFD.max_streamline_len=solveModelStreamline();
        solveModelBoundaryLayer();
        [CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff]=solveModelHypersonicViscid();
        CFD.max_heat_flux=solveModelHypersonicHeat();
    end
end

CFD.CD=CD;
CFD.CL=CL;
CFD.CSF=CSF;
CFD.CFx=CFx;
CFD.CFy=CFy;
CFD.CFz=CFz;
CFD.CMx=CMx;
CFD.CMy=CMy;
CFD.CMz=CMz;
CFD.CEff=CEff;

%% post model

postModel();
end