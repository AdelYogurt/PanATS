function [model,CATS]=PanATS(config)
% driver of PanATS
%

%% input model
if isa(config,'char')
    config_filestr=config;
    config=PanATSConfig(config_filestr);
    config=config.data_dict;
elseif isa(config,'struct')
elseif isa(config,'PanATSConfig')
    config=config.data_dict;
else
    error('PanATS: error input format');
end

model.config=config;
model.geometry=struct();
model.numerics=struct();
model.output=struct();

CATS=struct();

%% pre model

[model,geom_out]=preModelPanel(model);
CATS=mergeStruct(CATS,geom_out);

%% solve model

if strcmp(model.config.SOLVER,'PANEL_INVISCID') || strcmp(model.config.SOLVER,'PANEL_VISCID')
    % inviscid
    [model,CA_out]=solveModelHypersonicInviscid(model);
    CATS=mergeStruct(CATS,CA_out);

    if strcmp(model.config.SOLVER,'PANEL_VISCID')
        % viscid
        [model,SLE_out]=solveModelStreamlineGeom(model);
        [model,BL_out]=solveModelBoundaryLayer(model);
        [model,CA_out]=solveModelHypersonicViscid(model);
        [model,CT_out]=solveModelHypersonicHeat(model);

        CATS=mergeStruct(CATS,SLE_out);
        CATS=mergeStruct(CATS,BL_out);
        CATS=mergeStruct(CATS,CA_out);
        CATS=mergeStruct(CATS,CT_out);
    end
end

%% post model

model=postModel(model);
end

function T=mergeStruct(T,S)
if isempty(S),return;end
field_list=fieldnames(S);
for k=1:numel(field_list)
    field    =field_list{k};
    T.(field)=S.(field);
end
end