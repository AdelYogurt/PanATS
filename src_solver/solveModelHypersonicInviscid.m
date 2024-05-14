function [CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff]=solveModelHypersonicInviscid()
% base on Newton method and etc to calculate pressure of surface element
%
% copyright Adel 2023.03
%
global user_model

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

% reference value
ref_point=[config.REF_ORIGIN_MOMENT_X,config.REF_ORIGIN_MOMENT_Y,config.REF_ORIGIN_MOMENT_Z];
ref_area=config.REF_AREA;
ref_length=config.REF_LENGTH;
coord_vec=coordVecToOriSU2(config.AOA,config.SIDESLIP_ANGLE,dimension);

% calculate inflow vector
omega=[config.ANGULAR_VELOCITY_X,config.ANGULAR_VELOCITY_Y,config.ANGULAR_VELOCITY_Z]/180*pi;
free_flow_vector_list=coord_vec(:,1)'-cross(repmat(omega,length(element_list),1),center_point_list-ref_point,2);

% inflow condition
T_inf=config.FREESTREAM_TEMPERATURE;
P_inf=config.FREESTREAM_PRESSURE;
Ma_inf=config.MACH_NUMBER;
gamma=config.GAMMA_VALUE;
% Re_inf=config.REYNOLDS_NUMBER;

R=287.0955;
rho_inf=P_inf/R/T_inf;
a_inf=sqrt(gamma*R*T_inf);
V_inf=a_inf*Ma_inf;
V_inf_sq=V_inf*V_inf;
q_inf=rho_inf*V_inf_sq/2;

H_inf=(airEnthalpy(T_inf)+V_inf_sq/2); % Free-stream enthalpy J/kg equal to H_s

% solve prepare
Ma_1_sq=Ma_inf*Ma_inf;
gamma_plus=gamma+1;
gamma_sub=gamma-1;

% modified Dahlem Buck parameter
a=(6.0-0.3*Ma_inf)+sin((log(Ma_inf)-0.588)/1.2*pi);
n=-1.15-0.5*sin((log(Ma_inf)-0.916)/3.29*pi);

% Dejarnetle parameter
G=2.6054749-0.465998*Ma_inf+0.09309305*Ma_inf*Ma_inf+...
    -0.00817329*Ma_inf*Ma_inf*Ma_inf+0.00026447*Ma_inf*Ma_inf*Ma_inf*Ma_inf;
D=1.0081057-0.0132323*Ma_inf+0.00164956*Ma_inf*Ma_inf-0.00006797*Ma_inf*Ma_inf*Ma_inf;

% Prandtl Mayer parameter
theta_max=(sqrt(gamma_plus/gamma_sub)-1)*pi/2-...
    sqrt(gamma_plus/gamma_sub)*atan(sqrt(gamma_sub/gamma_plus*(Ma_1_sq-1)))+...
    atan(sqrt(Ma_1_sq-1));

% modified Newton Cp_max
Cp_stag=2/gamma/Ma_1_sq*(...
    (gamma_plus^2*Ma_1_sq/(4*gamma*Ma_1_sq-2*gamma_sub))^(gamma/gamma_sub)*((1-gamma+2*gamma*Ma_1_sq)/gamma_plus)...
    -1);

% if Ma_inf < 5
%     warning('solveModelHypersonicAerodynamic: Mach number less than 5');
%     HYPERSONIC_FLAG=0;
% elseif Ma_inf < 6.5
%     HYPERSONIC_FLAG=1;
% elseif Ma_inf < 8.5
%     HYPERSONIC_FLAG=2;
% else
%     HYPERSONIC_FLAG=3;
% end

if Ma_inf < 5
    warning('solveModelHypersonicAerodynamic: Mach number less than 5');
    HYPERSONIC_FLAG=0;
elseif Ma_inf < 10
    HYPERSONIC_FLAG=1;
else
    HYPERSONIC_FLAG=2;
end

% normal shock parameter
P_2__P_1=2*gamma/gamma_plus*Ma_1_sq-gamma_sub/gamma_plus;
rho_2__rho_1=gamma_plus*Ma_1_sq/(gamma_sub*Ma_1_sq+2);
P_2=P_2__P_1*P_inf;
rho_2=rho_2__rho_1*rho_inf;

% base on enthalpy calculate max pressure
P_max__rho_max=H_inf/(gamma/gamma_sub);
P_max__rho_max_p_gamma=P_2/rho_2^gamma;
P_max=(P_max__rho_max^gamma/P_max__rho_max_p_gamma)^(1/gamma_sub);

% calculate enthalpy correct factor
if HYPERSONIC_FLAG== 0
    CP_normal=calDahlemDuck(pi/2,1);
elseif HYPERSONIC_FLAG== 1
    CP_normal=transition(calDahlemDuck(pi/2,1),calModifyNewton(1),(Ma_inf-5)/5);
elseif HYPERSONIC_FLAG== 2
    CP_normal=calModifyNewton(1);
end
P_normal=CP_normal*q_inf+P_inf;

if P_normal > P_max, P_correct=P_max/P_normal;
else, P_correct=1;end

% fail is angle of normal vector and free flow vector
cos_fail_list=sum(normal_vector_list.*free_flow_vector_list,2);
fail_list=acos(cos_fail_list);
fail_list(cos_fail_list < -1)=pi;
fail_list(cos_fail_list > 1)=0;

% theta is attack angle, range is -pi/2 to pi/2
theta_list=fail_list-pi/2;
% notice cos(theta)^2== 1-cos_fail.^2
sin_theta_sq_list=cos_fail_list.*cos_fail_list;

% calculate surface Cp
CP_list=zeros(size(theta_list));
if HYPERSONIC_FLAG== 0
    Bool=theta_list < 0;
    CP_list(Bool)=calPrandtlMeyer(theta_list(Bool));

    Bool=theta_list >= 0 & theta_list < pi/6;
    CP_list(Bool)=calTangentWedge(theta_list(Bool),sin_theta_sq_list(Bool));

    Bool=theta_list >= pi/6 & theta_list < pi/3;
    CP_list(Bool)=calTangentCone(theta_list(Bool));
    
    Bool=theta_list >= pi/3;
    CP_list(Bool)=calDahlemDuck(theta_list(Bool),sin_theta_sq_list(Bool));
elseif HYPERSONIC_FLAG== 1
    Ma_t=(Ma_inf-5)/5;

    Bool=theta_list < 0;
    CP_list(Bool)=transition(...
        calPrandtlMeyer(theta_list(Bool)),...
        calACMEmpirical(theta_list(Bool)),...
        Ma_t);

    Bool=theta_list >= 0 & theta_list < pi/6;
    CP_list(Bool)=transition(...
        calTangentWedge(theta_list(Bool),sin_theta_sq_list(Bool)),...
        calModifyNewton(sin_theta_sq_list(Bool)),...
        Ma_t);

    Bool=theta_list >= pi/6 & theta_list < pi/3;
    CP_list(Bool)=transition(...
        calTangentCone(theta_list(Bool)),...
        calModifyNewton(sin_theta_sq_list(Bool)),...
        Ma_t);
    
    Bool=theta_list >= pi/3;
    CP_list(Bool)=transition(...
        calDahlemDuck(theta_list(Bool),sin_theta_sq_list(Bool)),...
        calModifyNewton(sin_theta_sq_list(Bool)),...
        Ma_t);
elseif HYPERSONIC_FLAG== 2
    Bool=theta_list < 0;
    CP_list(Bool)=calACMEmpirical(theta_list(Bool));
    CP_list(~Bool)=calModifyNewton(sin_theta_sq_list(~Bool));
end

% Pressure
P_list=q_inf*CP_list+P_inf;

% enthalpy correct
Bool_enthalpy_correct=theta_list > 0;
P_list(Bool_enthalpy_correct)=P_list(Bool_enthalpy_correct).*...
    (1-(1-P_correct)*sin(theta_list(Bool_enthalpy_correct)));
CP_list=(P_list-P_inf)/q_inf;

dFn_list=-normal_vector_list.*area_list.*P_list;
dMn_list=cross(center_point_list-ref_point,dFn_list);

% total force
force_inviscid=sum(dFn_list,1);
moment_inviscid=sum(dMn_list,1);

% calculate lift and drag coefficient
drag=force_inviscid*coord_vec(:,1);
slip=force_inviscid*coord_vec(:,2);
lift=force_inviscid*coord_vec(:,3);

% calculate velocity coefficient
CL=lift/ref_area/q_inf;
CD=drag/ref_area/q_inf;
CSF=slip/ref_area/q_inf;

% calculate force coefficient
CFx=force_inviscid(1)/ref_area/q_inf;
CFy=force_inviscid(2)/ref_area/q_inf;
CFz=force_inviscid(3)/ref_area/q_inf;

% calculate moment coefficient
CMx=moment_inviscid(1)/ref_area/ref_length/q_inf;
CMy=moment_inviscid(2)/ref_area/ref_length/q_inf;
CMz=moment_inviscid(3)/ref_area/ref_length/q_inf;

% calculate efficient coefficient
CEff=CL/(CD+eps);

% process SYMMETRY
if ~isempty(SYMMETRY)
    switch config.SYMMETRY
        case 'XOY'
            CFz=0;
            CMx=0;
            CMy=0;
        case 'YOZ'
            CFx=0;
            CMy=0;
            CMz=0;
        case 'ZOX'
            CFy=0;
            CMz=0;
            CMx=0;
        otherwise
            error('solveModelHypersonicInviscid: nuknown SYMMETRY type');
    end
end

output_inviscid.theta_list=theta_list;
output_inviscid.Cp_list=CP_list;
output_inviscid.P_list=P_list;
output_inviscid.dFn_list=dFn_list;
output_inviscid.dMn_list=dMn_list;
output_inviscid.force_inviscid=force_inviscid;
output_inviscid.moment_inviscid=moment_inviscid;

user_model.output_inviscid=output_inviscid;

if config.INFORMATION
    fprintf('solveModelHypersonicInviscid: hypersonic inviscid solve done!\n');
    fprintf('solveModelHypersonicInviscid: result\n');
    fprintf('CD:  %14f, CL:  %14f, CSF: %14f\nCFx:  %14f, CFy:  %14f, CFz:  %14f\nCMx: %14f, CMy: %14f, CMz: %14f\nCEff: %14f\n',...
        [CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff])
end

    function res=transition(A,B,t)
        % transition function between two condition
        %

        fc=-2*t^3+3*t^2;
        res=A*(1-fc)+B*fc;
        % res=A*cos(t*pi/2)^2+B*sin(t*pi/2)^2;
        % res=A*(1-t)+B*t;
    end

    function Cp=calDejarnetle(theta)
        % Dejarnetle modify Newton method
        %
        Cp=Cp_stag*(1-D*(cos(theta)).^G);
    end

    function Cp=calModifyNewton(sin_theta_sq)
        % modified Newton
        %
        Cp=Cp_stag*sin_theta_sq;
    end

    function Cp=calACMEmpirical(theta)
        % ACM empirical data
        %
        Cp=max((theta*3.8197)*(1/Ma_1_sq),-(1/Ma_1_sq)); % delta(deg)/15/Ma_1_sq
    end

    function Cp=calTangentCone(theta)
        % Each face element is regarded as a part of a cone whose half-top
        % Angle is equal to the Angle of the object surface of the face
        % element.
        % The conical surface pressure coefficient was obtained by
        % numerical fitting formula.
        %
        % reference: NASA TP 1539 (Appendix)
        %
        sin_theta=sin(theta);
        Ma_n=(0.87.*Ma_inf-0.544).*sin_theta + 0.53;
        Ma_n_sq=Ma_n.^2;
        Cp=48.*Ma_n_sq.*sin_theta.^2./(23.*Ma_n_sq-5);
    end

    function Cp=calTangentWedgeEmpirical(theta)
        % empirical tangent wedge method
        xx=Ma_inf.*sin(theta);
        Cp=((1.2.*xx+exp(-0.6.*xx)).^2-1.0)./(0.6.*Ma_inf.^2);
    end

    function Cp=calTangentWedge(theta,sin_theta_sq)
        % If theta is greater than 45.585 deg, then the shock is detached,
        % regardless of the value of mach. Use TangentWedgeEmpirical
        %
        Cp=theta;
        Ma_1_4=Ma_1_sq.^2;

        b=-(Ma_1_sq+2)./Ma_1_sq-gamma.*sin_theta_sq;
        c=(2.*Ma_1_sq+1)./Ma_1_4+sin_theta_sq.*((gamma_plus/2).^2+gamma_sub./Ma_1_sq);
        d=(sin_theta_sq-1)./Ma_1_4;
        q=(b.*b-3.*c)./9;
        r=(b.*(2.*b.^2-9.*c)+27.*d)./54;
        disc=q.*q.*q - r.*r;

        Bool_disc=disc > 0;
        Cp(~Bool_disc)=calTangentWedgeEmpirical(theta(~Bool_disc));

        if any(Bool_disc)
            costh=r(Bool_disc)./sqrt(q(Bool_disc).^3);
            acos_theta=acos(costh);
            c=-2.*sqrt(q(Bool_disc));
            d=b(Bool_disc)./3;
            root3=c.*cos((acos_theta+4.*pi)./3)-d;
            % root3 is supposed to be the middle root, but check it anyway
            % wave angle
            Cp(Bool_disc)=4.0.*(Ma_1_sq.*root3-1.0)./(Ma_1_sq.*(gamma+1.0));
        end

        Bool_theta=theta > 0.7956; % 45.585 deg
        Cp(Bool_theta)=calTangentWedgeEmpirical(theta(Bool_theta));

        % There can be numerical problems with very small wedge angles. Use the
        % equation on p. 92 of Liepmann and Roshko, Gasdynamics.
%         Bool_theta=theta < 0.035; % 45.585 deg
%         Cp(Bool_theta)=gamma.*Ma_1_sq.*theta(Bool_theta)./sqrt(Ma_1_sq.^2-1.0);
        
    end

    function Cp=calTangentWedgeInfiniteMach(theta)
        % PURPOSE - Calculate pressure using tangent wedge
        % infinite Mach method add theory
        %
        emns=0.5.*gamma_plus.*Ma_inf.*sin(theta)+exp(-0.25.*gamma_plus.*Ma_inf.*sin(theta));
        Cp=(4./gamma_plus).*(emns.^2-1)./Ma_inf.^2;
    end

    function Cp=calDahlemDuck(theta,sin_theta_sq)
        % modified Dahlem-Buck method
        %
        Cp=theta;
        Bool_neg=theta <= 0;
        Cp(Bool_neg)=0;

        if any(~Bool_neg)
            % first compute the original
            Cp_max=(1./sin(4*theta(~Bool_neg)).^0.75+1);
            Cp_max(theta(~Bool_neg) > 0.3927)=Cp_stag; % 22.5 deg
            Cp_max(Cp_max > 5.0)=5.0;
            Cp_max(Cp_max < Cp_stag)=Cp_stag;

            % modify for low mach number
            eta=a*(theta(~Bool_neg)/pi*180).^n+1.0;
            Cp(~Bool_neg)=eta.*Cp_max.*sin_theta_sq(~Bool_neg);
        end
    end

    function Cp=calPrandtlMayerFit(theta)
        theta(theta < theta_max)=calACMEmpirical(theta(theta < theta_max));
        Cp=-gamma_plus*theta.^2/2.*(sqrt(1+(4/gamma_plus/Ma_inf./theta).^2)-1);
        Cp(theta== 0)=0;
    end

    function Cp=calPrandtlMeyer(theta)
        % PURPOSE - Compute the pressure coefficient associated with an change of
        % freestream direction of deltar. Since deltar is negative in the calling
        % program, note that you subtract deltar from nuZero to get the nu on the
        % expansion surface.
        % radians  ( will be negative)
        %
        beta=sqrt(abs(Ma_1_sq-1.0));

        % free-stream nu
        nuzero=sqrt(6).*atan(beta./sqrt(6))-atan(beta);

        % nu on the expansion surface
        nu=nuzero-theta;

        Ma_local_sq=calInversePrandtlMeyer(nu).^2;
        bracket=(1.0+(gamma_sub/2).*Ma_local_sq)./(1.0+(gamma_sub/2).*Ma_1_sq);
        Cp=((2/gamma)./Ma_1_sq).*(bracket.^(-(gamma/gamma_sub))-1.0);

        % nu inf is 0.5*PI*(Sqrt(6)-1)
        % can't go any lower than this
        Cp(nu > 2.27685316)=-2.0./(gamma.*Ma_1_sq);
    end

    function Ma=calInversePrandtlMeyer(nu)
        % PURPOSE - Inverse Prandtl-Meyer Function. A simple rational polynomial
        % curve fit, good to 5 or 6 significant figures. Refer to the function
        % InversePrandtlMeyerPrecise if you need full doubleprecision accuracy.
        % 
        % reference: [1] Hall I M. Inversion of the Prandtl-Meyer
        % relation[J]. The Aeronautical Journal (1968), 1975, 79: 417 - 8.
        %
        y=(nu./2.27685316).^0.6666667;
        Ma=(1.0+y.*(1.3604+y.*(0.0962+y.*-0.5127)))./(1.0+y.*(-0.6722+y.*-0.3278));
    end

%     function Ma_1=calInversePrandtlMeyerprecise(nu)
%         % PURPOSE - Inverse Prandtl-Meyer function with high precision. Use Hall's
%         % approximation for a good first guess, then apply Newton's method to get
%         % greater accuracy. Instead of putting a test for convergence in the
%         % algorithm, I studied the function for various values of nu and found that
%         % four steps will give full convergence to doubleprecision (64 bits). One
%         % step would give adequate precision for single precision.
%         % Note the use of beta instead of Mach as the dependant variable until the
%         % very last step.
%         %
%         % first guess
%         beta=sqrt(abs((calInversePrandtlMeyer(nu)).^2 - one));
%         % use Newton MAX times
%         for i=1:max_ml
%             % error in nu
%             err=sqrt6.*atan(beta./sqrt6)-atan(beta)-nu;
%             betasq=beta.*beta;
%             % d(nu)/d(beta)
%             beta=beta-err.*(six+betasq).*(one+betasq)./(six.*term2.*betasq);
%         end
%         i=fix(max_ml+1);
%         Ma_1=sqrt(beta.*beta+one);
%     end

%     function Cp=calPrandtlMeyerxxx()
%         beta=sqrt(abs(Ma_1_sq-1.0));
%         %   free-stream nu }
%         nuzero=sqrt6.*atan(beta./sqrt6)-atan(beta);
%         %  should be > 0 }
%         del=-theta;
%         %      nu on the expansion surface }
%         nu=nuzero+del;
%         %   make a guess for beta - infinite Mach approx. }
%         beta=5.0./(numax-nu);
%         beta=sqrt(sqrt(nu./numax)).*5.0./(numax-nu);
%         %   error in nu at this beta}
%         err=sqrt6.*atan(beta./sqrt6)-atan(beta)-nu;
%         % use Newton-Raphson MAX times
%         for i=1:max_ml
%             betasq=beta.^2;
%             %  estimate
%             beta=beta-err.*(1+betasq./6.0).*(1+betasq)./(term2.*betasq);
%             %   new error in nu }
%             err=sqrt6.*atan(beta./sqrt6)-atan(beta)-nu;
%         end
%         i=fix(max_ml+1);
%         bracket=(1.0+term3.*(1+beta.*beta))./(1.0+term3.*Ma_1_sq);
%         Cp=(con7./Ma_1_sq).*(bracket.^(-expt) -1.0);
%     end
end
