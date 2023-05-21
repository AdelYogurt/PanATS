function [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicInviscid()
% base on Newton method and etc to calculate pressure of surface element
%
% point_list is coordinate of all node
% element_list contain []
% marker_list contain make{marker_name,marker_element_number,marker_element}
% marker_element contain HATSElement
%
% copyright Adel 2023.03
%
global user_model

dimension=user_model.geometry.dimension;

point_list=user_model.point_list;
marker_list=user_model.marker_list;

MARKER_MONITORING=user_model.MARKER_MONITORING;
SYMMETRY=user_model.SYMMETRY;

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
gamma=user_model.GAMMA_VALUE;
Re=user_model.REYNOLDS_NUMBER;

R=287.0955;
rho_1=P_1/R/T_1;
a_1=sqrt(gamma*R*T_1);
V_1=a_1*Ma_1;
q_1=rho_1*V_1*V_1/2;

% solve prepare
Ma_1_sq=Ma_1*Ma_1;
gamma_plus=gamma+1;
gamma_sub=gamma-1;

V_1_sq=V_1*V_1;
H_0=(airEnthalpy(T_1)+V_1_sq/2); % Free-stream enthalpy J/kg equal to H_s

% modified Dahlem Buck parameter
a=(6.0-0.3*Ma_1)+sin((log(Ma_1)-0.588)/1.2*pi);
n=-1.15-0.5*sin((log(Ma_1)-0.916)/3.29*pi);

% Dejarnetle parameter
G=2.6054749-0.465998*Ma_1+0.09309305*Ma_1*Ma_1+...
    -0.00817329*Ma_1*Ma_1*Ma_1+0.00026447*Ma_1*Ma_1*Ma_1*Ma_1;
D=1.0081057-0.0132323*Ma_1+0.00164956*Ma_1*Ma_1-0.00006797*Ma_1*Ma_1*Ma_1;

% Prandtl Mayer parameter
theta_max=(sqrt(gamma_plus/gamma_sub)-1)*pi/2-...
    sqrt(gamma_plus/gamma_sub)*atan(sqrt(gamma_sub/gamma_plus*(Ma_1_sq-1)))+...
    atan(sqrt(Ma_1_sq-1));

% modified Newton Cp_max
Cp_stag=2/gamma/Ma_1_sq*(...
    (gamma_plus^2*Ma_1_sq/(4*gamma*Ma_1_sq-2*gamma_sub))^(gamma/gamma_sub)*((1-gamma+2*gamma*Ma_1_sq)/gamma_plus)...
    -1);

if Ma_1 < 5
    warning('solveModelHypersonicAerodynamic: Mach number less than 5');
elseif Ma_1 < 8.5
    HIGH_HYPERSONIC_FLAG=false(1);
else
    HIGH_HYPERSONIC_FLAG=true(1);
end

% normal shock parameter
P_2__P_1=2*gamma/gamma_plus*Ma_1_sq-gamma_sub/gamma_plus;
rho_2__rho_1=gamma_plus*Ma_1_sq/(gamma_sub*Ma_1_sq+2);
P_2=P_2__P_1*P_1;
rho_2=rho_2__rho_1*rho_1;

% base on enthalpy calculate max pressure
P_max__rho_max=H_0/(gamma/gamma_sub);
P_max__rho_max_p_gamma=P_2/rho_2^gamma;
P_max=(P_max__rho_max^gamma/P_max__rho_max_p_gamma)^(1/gamma_sub);

% calculate enthalpy correct factor
theta=pi/2;
sin_theta_sq=1;
cos_theta_sq=0;
if HIGH_HYPERSONIC_FLAG
    Cp_normal=calDahlemDuck();
    P_normal=Cp_normal*q_1+P_1;
else
    Cp_normal=calModifyNewton();
    P_normal=Cp_normal*q_1+P_1;
end
if P_normal > P_max
    P_correct=P_max/P_normal;
else
    P_correct=1;
end

% initialize result sort array
% delta, Cp, P, dFn, dMn
theta_list=cell(length(marker_list),1);
Cp_list=cell(length(marker_list),1);
P_list=cell(length(marker_list),1);
dFn_list=cell(length(marker_list),1);
dMn_list=cell(length(marker_list),1);
force_inviscid=zeros(1,3);
moment_inviscid=zeros(1,3);

for monitor_index=1:length(MARKER_MONITORING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITORING(monitor_index),marker_list);

    theta_list_marker=zeros(marker_list(marker_index).element_number,1);
    Cp_list_marker=zeros(marker_list(marker_index).element_number,1);
    P_list_marker=zeros(marker_list(marker_index).element_number,1);
    dFn_list_marker=zeros(marker_list(marker_index).element_number,3);
    dMn_list_marker=zeros(marker_list(marker_index).element_number,3);

    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);

        normal_vector=element.normal_vector;
        area=element.area;

        % theta is angle of normal vector and free flow vector
        cos_fail=normal_vector*free_flow_vector;
        if cos_fail < -1
            fail=pi;
        elseif cos_fail > 1
            fail=0;
        else
            fail=acos(cos_fail);
        end
        % theta is attack angle, range is -pi/2 to pi/2
        theta=fail-pi/2;
        % notice cos(theta)^2 == 1-cos_fail.^2
        sin_theta_sq=cos_fail*cos_fail;
        cos_theta_sq=1-sin_theta_sq;

        % Cp
        if HIGH_HYPERSONIC_FLAG
            if theta < 0
                Cp=calDahlemDuck();
            elseif theta < pi/6
                Cp=calTangentWedge();
            elseif theta < pi/3
                Cp=calTangentCone();
            else
                Cp=calDahlemDuck();
            end
        else
            if theta < 0
                Cp=calACMEmpirical();
            else
                Cp=calModifyNewton();
            end
        end

        % Pressure
        P=q_1*Cp+P_1;

        % correct
        if theta > pi/3
            P=P*P_correct;
            Cp=(P-P_1)/q_1;
        end

        if P < 0
            disp('P < 0');
        end

        theta_list_marker(element_index,:)=theta;
        Cp_list_marker(element_index,:)=Cp;
        P_list_marker(element_index,:)=P;
        dFn_list_marker(element_index,:)=-normal_vector*area*P;
        dMn_list_marker(element_index,:)=cross(element.center_point-ref_point,dFn_list_marker(element_index,:));

        force_inviscid=force_inviscid+dFn_list_marker(element_index,:);
        moment_inviscid=moment_inviscid+dMn_list_marker(element_index,:);
    end

    theta_list{marker_index}=theta_list_marker;
    Cp_list{marker_index}=Cp_list_marker;
    P_list{marker_index}=P_list_marker;
    dFn_list{marker_index}=dFn_list_marker;
    dMn_list{marker_index}=dMn_list_marker;
end

% calculate lift and drag coefficient
drag=force_inviscid*free_flow_vector;
rotation_matrix=[0,0,1;
    0,1,0;
    -1,0,0]';
lift=force_inviscid*(rotation_matrix*free_flow_vector);
Cl=lift/ref_area/q_1;
Cd=drag/ref_area/q_1;
LDratio=Cl/Cd;

% calculate force coefficient
Cx=force_inviscid(1)/ref_area/q_1;
Cy=force_inviscid(2)/ref_area/q_1;
Cz=force_inviscid(3)/ref_area/q_1;

% calculate moment
Cmx=moment_inviscid(1)/ref_area/ref_length/q_1;
Cmy=moment_inviscid(2)/ref_area/ref_length/q_1;
Cmz=moment_inviscid(3)/ref_area/ref_length/q_1;

% process SYMMETRY
if ~isempty(SYMMETRY)
    switch user_model.SYMMETRY
        case 'XOY'
            Cz=0;
            Cmx=0;
            Cmy=0;
        case 'YOZ'
            Cx=0;
            Cmy=0;
            Cmz=0;
        case 'ZOX'
            Cy=0;
            Cmz=0;
            Cmx=0;
        otherwise
            error('solveModelHypersonicInviscid: nuknown SYMMETRY type');
    end
end

output_inviscid.theta_list=theta_list;
output_inviscid.Cp_list=Cp_list;
output_inviscid.P_list=P_list;
output_inviscid.dFn_list=dFn_list;
output_inviscid.dMn_list=dMn_list;
output_inviscid.force_inviscid=force_inviscid;
output_inviscid.moment_inviscid=moment_inviscid;

user_model.output_inviscid=output_inviscid;

if user_model.INFORMATION
    fprintf('solveModelHypersonicInviscid: hypersonic inviscid solve done!\n');
    fprintf('solveModelHypersonicInviscid: result\n');
    fprintf('Cl:  %14f, Cd:  %14f, L/D: %14f\nCx:  %14f, Cy:  %14f, Cz:  %14f\nCmx: %14f, Cmy: %14f, Cmz: %14f\n',...
        Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz)
end

    function Cp=calDejarnetle()
        % Dejarnetle modify Newton method
        %
        Cp=Cp_stag*(1-D*(cos(theta)).^G);
    end

    function Cp=calModifyNewton()
        % modified Newton
        %
        Cp=Cp_stag*sin_theta_sq;
    end

    function Cp=calACMEmpirical()
        % ACM empirical data
        %
        Cp=max((theta*3.8197)*(1/Ma_1_sq),-(1/Ma_1_sq)); % delta(deg)/15/Ma_1_sq
    end

    function Cp=calTangentCone()
        % Each face element is regarded as a part of a cone whose half-top
        % Angle is equal to the Angle of the object surface of the face
        % element.
        % The conical surface pressure coefficient was obtained by
        % numerical fitting formula.
        %
        % reference: NASA TP 1539 (Appendix)
        %
        sin_theta=sin(theta);
        Ma_n=(0.87.*Ma_1-0.544).*sin_theta + 0.53;
        Ma_n_sq=Ma_n.^2;
        Cp=48.*Ma_n_sq.*sin_theta.^2./(23.*Ma_n_sq-5);
    end

    function Cp=calTangentWedgeEmpirical()
        % empirical tangent wedge method
        xx=Ma_1.*sin(theta);
        Cp=((1.2.*xx+exp(-0.6.*xx)).^2-1.0)./(0.6.*Ma_1.^2);
    end

    function Cp=calTangentWedge()
        % If deltar is greater than 45.585 deg, then the shock is detached,
        % regardless of the value of mach. Use TangentWedgeEmpirical
        %
        if (theta > 0.7956) % 45.585 deg
            Cp=calTangentWedgeEmpirical();
            return;
        end

        Ma_1_4=Ma_1_sq.^2;
        sin_theta_sq=1-cos_theta_sq;

        % There can be numerical problems with very small wedge angles. Use the
        % equation on p. 92 of Liepmann and Roshko, Gasdynamics.
        if (theta < 0.035)
            Cp=gamma.*Ma_1_sq.*theta./sqrt(Ma_1_sq.^2-1.0);
        end

        b=-(Ma_1_sq+2)./Ma_1_sq-gamma.*sin_theta_sq;
        c=(2.*Ma_1_sq+1)./Ma_1_4+sin_theta_sq.*((gamma_plus/2).^2+gamma_sub./Ma_1_sq);
        d=(sin_theta_sq-1)./Ma_1_4;
        q=(b.*b-3.*c)./9;
        r=(b.*(2.*b.^2-9.*c)+27.*d)./54;
        disc=q.*q.*q - r.*r;
        if (disc < 0)
            Cp=calTangentWedgeEmpirical();
        else
            costh=r./sqrt(q.*q.*q);
            theta=acos(costh);
            c=-2.*sqrt(q);
            d=b./3;
            root1=c.*cos(theta./3)-d;
            root2=c.*cos((theta+2.*pi)./3)-d;
            root3=c.*cos((theta+4.*pi)./3)-d;
            % root3 is supposed to be the middle root, but check it anyway
            % wave angle
            theta=asin(sqrt(root3));
            Cp=4.0.*(Ma_1_sq.*root3-1.0)./(Ma_1_sq.*(gamma+1.0));
        end
    end

    function Cp=calTangentWedgeInfiniteMach()
        % PURPOSE - Calculate pressure using tangent wedge
        % infinite Mach method add theory
        %
        emns = 0.5.*gamma_plus.*Ma_1.*sin(theta)+exp(-0.25.*gamma_plus.*Ma_1.*sin(theta));
        Cp =(4./gamma_plus).*(emns.^2-1)./Ma_1.^2;
    end

    function Cp=calDahlemDuck()
        % modified Dahlem-Buck method
        %
        if (theta <= 0)
            Cp=0;
            return;
        end

        % first compute the original
        if (theta > 0.3927) % 22.5 deg
            Cp_max=2.0;
        else
            Cp_max=(1/sin(4*theta)^0.75+1);
            if(Cp_max > 5.0)
                Cp_max=5.0;
            end
            if(Cp_max < 2.0)
                Cp_max=2.0;
            end
        end

        % modify for low mach number
        eta=a*(theta/pi*180)^n+1.0;
        Cp=eta*Cp_max.*sin_theta_sq;
    end

    function Cp=calPrandtlMayerFit()
        if theta < theta_max
            theta = theta_max;
        end

        if theta == 0
            Cp=0;
        else
            Cp=-gamma_plus*theta*theta/2*(sqrt(1+(4/gamma_plus/Ma_1/theta)^2)-1);
        end
    end

    function Ma=calInversePrandtlMeyer(nu)
        % PURPOSE - Inverse Prandtl-Meyer Function. A simple rational polynomial
        % curve fit, good to 5 or 6 significant figures. Refer to the function
        % InversePrandtlMeyerPrecise if you need full doubleprecision accuracy.
        % 
        % reference: I.M. Hall Aeronautical Journal, Sept 1975, p.417
        %
        y=(nu./2.27685316).^0.6666667;
        Ma=(1.0+y.*(1.3604+y.*(0.0962+y.*-0.5127)))./(1.0+y.*(-0.6722+y.*-0.3278));
    end

    function Cp=calPrandtlMeyer()
        % PURPOSE - Compute the pressure coefficient associated with an change of
        % freestream direction of deltar. Since deltar is negative in the calling
        % program, note that you subtract deltar from nuZero to get the nu on the
        % expansion surface.
        % radians  ( will be negative)
        %
        beta=sqrt(abs(Ma_1_sq-1.0));

        % free-stream nu
        nuzero=sqrt6.*atan(beta./sqrt6)-atan(beta);

        % nu on the expansion surface
        nu=nuzero-theta;
        if(nu > 2.27685316) % NUMAX 0.5*PI*(Sqrt(6)-1)
            % can't go any lower than this
            Cp=-2.0./(gamma.*Ma_1_sq);
            return;
        end

        Ma_local_sq=calInversePrandtlMeyer(nu).^2;
        bracket=(1.0+(gamma_sub/2).*Ma_local_sq)./(1.0+(gamma_sub/2).*Ma_1_sq);
        Cp=((2/gamma)./Ma_1_sq).*(bracket.^(-(gamma/gamma_sub))-1.0);
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
