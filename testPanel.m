clc;
clear;
close all hidden;

theta_list=0:pi/180:pi/2;
% theta_list=-pi/2:pi/180:0;
Cp_list=zeros(1,length(theta_list));

for theta_index=1:length(theta_list)
    Cp_list(theta_index)=calCp(theta_list(theta_index),1);
end
line(theta_list/pi*180,Cp_list,'Color','r');

for theta_index=1:length(theta_list)
    Cp_list(theta_index)=calCp(theta_list(theta_index),2);
end
line(theta_list/pi*180,Cp_list,'Color','g');

for theta_index=1:length(theta_list)
    Cp_list(theta_index)=calCp(theta_list(theta_index),3);
end
line(theta_list/pi*180,Cp_list,'Color','b');

for theta_index=1:length(theta_list)
    Cp_list(theta_index)=calCp(theta_list(theta_index),4);
end
line(theta_list/pi*180,Cp_list,'Color','k');

function Cp=calCp(theta,type)

T_1=89.3;
P_1=951.5;
Ma_1=4;
gamma=1.4;

R=287.0955;
rou_1=P_1/R/T_1;
a_1=sqrt(gamma*R*T_1);
V_1=a_1*Ma_1;
q_1=rou_1*V_1*V_1/2;

% solve prepare
Ma_1_sq=Ma_1*Ma_1;
gamma_plus=gamma+1;
gamma_sub=gamma-1;

% modified Dahlem Buck parameter
a=(6.0-0.3*Ma_1)+sin((log(Ma_1)-0.588)/1.2*pi);
n=-1.15-0.5*sin((log(Ma_1)-0.916)/3.29*pi);

% Dejarnetle parameter
G=2.6054749-0.465998*Ma_1+0.09309305*Ma_1*Ma_1+...
    -0.00817329*Ma_1*Ma_1*Ma_1+0.00026447*Ma_1*Ma_1*Ma_1*Ma_1;
D=1.0081057-0.0132323*Ma_1+0.00164956*Ma_1*Ma_1-0.00006797*Ma_1*Ma_1*Ma_1;

% Prandtl Mayer parameter
theta_max=-(sqrt(gamma_plus/gamma_sub)-1)*pi/2-...
    sqrt(gamma_plus/gamma_sub)*atan(sqrt(gamma_sub/gamma_plus*(Ma_1_sq-1)))+...
    atan(sqrt(Ma_1_sq-1));

% modified Newton Cp_max
Cp_stag=2/gamma/Ma_1_sq*(...
    (gamma_plus^2*Ma_1_sq/(4*gamma*Ma_1_sq-2*gamma_sub))^(gamma/gamma_sub)*((1-gamma+2*gamma*Ma_1_sq)/gamma_plus)...
    -1);

cos_theta_sq=cos(theta).^2;
sin_theta_sq=1-cos_theta_sq;

switch type
    case 1
        Cp=calDahlemDuck();
    case 2
        Cp=calTangentCone();
    case 3
        Cp=calTangentWedge();
    case 4
        Cp=calTangentWedgeEmpirical();
end

% switch type
%     case 1
%         Cp=calDahlemDuck();
%     case 2
%         Cp=calACMEmpirical();
%     case 3
%         Cp=calPrandtlMayerFit();
%     case 4
%         Cp=calPrandtlMeyer();
% end

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
        nuzero=2.4494897.*atan(beta./2.4494897)-atan(beta);

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

end

