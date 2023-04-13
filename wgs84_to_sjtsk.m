function [Y, X, H] = wgs84_to_sjtsk(lat, lon, alt)
%WGS84_TO_SJTSK Script converting WGS84 to S-JTSK
%   Transform method described in [1]
%   lat = latitude  n x 1 vector in radians
%   lon = longitude n x 1 vector in radians
%   alt = altitude  n x 1 vector above WGS84 ellipsoid
%
%   Y = Y axis (direction west) in meters
%   X = X axis (diraction south) in meters
%   H = elevation axis in meters


    assert(isvector(lat), "Argument lat must be vector.");
    assert(isvector(lon), "Argument lon must be vector.");
    assert(isvector(alt), "Argument alt must be vector.");
    assert(all(size(lat) == size(lon)) && all(size(lat) == size(alt)),...
        "Argument vectors must have same length.");

    lat_wgs84 = reshape(lat, 1, 1, []);
    lon_wgs84 = reshape(lon, 1, 1, []);
    alt_wgs84 = reshape(alt, 1, 1, []);

    a_wgs84     = 6378137.0;
    f_wgs84_1   = 298.257223563; % f_wgs84^-1
    a_bessel    = 6377397.15508;
    f_bessel_1  = 299.152812853; % f_bessel^-1

    e_wgs84_2   = 1 - (1 - 1/f_wgs84_1)^2; % e_wgs84^2
    e_bessel_2  = 1 - (1 - 1/f_bessel_1)^2; % e_bessel^2
    %e_wgs84     = sqrt(e_wgs84_2);
    e_bessel    = sqrt(e_bessel_2);

    rho_wgs84 = a_wgs84/((1-e_wgs84_2*sin(lat_wgs84).^2).^(1/2));
    
    r_wgs84 = [
        (rho_wgs84 + alt_wgs84).*cos(lat_wgs84).*cos(lon_wgs84);
        (rho_wgs84 + alt_wgs84).*cos(lat_wgs84).*sin(lon_wgs84);
        ((1-e_wgs84_2).*rho_wgs84 + alt_wgs84).*sin(lat_wgs84);
    ];

    R_0 = [-533.230; -75.375; -452.045];
    m = -8.750e-6;
    omega = [
        deg2rad(5.514/(60^2));
        deg2rad(2.471/(60^2));
        deg2rad(6.115/(60^2));
    ];
    cos_o = cos(omega);
    sin_o = sin(omega);
    R = [
        cos_o(2)*cos_o(3), cos_o(2)*sin_o(3), -sin_o(2);
        sin_o(1)*sin_o(2)*cos_o(3)-cos_o(1)*sin_o(3), cos_o(1)*cos_o(3)+sin_o(1)*sin_o(2)*sin_o(3), sin_o(1)*cos_o(2);
        cos_o(1)*sin_o(2)*cos_o(3)+sin_o(1)*sin_o(3), cos_o(1)*sin_o(2)*sin_o(3)-sin_o(1)*cos_o(3), cos_o(1)*cos_o(2);
    ];
    
    r_sjtsk = R_0 + (1+m)*pagemtimes(R, r_wgs84);

    p_sjtsk = (r_sjtsk(1, :, :).^2 + r_sjtsk(2, :, :).^2).^(1/2);

    lon_sjtsk = 2*atan(r_sjtsk(2, :, :)./(r_sjtsk(1, :, :) + p_sjtsk));

    t_prew = r_sjtsk(3, :, :)./((1-e_bessel_2)*p_sjtsk);
    t = zeros(size(t_prew));
    while true
        t = r_sjtsk(3, :, :)./(p_sjtsk - (a_bessel*e_bessel_2)./( (1+(1-e_bessel_2)*t_prew.^2).^(1/2) ));
        if all(abs(t-t_prew) < 1e-9)
            break
        end
        t_prew = t;
    end

    lat_sjtsk = atan(t);

    alt_sjtsk = ( (1+t.^2).^(1/2) ) .* (p_sjtsk - a_bessel/( (1+(1-e_bessel_2).*t.^2)).^(1/2) );


    lat_0 = deg2rad(49.5);
    alpha = sqrt(1 + e_bessel_2/(1-e_bessel_2) * cos(lat_0)^4);
    U_0 = asin(sin(lat_0)/alpha);

    k = tan(U_0/2 + pi/4) * ( 1/tan(lat_0/2 + pi/4) * (( 1+e_bessel*sin(lat_0) )/( 1-e_bessel*sin(lat_0) )).^(e_bessel/2) )^alpha;
    
    U = 2*atan(k*( ( ( 1-e_bessel*sin(lat_sjtsk) )./( 1+e_bessel*sin(lat_sjtsk) ) ).^(e_bessel/2) .* tan(lat_sjtsk/2 + pi/4) ).^alpha)-pi/2;
    V = alpha*(lon_sjtsk + deg2rad(53/3));

    lat_Q = deg2rad(48.25);
    lon_Q = deg2rad(42.5);

    U_Q = 2*atan(k*( ( ( 1-e_bessel*sin(lat_Q) )./( 1+e_bessel*sin(lat_Q) ) ).^(e_bessel/2) .* tan(lat_Q/2 + pi/4) ).^alpha)-pi/2 + deg2rad(11.5);
    delta_V = alpha*lon_Q - V;

    S = asin(sin(U_Q).*sin(U) + cos(U_Q).*cos(U).*cos(delta_V));
    D = asin(sin(delta_V) .* (cos(U)./cos(S)));

    S_0 = deg2rad(78.5);
    n = sin(S_0);
    
    rho_sjtsk_0 = 0.9999 * a_bessel * sqrt(1-e_bessel_2)/(1-e_bessel_2*(sin(lat_0)^2)) * 1/tan(S_0);
    
    rho_sjtsk = rho_sjtsk_0*( tan(S_0/2 + pi/4)./tan(S/2 + pi/4)).^n;
    epsilon = n*D;

    % additional tranform for correction sourced from [2]
    A = [
        0.2946529277d-01;
        0.2515965696d-01;
        0.1193845912d-06;
        -0.4668270147d-06;
        0.9233980362d-11;
        0.1523735715d-11;
        0.1696780024d-17;
        0.4408314235d-17;
        -0.8331083518d-23;
        -0.3689471323d-23;
    ];

    Y_pre = rho_sjtsk.*sin(epsilon);
    X_pre = rho_sjtsk.*cos(epsilon);
    
    Y_red = Y_pre - 654000.0;
    X_red = X_pre - 1089000.0;
    delta_Y = A(2) + A(3)*Y_red + A(4)*X_red + 2*A(5)*Y_red.*X_red + A(6)*(X_red.^2 + Y_red.^2) + A(8)*X_red.*(X_red.^2 - 3*Y_red.^2) + A(7)*Y_red.*(3*X_red.^2-Y_red.^2)-4*A(10)*Y_red.*X_red.*(X_red.^2 - Y_red.^2) + A(9)*(X_red.^4 + Y_red.^4 - 6*(X_red.^2) .* (Y_red.^2));
    delta_X = A(1) + A(3)*X_red - A(4)*Y_red - 2*A(6)*Y_red.*X_red + A(5)*(X_red.^2 + Y_red.^2) + A(7)*X_red.*(X_red.^2 - 3*Y_red.^2) - A(8)*(3*X_red.^2 - Y_red.^2) + 4*A(9)*Y_red.*X_red.*(X_red.^2 - Y_red.^2) + A(10)*(X_red.^4 + Y_red.^4 - 6*(X_red.^2) .* (Y_red.^2));

    Y = reshape(Y_pre + delta_Y, [], 1, 1);
    X = reshape(X_pre + delta_X, [], 1, 1);
    H = reshape(alt_sjtsk, [], 1, 1);
end

% [1] https://www.ibot.cas.cz/personal/wild/data/WGS_JTSK.pdf
% [2] https://www.cuzk.cz/Zememerictvi/Geodeticke-zaklady-na-uzemi-CR/GNSS/Nova-realizace-systemu-ETRS89-v-CR/Metodika-prevodu-ETRF2000-vs-S-JTSK-var2(101208).aspx

