%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function Kep2States = Kep2States(a, e, i, RAAN, ArgP, L)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function converts Keplerian elements into state vectors through a 
% series of calculations. The process involves determining the mean and 
% eccentric anomalies to find the true anomaly, followed by computing the
% orbit in the perifocal reference frame. By applying a transformation 
% matrix, we can determine the state vectors in inertial coordinates.

%--------------------------------------------------------------------------
% Input     - Keplerian elements of planets
% Output    - State vectors in inertial coordinates
%--------------------------------------------------------------------------

    % Constants.
    mu = 132712000000; % Standard gravitational parameter for the Sun.
    
    % Compute the argument of periapsis in the perifocal coordinate system.
    ArgP_P = deg2rad(ArgP - RAAN);
    % Calculate the mean anomaly.
    M = (L - ArgP);
    
     % Ensure mean anomaly is within [0, 360) degrees.
    if M < 0
       M = deg2rad(M + 360);
    else
       M = deg2rad(M);
    end
    
    % Calculate eccentric anomaly using iterative method.
    if M < 3.1415
       E = M + (e / 2);
    elseif M > 3.1415
       E = M - (e / 2);
    end
    
    % Declaring calculation tolerence value
    ratio = 1e-4;
    while ratio > 1e-6
        f = E - e * sin(E) - M;
        fdash = 1 - e * cos(E);
        ratio = f / fdash;
        E = E - ratio;
    end
    E = mod((rad2deg(E)), 360);
    
    % Calculate true anomaly
    theta = 2 * atan((sqrt((1 + e) / (1 - e)) * tan(deg2rad(E / 2))));
    if rad2deg(theta) < 0
        theta = rad2deg(theta) + 360;
    else
        theta = rad2deg(theta);
    end
    
    % Calculate specific angular momentum and inclination & RAAN in radians.
    h = sqrt(mu * a * (1 - e^2));
    RAAN = deg2rad(RAAN);
    i = deg2rad(i);
        
    % Transform perifocal coordinates to inertial coordinates.
    r_xbar = [((h^2 / mu) * (1 / (1 + e * cos(deg2rad(theta))))) * cos(deg2rad(theta)), ((h^2 / mu) * (1 / (1 + e * cos(deg2rad(theta))))) * sin(deg2rad(theta)), 0];
    v_xbar = (mu / h) * [-sin(deg2rad(theta)), e + cos(deg2rad(theta)), 0];
    
    % Transformation matrix from perifocal to inertial coordinates.
    Q = [-sin(RAAN) * cos(i) * sin(ArgP_P) + cos(RAAN) * cos(ArgP_P), cos(RAAN) * cos(i) * sin(ArgP_P) + sin(RAAN) * cos(ArgP_P), sin(i) * sin(ArgP_P);
         -sin(RAAN) * cos(i) * cos(ArgP_P) - cos(RAAN) * sin(ArgP_P), cos(RAAN) * cos(i) * cos(ArgP_P) - sin(RAAN) * sin(ArgP_P), sin(i) * cos(ArgP_P);
          sin(RAAN) * sin(i), -cos(RAAN) * sin(i), cos(i)];
    
    % Transpose of the transformation matrix.
    Qt = Q.';
    
    % Compute position and velocity vectors in inertial coordinates.
    r = Qt * r_xbar.';
    v = Qt * v_xbar.';
        
    Kep2States = [r(1) r(2) r(3) v(1) v(2) v(3)];
end