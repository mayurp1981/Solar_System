%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function PlanetStateVec = PlanetStateVec(T0)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function determines the initial Keplerian elements of planets based 
% on the input Julian date. The following table consists of Planetary 
% Orbital Elements and their centennial rates from the book "Orbital 
% Mechanics for Engineering Students" by H. Curtis. % This function is
% explained in more detail in Algorithm 8.1.
% Subsequent conversion of state vectors is performed for the purpose 
% of orbit propagation using Kep2State function.

%-------------------------------------------------------------------------
% Input     - Initial Julian Date
% Output    - State vectors of Planets
%-------------------------------------------------------------------------

    % Constants
    AU = 1.496e+8;  % To convert astronomical unit in km
    
    % Planetary Orbital Elements and their centennial rates
    PlOE = [1, 0.38709927, 0.00000037, 0.20563593, 0.00001906, 7.00497902, -0.00594749, 48.33076593, -0.12534081, 77.45779628, 0.16047689, 252.25032350, 149472.67411175;
            2, 0.72333566, 0.00000390, 0.00677672, -0.00004107, 3.39467605, -0.00078890, 76.67984255, -0.27769418, 131.60246718, 0.00268329, 181.97909950, 58517.81538729;
            3, 1.00000261, 0.00000562, 0.01671123, -0.00004392, -0.00001531, -0.01294668, 0, 0, 102.93768193, 0.032327364, 100.46457166, 35999.37244981;
            4, 1.52371034, 0.0001847, 0.093339410, 0.00007882, 1.84969142, -0.00813131, 49.55953891, -0.29257343, -23.94362959, 0.44441088, -4.55343205, 19140.30268499;
            5, 5.20288700, -0.00011607, 0.04838624, 0.00013253, 1.30439695, -0.00183714, 100.47390909, 0.20469106, 14.72847983, 0.21252668, 34.39644501, 3034.74612775;
            6, 9.53667594, -0.00125060, 0.05386179, -0.00050991, 2.48599187, 0.00193609, 113.66242448, -0.28867794, 92.59887831, -0.41897216, 49.95424423, 1222.49362201;
            7, 19.18916464, -0.00196176, 0.04725744, -0.00004397, 0.77263783, -0.00242939, 74.01692503, 0.04240589, 170.95427630, 0.40805281, 313.23810451, 428.48202785;
            8, 30.06992276, 0.00026291, 0.00859048, 0.00005105, 1.77004347, 0.00035372, 131.78422574, -0.00508664, 44.96476227, -0.32241464, -55.12002969, 218.45945325];
    
    % Determine Keplerian elements for Mercury
    aMc = AU * (PlOE(1, 2) + PlOE(1, 3) * T0);   
    eMc = PlOE(1, 4) + PlOE(1, 5) * T0;
    iMc = (PlOE(1, 6) + PlOE(1, 7) * T0);
    if iMc < 0
        iMc = iMc + 360;
    else
        iMc = iMc;
    end
    RAANMc = (PlOE(1, 8) + PlOE(1, 9) * T0);
    if RAANMc < 0
        RAANMc = RAANMc + 360;
    else
        RAANMc = RAANMc;
    end
    ArgPMc = (PlOE(1, 10) + PlOE(1, 11) * T0);
    if ArgPMc < 0
        ArgPMc = ArgPMc + 360;
    else
        ArgPMc = ArgPMc;
    end
    LMc = mod((PlOE(1, 12) + PlOE(1, 13) * T0), 360);
    Mercury_posnvelo = Kep2States(aMc, eMc, iMc, RAANMc, ArgPMc, LMc);
    rMc = [Mercury_posnvelo(1), Mercury_posnvelo(2), Mercury_posnvelo(3)];
    vMc = [Mercury_posnvelo(4), Mercury_posnvelo(5), Mercury_posnvelo(6)];
    
    % Determine Keplerian elements for Venus
    aV = AU * (PlOE(2, 2) + PlOE(2, 3) * T0);
    eV = PlOE(2, 4) + PlOE(2, 5) * T0;
    iV = (PlOE(2, 6) + PlOE(2, 7) * T0);
    if iV < 0
        iV = iV + 360;
    else
        iV = iV;
    end
    RAANV = (PlOE(2, 8) + PlOE(2, 9) * T0);
    if RAANV < 0
        RAANV = RAANV + 360;
    else
        RAANV = RAANV;
    end
    ArgPV = (PlOE(2, 10) + PlOE(2, 11) * T0);
    if ArgPV < 0
        ArgPV = ArgPV + 360;
    else
        ArgPV = ArgPV;
    end
    LV = mod((PlOE(2, 12) + PlOE(2, 13) * T0), 360);
    Venus_posnvelo = Kep2States(aV, eV, iV, RAANV, ArgPV, LV);
    rV = [Venus_posnvelo(1), Venus_posnvelo(2), Venus_posnvelo(3)];
    vV = [Venus_posnvelo(4), Venus_posnvelo(5), Venus_posnvelo(6)];
    
    % Determine Keplerian elements for Earth
    aE = AU * (PlOE(3, 2) + PlOE(3, 3) * T0);
    eE = PlOE(3, 4) + PlOE(3, 5) * T0;
    iE = PlOE(3, 6) + PlOE(3, 7) * T0;
    RAANE = PlOE(3, 8) + PlOE(3, 9) * T0;
    ArgPE = PlOE(3, 10) + PlOE(3, 11) * T0;
    LE = mod((PlOE(3, 12) + PlOE(3, 13) * T0), 360);
    Earth_posnvelo = Kep2States(aE, eE, iE, RAANE, ArgPE, LE);
    rE = [Earth_posnvelo(1) Earth_posnvelo(2) Earth_posnvelo(3)];
    vE = [Earth_posnvelo(4) Earth_posnvelo(5) Earth_posnvelo(6)];
    
    % Determine Keplerian elements for Mars
    aM = AU * (PlOE(4, 2) + PlOE(4, 3) * T0);
    eM = PlOE(4, 4) + PlOE(4, 5) * T0;
    iM = (PlOE(4, 6) + PlOE(4, 7) * T0);
    if iM < 0
        iM = iM + 360;
    else
        iM = iM;
    end
    RAANM = (PlOE(4, 8) + PlOE(4, 9) * T0);
    if RAANM < 0
        RAANM = RAANM + 360;
    else
        RAANM = RAANM;
    end
    ArgPM = (PlOE(4, 10) + PlOE(4, 11) * T0);
    if ArgPM < 0
        ArgPM = ArgPM + 360;
    else
        ArgPM = ArgPM;
    end
    LM = mod((PlOE(4, 12) + PlOE(4, 13) * T0), 360);
    Mars_posnvelo = Kep2States(aM, eM, iM, RAANM, ArgPM, LM);
    rM = [Mars_posnvelo(1), Mars_posnvelo(2), Mars_posnvelo(3)];
    vM = [Mars_posnvelo(4), Mars_posnvelo(5), Mars_posnvelo(6)];
    
    % Determine Keplerian elements for Jupiter
    aJ = AU * (PlOE(5, 2) + PlOE(5, 3) * T0);
    eJ = PlOE(5, 4) + PlOE(5, 5) * T0;
    iJ = (PlOE(5, 6) + PlOE(5, 7) * T0);
    if iJ < 0
        iJ = iJ + 360;
    else
        iJ = iJ;
    end
    RAANJ = (PlOE(5, 8) + PlOE(5, 9) * T0);
    if RAANJ < 0
        RAANJ = RAANJ + 360;
    else
        RAANJ = RAANJ;
    end
    ArgPJ = (PlOE(5, 10) + PlOE(5, 11) * T0);
    if ArgPJ < 0
        ArgPJ = ArgPJ + 360;
    else
        ArgPJ = ArgPJ;
    end
    LJ = mod((PlOE(5, 12) + PlOE(5, 13) * T0), 360);
    Jupiter_posnvelo = Kep2States(aJ, eJ, iJ, RAANJ, ArgPJ, LJ);
    rJ = [Jupiter_posnvelo(1), Jupiter_posnvelo(2), Jupiter_posnvelo(3)];
    vJ = [Jupiter_posnvelo(4), Jupiter_posnvelo(5), Jupiter_posnvelo(6)];
    
    % Determine Keplerian elements for Saturn
    aSa = AU * (PlOE(6, 2) + PlOE(6, 3) * T0);
    eSa = PlOE(6, 4) + PlOE(6, 5) * T0;
    iSa = (PlOE(6, 6) + PlOE(6, 7) * T0);
    if iSa < 0
        iSa = iSa + 360;
    else
        iSa = iSa;
    end
    RAANSa = (PlOE(6, 8) + PlOE(6, 9) * T0);
    if RAANSa < 0
        RAANSa = RAANSa + 360;
    else
        RAANSa = RAANSa;
    end
    ArgPSa = (PlOE(6, 10) + PlOE(6, 11) * T0);
    if ArgPSa < 0
        ArgPSa = ArgPSa + 360;
    else
        ArgPSa = ArgPSa;
    end
    LSa = mod((PlOE(6, 12) + PlOE(6, 13) * T0), 360);
    Saturn_posnvelo = Kep2States(aSa, eSa, iSa, RAANSa, ArgPSa, LSa);
    rSa = [Saturn_posnvelo(1), Saturn_posnvelo(2), Saturn_posnvelo(3)];
    vSa = [Saturn_posnvelo(4), Saturn_posnvelo(5), Saturn_posnvelo(6)];
    
    % Determine Keplerian elements for Uranus
    aU = AU * (PlOE(7, 2) + PlOE(7, 3) * T0);
    eU = PlOE(7, 4) + PlOE(7, 5) * T0;
    iU = (PlOE(7, 6) + PlOE(7, 7) * T0);
    if iU < 0
        iU = iU + 360;
    else
        iU = iU;
    end
    RAANU = (PlOE(7, 8) + PlOE(7, 9) * T0);
    
    if RAANU < 0
        RAANU = RAANU + 360;
    else
        RAANU = RAANU;
    end
    ArgPU = (PlOE(7, 10) + PlOE(7, 11) * T0);
    if ArgPU < 0
        ArgPU = ArgPU + 360;
    else
        ArgPU = ArgPU;
    end
    LU = mod((PlOE(7, 12) + PlOE(7, 13) * T0), 360);
    Uranus_posnvelo = Kep2States(aU, eU, iU, RAANU, ArgPU, LU);
    rU = [Uranus_posnvelo(1), Uranus_posnvelo(2), Uranus_posnvelo(3)];
    vU = [Uranus_posnvelo(4), Uranus_posnvelo(5), Uranus_posnvelo(6)];
    
    % Determine Keplerian elements for Neptune
    aN = AU * (PlOE(8, 2) + PlOE(8, 3) * T0);
    eN = PlOE(8, 4) + PlOE(8, 5) * T0;
    iN = (PlOE(8, 6) + PlOE(8, 7) * T0);
    if iN < 0
        iN = iN + 360;
    else
        iN = iN;
    end
    RAANN = (PlOE(8, 8) + PlOE(8, 9) * T0);
    if RAANN < 0
        RAANN = RAANN + 360;
    else
        RAANN = RAANN;
    end
    ArgPN = (PlOE(8, 10) + PlOE(8, 11) * T0);
    if ArgPN < 0
        ArgPN = ArgPN + 360;
    else
        ArgPN = ArgPN;
    end
    LN = mod((PlOE(8, 12) + PlOE(8, 13) * T0), 360);
    Neptune_posnvelo = Kep2States(aN, eN, iN, RAANN, ArgPN, LN);
    rN = [Neptune_posnvelo(1), Neptune_posnvelo(2), Neptune_posnvelo(3)];
    vN = [Neptune_posnvelo(4), Neptune_posnvelo(5), Neptune_posnvelo(6)];
    
    % Combining state vectors of each planet into a single vector for simplicity.
    PlanetStateVec = [rMc vMc rV vV rE vE rM vM rJ vJ rSa vSa rU vU rN vN];
end