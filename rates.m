function dydt2D = rates2D(t, y)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function computes the acceleration of each planet in the solar 
% system at certain time step (t) based on their positions and velocities. 
% It takes into account the gravitational force not only from the Sun but 
% also from every other planet, contributing to the n-body problem in 
% Astrodynamics.
%--------------------------------------------------------------------------
% Input     - time step (t) and position & velocity vectors at t (y[x, v])
% Output    - Gravitational acceleration and updated velocity of 
%             each planet (dydt[V, A])
%--------------------------------------------------------------------------
    % Constants
    G = 6.6743e-20;     % Universal Gravitational Constant
    mS = 1.989e30;      % Mass of Sun
    mMc = 330.2e21;     % Mass of Mercury
    mV = 4.869e24;      % Mass of Venus
    mE = 5.972e24;      % Mass of Earth
    mM = 641.9e21;      % Mass of Mars
    mJ = 1.899e27;      % Mass of Jupiter
    mSa = 568.5e24;     % Mass of Saturn
    mU = 86.83e24;      % Mass of Uranus
    mN = 102.4e24;      % Mass of Neptune
    
    % Extracting position components from (y) vector
    X1 = y( 1);
    Y1 = y( 2);
    
    X2 = y( 3);
    Y2 = y( 4);
    
    X3 = y( 5);
    Y3 = y( 6);
    
    X4 = y( 7);
    Y4 = y( 8);
    
    X5 = y( 9);
    Y5 = y( 10);
    
    X6 = y( 11);
    Y6 = y( 12);
    
    X7 = y( 13);
    Y7 = y( 14);
    
    X8 = y( 15);
    Y8 = y( 16);
    
    X9 = y( 17);
    Y9 = y( 18);
    
    % Extracting velocity components from (y) vector
    VX1 = y( 19);
    VY1 = y( 20);
    
    VX2 = y( 21);
    VY2 = y( 22);
    
    VX3 = y( 23);
    VY3 = y( 24);
    
    VX4 = y( 25);
    VY4 = y( 26);
    
    VX5 = y( 27);
    VY5 = y( 28);
    
    VX6 = y( 29);
    VY6 = y( 30);
    
    VX7 = y( 31);
    VY7 = y( 32);
    
    VX8 = y( 33);
    VY8 = y( 34);
    
    VX9 = y( 35);
    VY9 = y( 36);
    
    % Calculating distance between planets and sun
    R12 = norm([X2 - X1, Y2 - Y1])^3;
    R13 = norm([X3 - X1, Y3 - Y1])^3;
    R14 = norm([X4 - X1, Y4 - Y1])^3;
    R15 = norm([X5 - X1, Y5 - Y1])^3;
    R16 = norm([X6 - X1, Y6 - Y1])^3;
    R17 = norm([X7 - X1, Y7 - Y1])^3;
    R18 = norm([X8 - X1, Y8 - Y1])^3;
    R19 = norm([X9 - X1, Y9 - Y1])^3;
    R23 = norm([X3 - X2, Y3 - Y2])^3;
    R24 = norm([X4 - X2, Y4 - Y2])^3;
    R25 = norm([X5 - X2, Y5 - Y2])^3;
    R26 = norm([X6 - X2, Y6 - Y2])^3;
    R27 = norm([X7 - X2, Y7 - Y2])^3;
    R28 = norm([X8 - X2, Y8 - Y2])^3;
    R29 = norm([X9 - X2, Y9 - Y2])^3;
    R34 = norm([X4 - X3, Y4 - Y3])^3;
    R35 = norm([X5 - X3, Y5 - Y3])^3;
    R36 = norm([X6 - X3, Y6 - Y3])^3;
    R37 = norm([X7 - X3, Y7 - Y3])^3;
    R38 = norm([X8 - X3, Y8 - Y3])^3;
    R39 = norm([X9 - X3, Y9 - Y3])^3;
    R45 = norm([X5 - X4, Y5 - Y4])^3;
    R46 = norm([X6 - X4, Y6 - Y4])^3;
    R47 = norm([X7 - X4, Y7 - Y4])^3;
    R48 = norm([X8 - X4, Y8 - Y4])^3;
    R49 = norm([X9 - X4, Y9 - Y4])^3;
    R56 = norm([X6 - X5, Y6 - Y5])^3;
    R57 = norm([X7 - X5, Y7 - Y5])^3;
    R58 = norm([X8 - X5, Y8 - Y5])^3;
    R59 = norm([X9 - X5, Y9 - Y5])^3;
    R67 = norm([X7 - X6, Y7 - Y6])^3;
    R68 = norm([X8 - X6, Y8 - Y6])^3;
    R69 = norm([X9 - X6, Y9 - Y6])^3;
    R78 = norm([X8 - X7, Y8 - Y7])^3;
    R79 = norm([X9 - X7, Y9 - Y7])^3;
    R89 = norm([X9 - X8, Y9 - Y8])^3;
    
    % Calculating acceleration components of Sun
    AX1 = G * mMc * (X2 - X1) / R12 + G * mV * (X3 - X1) / R13 + G * mE * (X4 - X1) / R14 + G * mM * (X5 - X1) / R15 + G * mJ * (X6 - X1) / R16 + G * mSa * (X7 - X1) / R17 + G * mU * (X8 - X1) / R18 + G * mN * (X9 - X1) / R19;
    AY1 = G * mMc * (Y2 - Y1) / R12 + G * mV * (Y3 - Y1) / R13 + G * mE * (Y3 - Y1) / R13 + G * mM * (Y4 - Y1) / R14 + G * mJ * (Y6 - Y1) / R16 + G * mSa * (Y7 - Y1) / R17 + G * mU * (Y8 - Y1) / R18 + G * mN * (Y9 - Y1) / R19;
    
    % Calculating acceleration components of Mercury
    AX2 = G * mS * (X1 - X2) / R12 + G * mV * (X3 - X2) / R23 + G * mE * (X4 - X2) / R24 + G * mM * (X5 - X2) / R25 + G * mJ * (X6 - X2) / R26 + G * mSa * (X7 - X2) / R27 + G * mU * (X8 - X2) / R28 + G * mN * (X9 - X2) / R29;
    AY2 = G * mS * (Y1 - Y2) / R12 + G * mV * (Y3 - Y2) / R23 + G * mE * (Y4 - Y2) / R24 + G * mM * (Y5 - Y2) / R25 + G * mJ * (Y6 - Y2) / R26 + G * mSa * (Y7 - Y2) / R27 + G * mU * (Y8 - Y2) / R28 + G * mN * (Y9 - Y2) / R29;
    
    % Calculating acceleration components of Venus
    AX3 = G * mS * (X1 - X3) / R13 + G * mMc * (X2 - X3) / R23 + G * mE * (X4 - X3) / R34 + G * mM * (X5 - X3) / R35 + G * mJ * (X6 - X3) / R36 + G * mSa * (X7 - X3) / R37 + G * mU * (X8 - X3) / R38 + G * mN * (X9 - X3) / R39;
    AY3 = G * mS * (Y1 - Y3) / R13 + G * mMc * (Y2 - Y3) / R23 + G * mE * (Y4 - Y3) / R34 + G * mM * (Y5 - Y3) / R35 + G * mJ * (Y6 - Y3) / R36 + G * mSa * (Y7 - Y3) / R37 + G * mU * (Y8 - Y3) / R38 + G * mN * (Y9 - Y3) / R39;
    
    % Calculating acceleration components of Earth
    AX4 = G * mS * (X1 - X4) / R14 + G * mV * (X3 - X4) / R34 + G * mMc * (X2 - X4) / R24 + G * mM * (X5 - X4) / R45 + G * mJ * (X6 - X4) / R46 + G * mSa * (X7 - X4) / R47 + G * mU * (X8 - X4) / R48 + G * mN * (X9 - X4) / R49;
    AY4 = G * mS * (Y1 - Y4) / R14 + G * mV * (Y3 - Y4) / R34 + G * mMc * (Y2 - Y4) / R24 + G * mM * (Y5 - Y4) / R45 + G * mJ * (Y6 - Y4) / R46 + G * mSa * (Y7 - Y4) / R47 + G * mU * (Y8 - Y4) / R48 + G * mN * (Y9 - Y4) / R49;
    
    % Calculating acceleration components of Mars
    AX5 = G * mS * (X1 - X5) / R15 + G * mMc * (X2 - X5) / R25 + G * mV * (X3 - X5) / R35 + G * mE * (X4 - X5) / R45 + G * mJ * (X6 - X5) / R56 + G * mSa * (X7 - X5) / R57 + G * mU * (X8 - X5) / R58 + G * mN * (X9 - X5) / R59;
    AY5 = G * mS * (Y1 - Y5) / R15 + G * mMc * (Y2 - Y5) / R25 + G * mV * (Y3 - Y5) / R35 + G * mE * (Y4 - Y5) / R45 + G * mJ * (Y6 - Y5) / R56 + G * mSa * (Y7 - Y5) / R57 + G * mU * (Y8 - Y5) / R58 + G * mN * (Y9 - Y5) / R59;
    
    % Calculating acceleration components of Jupiter
    AX6 = G * mS * (X1 - X6) / R16 + G * mMc * (X2 - X6) / R26 + G * mV * (X3 - X6) / R36 + G * mE * (X4 - X6) / R46 + G * mM * (X5 - X6) / R56 + G * mSa * (X7 - X6) / R67 + G * mU * (X8 - X6) / R68 + G * mN * (X9 - X6) / R69;
    AY6 = G * mS * (Y1 - Y6) / R16 + G * mMc * (Y2 - Y6) / R26 + G * mV * (Y3 - Y6) / R36 + G * mE * (Y4 - Y6) / R46 + G * mM * (Y5 - Y6) / R56 + G * mSa * (Y7 - Y6) / R67 + G * mU * (Y8 - Y6) / R68 + G * mN * (Y9 - Y6) / R69;
    
    % Calculating acceleration components of Saturn
    AX7 = G * mS * (X1 - X7) / R17 + G * mMc * (X2 - X7) / R27 + G * mV * (X3 - X7) / R37 + G * mE * (X4 - X7) / R47 + G * mM * (X5 - X7) / R57 + G * mJ * (X6 - X7) / R67 + G * mU * (X8 - X7) / R78 + G * mN * (X9 - X7) / R79;
    AY7 = G * mS * (Y1 - Y7) / R17 + G * mMc * (Y2 - Y7) / R27 + G * mV * (Y3 - Y7) / R37 + G * mE * (Y4 - Y7) / R47 + G * mM * (Y5 - Y7) / R57 + G * mJ * (Y6 - Y7) / R67 + G * mU * (Y8 - Y7) / R78 + G * mN * (Y9 - Y7) / R79;
    
    % Calculating acceleration components of Uranus
    AX8 = G * mS * (X1 - X8) / R18 + G * mMc * (X2 - X8) / R28 + G * mV * (X3 - X8) / R38 + G * mE * (X4 - X8) / R48 + G * mM * (X5 - X8) / R58 + G * mJ * (X6 - X8) / R68 + G * mSa * (X7 - X8) / R78 + G * mN * (X8 - X9) / R89;
    AY8 = G * mS * (Y1 - Y8) / R18 + G * mMc * (Y2 - Y8) / R28 + G * mV * (Y3 - Y8) / R38 + G * mE * (Y4 - Y8) / R48 + G * mM * (Y5 - Y8) / R58 + G * mJ * (Y6 - Y8) / R68 + G * mSa * (Y7 - Y8) / R78 + G * mN * (Y8 - Y9) / R89;
    
    % Calculating acceleration components of Neptune
    AX9 = G * mS * (X1 - X9) / R19 + G * mMc * (X2 - X9) / R29 + G * mV * (X3 - X9) / R39 + G * mE * (X4 - X9) / R49 + G * mM * (X5 - X9) / R59 + G * mJ * (X6 - X9) / R69 + G * mSa * (X7 - X9) / R79 + G * mU * (X9 - X8) / R89;
    AY9 = G * mS * (Y1 - Y9) / R19 + G * mMc * (Y2 - Y9) / R29 + G * mV * (Y3 - Y9) / R39 + G * mE * (Y4 - Y9) / R49 + G * mM * (Y5 - Y9) / R59 + G * mJ * (Y6 - Y9) / R69 + G * mSa * (Y7 - Y9) / R79 + G * mU * (Y9 - Y8) / R89;
    
    % Combining velocities and acceleration components of all planet and sun
    dydt2D = [VX1 VY1 VX2 VY2 VX3 VY3 VX4 VY4 VX5 VY5 VX6 VY6 VX7 VY7 VX8 VY8 VX9 VY9 AX1 AY1 AX2 AY2 AX3 AY3 AX4 AY4 AX5 AY5 AX6 AY6 AX7 AY7 AX8 AY8 AX9 AY9]';
end