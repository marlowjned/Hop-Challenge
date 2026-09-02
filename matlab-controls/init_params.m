% Centralized script for initializing all parameters relating to 
% SEB's dispersion analysis model
% Author: Marlow Nedelchev (add your name as you contribute)

% NOTE: include units for all parameters!

% Launch location parameters
initLatitude = 0;   % degrees
initLongitude = 0;  % degrees
initAltitude = 627; % m

initHeading = 45; % degrees

launchRail = 10; % m

%% ROCKET PARAMETERS
rocketAref = pi * (4 * 0.0254)^2; % m^2

%% WIND DATA/PARAMETERS
%  Wind shear data lookup
%  Preprocessed in wind_parser.py using CDS
windData = readtable('wind_combined.csv');
% TODO: in knots? ************
windPressures = flipud(windData{:, 1});
windEast = flipud(windData{:, 2});
windNorth = flipud(windData{:, 3});

gustStart = 0.5; % starts half a second into the simulation, 
                 % recommended to set to a non-zero value


%% RECOVERY PARAMS
%  Drogue
drogueCd = 0.8;
drogueAref = pi * (4 * 0.0254)^2; % m^2

%  Main
mainCd = 0.8;
mainAref = pi * (4 * 0.0254)^2; % m^2
hOpening = 2000;

earlyDeployThreshold = 1000; % prevents parachuts from deploying before 
                             % the rocket passes this altitude


%%  AERODYNAMIC COEFFICIENT LOOKUP (RASAero Export)
%%% Edited by Kelvin
% Load RASAero Aero Table
aeroData = readtable('Eureka3_Machnum_RASAero2_sim.CSV');

% Print available column names (use this to verify spelling)
%disp('Available CSV Columns:')
%disp(aeroData.Properties.VariableNames)

% Breakpoints (Unique Mach & Alpha values)
machVec  = unique(aeroData.Mach);
alphaVec = unique(aeroData.Alpha);

%  SELECT WHICH Cd COLUMN TO USE
%
%  Change ONLY the string inside quotes below
%  depending on which aerodynamic model you want.
%
%  Examples:
%     aeroData.("CDPower_On")
%     aeroData.("CDPower_Off")
%     aeroData.("CD")
%
%  ⚠ IMPORTANT:
%  Use EXACT spelling from the CSV header.
%  No spaces. Case sensitive.

CA_col = aeroData.("CAPower_Off");
CN_col = aeroData.("CN");
CP_col = aeroData.("CP");

% Build Cd Grid: size = [length(machVec) x length(alphaVec)]

[~, iM] = ismember(aeroData.Mach, machVec);
[~, iA] = ismember(aeroData.Alpha, alphaVec);

CAGrid = nan(numel(machVec), numel(alphaVec));
CNGrid = nan(numel(machVec), numel(alphaVec));
CPGrid = nan(numel(machVec), numel(alphaVec));

for i = 1:height(aeroData)
    CAGrid(iM(i), iA(i)) = CA_col(i);
    CNGrid(iM(i), iA(i)) = CN_col(i);
    CPGrid(iM(i), iA(i)) = CN_col(i);
end





%% DATA TABLE
% data from one simulation of a given rocket; set below
% rocketData = readtable('p2simulinkrefNEW.csv'); % PROSPECT 2
rocketData = readtable('p1simulinkrefNEW.csv'); % PROSPECT 1
% rocketData = readtable('ladhadsimulinkrefNEW.csv'); % LADHAD

% longitudinal dI/dt
moiderivativeL = zeros(height(rocketData),1);
for i = 1:(height(rocketData)-1)
    I = rocketData{i:i+1,"LongitudinalMomentOfInertia_kg_m__"};
    dI = I(2,1) - I(1,1);
    t = rocketData{i:i+1,"x_Time_s_"};
    dt = t(2,1) - t(1,1);
    dIdt = dI/dt;
    moiderivativeL(i,1) = dIdt;
end
moiderivativeL(height(rocketData)-1,1) = moiderivativeL(height(rocketData)-1); % make column heights match
rocketData.("LongitudialMomentOfInertiaOverTime_kg_m_s") = moiderivativeL;

% rotational dI/dt
moiderivativeR = zeros(height(rocketData),1);
for i = 1:(height(rocketData)-1)
    I = rocketData{i:i+1,"RotationalMomentOfInertia_kg_m__"};
    dI = I(2,1) - I(1,1);
    t = rocketData{i:i+1,"x_Time_s_"};
    dt = t(2,1) - t(1,1);
    dIdt = dI/dt;
    moiderivativeR(i,1) = dIdt;
end
moiderivativeR(height(rocketData)-1,1) = moiderivativeR(height(rocketData)-1); % make column heights match
rocketData.("RotationalMomentOfInertiaOverTime_kg_m_s") = moiderivativeR;

% set data to arrays
rocketTime = rocketData.("x_Time_s_");
rocketMass = rocketData.("Mass_g_");
thrust = rocketData.("Thrust_N_");
CG = rocketData.("CGLocation_cm_");
I_L = rocketData.("LongitudinalMomentOfInertia_kg_m__");
I_R = rocketData.("RotationalMomentOfInertia_kg_m__"); % oscillates kinda weirdly
I_Ldt = moiderivativeL;
I_Rdt = moiderivativeR;




%% GET W20
% wind parser results
% parser_east_wind = readtable("east_wind_means_ladhad.csv").("Var1"); % LADHAD
% parser_north_wind = readtable("north_wind_means_ladhad.csv").("Var1"); % LADHAD
% parser_pressure_levels = readtable("pressure_levels_ladhad.csv").("Var1"); % LADHAD
parser_east_wind = readtable("east_wind_means_p1.csv").("Var1"); % PROSPECT 1
parser_north_wind = readtable("north_wind_means_p1.csv").("Var1"); % PROSPECT 1
parser_pressure_levels = readtable("pressure_levels_p1.csv").("Var1"); % PROSPECT 1

% COESA model to get pressure at 20ft (6m)
[T,a,P_20ft,rho] = atmoscoesa(6);
P_20ft = P_20ft * 0.01; % hPa

% interpolate; get wind speeds at 20ft pressure
east_wind_20ft = interp1(parser_pressure_levels,parser_east_wind,P_20ft,'linear','extrap');
north_wind_20ft = interp1(parser_pressure_levels,parser_north_wind,P_20ft,'linear','extrap');

% calculate results
W20 = norm([east_wind_20ft north_wind_20ft]);
W20_angle = atan(east_wind_20ft/north_wind_20ft) * 180 / pi;





% TODO: Add custom .mat file generation/replacement when the file read
% changes

