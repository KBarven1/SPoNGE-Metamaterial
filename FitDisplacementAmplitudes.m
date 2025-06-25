clear;
clc;

energy = readmatrix("NA1e-6_TS100_St50_2x2_L30_W556_R100_wt30_tt15_ef0_per0e-1_EnergyResults.txt");
displacement = readmatrix("NA1e-6_TS100_St50_2x2_L30_W556_R100_wt30_tt15_ef0_per0e-1_inf_eigenfreq_st_DispRes.txt");
 
% Separate by unique frame
frame = displacement(:, 1);
[unFrame, ~, unInds] = unique(frame);

% For printing
disptemp = displacement(:, [1 3 4 5 6 7 8]);
%disptemp(:, 1) = unInds;
%

%%
f = figure(1);
clf;
hold on
daspect([1 1 1])

for iFr = 1:length(unFrame)-1
    inds = frame == unFrame(iFr);
    displacementSubset = displacement(inds, :);
    
    % Split into cylinders and plates
    cylInds = displacementSubset(:, 2) == 0;
    plateInds = displacementSubset(:, 2) == 1;
    cylDisp = displacementSubset(cylInds, :);
    plateDisp = displacementSubset(plateInds, :);

    % Split into bottom and top plates
    maxPlateZ = max(plateDisp(:, 5));
    minPlateZ = min(plateDisp(:, 5));
    % Find min and max plate
    bottomPlateDisp = plateDisp(plateDisp(:, 5) == 0, :);

    % Loop through plates
    plateDispShape = size(bottomPlateDisp);

    xRec = bottomPlateDisp(:, 3) - 30.0;
    yRec = bottomPlateDisp(:, 4) - 30.0;
    zRec = bottomPlateDisp(:, 5);

    % All fitting terms
    sx = sin(pi*xRec/30.0);
    cx = cos(pi*xRec/30.0);
    sy = sin(pi*yRec/30.0);
    cy = cos(pi*yRec/30.0);
    s2x = sin(2*pi*xRec/30.0);
    c2x = cos(2*pi*xRec/30.0);
    s2y = sin(2*pi*yRec/30.0);
    c2y = cos(2*pi*yRec/30.0);

    % Order:
    % Sin(pi x / b) * Sin(pi y / b)
    % Sin(pi x / b) * Cos(pi y / b)
    % Cos(pi x / b) * Sin(pi y / b)
    % Cos(pi x / b) * Cos(pi y / b)
    % Sin(2 pi x / b) * Sin(2 pi y / b)
    % Sin(2 pi x / b) * Cos(2 pi y / b)
    % Cos(2 pi x / b) * Sin(2 pi y / b)
    % Cos(2 pi x / b) * Cos(2 pi y / b)
    % x
    % y
    % C

    % Fit X Displacement
    objDX = @(x) sum((bottomPlateDisp(:, 6) - (x(1).*sx.*sy + x(2).*sx.*cy + x(3).*cx.*sy + x(4).*cx.*cy + ...
        + x(5).*s2x.*s2y + x(6).*s2x.*c2y + x(7).*c2x.*s2y + x(8).*c2x.*c2y + ...
        x(9).*xRec + x(10).*yRec + x(11))).^2);
    solnDX = fmincon(objDX, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    fittedDX = (solnDX(1).*sx.*sy + solnDX(2).*sx.*cy + solnDX(3).*cx.*sy + solnDX(4).*cx.*cy + ...
        solnDX(5).*s2x.*s2y + solnDX(6).*s2x.*c2y + solnDX(7).*c2x.*s2y + solnDX(8).*c2x.*c2y + ...
        solnDX(9).*xRec + solnDX(10).*yRec + solnDX(11));
    
    % Fit Y Displacement
    objDY = @(x) sum((bottomPlateDisp(:, 7) - (x(1).*sx.*sy + x(2).*sx.*cy + x(3).*cx.*sy + x(4).*cx.*cy + ...
        + x(5).*s2x.*s2y + x(6).*s2x.*c2y + x(7).*c2x.*s2y + x(8).*c2x.*c2y + ...
        x(9).*xRec + x(10).*yRec + x(11))).^2);
    solnDY = fmincon(objDY, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    fittedDY = (solnDY(1).*sx.*sy + solnDY(2).*sx.*cy + solnDY(3).*cx.*sy + solnDY(4).*cx.*cy + ...
        solnDY(5).*s2x.*s2y + solnDY(6).*s2x.*c2y + solnDY(7).*c2x.*s2y + solnDY(8).*c2x.*c2y + ...
        solnDY(9).*xRec + solnDY(10).*yRec + solnDY(11));
    
    % Fit Z Displacement
    objDZ = @(x) sum((bottomPlateDisp(:, 8) - (x(1).*sx.*sy + x(2).*sx.*cy + x(3).*cx.*sy + x(4).*cx.*cy + ...
        + x(5).*s2x.*s2y + x(6).*s2x.*c2y + x(7).*c2x.*s2y + x(8).*c2x.*c2y + ...
        x(9).*xRec + x(10).*yRec + x(11))).^2);
    solnDZ = fmincon(objDZ, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    fittedDZ = (solnDZ(1).*sx.*sy + solnDZ(2).*sx.*cy + solnDZ(3).*cx.*sy + solnDZ(4).*cx.*cy + ...
        solnDZ(5).*s2x.*s2y + solnDZ(6).*s2x.*c2y + solnDZ(7).*c2x.*s2y + solnDZ(8).*c2x.*c2y + ...
        solnDZ(9).*xRec + solnDZ(10).*yRec + solnDZ(11));

    % Scatter full displacement
     scatter3(xRec+bottomPlateDisp(:,6), yRec+bottomPlateDisp(:,7), ...
         bottomPlateDisp(:, 5)+bottomPlateDisp(:,8), 'b.')
     scatter3(xRec+fittedDX(:), yRec+fittedDY(:), ...
         bottomPlateDisp(:, 5)+fittedDZ(:), 'r.')

   
    saverows(iFr, :) = [energy(unFrame(iFr), :) solnDX solnDY solnDZ];
end

%%

save W370RecBot saverows

%% Find the point that best matches our current criteria for 370

% Cylinder bending energy closest to 4.5
load("W370Rec.mat")
%cylBending = saverows(:, 3);
%[~, ids] = min(abs(cylBending - 4.5));
totalEn = saverows(:, 1)+saverows(:, 2)+saverows(:, 3)+saverows(:, 4);
[~, ids] = min(abs(totalEn/946.94 - 0.0125)); % 946.94 is the buckling pressure of alternating ellipse with these dimensions in pa

% Output matrices to csv
outputLineTop = saverows(ids, :);
load("W370RecBot.mat")
outputLineBot = saverows(ids, :);

%%

writematrix(outputLineTop, 'w370_Top_EoP.csv')
writematrix(outputLineBot, 'w370_Bot_EoP.csv')

%% Find the point that best matches our current criteria for 556

% Cylinder bending energy closest to 4.5
load("W556Rec.mat")
%cylBending = saverows(:, 3);
%[~, ids] = min(abs(cylBending - 4.5));
totalEn = saverows(:, 1)+saverows(:, 2)+saverows(:, 3)+saverows(:, 4);
[~, ids] = min(abs(totalEn/442.71 - 0.0125)); % 442.71 is the buckling pressure of diamond with these dimensions in pa

% Output matrices to csv
outputLineTop = saverows(ids, :);
load("W556RecBot.mat")
outputLineBot = saverows(ids, :);

%%

writematrix(outputLineTop, 'w556_Top_EoP.csv')
writematrix(outputLineBot, 'w556_Bot_EoP.csv')

%%

clear;
clc;
load W556Rec.mat

