%ELEC 4700 
%Assignment 3
%Tariq Aboushaer
%101064544

close all;
clear;
clc;
set(0, 'DefaultFigureWindowStyle', 'docked')

%%
    global C X Y
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons
    
%% Part 1

TimeStep = 1e-14; %10fs
RegionSizeX = 200e-9; %meter
RegionSizeY = 100e-9; %meter
NumberOfParticles = 10000;
SimLength = 1000; %number of iterations of the simulation
Temperature = 300; %Kelvin


% Question 1a
% Equation used:
%  E = V/d

Voltage = 0.1; %Volts

Xregion = 2e-7;%Meters
Yregion = 1e-7;%Meters

EField = Voltage / Xregion;
disp(['Question 1a) Electric field is: ', num2str(EField), ' V/m']);


% Question 1b
% Equation used:
% F = qE

ElectronForce = C.q_0 * EField;
disp(['Question 1b) Force on each electron: ', num2str(ElectronForce), ' N']);


% Question 1c
% Equation used:
% A = F/m 

AccelerationElectron = ElectronForce/C.m_n;

disp(['Question 1c) Acceleration on each electron due to E-field: ', num2str(AccelerationElectron), ' m/s^2']);

%%

velocity_thermal = sqrt(C.kb*Temperature/C.m_n); %thermal velocity
collisionTime = 0.2e-12; %time between particle collisions


ParticleArray = [];

ParticleXPosArray = []; %array to store 7 X positions history
ParticleYPosArray = []; %array to store 7 Y positions history
ScatterPathSum = 0; %store the summed path between scatters
ScatterTimeSum = 0; %store the summed time between scatters
NumScatters = 0; %total number of particles scattered
xCurrent_Array = []; %Array to store the currents in X direction for every iteration
electron_Concentration = 10e15*10000; %m^-2
densityMap_Array = [];
temperatureMap_Array = [];
yResolution = 20;
xResolution = yResolution * 2;

scatterProbability = 1 - exp(-TimeStep/collisionTime); %probability of particle being scattered


for i = 1:NumberOfParticles
    %init positions
    ParticleArray(i,1) = rand * RegionSizeX;
    ParticleArray(i,2) = rand * RegionSizeY;
    vth = randn*velocity_thermal + velocity_thermal;
    
    %init velocity
    ParticleArray(i,3) = vth; %thermal velocity magnitude
    ParticleArray(i,4) = ((rand * 2) - 1)*vth; %xVelocity
    ParticleArray(i,5) = sqrt(ParticleArray(i,3)^2 - ParticleArray(i,4)^2); %yVelocity
    if(rand > 0.5)
        ParticleArray(i,5) = ParticleArray(i,5) * -1;
    end
    
    ParticleArray(i,6) = 0; %set time since last scatter to 0
        
end


for simCount = 1:SimLength
    currTime = simCount * TimeStep;
        
 
    for i = 1:NumberOfParticles

        
        ParticleArray(i,6) = ParticleArray(i,6) + TimeStep;
        
        if(rand <= scatterProbability) 
            ScatterPathSum = ScatterPathSum + (ParticleArray(i,6)*ParticleArray(i,3)); 
            ScatterTimeSum = ScatterTimeSum + ParticleArray(i,6); 
            ParticleArray(i,6) = 0;             
            NumScatters = NumScatters + 1;
            ParticleArray(i,3) = randn*velocity_thermal + velocity_thermal; 
            ParticleArray(i,4) = ((rand * 2) - 1)*ParticleArray(i,3); 
            ParticleArray(i,5) = sqrt(ParticleArray(i,3)^2 - ParticleArray(i,4)^2); 
            if(rand > 0.5)
                ParticleArray(i,5) = ParticleArray(i,5) * -1;
            end
        end
    end
    

    ParticleArray(:,4) = ParticleArray(:,4) + TimeStep * AccelerationElectron;
    

    new_xPos = ParticleArray(:,1) + TimeStep * ParticleArray(:,4);
    new_yPos = ParticleArray(:,2) + TimeStep * ParticleArray(:,5);    
    
    
    %checking boundary conditions
    for i = 1:NumberOfParticles
        if(new_xPos(i) < 0) %pass through (x-dir)
            new_xPos(i) = new_xPos(i)+RegionSizeX;
        elseif(new_xPos(i) > RegionSizeX) %bounce off boundary
            new_xPos(i) = new_xPos(i) - RegionSizeX;
        end
        
        if(new_yPos(i) < 0) %bounce off boundary (y-dir)
            new_yPos(i) = abs(new_yPos(i));
            ParticleArray(i,5) = ParticleArray(i,5) * -1; %swap direction
        elseif(new_yPos(i) > RegionSizeY) %bounce off boundary
            new_yPos(i) = 2*RegionSizeY - new_yPos(i);
            ParticleArray(i,5) = ParticleArray(i,5) * -1; %swap direction
        end
       
    end
    
    ParticleArray(:,1) = new_xPos;
    ParticleArray(:,2) = new_yPos;

    ParticleXPosArray(1:7,simCount) = ParticleArray(1:7,1);
    ParticleYPosArray(1:7,simCount) = ParticleArray(1:7,2);

    xCurrent_Array(simCount) = C.q_0 * electron_Concentration * mean(ParticleArray(:,4)) * RegionSizeX;


    
    if(simCount == SimLength) %last iteration, populate Maps
        
        xDiv = RegionSizeX / xResolution;
        yDiv = RegionSizeY / yResolution;
        
        for yCount = 1:yResolution
            for xCount = 1:xResolution
                densityMap_Array(xCount,yCount) = 0;
                temperatureMap_Array(xCount,yCount) = 0;
                
                for eCount = 1:NumberOfParticles
                    if(ParticleArray(eCount,1) >= (xCount*xDiv-xDiv) && ParticleArray(eCount,1) < (xCount*xDiv))
                        if(ParticleArray(eCount,2) >= (yCount*yDiv-yDiv) && ParticleArray(eCount,2) < (yCount*yDiv))
                            densityMap_Array(xCount,yCount) = densityMap_Array(xCount,yCount) + 1;
                            temptemp = sqrt(ParticleArray(eCount,4).^2+ParticleArray(eCount,5).^2)*C.m_n/(3*C.kb);
                            temperatureMap_Array(xCount,yCount) = temperatureMap_Array(xCount,yCount) + temptemp;
                        end                        
                    end
                end                
            end
        end   
        temperatureMap_Array = temperatureMap_Array ./ densityMap_Array;
    end
end


disp('Simulation Specs:');
disp(['Timestep Size: ', num2str(TimeStep), ' seconds']);
disp(['Number of Simulation Iterations: ', num2str(SimLength)]);
disp(['Number of Particles: ', num2str(NumberOfParticles)]);
disp('----------');
disp('Simulation Results:');
disp(['Mean Free Path (MFP): ', num2str(ScatterPathSum/NumScatters), ' meters']);
disp(['Time between Collisions: ', num2str(ScatterTimeSum/NumScatters), ' seconds']);

figure(13);
hold on;
for i = 1:7
    plot(ParticleXPosArray(i,:),ParticleYPosArray(i,:));
end
xlim([0,RegionSizeX]);
ylim([0,RegionSizeY]);
title(['Particle Trajectory, Simulation Count: ', num2str(simCount), ', Timestep: ', num2str(TimeStep), '(TA 101064544)']);
xlabel('X (m)');
ylabel('Y (m)');

figure(14);
x = linspace(TimeStep,SimLength*TimeStep,SimLength);
plot(x,xCurrent_Array);
title('X Current over Time (TA 101064544)');
xlabel('Time (s)');
ylabel('Current (A)');

figure(15);
[xMesh, yMesh] = meshgrid(linspace(0,RegionSizeX,xResolution),linspace(0,RegionSizeY,yResolution));
surf(xMesh,yMesh,transpose(densityMap_Array));
title('Density Map (TA 101064544)');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Count');

figure(152);
[xMesh, yMesh] = meshgrid(linspace(0,RegionSizeX,xResolution),linspace(0,RegionSizeY,yResolution));
surf(xMesh,yMesh,transpose(temperatureMap_Array));
title('Temperature Map (TA 101064544)');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Temperature (K)');
pause(0.1);

%% Part 2

disp('Question 2:');
SizeY = 100; %size of the matrix
SizeX = SizeY*2; %size of the matrix
xDel = RegionSizeX/SizeX; %delta X
yDel = RegionSizeY/SizeY; %delta Y
maxIterations = 10000; %maximum number of iterations

vMatrix = zeros(SizeY,SizeX); %create our matrix of values
oldvMatrix = zeros(SizeY,SizeX); %another matrix to be used to clone the new one
vecMatrix = zeros(SizeY,SizeX); %matrix of the vectors
[xMesh, yMesh] = meshgrid(linspace(0,RegionSizeX,SizeX),linspace(0,RegionSizeY,SizeY));

%set boundary values
TBC_Left = 0.1;
TBC_Right = 0;
TBC_Top = 0; %this will be unbounded
TBC_Bottom = 0; %this will be unbounded

barrier_Height = 30;
barrier_Width = 40;

for i = 1:maxIterations %keep running until reached max iterations
    for m = 1:SizeX %loop through columns
        for n = 1:SizeY %loop through rows
            
            if(m == 1) %left boundary
                vMatrix(n,m) = TBC_Left;
            elseif(m == SizeX) %right boundary
                vMatrix(n,m) = TBC_Right;
            elseif(n == 1) %top boundary
                %vMatrix(n,m) = TBC_Top;
                vMatrix(n,m) = (oldvMatrix(n,m-1) + oldvMatrix(n,m+1) + oldvMatrix(n+1,m))/3;
            elseif(n == SizeY) %bottom boundary
                %vMatrix(n,m) = TBC_Bottom;
                vMatrix(n,m) = (oldvMatrix(n,m-1) + oldvMatrix(n,m+1) + oldvMatrix(n-1,m))/3;
            else
                vMatrix(n,m) = (oldvMatrix(n,m-1) + oldvMatrix(n,m+1) + oldvMatrix(n-1,m) + oldvMatrix(n+1,m))/4;
            end   
            
            if(m>=SizeX/2-barrier_Width/2 && m<=SizeX/2+barrier_Width/2 && (n<=barrier_Height || n>=SizeY-barrier_Height)) %inside one of the boxes
                vMatrix(n,m) = vMatrix(n,m)*10^-2;
            end
        end
    end        
    [xVectors, yVectors] = gradient(oldvMatrix(:,:));   

    oldvMatrix = vMatrix;
end

figure(21)
surf(xMesh, yMesh, vMatrix);
shading interp
title('V(x,y) for Bottle-neck (TA 101064544)');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Voltage (V)');

figure(22)
quiver(xMesh, yMesh); % -xVectors, -yVectors
title('Electric Field (TA 101064544)');
xlabel('X (m)');
ylabel('Y (m)');


particleArray = [];


particleXPos7Array = []; %array to store 7 X positions history
particleYPos7Array = []; %array to store 7 Y positions history
scatterPathSum = 0; %store the summed path between scatters
scatterTimeSum = 0; %store the summed time between scatters
numScatters = 0; %total number of particles scattered
xCurrent_Array = []; %Array to store the currents in X direction for every iteration
electron_Concentration = 10e15*10000; %m^-2
densityMap_Array = [];
temperatureMap_Array = [];
yResolution = 20;
xResolution = yResolution * 2;

scatterProbability = 1 - exp(-TimeStep/collisionTime); %probability of particle being scattered


%initialize the particles
numParticles = 1000;
for i = 1:numParticles
    %init positions
    particleArray(i,1) = rand * RegionSizeX;
    particleArray(i,2) = rand * RegionSizeY;
    vth = randn*velocity_thermal + velocity_thermal;
    
    while(particleArray(i,1)>=(RegionSizeX/2-barrier_Width/2*1e-9) && particleArray(i,1)<=(RegionSizeX/2+barrier_Width/2*1e-9) && (particleArray(i,2)<=(barrier_Height*1e-9) || particleArray(i,2)>=(RegionSizeY-barrier_Height*1e-9)))
        %pick new location, you're in the boxes!
        particleArray(i,1) = rand * RegionSizeX;
        particleArray(i,2) = rand * RegionSizeY;
    end
    
    %init velocity
    particleArray(i,3) = vth; %thermal velocity magnitude
    particleArray(i,4) = ((rand * 2) - 1)*vth; %xVelocity
    particleArray(i,5) = sqrt(particleArray(i,3)^2 - particleArray(i,4)^2); %yVelocity
    if(rand > 0.5)
        particleArray(i,5) = particleArray(i,5) * -1;
    end
    
    particleArray(i,6) = 0; %set time since last scatter to 0
        
end

%run simulation
simLength = 1000;
for simCount = 1:simLength
    currTime = simCount * TimeStep;
        
    %scatter particles
    for i = 1:numParticles
        %update time since last scatter
        particleArray(i,6) = particleArray(i,6) + TimeStep;
        
        if(rand <= scatterProbability) %scatter the particle
            scatterPathSum = scatterPathSum + (particleArray(i,6)*particleArray(i,3)); %store path between scatters
            scatterTimeSum = scatterTimeSum + particleArray(i,6); %store time between scatters
            particleArray(i,6) = 0; %reset time since last scatter            
            numScatters = numScatters + 1;
            particleArray(i,3) = randn*velocity_thermal + velocity_thermal; %randomize velocity
            particleArray(i,4) = ((rand * 2) - 1)*particleArray(i,3); %xVelocity
            particleArray(i,5) = sqrt(particleArray(i,3)^2 - particleArray(i,4)^2); %yVelocity
            if(rand > 0.5)
                particleArray(i,5) = particleArray(i,5) * -1;
            end
        end
    end
    
    %update x velocity due to E-field    
    for eCount = 1:numParticles
        for yCount = 1:SizeY
            for xCount = 1:SizeX            
                if(particleArray(eCount,1) >= (xCount*xDel-xDel) && particleArray(eCount,1) < (xCount*xDel))
                    if(particleArray(eCount,2) >= (yCount*yDel-yDel) && particleArray(eCount,2) < (yCount*yDel))
                        particleArray(eCount,4) = particleArray(eCount,4) + TimeStep * -xVectors(yCount,xCount)*C.q_0/C.m_n/xDel;
                        particleArray(eCount,5) = particleArray(eCount,5) + TimeStep * -yVectors(yCount,xCount)*C.q_0/C.m_n/xDel; 
                    end                        
                end
            end                
        end
    end  
    
    %update particle positions
    new_xPos = particleArray(:,1) + TimeStep * particleArray(:,4);
    new_yPos = particleArray(:,2) + TimeStep * particleArray(:,5);    
    
    
    %check boundary conditions
    for i = 1:numParticles
        if(new_xPos(i) < 0) %pass through (x-dir)
            new_xPos(i) = new_xPos(i)+RegionSizeX;
        elseif(new_xPos(i) > RegionSizeX) %bounce off boundary
            new_xPos(i) = new_xPos(i) - RegionSizeX;
        end
        
        if(new_yPos(i) < 0) %bounce off boundary (y-dir)
            new_yPos(i) = abs(new_yPos(i));
            particleArray(i,5) = particleArray(i,5) * -1; %swap direction
        elseif(new_yPos(i) > RegionSizeY) %bounce off boundary
            new_yPos(i) = 2*RegionSizeY - new_yPos(i);
            particleArray(i,5) = particleArray(i,5) * -1; %swap direction
        end
       
        %You're in the boxes
        box_Left = RegionSizeX/2-barrier_Width/2*1e-9;
        box_Right = RegionSizeX/2+barrier_Width/2*1e-9;
        box_Bottom = barrier_Height*1e-9;
        box_Top = RegionSizeY-barrier_Height*1e-9;
        if(new_xPos(i) >= box_Left && new_xPos(i) <= box_Right && (new_yPos(i) <= box_Bottom || new_yPos(i) >= box_Top))
            if(particleArray(i,1) < box_Left)
                tempx = box_Left - abs(new_xPos(i) - particleArray(i,1));
                if(~(tempx >= box_Left && tempx <= box_Right && (new_yPos(i) <= box_Bottom || new_yPos(i) >= box_Top))) %not in box
                    new_xPos(i) = tempx;
                    particleArray(i,4) = particleArray(i,4) * -1; %swap direction
                end
            elseif(particleArray(i,1) > box_Right)
                tempx = box_Right + abs(new_xPos(i) - particleArray(i,1));
                if(~(tempx >= box_Left && tempx <= box_Right && (new_yPos(i) <= box_Bottom || new_yPos(i) >= box_Top))) %not in box
                    new_xPos(i) = tempx;
                    particleArray(i,4) = particleArray(i,4) * -1; %swap direction
                end
            elseif(particleArray(i,2) < box_Top)
                tempy = box_Top - abs(new_yPos(i) - particleArray(i,2));
                if(~(new_xPos(i) >= box_Left && new_xPos(i) <= box_Right && (tempy <= box_Bottom || tempy >= box_Top))) %not in box
                    new_yPos(i) = tempy;
                    particleArray(i,5) = particleArray(i,5) * -1; %swap direction
                end
            elseif(particleArray(i,2) > box_Bottom)
                tempy = box_Bottom + abs(new_yPos(i) - particleArray(i,2));
                if(~(new_xPos(i) >= box_Left && new_xPos(i) <= box_Right && (tempy <= box_Bottom || tempy >= box_Top))) %not in box
                    new_yPos(i) = tempy;
                    particleArray(i,5) = particleArray(i,5) * -1; %swap direction
                end
            end
        end
        
    end
    
    particleArray(:,1) = new_xPos;
    particleArray(:,2) = new_yPos;
    
    
    particleXPos7Array(1:7,simCount) = particleArray(1:7,1);
    particleYPos7Array(1:7,simCount) = particleArray(1:7,2);

    xCurrent_Array(simCount) = C.q_0 * electron_Concentration * mean(particleArray(:,4)) * RegionSizeX;


    if(simCount == simLength) %last iteration, populate Maps
        
        xDiv = RegionSizeX / xResolution;
        yDiv = RegionSizeY / yResolution;
        
        for yCount = 1:yResolution
            for xCount = 1:xResolution
                densityMap_Array(xCount,yCount) = 0;
                temperatureMap_Array(xCount,yCount) = 0;
                
                for eCount = 1:numParticles
                    if(particleArray(eCount,1) >= (xCount*xDiv-xDiv) && particleArray(eCount,1) < (xCount*xDiv))
                        if(particleArray(eCount,2) >= (yCount*yDiv-yDiv) && particleArray(eCount,2) < (yCount*yDiv))
                            densityMap_Array(xCount,yCount) = densityMap_Array(xCount,yCount) + 1;
                            temptemp = sqrt(particleArray(eCount,4).^2+particleArray(eCount,5).^2)*C.m_n/(3*C.kb);
                            temperatureMap_Array(xCount,yCount) = temperatureMap_Array(xCount,yCount) + temptemp;
                        end                        
                    end
                end                
            end
        end   
        temperatureMap_Array = temperatureMap_Array ./ densityMap_Array;
        temperatureMap_Array(isnan(temperatureMap_Array)) = 0;
    end
end


disp('Simulation Specs:');
disp(['Timestep Size: ', num2str(TimeStep), ' seconds']);
disp(['Number of Simulation Iterations: ', num2str(simLength)]);
disp(['Number of Particles: ', num2str(numParticles)]);
disp('----------');
disp('Simulation Results:');
disp(['Mean Free Path (MFP): ', num2str(scatterPathSum/numScatters), ' meters']);
disp(['Time between Collisions: ', num2str(scatterTimeSum/numScatters), ' seconds']);



figure(23);
hold on;
for i = 1:7
    plot(particleXPos7Array(i,:),particleYPos7Array(i,:));
end
xlim([0,RegionSizeX]);
ylim([0,RegionSizeY]);
title(['Particle Trajectory, Simulation Count: ', num2str(simCount), ', Timestep: ', num2str(TimeStep), '(TA 101064544)']);
xlabel('X (m)');
ylabel('Y (m)');

%plot the density map
figure(31);
[xMesh, yMesh] = meshgrid(linspace(0,RegionSizeX,xResolution),linspace(0,RegionSizeY,yResolution));
surf(xMesh,yMesh,transpose(densityMap_Array));
title('Density Map (TA 101064544)');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Count');

%Plot the X current over time
figure(32);
x = linspace(TimeStep,simLength*TimeStep,simLength);
plot(x,xCurrent_Array);
title('X Current over Time (TA 101064544)');
xlabel('Time (s)');
ylabel('Current (A)');


%plot the avg current vs bottleneck width
avgCurArr = [0.053111, 0.069134, 0.088485, 0.096499, 0.10608, 0.11388, 0.1217, 0.12317, 0.12871];
bneckWArr = [10, 20, 30, 40, 50, 60, 70, 80, 90];
figure(32);
plot(bneckWArr,avgCurArr);
title('Average Current vs Bottle-neck Width (TA 101064544)');
xlabel('Bottle-neck Width (nm)');
ylabel('Current (A)');


%plot the temperature map
figure(33);
[xMesh, yMesh] = meshgrid(linspace(0,RegionSizeX,xResolution),linspace(0,RegionSizeY,yResolution));
surf(xMesh,yMesh,transpose(temperatureMap_Array));
title('Temperature Map (TA 101064544)');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Temperature (K)');

avgCurrent = mean(xCurrent_Array(100:end));
BottleneckWidth = 100-barrier_Height*2;
disp(['Avg Curr: ', num2str(avgCurrent)]);
disp(['B-n width: ', num2str(BottleneckWidth)]);




