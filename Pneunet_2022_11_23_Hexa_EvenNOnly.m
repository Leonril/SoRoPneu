%% Replicating previous models using Hexa-hedral elemetents
% this allows:
%
% * More complex geometry
% * Simpler material handling 
%
% This route was not used previously due to
%
% * Awkward shape generation


%% Keywords
%
% * febio_spec version 3.0
% * febio, FEBio
% * pressure loading
% * hexahedral elements, hex8
% * pneunet actuator
% * soft robotic
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
markerSize2=20;
lineWidth=3;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force


% FEA control settings
numTimeSteps=25; %Number of time steps desired
opt_iter=25; %Optimum number of iterations
max_refs=opt_iter*2; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps)*4; %Maximum time step size

runMode='internal';%'internal';

%% Load Inputs
%Load
appliedPressure1=0.15; 
appliedPressure2=appliedPressure1/30; 

%Define applied displacement perturbation
prescribedDisplacement_X=2;

%% Material Properties
%Material parameter set
c1=1; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=100; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

c2=c1*2; %Shear-modulus-like parameter
m2=2; %Material parameter setting degree of non-linearity
k2=c2*k_factor; %Bulk modulus

%% Geometry Inputs:

n=4; %no. of chambers 

% X direction 

Chamber_length=8; % internal chamber length (mm)
Chamber_wt=0.5; % wall thickness expanding (mm)
Gap_length=1; %(mm)

% y direction 

Side_wall_thickness=1;
Chamber_width=8; % internal chamber width (mm)
Channel_width=2; % width of the intenal channel (mm)

% Z direction 

SLL_thickness=1; %Strain Limiting Layer (mm)
Mat1_base_thickness=0.5; %mat 1 basse layer thickness (mm)
Channel_height=1; %height of internal channel between chambers (mm)
Channel_roof_thickness=1;% (mm)
Chamber_roof_thickness=2.5; % (mm)
Chamber_height=10; %internal height of the chamber (mm)

%% Geometry calculations
% X direction
Total_length=((n)*((2*Chamber_wt)+Chamber_length))+((n-1)*(Gap_length));
StepX=0.5;

% Y direction
Total_width=(2*Side_wall_thickness)+Chamber_width;
StepY=1;

% Z direction
Total_height=SLL_thickness+Mat1_base_thickness+Chamber_height+Chamber_roof_thickness;
StepZ=0.5;

numElementsLength=Total_length/StepX;
numElementsWidth=Total_width/StepY;
numElementsHeight=Total_height/StepZ;


boxDim=[Total_length Total_width Total_height];
boxE1=[numElementsLength numElementsWidth numElementsHeight];

[meshStruct]=hexMeshBox(boxDim,boxE1);

E_bar=meshStruct.E;
V_bar=meshStruct.V;
F_bar=meshStruct.F;
Fb_bar=meshStruct.Fb;
Cb_bar=meshStruct.faceBoundaryMarker;


cFigure;
gpatch(Fb_bar,V_bar,Cb_bar);
axisGeom;

VE_bar=patchCentre(E_bar,V_bar); 


%% HEAD MELTER: OUTER GEOMETRY
CX1=VE_bar(:,1);
CX=abs(CX1);
CZCh=VE_bar(:,3);
CZCh=CZCh-min(CZCh)+(StepZ/2);
imax=size(CX,1);
upperlimZCh=SLL_thickness+Mat1_base_thickness+Channel_height+Channel_roof_thickness;

evenTest=mod(n,2);
if evenTest == 0
    jmax=n+1;
    for i=1:1:imax
    CXtemp=0;
    CZtemp=0;
    lowerlimX=0;
    upperlimX=Gap_length/2;
        for j=1:1:jmax
            if CX(i) >= lowerlimX && CX(i) <= upperlimX
                CXtemp=1;
            elseif CZCh(i)<= upperlimZCh
                CZtemp=1;
            end
            lowerlimX=(0.5*Gap_length)+((j-1)*(Gap_length))+(j*((2*Chamber_wt)+Chamber_length));
            upperlimX=lowerlimX+Gap_length;
        end
        if CXtemp == 1
            CX(i)=0;
        end
        if CZtemp==1
            CZCh(i)=0;
        end
    end
end

%% HEAD MELTER: INTERNAL SURFACE
CZI=VE_bar(:,3);
CZI=CZI-min(CZI)+(StepZ/2);
CZIChannel=CZI;
CXI=VE_bar(:,1);
CXI=abs(CXI);
CXIChannel=CXI;
CYI=VE_bar(:,2);
CYI=abs(CYI);
CYIChannel=CYI;

%Static Limits
UpperLimZI_Chamb=SLL_thickness+Mat1_base_thickness+Chamber_height;
UpperLimZI_Chann=SLL_thickness+Mat1_base_thickness+Channel_height;
LowerLimZI=SLL_thickness+Mat1_base_thickness;
LimY=(Total_width/2)-Side_wall_thickness;
lowerlimX_Chann=0;
upperlimX_Chann=(Total_length/2)-Chamber_wt;


if evenTest == 0
    for i=1:1:imax
    CXtemp=0;
    CYtemp=0;
    CZtemp=0;
    C_Channel_temp=0;
    lowerlimX=(Gap_length/2)+Chamber_wt;
    upperlimX=(Gap_length/2)+Chamber_wt+Chamber_length;
        for j=1:1:jmax
            if CXI(i) >= lowerlimX && CXI(i) <= upperlimX
                CXtemp=1;
            elseif CXI(i) >= lowerlimX_Chann && CXI(i) <= upperlimX_Chann
                C_Channel_temp=1;
            end
            lowerlimX=lowerlimX+(j*(Gap_length+(2*Chamber_wt)+Chamber_length));
            upperlimX=upperlimX+(j*(Gap_length+(2*Chamber_wt)+Chamber_length));
        end
        if CZI(i) >= LowerLimZI && CZI(i) <= UpperLimZI_Chamb
            CZI(i)=0;
        end
        if CZIChannel(i) >= LowerLimZI && CZIChannel(i) <= UpperLimZI_Chann
            CZIChannel(i)=0;
        end
        if CYI(i) >= 0 && CYI(i) <= LimY
            CYI(i)=0;
        end
        if CYIChannel(i) >= 0 && CYIChannel(i) <= Channel_width
            CYIChannel(i)=0;
        end
        if CXtemp == 1
            CXI(i)=0;
        end
        if C_Channel_temp == 1
            CXIChannel(i)=0;
        end
    end
end




%%

CW=VE_bar(:,3);
CW=CW-min(CW); 
CW=CW./max(CW); 
CW=round((CW.*(numElementsHeight-1)))+1;


% numElementsPeriod=(Chamber_length+(2*Chamber_wt))/StepX;
% CD=rem(CZ,numElementsPeriod); 

logicKeep1=~(CX==0 & CZCh>0);


logicKeep3=~(CXI ==0 & CYI==0 & CZI==0  ); %Chamber logic
logicKeep4=~(CXIChannel==0 & CYIChannel==0 & CZIChannel==0); %Channel logic

logicKeep2=logicKeep3.*logicKeep4;
logicKeep2=logicKeep2(logicKeep1,:);
logicKeep2=logical(logicKeep2);


E1=E_bar(logicKeep1,:);
F1=element2patch(E1);
[indBoundary1]=tesBoundary(F1);

F2=element2patch(E1(logicKeep2,:));
[indBoundary2]=tesBoundary(F2);

cFigure
gpatch(E1,V_bar)
axisGeom


% logicKeep2=logicKeep3*logicKeep4;
% F2=element2patch(E1(logicKeep2,:));
% [indBoundary2]=tesBoundary(F2);

Fb=F2(indBoundary2,:);
Cb=7*ones(size(Fb,1),1);

for q=1:1:6
    F_Cb1=Fb_bar(Cb_bar==q,:);
    logicNow=all(ismember(Fb,F_Cb1),2);
    Cb(logicNow)=q;
end

Cb(~any(ismember(Fb,F1(indBoundary1,:)),2))=0;

% Remove unused nodes and clean up index matrices
 
[E,V,indFix2]=patchCleanUnused(E1(logicKeep2,:),V_bar);
Fb=indFix2(Fb);

F=indFix2(F2);

%%

cFigure; 
gpatch(Fb,V,Cb,'k',0.5);
axisGeom; 
colormap(turbo(250)); icolorbar; 
camlight headlight; 
gdrawnow; 



