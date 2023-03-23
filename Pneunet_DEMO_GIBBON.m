%% DEMO_febio_00XX_pneunet_actuator_multi
% Below is a demonstration for:
% 
% * Building geometry for a parameterised pneunet actuator
% * Mesh based geomerty, followed by a scaling process
% * Three methods of modelling the strain limiting layer of the actuator
% * Defining the boundary conditions
% * Pneunet tip force analysis
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement,force and stress results

%% Keywords
%
% * febio_spec version 3.0
% * febio, FEBio
% * pressure loading
% * hexahedral elements, hex8
% * quad elements, quad4
% * hexahedral elements, hex8
% * pneunet actuator
% * soft robotic
% * static, solid
% * hyperelastic, Ogden
% * neo-Hookean
% * multiple steps
% * displacement logfile
% * stress logfile
% * tip force
% * automated geometry .stl output
% * automated mould .stl output
%%

clear; close all; clc;

%% Inputs:
%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
markerSize2=20;
linewidth=3;

%% Control parameters

%Mesh Tool
meshCreationTool=2;%
% == 1: Manual element seeded for every dimension of the pneunet
% == 2: Global desired element size dominates mesh control

%SLL method
SLL_control=2;
%==1: Adds a second material layer on the bottom of the pneunet of
%specified thickness (hexa8)
%==2: Adds a thin shell element the base of the structure (quad4). Second
%material.
%==3: Adds thin shell SLL elements at a the specified SLL thickness used in
%Method 1(quad4). Bottom layer is now the same material as the top and 
%material 2 is used for the shell.

%Chamber Wall Contact
contact_control=2;
%==1: Contact modelling inactive
%==2: Contact modelling active

%Fabrication file outputs
stl_control=2;
%==1: Body .stls (for full print of body and/or SLL) - dependant on chosen
%SLL method
%==2: Moulds (for moulding pneunet, as per SLL method 3)


%% File saving locations
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

stlPnuenetBody=fullfile(savePath,[febioFebFileNamePart,'_PneunetBody.stl']); %Log file name for exporting .stl
stlPnuenetSLL=fullfile(savePath,[febioFebFileNamePart,'_PneunetSLL.stl']); %Log file name for exporting .stl

stlTopAMould=fullfile(savePath,[febioFebFileNamePart,'_TopAMould.stl']); %Log file name for exporting .stl
stlTopBMould=fullfile(savePath,[febioFebFileNamePart,'_TopBMould.stl']); %Log file name for exporting .stl
stlBottomMould=fullfile(savePath,[febioFebFileNamePart,'_BottomMould.stl']); %Log file name for exporting .stl


%% FEA control settings
numTimeSteps=50; %Number of time steps desired
opt_iter=7; %Optimum number of iterations
max_refs=opt_iter*4; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
max_retries=10; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps)*2;%*4; %Maximum time step size

runMode='external';%'internal';


%% Contact Parameters 
contactPenalty=5;%5
laugon=0;
minaug=1;
maxaug=15;
fric_coeff=0;


%% Load Inputs
%Load
designPressureAngle=0.07; %(MPa)-target pressure for required bend
designPressureForce=0.1; % Maximum pressure - end of step 2


%% Material Properties

%Pneunet body
c1=1;%Shear-modulus-like parameter
m1=2;%Material parameter setting degree of non-linearity
k_factor=100; %Bulk modulus factor 
k1=c1*k_factor; %Bulk modulus


%Thick SLL material
c2=c1*50; %Shear-modulus-like parameter
m2=2; %Material parameter setting degree of non-linearity
k2=c2*k_factor; %Bulk modulus

%Paper - used for shell strain limiting layers
E_material2=6.5*(10^3);
Poissons_material2=0.2;



%% Pneunet Geometry Inputs:

n=7; %no. of chambers 
pointSpacing=2;%Only active if meshCreationTool == 2, designated aproximate size of desired cube element length and width if possible (mm) 

% X direction 

chamberlength=8; % internal chamber length (mm)
chamberWallThickness=1.2; % wall thickness expanding (mm)
gaplength=2.5; %(mm)

% y direction 

sidewallThickness=2;
chamberwidth=15; % internal chamber width (mm)
channelwidth=2; % width of the intenal channel (mm)

% Z direction 

SLLThickness=1; %Strain Limiting Layer (mm)% for Method 2, this could be made zero and Mat1 increased
SLLThicknessShell=0.1;

baseThickness=1; %material 1 base layer thickness (mm)
channelheight=3; %height of internal channel between chambers (mm)
channelRoofThickness=1;% (mm)
chamberRoofThickness=4; % (mm)
chamberheight=25; %internal height of the chamber (mm)


%% Desired Mesh Inputs
% If meshCreationTool == 1 input desired mesh geometry manually below

if meshCreationTool == 1
% X direction

chamberlengthElements=4; % no of mesh elements on this length
chamberWallThicknessElements=2; % no of mesh elements on this length
gaplengthElements=3; % no of mesh elements on this length

% y direction 

sidewallThicknessElements=1; % no of mesh elements on this length
chamberwidthElements=6; % no of mesh elements on this length (Inclusive of channel width)
channelwidthElements=2; % no of mesh elements on this length must even if two Chamber width is even, off if chamber width is odd

% Z direction 

SLLThicknessElements=1; %Elements used in SLL layer if method 1 applied

baseThicknessElements=1; %mat 1 basse layer thickness (mm)
channelheightElements=1; %height of internal channel between chambers (mm)
channelRoofThicknessElements=2;% (mm)
chamberRoofThicknessElements=2; % (mm)
chamberheightElements=8; % no of mesh elements on this length (Inclusive of channel parameters)

% Generalised mesh spacing calculations
elseif meshCreationTool == 2 % will auto compute based on spacing
    
chamberlengthElements=ceil(chamberlength/pointSpacing); % no of mesh elements on this length
chamberWallThicknessElements=ceil(chamberWallThickness/pointSpacing); % no of mesh elements on this length
gaplengthElements=ceil(gaplength/pointSpacing); % no of mesh elements on this length

% y direction 

sidewallThicknessElements=ceil(sidewallThickness/pointSpacing); % no of mesh elements on this length
chamberwidthElements=ceil(chamberwidth/pointSpacing); % no of mesh elements on this length (Inclusive of channel width)
channelwidthElements=ceil(channelwidth/pointSpacing); % no of mesh elements on this length must even if two Chamber width is even, off if chamber width is odd

if mod(chamberwidthElements,2) ~= mod(channelwidthElements,2)
    chamberwidthElements=chamberwidthElements+1;
end

% Z direction 

SLLThicknessElements=ceil(SLLThickness/pointSpacing); %Elements used in SLL layer if method 1 applied

baseThicknessElements=ceil(baseThickness/pointSpacing); %mat 1 basse layer thickness (mm)
channelheightElements=ceil(channelheight/pointSpacing); %height of internal channel between chambers (mm)
channelRoofThicknessElements=ceil(channelRoofThickness/pointSpacing);% (mm)
chamberRoofThicknessElements=ceil(chamberRoofThickness/pointSpacing); % (mm)
chamberheightElements=ceil(chamberheight/pointSpacing); % no of mesh elements on this length (Inclusive of channel parameters)
end 

if SLL_control==2
    SLLThickness=0;
    SLLThicknessElements=0;
end

%% Simplfied geometry creation
% Simple geometry is defined using mesh elements only

% X direction elements needed
length=(n*((2*chamberWallThicknessElements)+chamberlengthElements))+((n-1)*gaplengthElements);

% Y direction elements needed
width=(2*sidewallThicknessElements)+chamberwidthElements;

% Z direction elements needed
height=SLLThicknessElements+baseThicknessElements+chamberheightElements+chamberRoofThicknessElements;

boxDim=[length width height];
boxE1=[length width height];

[meshStruct]=hexMeshBox(boxDim,boxE1);

E_bar=meshStruct.E;
V_bar=meshStruct.V;
F_bar=meshStruct.F;
Fb_bar=meshStruct.Fb;
Cb_bar=meshStruct.faceBoundaryMarker;


V_bar(:,1)=V_bar(:,1)-min(V_bar(:,1)); %undoes centre alignment on X axis
V_bar(:,3)=V_bar(:,3)-min(V_bar(:,3)); %undoes centre alignment on Z axis Y alignment not changed yet as it is used to simplify scaling


VE_bar=patchCentre(E_bar,V_bar);%getting element centre for one value control of elements - better than 8 values for nodes
VE_bar=abs(VE_bar);% Abs values for easier control


%X axis logic - Inner & Outer 
lowerLimChannel=chamberWallThicknessElements;
upperLimChannel=((length)-chamberWallThicknessElements);

VE_x=VE_bar(:,1);
logicDeleteGapLengthX=zeros(1,size(VE_bar,1));
logicDeleteChamberX=zeros(1,size(VE_bar,1));
logicDeleteChannelX=zeros(1,size(VE_bar,1));

for i=1:1:size(VE_bar,1)
    
    for j=1:1:n
        %Setting Changing Limits
        LowerLimChamber=chamberWallThicknessElements + ((j-1)*((2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements));
        UpperLimChamber=chamberWallThicknessElements+chamberlengthElements  + ((j-1)*((2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements));
        
        LowerLimgaplength=((2*chamberWallThicknessElements)+chamberlengthElements)+((j-1)*((2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements));
        UpperLimgaplength=((j)*((2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements));
        
        %Outer geometry -  if needs removal -> mark as 1
        if VE_x(i) > LowerLimgaplength && VE_x(i) < UpperLimgaplength
            logicDeleteGapLengthX(i)=1;
        end
        
        if VE_x(i) > LowerLimChamber && VE_x(i) < UpperLimChamber
           logicDeleteChamberX(i)=1;
        end
        
        if VE_x(i) > lowerLimChannel && VE_x(i) < upperLimChannel
            logicDeleteChannelX(i)=1;
        end
        
    end
    
end
logicDeleteOuterX=logical(logicDeleteGapLengthX)';%converting to logical
logicDeleteInnerX1=logical(logicDeleteChamberX)';
logicDeleteInnerX2=logical(logicDeleteChannelX)';

% Y and Z logic -  Inner and Outer _ simpler as not effected by no. of
% chambers
logicDeleteOuterZ=VE_bar(:,3)>(SLLThicknessElements+baseThicknessElements+channelheightElements+channelRoofThicknessElements);
logicKeepOuter=~(logicDeleteOuterX.*logicDeleteOuterZ);
logicKeepOuter=logical(logicKeepOuter);

logicDeleteInnerY1=~(VE_bar(:,2)>(chamberwidthElements/2));
logicDeleteInnerY2=~(VE_bar(:,2)>(channelwidthElements/2));
logicDeleteInnerZ1=(VE_bar(:,3)>(SLLThicknessElements+baseThicknessElements) & VE_bar(:,3)<(height-chamberRoofThicknessElements));
logicDeleteInnerZ2=(VE_bar(:,3)>(SLLThicknessElements+baseThicknessElements) & VE_bar(:,3)<(SLLThicknessElements+baseThicknessElements+channelheightElements));

% uses logic vectors to choose which elements to keep
logicKeepChamber=~(logicDeleteInnerX1 == 1 & logicDeleteInnerY1 == 1 & logicDeleteInnerZ1 == 1);
logicKeepChannel=~(logicDeleteInnerX2 == 1 & logicDeleteInnerY2 == 1 & logicDeleteInnerZ2 == 1);
logicKeepInner=logicKeepChamber.*logicKeepChannel;

logicKeepInner=logicKeepInner(logicKeepOuter,:);
logicKeepInner=logical(logicKeepInner);

E_Body=E_bar(logicKeepOuter,:); %removes outer elements

F1=element2patch(E_Body);
[indBoundary1]=tesBoundary(F1);

F2=element2patch(E_Body(logicKeepInner,:));%removes internal elements
[indBoundary2]=tesBoundary(F2);


Fb=F2(indBoundary2,:);
Cb=7*ones(size(Fb,1),1);% Sample C matrix

for q=1:1:6 %Filling C matrix
    F_Cb1=Fb_bar(Cb_bar==q,:);
    logicNow=all(ismember(Fb,F_Cb1),2);
    Cb(logicNow)=q;
end

Cb(~any(ismember(Fb,F1(indBoundary1,:)),2))=0;%sets pressure faces to zero

% Remove unused nodes and clean up index matrices
 
[E,V,indFix2]=patchCleanUnused(E_Body(logicKeepInner,:),V_bar);
Fb=indFix2(Fb);

F=indFix2(F2);


cFigure; 
title('Unscaled simplified geometry','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',1);%0.5 transperancy normally hold on
hold on;% plotV(V,'k.','MarkerSize',markerSize/2);
axisGeom; 
colormap(turbo(250)); icolorbar; 
camlight headlight; 
gdrawnow; 

V_element=V;%copy of V based on element data only - won't be scaled
V_element(:,2)=V_element(:,2)-min(V_element(:,2));%undoing centre alignment in Z axis


%% Defining the boundary conditions
% The visualization of the model boundary shows colours. These labels can be used to define boundary conditions. 

%Define supported node sets
bcSupportList=unique(Fb(Cb==1,:)); %Node set part of selected face

%Get pressure faces
F_pressure=Fb(Cb==0,:); 

%Get end face
F_end=Fb(Cb==2,:);

%Get end face corner nodes
Spine_list=unique(Fb(Cb==5,:));%List of nodes along bottom surface of Pneunet

[LIA,indForceNodes]=ismember(unique(F_end),Spine_list);%unique index when both end conditions are met
indForceNodes=indForceNodes(indForceNodes~=0);%index list of Nodes where reaction force will be investigated

%% Defining Contact surfaces
logicContactSurf=Cb==7;%chamber wall faces
Normals=patchNormal(Fb,V);%gets normal vector of all facets

logicContactPosX=Normals(:,1)==1;%where norm v is pos x 
logicContactNegX=Normals(:,1)==-1;%where norm v is neg x 

logicContactAllPrimarySets=logical(logicContactSurf.*logicContactPosX);
logicContactAllSecondarySets=logical(logicContactSurf.*logicContactNegX);

% F_contactPrimary=Fb(logicContactPrimary,:);
% F_contactSecondary=Fb(logicContactSecondary,:);


Vm=patchCentre(Fb,V);% centre location of faces, where used here should not affect x axis as contact faces are parrallel, this is used to find the faces on the x ccordinate of interest
chamberAngleLeft={1 n-1};%empty for analysis of bending angle
chamberAngleRight={1 n-2};

if n>1
for q=1:1:n-1
%where the primary and secondary for each contact pair should be found
desiredPrimary=(2*chamberWallThicknessElements)+chamberlengthElements+((q-1)*((2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements));
desiredSecondary=(2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements+((q-1)*((2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements));

logicContactCoordinatePrimary=Vm(:,1)==desiredPrimary;
logicContactCoordinateSecondary=Vm(:,1)==desiredSecondary;

logicContactPrimary=logical(logicContactCoordinatePrimary.*logicContactAllPrimarySets);
logicContactSecondary=logical(logicContactCoordinateSecondary.*logicContactAllSecondarySets);


if q<n
chamberAngleLeft{q}=F(logicContactSecondary,:);
end
if q>1
chamberAngleRight{q-1}=F(logicContactPrimary,:);
end

ContactPair.Primary{q}=Fb(logicContactPrimary,:);
ContactPair.Secondary{q}=Fb(logicContactSecondary,:);

end

%display contact pairs

cFigure;
hold on;
title('Contact Pair Surfaces')
gpatch(Fb,V,'bw','k',0.25);
for q=1:1:n-1
    gpatch(ContactPair.Primary{q},V,'rw','k');
    gpatch(ContactPair.Secondary{q},V,'y','k');
end
axisGeom;

end


%% SLL Layer
if SLL_control == 1 % Thick SLL material

VE=patchCentre(E,V); %
logicSLL=VE(:,3)<SLLThicknessElements;
E_Body=E(~logicSLL,:);%main body
E_SLL=E(logicSLL,:);%SLL layer
E=[E_Body;E_SLL];

elseif SLL_control == 2 % Shell SLL material

logic_FSLL_height=Cb==5; % bottom face
FSLL=Fb(logic_FSLL_height,:);

E_Body=E;
E_SLL=FSLL;
E={};
E{1}=E_Body;
E{2}=E_SLL;

cFigure;
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; gpatch(E{2},V,'g'); % SLL shell elements
title('SLL shell elements')
axisGeom;

elseif SLL_control == 3 
FE=patchCentre(F,V); %   
logic_FSLL_height=FE(:,3)==SLLThicknessElements;
Normals_internal=patchNormal(F,V); % use normals to avoid takeing to face sets
logic_FSLL_direction=Normals_internal(:,3)==1;
logic_FSLL=logical(logic_FSLL_direction.*logic_FSLL_height);
FSLL=F(logic_FSLL,:);

E_Body=E;
E_SLL=FSLL;
E={};
E{1}=E_Body;
E{2}=E_SLL;

cFigure;
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; gpatch(E{2},V,'g');%SLL shell elements
title('SLL shell elements')
axisGeom;

end

%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,Cb,'none',0.35);%normally set transp to 0.5, 0.05 shows pressure face better
hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',15);
hl(2)=gpatch(F_pressure,V,'r','k',1);
hl(3)=plotV(V(indForceNodes,:),'r.','MarkerSize',15);
% hl(3)=gpatch(F2,V,'c','k',1);

patchNormPlot(F_pressure,V);
legend(hl,{'BC full support','Pressure surface', 'End Points'});%'SLLsurface'
legend('Location','northeastoutside')
axisGeom(gca,fontSize);
camlight headlight; %icolorbar;
gdrawnow; 


%% Scaling to match desired geometry
%Creating vectors of needed Coordinates
V_height=zeros(1,height); %creating empty vectors to reduce memory
V_length=zeros(1,length);
V_width=zeros(1,width);

%height scaling vector creation
current_height=0;
for i=1:1:height
    
    if i <= SLLThicknessElements
        V_height(i)=current_height+(SLLThickness*(1/SLLThicknessElements));
        
    elseif i <= baseThicknessElements+SLLThicknessElements
        V_height(i)=current_height+(baseThickness*(1/baseThicknessElements));
        
    elseif i <= channelheightElements+baseThicknessElements+SLLThicknessElements
        V_height(i)=current_height+(channelheight*(1/channelheightElements));
        
    elseif i <= channelRoofThicknessElements+channelheightElements+baseThicknessElements+SLLThicknessElements 
        V_height(i)=current_height+(channelRoofThickness*(1/channelRoofThicknessElements));
        
    elseif i <= chamberheightElements+baseThicknessElements+SLLThicknessElements
        V_height(i)=current_height+((chamberheight-channelheight-channelRoofThickness)*(1/(chamberheightElements-channelheightElements-channelRoofThicknessElements)));
        
    elseif  i <= chamberRoofThicknessElements+chamberheightElements+baseThicknessElements+SLLThicknessElements
        V_height(i)=current_height+(chamberRoofThickness*(1/chamberRoofThicknessElements));
        
    end
    current_height=V_height(i);
    
end
V_height=[0 V_height];

%width scaling vector creation
V(:,2)=V(:,2)-min(V(:,2));%undoes centre alignment on Y axis
current_width=0;
for i=1:1:width
    
    if i <= sidewallThicknessElements
        V_width(i)=current_width+(sidewallThickness*(1/sidewallThicknessElements));
    
    elseif i <= sidewallThicknessElements+((chamberwidthElements/2)-(channelwidthElements/2))
        V_width(i)=current_width+((chamberwidth-channelwidth)*(1/(chamberwidthElements-channelwidthElements)));
        
    elseif i <= sidewallThicknessElements+((chamberwidthElements/2)+(channelwidthElements/2))
        V_width(i)=current_width+((channelwidth)*(1/(channelwidthElements)));
        
    elseif i <= sidewallThicknessElements+chamberwidthElements
        V_width(i)=current_width+((chamberwidth-channelwidth)*(1/(chamberwidthElements-channelwidthElements)));
    
    elseif i<= (2*sidewallThicknessElements)+chamberwidthElements
        V_width(i)=current_width+(sidewallThickness*(1/sidewallThicknessElements));
    end
    current_width=V_width(i);
    
end
V_width=[0 V_width];

%length scaling vector creation
current_length=0;
for i=1:1:length
    itemp=i;
    total_chamberlengthElements=((2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements);
    chamber_count=floor(itemp/total_chamberlengthElements);
    itemp=itemp-(chamber_count*((2*chamberWallThicknessElements)+chamberlengthElements+gaplengthElements));
    
    
        if itemp == 0
            V_length(i)=current_length+(gaplength*(1/gaplengthElements));
            
        elseif itemp <= chamberWallThicknessElements
            V_length(i)=current_length+(chamberWallThickness*(1/chamberWallThicknessElements));

        elseif itemp <= chamberWallThicknessElements+chamberlengthElements
            V_length(i)=current_length+(chamberlength*(1/chamberlengthElements));

        elseif itemp <= (2*chamberWallThicknessElements)+chamberlengthElements
            V_length(i)=current_length+(chamberWallThickness*(1/chamberWallThicknessElements));

        elseif itemp <= gaplengthElements+(2*chamberWallThicknessElements)+chamberlengthElements
            V_length(i)=current_length+(gaplength*(1/gaplengthElements));

        end
    
    
    current_length=V_length(i);
end
V_length=[0 V_length];


for i=1:1:size(V,1) %applies true coordinate to element based system
    v_indexX=V(i,1)+1;%takes the element coordinate as an index to true value vectors
    v_indexY=V(i,2)+1;
    v_indexZ=V(i,3)+1;
    
    VXreal=V_length(v_indexX);
    VYreal=V_width(v_indexY);
    VZreal=V_height(v_indexZ);
    
    V(i,1)=VXreal;
    V(i,2)=VYreal;
    V(i,3)=VZreal;
    
end

cFigure; %Displays scaled geometry of two cells
title('Scaled geometry','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; plotV(V,'k.','MarkerSize',markerSize/2);
axisGeom; 
colormap(turbo(250)); icolorbar; 
% camlight headlight; 
gdrawnow; 


%% .stl file creation
% Pneunet body .stl files
if stl_control==1% if 3D printing .stls are needed
    [F_PneunetBody,~,~]=element2patch(E_Body,[],'hex8');
    Fb_PneunetBody=F_PneunetBody(tesBoundary(F_PneunetBody),:);%getting boundary of Pneunet(no SLL)
    Fb_PneunetBodyTri=quad2tri(Fb_PneunetBody,V);%changing to tri for .stl
    [Fb_PneunetBodyTri,V_stl_temp]=patchCleanUnused(Fb_PneunetBodyTri,V);%removing unused V to stop warning in command window
    
    TR_PnuenetBody=triangulation(Fb_PneunetBodyTri,V_stl_temp);
    stlwrite(TR_PnuenetBody,stlPnuenetBody);%writing .stl file
    
   if SLL_control==1 % create SLL .stl if required
       [F_PneunetSLL,~,~]=element2patch(E_SLL,[],'hex8');
        Fb_PneunetSLL=F_PneunetBody(tesBoundary(F_PneunetSLL),:);%getting boundary of SLL
        Fb_PneunetSLLTri=quad2tri(Fb_PneunetSLL,V);%changing to tri for .stl
        [Fb_PneunetSLLTri,V_stl_temp]=patchCleanUnused(Fb_PneunetSLLTri,V);%removing unused V to stop warning in command window
    
        TR_PnuenetSLL=triangulation(Fb_PneunetSLLTri,V_stl_temp);
        stlwrite(TR_PnuenetSLL,stlPnuenetSLL);%writing .stl file
       
       
   end

% Mould .stl files
elseif stl_control==2
    
%%% Bottom mould part %%%

%sample model to cut into
bottomDimEl=[3 3 2];
stlThickness=3;
[meshStruct]=hexMeshBox(bottomDimEl,bottomDimEl);

E_STLBase=meshStruct.E;
V_STLBase=meshStruct.V;
F_STLBase=meshStruct.F;
Fb_STLBase=meshStruct.Fb;
Cb_STLBase=meshStruct.faceBoundaryMarker;

%Getting centre for easier control of elements
VE_STLBase=patchCentre(E_STLBase,V_STLBase);

%Shifting data for easier deletion of sink elements
VE_STLBase(:,1)=abs(VE_STLBase(:,1));
VE_STLBase(:,2)=abs(VE_STLBase(:,2));
VE_STLBase(:,3)=VE_STLBase(:,3)-min(VE_STLBase(:,3))+0.5;

%Deleting sink elements for mould pouring
logicDeleteSTLBase= VE_STLBase(:,3)>1 & VE_STLBase(:,2) ==0 & VE_STLBase(:,1) == 0;
E_bottomMould=E_STLBase(~logicDeleteSTLBase,:);
F_bottomMould=element2patch(E_bottomMould,V_STLBase);
Fb_bottomMould=F_bottomMould(tesBoundary(F_bottomMould),:);

%scaling factors for the bottom mould
bottomMouldDimX=[0 stlThickness stlThickness+max(V(:,1)) 2*stlThickness+max(V(:,1))];
bottomMouldDimY=[0 stlThickness stlThickness+max(V(:,2)) 2*stlThickness+max(V(:,2))];
bottomMouldDimZ=[0 stlThickness stlThickness+SLLThickness+baseThickness];

%Moving axes
V_STLBase(:,1)=V_STLBase(:,1)-min(V_STLBase(:,1));
V_STLBase(:,2)=V_STLBase(:,2)-min(V_STLBase(:,2));
V_STLBase(:,3)=V_STLBase(:,3)-min(V_STLBase(:,3));

%scaling
for pp=1:1:size(V_STLBase,1)
    V_STLBase(pp,1)=bottomMouldDimX(V_STLBase(pp,1)+1);
    V_STLBase(pp,2)=bottomMouldDimY(V_STLBase(pp,2)+1);
    V_STLBase(pp,3)=bottomMouldDimZ(V_STLBase(pp,3)+1);
end

%plotting
cFigure;
gpatch(Fb_bottomMould,V_STLBase,'rw')
axisGeom;


%%% Top mould part A %%%

%sample model to cut into
bottomDimEl=[(4+n+n-1) 7 3];

[meshStruct]=hexMeshBox(bottomDimEl,bottomDimEl);

E_stlTopA=meshStruct.E;
V_stlTopA=meshStruct.V;
F_stlTopA=meshStruct.F;
Cb_stlTopA=meshStruct.faceBoundaryMarker;

%Getting centre for easier control of elements
VE_stlTopA=patchCentre(E_stlTopA,V_stlTopA);

%Shifting data for easier deletion of sink elements
V_stlTopA(:,1)=V_stlTopA(:,1)-min(V_stlTopA(:,1));
VE_stlTopA_temp=VE_stlTopA(:,1)-min(VE_stlTopA(:,1))+0.5;
VE_stlTopA(:,1)=abs(VE_stlTopA(:,1)); % only direction where changes are symmetric
VE_stlTopA(:,2)=abs(VE_stlTopA(:,2)); % only direction where changes are symmetric

VE_stlTopA(:,3)=VE_stlTopA(:,3)-min(VE_stlTopA(:,3))+0.5;

%Deleting general area  elements for mould pouring
logicDeletestlTopAXY1=VE_stlTopA(:,3)>1 & VE_stlTopA(:,2) ==2 & VE_stlTopA(:,1)<=max(VE_stlTopA(:,1))-1;
logicDeletestlTopAXY2=VE_stlTopA(:,3)>1 & VE_stlTopA(:,2) <=2 & VE_stlTopA(:,1)<=max(VE_stlTopA(:,1))-1 & VE_stlTopA(:,1)>max(VE_stlTopA(:,1))-2;
logicKeepstlTopASink=logical(~logicDeletestlTopAXY1.*~logicDeletestlTopAXY2);%general sink shape

%Deleting the elements required for channel
VE_stlTopA(:,1)=VE_stlTopA_temp;

DeleteElTopAXCoord=3.5:2:(max(V_stlTopA(:,1))-3.5); %element centres (in x) which must be deleted to make general shape

logicDeletestlTopChannelX=ismember(VE_stlTopA(:,1),DeleteElTopAXCoord);
logicDeletestlTopChannelYZ=VE_stlTopA(:,3)>1 & VE_stlTopA(:,2) <=2 & VE_stlTopA(:,2) > 0;
logicDeletestlTopChannelYZ2=VE_stlTopA(:,3)>2 & VE_stlTopA(:,2) ==0;%second logic for vertical grooves

logicDeletestlTopChannel=logical(logicDeletestlTopChannelX.*logicDeletestlTopChannelYZ);
logicDeletestlTopChannel2=logical(logicDeletestlTopChannelX.*logicDeletestlTopChannelYZ2);
logicKeepstlTopChannel=logical(~logicDeletestlTopChannel.*~logicDeletestlTopChannel2);
logicKeepstlTopA=logical(logicKeepstlTopASink.*logicKeepstlTopChannel);%combining all logic

%contolling faces and element identification
E_TopA=E_stlTopA(logicKeepstlTopA,:);
F_TopA=element2patch(E_TopA,V_stlTopA);
Fb_stlTopA=F_TopA(tesBoundary(F_TopA),:);

%scaling factors for the bottom mould
% bottomMouldDimX=[0 stl_gen_thickness stl_gen_thickness+max(V(:,1)) 2*stl_gen_thickness+max(V(:,1))];
TopADimY=[0 stlThickness ((2*stlThickness)+sidewallThickness) ((2*stlThickness)+sidewallThickness+((chamberwidth/2)-(channelwidth/2)))...
    ((2*stlThickness)+sidewallThickness+((chamberwidth/2)+(channelwidth/2))) ((2*stlThickness)+(sidewallThickness)+((chamberwidth))) ...
    ((3*stlThickness)+(2*sidewallThickness)+((chamberwidth))) ((4*stlThickness)+(2*sidewallThickness)+((chamberwidth)))];
TopADimZEdge=[0 stlThickness stlThickness*1.5 stlThickness*2];
TopADimZCentre=[0 stlThickness stlThickness+channelheight stlThickness+chamberheight];

%lenght scaling factor
lengthTopA=1:1:n+n-1; %general length vector - focused on the chamber geomertry
current_length=0;
for i=1:1:n+n-1 %scaling general length vector - focused on the chamber geomertry
   
    if iseven(i)==0
       
        lengthTopA(i)=current_length+chamberlength;
        
    elseif iseven(i)==1
        
        lengthTopA(i)=current_length+((2*chamberWallThickness)+gaplength);
        
    end
    current_length=lengthTopA(i);
end
lengthTopA=lengthTopA+(2*stlThickness)+chamberWallThickness;%shifting along to account for the lips that are added next
lengthTopA=[0 stlThickness (2*stlThickness)+chamberWallThickness lengthTopA (max(lengthTopA)+chamberWallThickness+stlThickness) (max(lengthTopA)+chamberWallThickness+(2*stlThickness))];

%Shifting V nodes
V_stlTopA(:,2)=V_stlTopA(:,2)-min(V_stlTopA(:,2));
V_stlTopA(:,3)=V_stlTopA(:,3)-min(V_stlTopA(:,3));

%applying scaling 
indEdgeTopA=find(V_stlTopA(:,2)<=1 | V_stlTopA(:,2)>=6 | V_stlTopA(:,1)<=1 | V_stlTopA(:,1)>=max(V_stlTopA(:,1)-1));

for pp=1:1:size(V_stlTopA,1)
    %  height (z)
    if ismember(pp,indEdgeTopA)==1
       
        V_stlTopA(pp,3)=TopADimZEdge(V_stlTopA(pp,3)+1);
        
    else 
        
        V_stlTopA(pp,3)=TopADimZCentre(V_stlTopA(pp,3)+1);
        
    end
    
    % x direction
    V_stlTopA(pp,1)=lengthTopA(V_stlTopA(pp,1)+1);
    
    % y direction
    
    V_stlTopA(pp,2)=TopADimY(V_stlTopA(pp,2)+1);
    
end

%plotting
cFigure;
gpatch(Fb_stlTopA,V_stlTopA,'rw')
axisGeom;

%%% Top Mould part B

%sample model to cut into
bottomDimEl=[(2+n+n-1) 3 2];

[meshStruct]=hexMeshBox(bottomDimEl,bottomDimEl);

E_stlTopB=meshStruct.E;
V_stlTopB=meshStruct.V;
Cb_stlTopB=meshStruct.faceBoundaryMarker;

%Getting centre for easier control of elements
VE_stlTopB=patchCentre(E_stlTopB,V_stlTopB);

%Shifting data for easier deletion of sink elements
V_stlTopB(:,1)=V_stlTopB(:,1)-min(V_stlTopB(:,1));

VE_stlTopB(:,1)=VE_stlTopB(:,1)-min(VE_stlTopB(:,1))+0.5;
VE_stlTopB(:,2)=abs(VE_stlTopB(:,2)); % 
VE_stlTopB(:,3)=VE_stlTopB(:,3)-min(VE_stlTopB(:,3))+0.5;

%Deleting general area  elements for mould pouring
% DeleteElTopBFill=1.5:2:(max(V_stlTopB(:,1))-1.5); %element centres (in x) which must be deleted to make general shape
DeleteElTopBWall=1.5:2:(max(V_stlTopB(:,1))-1.5); %element centres (in x) which must be deleted to make general shape

logicDeleteTopBTrench=VE_stlTopB(:,3)<1 & VE_stlTopB(:,2) <1 & VE_stlTopB(:,1) > 1 & VE_stlTopB(:,1) < max(V_stlTopB(:,1))-1; % creating trench
logicDeleteTopBWallX=ismember(VE_stlTopB(:,1),DeleteElTopBWall);%x direction for deleting scraper elements
logicDeleteTopBWallYZ=VE_stlTopB(:,3)>1 & VE_stlTopB(:,2) <1;%y z direction for deleting scraper elements
logicDeleteTopBWall=logical(logicDeleteTopBWallX.*logicDeleteTopBWallYZ); %combining
logicDeleteTopBWall=logical(logicDeleteTopBWall(~logicDeleteTopBTrench,:)); %accounting for other logic

%controlling faces and element identification & applying logic
E_stlTopB=E_stlTopB(~logicDeleteTopBTrench,:);
E_stlTopB=E_stlTopB(~logicDeleteTopBWall,:);

F_stlTopB=element2patch(E_stlTopB,V_stlTopB);
Fb_stlTopB=F_stlTopB(tesBoundary(F_stlTopB),:);

%scaling
TopBDimY=[0 stlThickness (stlThickness+(2*sidewallThickness)+chamberwidth) ((2*stlThickness)+(2*sidewallThickness)+chamberwidth)];
TopBDimZ=[0 (channelheight+channelRoofThickness) (chamberheight+chamberRoofThickness)];

lengthTopB=1:1:n+n-1;
current_length=0;
for i=1:1:n+n-1 %scaling general length vector - focused on the chamber geomertry
   
    if iseven(i)==0
       
        lengthTopB(i)=current_length+(2*chamberWallThickness)+chamberlength;
        
    elseif iseven(i)==1
        
        lengthTopB(i)=current_length+gaplength;
        
    end
    current_length=lengthTopB(i);
end
lengthTopB=lengthTopB+stlThickness;%shifting along to account for the lips that are added next
lengthTopB=[0 stlThickness lengthTopB (max(lengthTopB)+stlThickness)];

%shifting
V_stlTopB(:,2)=V_stlTopB(:,2)-min(V_stlTopB(:,2));
V_stlTopB(:,3)=V_stlTopB(:,3)-min(V_stlTopB(:,3));


%assigning scaled values
for pp=1:1:size(V_stlTopB,1)
    %  height (z)
    V_stlTopB(pp,3)=TopBDimZ(V_stlTopB(pp,3)+1);
    
    % x direction
    V_stlTopB(pp,1)=lengthTopB(V_stlTopB(pp,1)+1);
    
    % y direction
    
    V_stlTopB(pp,2)=TopBDimY(V_stlTopB(pp,2)+1);
    
end

%plotting
cFigure;
gpatch(Fb_stlTopB,V_stlTopB,'rw')
axisGeom;


%%% .stl file creation
% Bottom
Fb_bottomMould=quad2tri(Fb_bottomMould,V_STLBase);%changing to tri for .stl
[Fb_bottomMould,V_STLBase]=patchCleanUnused(Fb_bottomMould,V_STLBase);%removing unused V to stop warning in command window
TR_TopA=triangulation(Fb_bottomMould,V_STLBase);
stlwrite(TR_TopA,stlBottomMould);%writing .stl file

% Top A
Fb_stlTopA=quad2tri(Fb_stlTopA,V_stlTopA);%changing to tri for .stl
[Fb_stlTopA,V_stlTopA]=patchCleanUnused(Fb_stlTopA,V_stlTopA);%removing unused V to stop warning in command window
TR_TopA=triangulation(Fb_stlTopA,V_stlTopA);
stlwrite(TR_TopA,stlTopAMould);%writing .stl file

% Top B
Fb_stlTopB=quad2tri(Fb_stlTopB,V_stlTopB);%changing to tri for .stl
[Fb_stlTopB,V_stlTopB]=patchCleanUnused(Fb_stlTopB,V_stlTopB);%removing unused V to stop warning in command window
TR_TopB=triangulation(Fb_stlTopB,V_stlTopB);
stlwrite(TR_TopB,stlTopBMould);%writing .stl file
       
end


%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='3.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Create control structure for use by all steps
stepStruct.Control.analysis='STATIC';
stepStruct.Control.time_steps=numTimeSteps;
stepStruct.Control.step_size=1/numTimeSteps;
stepStruct.Control.solver.max_refs=max_refs;
stepStruct.Control.solver.max_ups=max_ups;
stepStruct.Control.solver.symmetric_stiffness=0;
stepStruct.Control.time_stepper.dtmin=dtmin;
stepStruct.Control.time_stepper.dtmax=dtmax; 
stepStruct.Control.time_stepper.max_retries=max_retries;
stepStruct.Control.time_stepper.opt_iter=opt_iter;

%Add template based default settings to proposed control section
[stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

febio_spec.Step.step{1}.Control=stepStruct.Control;
febio_spec.Step.step{1}.ATTR.id=1; %inflate to bending pressure
febio_spec.Step.step{2}.Control=stepStruct.Control;
febio_spec.Step.step{2}.ATTR.id=2; %inflate further (used to fix tip)

%Material section
materialName1='Material1'; % pneunet body
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k1;

if SLL_control == 1

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.c1=c2;
febio_spec.Material.material{2}.m1=m2;
febio_spec.Material.material{2}.c2=c2;
febio_spec.Material.material{2}.m2=-m2;
febio_spec.Material.material{2}.k=k2;

elseif SLL_control >= 2
    
materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='neo-Hookean';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.E=E_material2;
febio_spec.Material.material{2}.v=Poissons_material2;
end

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates 

% -> Elements

partBodyPneunet='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partBodyPneunet; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E_Body,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E_Body; %The element matrix

if SLL_control == 1 %Thick SLL -> hexa8

    partSLL='Part2';
    febio_spec.Mesh.Elements{2}.ATTR.name=partSLL; %Name of this part
    febio_spec.Mesh.Elements{2}.ATTR.type='hex8'; %Element type 
    febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E_Body,1)+(1:1:size(E_SLL,1))'; %Element id's
    febio_spec.Mesh.Elements{2}.elem.VAL=E_SLL; %The element matrix

elseif SLL_control >= 2 %Thin SLL -> quad4

    partSLL='Part2';
    febio_spec.Mesh.Elements{2}.ATTR.name=partSLL; %Name of this part
    febio_spec.Mesh.Elements{2}.ATTR.type='quad4'; %Element type 
    febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E_Body,1)+(1:1:size(E_SLL,1))'; %Element id's
    febio_spec.Mesh.Elements{2}.elem.VAL=E_SLL; %The element matrix

end


% -> Surfaces
surfaceName1='LoadedSurface';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure;

leaderStringSurf='contactSurface';
leaderStringContact='Contact';

if contact_control == 2 % contact surfaces
    contactSurfacesValPrimary=2:2:(2*(n-1));% attribute number of contact faces
    contactSurfacesValSecondary=3:2:((2*(n-1))+1);% attribute number of contact faces

    if n>1 %only need face contact for more than one chamber
        for q=1:1:n-1
            numStringPrimary=num2str(contactSurfacesValPrimary(q));
            numStringSecondary=num2str(contactSurfacesValSecondary(q));
            numStringPair=num2str(q);

            contactSurfacesStringPrimary=strcat(leaderStringSurf,numStringPrimary);%creates string 'ContactSurfaceINT'
            febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.ATTR.name=contactSurfacesStringPrimary;
            febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.quad4.ATTR.id=(1:1:size(ContactPair.Primary{q},1))';%size(F_pressure,1)+
            febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.quad4.VAL=ContactPair.Primary{q};

            contactSurfacesStringSecondary=strcat(leaderStringSurf,numStringSecondary);%creates string 'ContactSurfaceINT'
            febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.ATTR.name=contactSurfacesStringSecondary;
            febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.quad4.ATTR.id=(1:1:size(ContactPair.Secondary{q},1))';%size(F_pressure,1)+size(F_contactPrimary,1)+
            febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.quad4.VAL=ContactPair.Secondary{q};

            % -> Surface pairs
            febio_spec.Mesh.SurfacePair{q}.ATTR.name=strcat(leaderStringContact,numStringPair);
            febio_spec.Mesh.SurfacePair{q}.primary=contactSurfacesStringPrimary;
            febio_spec.Mesh.SurfacePair{q}.secondary=contactSurfacesStringSecondary;
        end
    end

end


% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

nodeSetName2='bcTipList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.node.ATTR.id=indForceNodes(:);


%MeshDomains section

febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partBodyPneunet;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

if SLL_control == 1 %SLL domain

    febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partSLL;
    febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;
elseif SLL_control >= 2

    febio_spec.MeshDomains.ShellDomain{1}.ATTR.name=partSLL;
    febio_spec.MeshDomains.ShellDomain{1}.ATTR.mat=materialName2;
end


%MeshData secion
%-> Element data
if SLL_control >= 2
    febio_spec.MeshData.ElementData{1}.ATTR.var='shell thickness';
    febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partSLL;
    febio_spec.MeshData.ElementData{1}.elem.ATTR.lid=(1:1:size(E_SLL,1))';
    febio_spec.MeshData.ElementData{1}.elem.VAL=SLLThicknessShell*ones(size(E_SLL,1),size(E_SLL,2));
end

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.type='fix';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dofs='x,y,z';

febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.type='prescribe';
febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.node_set=nodeSetName2;
febio_spec.Step.step{2}.Boundary.bc{1}.dof='z';
febio_spec.Step.step{2}.Boundary.bc{1}.scale.ATTR.lc=2;
febio_spec.Step.step{2}.Boundary.bc{1}.scale.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{1}.relative=1;

%Loads section
% -> Surface load
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=designPressureForce;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;


%Contact section
q=0;
if contact_control == 2
    for q=1:1:n-1
        febio_spec.Contact.contact{q}.ATTR.type='sliding-elastic';
        febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{q}.ATTR.name;
        febio_spec.Contact.contact{q}.two_pass=1;
        febio_spec.Contact.contact{q}.laugon=laugon;
        febio_spec.Contact.contact{q}.tolerance=0.2;
        febio_spec.Contact.contact{q}.gaptol=0;
        febio_spec.Contact.contact{q}.minaug=minaug;
        febio_spec.Contact.contact{q}.maxaug=maxaug;
        febio_spec.Contact.contact{q}.search_tol=0.01;
        febio_spec.Contact.contact{q}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));
        febio_spec.Contact.contact{q}.symmetric_stiffness=0;
        febio_spec.Contact.contact{q}.auto_penalty=1;
        febio_spec.Contact.contact{q}.penalty=contactPenalty;
        febio_spec.Contact.contact{q}.fric_coeff=fric_coeff;
    end

end


%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 (designPressureAngle/designPressureForce); 2 1];

febio_spec.LoadData.load_controller{2}.ATTR.id=2;%loadcurve ID no.
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='STEP'; 
febio_spec.LoadData.load_controller{2}.points.point.VAL=[1 0; 2 1];


%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
if SLL_control==1
    febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);
elseif SLL_control>=2
    febio_spec.Output.logfile.element_data{1}.VAL=1:(size(E{1},1)+size(E{2},1));        
end    
%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% |febView(febio_spec); %Viewing the febio file|

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
% febView(febioFebFileName); 

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.runMode=runMode;
febioAnalysis.maxLogCheckTime=100; %Max log file checking time - EDITED 

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
     %% 
    % Importing nodal displacements from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    
    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    timeVec=dataStruct.time; %Time
    
    %Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
               
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
        
    % Create basic view and store graphics handle to initiate animation
    
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
%     hp.Marker='.';
%     hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude
                
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;      
    %%
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
    
    %Access data
    E_stress_mat=dataStruct.data;
    
    E_stress_mat(isnan(E_stress_mat))=0;
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
     if SLL_control==1
    [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,end));
     else
    [CV]=faceToVertexMeasure(E_Body,V,E_stress_mat(:,:,end));
     end
    % Create basic view and store graphics handle to initiate animation

    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{1}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
%     hp.Marker='.';
%     hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([min(E_stress_mat(:)) max(E_stress_mat(:))]);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        if SLL_control==1
            [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,qt));
        else
            [CV]=faceToVertexMeasure(E_Body,V,E_stress_mat(:,:,qt));
        end    
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
    


%% Bending Angle of Each Chamber
theta=zeros(1,n-2);


if n>4 % results only comparable with 5+ chambers as there are no chambers far enough from an end of the Pneunet
    
    for i=1:1:n-2 
    normLeft=patchNormal(ContactPair.Secondary{i},V_DEF(:,:,size(V_DEF,3)));%normals of left face of chamber wall
    normRight=patchNormal(ContactPair.Primary{i+1},V_DEF(:,:,size(V_DEF,3)));%normals of right face of chamber wall
    
    normLeft=mean(normLeft,1);%average normal for left side
    normRight=mean(normRight,1);%average normal for right side
    
    theta(i)=rad2deg(acos((dot(normLeft,normRight,2))/((norm(normLeft)*(norm(normRight))))));% angle between chamber wall normals
    
    
    
    end
    
    cFigure;
    gpatch(Fb,V_DEF(:,:,size(V_DEF,3)),Cb,'k',0.2); hold on;
    title('Contact surface normal vectors')
    for i=1:1:n-2 
        patchNormPlot(ContactPair.Primary{i+1},V_DEF(:,:,size(V_DEF,3)));
        patchNormPlot(ContactPair.Secondary{i},V_DEF(:,:,size(V_DEF,3))); 
    end
    axisGeom;
    %disp(theta)%display chamber angles to command window
    
    
end


%% Bending Angle using Tip Nodes
bendingAngle=zeros(size(V_DEF,3),3);%empty vector of bending angles over time
bendingAngle(:,3)=timeVec;%assigning time values
YZNormVec=[1 0 0];

bendingAngleRefPoint=mean(V(indForceNodes,:),1);

for i=1:1:size(V_DEF,3)
    
    normEnd=mean(patchNormal(F_end,V_DEF(:,:,i)),1);%avergage normal vector of the last face
    nodeEnd=mean(V_DEF(indForceNodes,:,i),1);%averaged end Node tracking over time
    
    bendingAngle(i,1)=rad2deg(acos((dot(normEnd,YZNormVec,2))/((norm(normEnd)*(norm(YZNormVec))))));%angle of end face normal vector using 2D dot product
    bendingAngle(i,2)=rad2deg(atan((abs(nodeEnd(3)))/(nodeEnd(1))));%saving angle of end node
    
end

cFigure;
ha(1)=plot(bendingAngle(:,3),bendingAngle(:,1),'r');
hold on;
ha(2)=plot(bendingAngle(:,3),bendingAngle(:,2),'b');
legend(ha,{'Normal Vector Angle', 'End Point Angle'});%'SLLsurface'
title('Bending Angle');

%% Force

% Importing nodal displacements from a log file
dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),1,1);

%Access data
forceData=dataStruct.data; %All force values
timeVec=dataStruct.time; %Time

% Extract Forces at Restricted nodes
forceMean=zeros(1,size(forceData,3));
forceTotal=zeros(1,size(forceData,3));

for i=1:1:size(forceData,3)
    Force_nodal=forceData(indForceNodes,3,i); %Force only at restricted nodes
    forceMean(i)=mean(Force_nodal);%mean force at the end nodes due to Z restriction
    forceTotal(i)=sum(Force_nodal);% sum of all the end node forces due to Z restriction
end

%Plotting
cFigure;
hforce=plot(timeVec,forceTotal,'b');
title('Nodal Force (N)');
xlabel('Time (s)'); ylabel('Average End-Node Force ()');

end % successful run check


%% 
cFigure;
gpatch(Fb,V,Cb)
axisGeom;

