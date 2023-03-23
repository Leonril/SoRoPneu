%% Revisions
%% Rev I
% * 
% * Notes:
% * Adding stl output
% * 
% * 

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

%% Inputs:
%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
markerSize2=20;
lineWidth=3;

%% Control parameters

%Mesh Tool
MeshCreationTool=2;%
% == 1: Uses 'must points' at key variable locations
% == 2: Used designated 'mesh size' value to add more homogentity to the
% mesh - this should reduce some convergence issues

%SLL method
SLL_control=3;
%==1: Adds a second material layer on the bottom of the pneunet of
%specified thickness
%==2: Adds a thin shell element the base of the structure
%==3: Adds thin shell SLL elements at in between the top and bottom layers
%created in method one. Bottom layer is now the sane material as the top
%and material 2 is the shell.

%Object type
object_control=0;
%==1: Pneunet fixed in Z axis at time=1. Least computationally expensive
%==2: Rigid body. At a specificed Z coordinate or plane, a ridid face is added to restrict bending and simulate a rigid object.
%==3: Soft object: A 3D soft body is modelled. Most computationally expensive.
softObject_shape=2;
%==1: Half of a sphere is modelled at a specified height.
%==2: Half of a cube is modelled, at specified height. 

%Chamber Wall Contact
contact_control=2;
%==1: Contact modelling inactive
%==2: Contact modelling active

%Plots
plot_control=2;
%==1: Plots on
%==2: Plots off
plot_vector=[7 8 10 11];% list of plots to be active regardless of plot control - see each plot's code for identifiers

%Fabrication file outputs
stl_control=0;
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

STL_Pnuenet_Body=fullfile(savePath,[febioFebFileNamePart,'_PneunetBody.stl']); %Log file name for exporting force
STL_Pnuenet_SLL=fullfile(savePath,[febioFebFileNamePart,'_PneunetSLL.stl']); %Log file name for exporting force

STL_TopMould=fullfile(savePath,[febioFebFileNamePart,'_TopMould.stl']); %Log file name for exporting force
STL_BottomMould=fullfile(savePath,[febioFebFileNamePart,'_BottomMould.stl']); %Log file name for exporting force
%% FEA control settings
numTimeSteps=50; %Number of time steps desired
opt_iter=25; %Optimum number of iterations
max_refs=opt_iter*4; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
max_retries=10; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps)*4; %Maximum time step size

runMode='external';%'internal';

%% Contact Parameters 
contactPenalty=5;%5
laugon=0;
minaug=1;
maxaug=15;
fric_coeff=0;


%% Load Inputs
%Load
appliedPressure=0.1; %0.15

%% Pneunet Material Properties
%Material_Bank - takes from bank of previous materials
materialBank=Material_Bank;


%Material parameter set

% %Elastosil
% c1=75*(10^(-3)); %Shear-modulus-like parameter
% m1=2.749; %Material parameter setting degree of non-linearity

c1=1;
m1=2;
k_factor=100; %Bulk modulus factor 
k1=c1*k_factor; %Bulk modulus

%Paper
E_material2=6.5*(10^3);
Poissons_material2=0.2;

%Backup sample material 2
c2=c1*50; %Shear-modulus-like parameter
m2=2; %Material parameter setting degree of non-linearity
k2=c2*k_factor; %Bulk modulus



%% Pneunet Geometry Inputs:

n=4; %no. of chambers 
pointSpacing=1;%Only active if MeshCreationTool == 2, designated aproximate size of desired cube element length and width if possible (mm) 

% X direction 

Chamber_length=8; % internal chamber length (mm)
Chamber_wt=1.2; % wall thickness expanding (mm)
Gap_length=2.5; %(mm)

% y direction 

Side_wall_thickness=2;
Chamber_width=15; % internal chamber width (mm)
Channel_width=2; % width of the intenal channel (mm)

% Z direction 

SLL_thickness=1; %Strain Limiting Layer (mm)% for Method 2, this could be made zero and Mat1 increased
SLLThicknessShell=0.1;

Mat1_base_thickness=1; %mat 1 basse layer thickness (mm)
Channel_height=3; %height of internal channel between chambers (mm)
Channel_roof_thickness=1;% (mm)
Chamber_roof_thickness=4; % (mm)
Chamber_height=25; %internal height of the chamber (mm)

%% Object controls

% Soft object controls
ObjectCentreDist=[0 0 -10]; % distance from Pneuent tip to object centre (mm)

% Rigid plane distance from Pneunet bottom (Z) (for object_control==1)
rigidLength=30;
rigidWidth=20;

% Soft object material properties
Object_properties=materialBank.SoftMaterialSample1;


% Hemisphere shape controls
sphereRadius=15; % (mm)


% Half cuboid shape controls
SoftObject_PointSpacing=1; %mesh refinement parameter
SoftObjectLength=15;
SoftObjectWidth=30;
SoftObjectHeight=15; %from midplane - so half of true body height  

%% Desired Mesh Inputs
% If MeshCreationTool == 1 input desired mesh geometry manually below

if MeshCreationTool == 1
% X direction

Chamber_length_elements=4; % no of mesh elements on this length
Chamber_wt_elements=2; % no of mesh elements on this length
Gap_length_elements=3; % no of mesh elements on this length

% y direction 

Side_wall_thickness_elements=1; % no of mesh elements on this length
Chamber_width_elements=6; % no of mesh elements on this length (Inclusive of channel width)
Channel_width_elements=2; % no of mesh elements on this length must even if two Chamber width is even, off if chamber width is odd

% Z direction 

SLL_thickness_elements=1; %Elements used in SLL layer if method 1 applied

Mat1_base_thickness_elements=1; %mat 1 basse layer thickness (mm)
Channel_height_elements=1; %height of internal channel between chambers (mm)
Channel_roof_thickness_elements=2;% (mm)
Chamber_roof_thickness_elements=2; % (mm)
Chamber_height_elements=8; % no of mesh elements on this length (Inclusive of channel parameters)

% Generalised mesh spacing calculations
elseif MeshCreationTool == 2
    
Chamber_length_elements=ceil(Chamber_length/pointSpacing); % no of mesh elements on this length
Chamber_wt_elements=ceil(Chamber_wt/pointSpacing); % no of mesh elements on this length
Gap_length_elements=ceil(Gap_length/pointSpacing); % no of mesh elements on this length

% y direction 

Side_wall_thickness_elements=ceil(Side_wall_thickness/pointSpacing); % no of mesh elements on this length
Chamber_width_elements=ceil(Chamber_width/pointSpacing); % no of mesh elements on this length (Inclusive of channel width)
Channel_width_elements=ceil(Channel_width/pointSpacing); % no of mesh elements on this length must even if two Chamber width is even, off if chamber width is odd

if mod(Chamber_width_elements,2) ~= mod(Channel_width_elements,2)
    Chamber_width_elements=Chamber_width_elements+1;
end

% Z direction 

SLL_thickness_elements=ceil(SLL_thickness/pointSpacing); %Elements used in SLL layer if method 1 applied

Mat1_base_thickness_elements=ceil(Mat1_base_thickness/pointSpacing); %mat 1 basse layer thickness (mm)
Channel_height_elements=ceil(Channel_height/pointSpacing); %height of internal channel between chambers (mm)
Channel_roof_thickness_elements=ceil(Channel_roof_thickness/pointSpacing);% (mm)
Chamber_roof_thickness_elements=ceil(Chamber_roof_thickness/pointSpacing); % (mm)
Chamber_height_elements=ceil(Chamber_height/pointSpacing); % no of mesh elements on this length (Inclusive of channel parameters)
end 

if SLL_control==2
    SLL_thickness=0;
    SLL_thickness_elements=0;
end

%% Simplfied geometry creation
% Simple geometry is defined using mesh elements only

% X direction elements needed
Length=(n*((2*Chamber_wt_elements)+Chamber_length_elements))+((n-1)*Gap_length_elements);

% Y direction elements needed
Width=(2*Side_wall_thickness_elements)+Chamber_width_elements;

% Z direction elements needed
Height=SLL_thickness_elements+Mat1_base_thickness_elements+Chamber_height_elements+Chamber_roof_thickness_elements;

boxDim=[Length Width Height];
boxE1=[Length Width Height];

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
LowerLimChannel=Chamber_wt_elements;
UpperLimChannel=((Length)-Chamber_wt_elements);

CX=VE_bar(:,1);
CXGap_length=zeros(1,size(VE_bar,1));
CXChamber=zeros(1,size(VE_bar,1));
CXChannel=zeros(1,size(VE_bar,1));

for i=1:1:size(VE_bar,1)
    
    for j=1:1:n
        %Setting Changing Limits
        LowerLimChamber=Chamber_wt_elements + ((j-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
        UpperLimChamber=Chamber_wt_elements+Chamber_length_elements  + ((j-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
        
        LowerLimGap_length=((2*Chamber_wt_elements)+Chamber_length_elements)+((j-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
        UpperLimGap_length=((j)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
        
        %Outer geometry -  if needs removal -> mark as 1
        if CX(i) > LowerLimGap_length && CX(i) < UpperLimGap_length
            CXGap_length(i)=1;
        end
        
        if CX(i) > LowerLimChamber && CX(i) < UpperLimChamber
           CXChamber(i)=1;
        end
        
        if CX(i) > LowerLimChannel && CX(i) < UpperLimChannel
            CXChannel(i)=1;
        end
        
    end
    
end
logicDeleteOuterX=logical(CXGap_length)';%converting to logical
logicDeleteInnerX1=logical(CXChamber)';
logicDeleteInnerX2=logical(CXChannel)';

% Y and Z logic -  Inner and Outer _ simpler as not effected by no. of
% chambers
logicDeleteOuterZ=VE_bar(:,3)>(SLL_thickness_elements+Mat1_base_thickness_elements+Channel_height_elements+Channel_roof_thickness_elements);
logicKeepOuter=~(logicDeleteOuterX.*logicDeleteOuterZ);
logicKeepOuter=logical(logicKeepOuter);

logicDeleteInnerY1=~(VE_bar(:,2)>(Chamber_width_elements/2));
logicDeleteInnerY2=~(VE_bar(:,2)>(Channel_width_elements/2));
logicDeleteInnerZ1=(VE_bar(:,3)>(SLL_thickness_elements+Mat1_base_thickness_elements) & VE_bar(:,3)<(Height-Chamber_roof_thickness_elements));
logicDeleteInnerZ2=(VE_bar(:,3)>(SLL_thickness_elements+Mat1_base_thickness_elements) & VE_bar(:,3)<(SLL_thickness_elements+Mat1_base_thickness_elements+Channel_height_elements));



% uses logic vectors to choose which elements to keep
logicKeepChamber=~(logicDeleteInnerX1 == 1 & logicDeleteInnerY1 == 1 & logicDeleteInnerZ1 == 1);
logicKeepChannel=~(logicDeleteInnerX2 == 1 & logicDeleteInnerY2 == 1 & logicDeleteInnerZ2 == 1);
logicKeepInner=logicKeepChamber.*logicKeepChannel;

logicKeepInner=logicKeepInner(logicKeepOuter,:);
logicKeepInner=logical(logicKeepInner);

E1=E_bar(logicKeepOuter,:); %removes outer elements

F1=element2patch(E1);
[indBoundary1]=tesBoundary(F1);

F2=element2patch(E1(logicKeepInner,:));%removes internal elements
[indBoundary2]=tesBoundary(F2);


Fb=F2(indBoundary2,:);
Cb=7*ones(size(Fb,1),1);% Sample C matrix

for q=1:1:6
    F_Cb1=Fb_bar(Cb_bar==q,:);
    logicNow=all(ismember(Fb,F_Cb1),2);
    Cb(logicNow)=q;
end

Cb(~any(ismember(Fb,F1(indBoundary1,:)),2))=0;%sets pressure faces to zero

% Remove unused nodes and clean up index matrices
 
[E,V,indFix2]=patchCleanUnused(E1(logicKeepInner,:),V_bar);
Fb=indFix2(Fb);

F=indFix2(F2);

plot_number=1;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure; 
title('Unscaled simplified geometry','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',1);%0.5 transperancy normally hold on
hold on;% plotV(V,'k.','MarkerSize',markerSize/2);
axisGeom; 
colormap(turbo(250)); icolorbar; 
camlight headlight; 
gdrawnow; 
end

V_element=V;%copy of V based on element data only - won't be scaled
V_element(:,2)=V_element(:,2)-min(V_element(:,2));%undoing centre alignment in Z axis


%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% disc. These labels can be used to define boundary conditions. 

%Define supported node sets
bcSupportList=unique(Fb(Cb==1,:)); %Node set part of selected face

%Get pressure faces
F_pressure=Fb(Cb==0,:); 

%Get end face
F_end=Fb(Cb==2,:);

bcTipList=unique(F_end);

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
ChamberAngleLeft={1 n-1};%empty for analysis of bending angle
ChamberAngleRight={1 n-2};

if n>1
for q=1:1:n-1
%where the primary and secondary for each contact pair should be found
desiredPrimary=(2*Chamber_wt_elements)+Chamber_length_elements+((q-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
desiredSecondary=(2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements+((q-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));

logicContactCoordinatePrimary=Vm(:,1)==desiredPrimary;
logicContactCoordinateSecondary=Vm(:,1)==desiredSecondary;

logicContactPrimary=logical(logicContactCoordinatePrimary.*logicContactAllPrimarySets);
logicContactSecondary=logical(logicContactCoordinateSecondary.*logicContactAllSecondarySets);


if q<n
ChamberAngleLeft{q}=F(logicContactSecondary,:);
end
if q>1
ChamberAngleRight{q-1}=F(logicContactPrimary,:);
end

ContactPair.Primary{q}=Fb(logicContactPrimary,:);
ContactPair.Secondary{q}=Fb(logicContactSecondary,:);

end

%display contact pairs
plot_number=2;
if plot_control ==1 || ismember(plot_number,plot_vector)==1
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


end


%% SLL Layer
if SLL_control == 1 % Thick SLL material

VE=patchCentre(E,V); %
logicSLL=VE(:,3)<SLL_thickness_elements;
E1=E(~logicSLL,:);%main body
E2=E(logicSLL,:);%SLL layer
E=[E1;E2];

elseif SLL_control == 2 % Shell SLL material

logic_FSLL_height=Cb==5;
FSLL=Fb(logic_FSLL_height,:);

E1=E;
E2=FSLL;
E={};
E{1}=E1;
E{2}=E2;

plot_number=3;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure;
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; gpatch(E{2},V,'g'); % SLL shell elements
title('SLL shell elements')
axisGeom;
end

elseif SLL_control == 3 
FE=patchCentre(F,V); %   
logic_FSLL_height=FE(:,3)==SLL_thickness_elements;
Normals_internal=patchNormal(F,V);
logic_FSLL_direction=Normals_internal(:,3)==1;
logic_FSLL=logical(logic_FSLL_direction.*logic_FSLL_height);
FSLL=F(logic_FSLL,:);

E1=E;
E2=FSLL;
E={};
E{1}=E1;
E{2}=E2;

plot_number=4;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure;
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; gpatch(E{2},V,'g');%SLL shell elements
title('SLL shell elements')
axisGeom;
end

end



%% Scaling to match desired geometry
%Creating vectors of needed Coordinates
V_height=zeros(1,Height); %creating empty vectors to reduce memory
V_length=zeros(1,Length);
V_width=zeros(1,Width);

%Height scaling vector creation
current_height=0;
for i=1:1:Height
    
    if i <= SLL_thickness_elements
        V_height(i)=current_height+(SLL_thickness*(1/SLL_thickness_elements));
        
    elseif i <= Mat1_base_thickness_elements+SLL_thickness_elements
        V_height(i)=current_height+(Mat1_base_thickness*(1/Mat1_base_thickness_elements));
        
    elseif i <= Channel_height_elements+Mat1_base_thickness_elements+SLL_thickness_elements
        V_height(i)=current_height+(Channel_height*(1/Channel_height_elements));
        
    elseif i <= Channel_roof_thickness_elements+Channel_height_elements+Mat1_base_thickness_elements+SLL_thickness_elements 
        V_height(i)=current_height+(Channel_roof_thickness*(1/Channel_roof_thickness_elements));
        
    elseif i <= Chamber_height_elements+Mat1_base_thickness_elements+SLL_thickness_elements
        V_height(i)=current_height+((Chamber_height-Channel_height-Channel_roof_thickness)*(1/(Chamber_height_elements-Channel_height_elements-Channel_roof_thickness_elements)));
        
    elseif  i <= Chamber_roof_thickness_elements+Chamber_height_elements+Mat1_base_thickness_elements+SLL_thickness_elements
        V_height(i)=current_height+(Chamber_roof_thickness*(1/Chamber_roof_thickness_elements));
        
    end
    current_height=V_height(i);
    
end
V_height=[0 V_height];

%Width scaling vector creation
V(:,2)=V(:,2)-min(V(:,2));%undoes centre alignment on Y axis
current_width=0;
for i=1:1:Width
    
    if i <= Side_wall_thickness_elements
        V_width(i)=current_width+(Side_wall_thickness*(1/Side_wall_thickness_elements));
    
    elseif i <= Side_wall_thickness_elements+((Chamber_width_elements/2)-(Channel_width_elements/2))
        V_width(i)=current_width+((Chamber_width-Channel_width)*(1/(Chamber_width_elements-Channel_width_elements)));
        
    elseif i <= Side_wall_thickness_elements+((Chamber_width_elements/2)+(Channel_width_elements/2))
        V_width(i)=current_width+((Channel_width)*(1/(Channel_width_elements)));
        
    elseif i <= Side_wall_thickness_elements+Chamber_width_elements
        V_width(i)=current_width+((Chamber_width-Channel_width)*(1/(Chamber_width_elements-Channel_width_elements)));
    
    elseif i<= (2*Side_wall_thickness_elements)+Chamber_width_elements
        V_width(i)=current_width+(Side_wall_thickness*(1/Side_wall_thickness_elements));
    end
    current_width=V_width(i);
    
end
V_width=[0 V_width];

%Length scaling vector creation
current_length=0;
for i=1:1:Length
    itemp=i;
    total_chamber_length_elements=((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements);
    chamber_count=floor(itemp/total_chamber_length_elements);
    itemp=itemp-(chamber_count*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
    
    
        if itemp == 0
            V_length(i)=current_length+(Gap_length*(1/Gap_length_elements));
            
        elseif itemp <= Chamber_wt_elements
            V_length(i)=current_length+(Chamber_wt*(1/Chamber_wt_elements));

        elseif itemp <= Chamber_wt_elements+Chamber_length_elements
            V_length(i)=current_length+(Chamber_length*(1/Chamber_length_elements));

        elseif itemp <= (2*Chamber_wt_elements)+Chamber_length_elements
            V_length(i)=current_length+(Chamber_wt*(1/Chamber_wt_elements));

        elseif itemp <= Gap_length_elements+(2*Chamber_wt_elements)+Chamber_length_elements
            V_length(i)=current_length+(Gap_length*(1/Gap_length_elements));

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

plot_number=5;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure; %Displays scaled geometry of two cells
title('Scaled geometry','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; plotV(V,'k.','MarkerSize',markerSize/2);
axisGeom; 
colormap(turbo(250)); icolorbar; 
% camlight headlight; 
gdrawnow; 
end




%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 
plot_number=6;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,Cb,'none',0.5);%normally set transp to 0.5, 0.05 shows pressure face better
hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',15);
hl(2)=gpatch(F_pressure,V,'r','k',1);
hl(3)=plotV(V(indForceNodes,:),'r.','MarkerSize',15);
% hl(3)=gpatch(F2,V,'c','k',1);

patchNormPlot(F_pressure,V);
legend(hl,{'BC full support','Pressure surface', 'End Points'});%'SLLsurface'

axisGeom(gca,fontSize);
camlight headlight; icolorbar;
gdrawnow; 
end


%% Object
if object_control==2
    
    % Drawing the rigid body
    V_rigid=[-rigidLength/2 rigidWidth/2; -rigidLength/2 -rigidWidth/2; rigidLength/2 -rigidWidth/2; rigidLength/2 rigidWidth/2]; %unshifted shape of rigid face
    
    regionCellRigid={V_rigid}; %A region between V1 and V2 (V2 forms a hole inside V1)
    plotOnRigid=0; %This turns on/off plotting
    pointSpacingRigid=1; %Desired point spacing
    resampleCurveOpt=1; %Option to turn on/off resampling of input boundary curves

    [F_rigid,V_rigid]=regionTriMesh2D(regionCellRigid,pointSpacingRigid,resampleCurveOpt,plotOnRigid); %creates tri mesh of rigid body
    
    %shifting position of the rigid body
    V_rigidZ=zeros(size(V_rigid,1),1); %Vector of Z heights of rigid plate
    V_rigid=[V_rigid(:,1) V_rigid(:,2) V_rigidZ]; % Adding the Z coords to the plate
    
    V_rigid(:,1)=V_rigid(:,1)+max(V(:,1))++ObjectCentreDist(1);% X direction shift
    V_rigid(:,2)=V_rigid(:,2)+(max(V(:,2))/2)+ObjectCentreDist(2);% Y direction shift - centre of Y
    V_rigid(:,3)=V_rigid(:,3)+min(V(:,3))+ObjectCentreDist(3);% Z direction shift
    
    % Merging F,V,
    F_rigid=F_rigid+size(V,1);
    Cb_rigid=ones(size(F_rigid,1),1)*(max(Cb)+1); %adding a boundary colour for plotting
    
    V=[V; V_rigid]; %updated V to include rigid body
    center_of_mass_rigid=mean(V_rigid);
    
    %contact surfaces
    contactObjectPrimary=Fb(Cb==2 | Cb==5,:); % object contact faces of the Pneunet
    contactObjectSecondary=F_rigid;
    
    
    cFigure;
    gpatch(Fb,V,Cb);
    hold on;
    gpatch(F_rigid, V, Cb_rigid)
    axisGeom; icolorbar;
    
    
elseif object_control == 3
    
    %Sphere
    if softObject_shape == 1 
        
        SoftObjectStruct.sphereRadius=sphereRadius;
        SoftObjectStruct.coreRadius=SoftObjectStruct.sphereRadius/2;
        SoftObjectStruct.numElementsMantel=3;
        SoftObjectStruct.numElementsCore=SoftObjectStruct.numElementsMantel*2; 
        
        [SoftObjectStruct]=hexMeshHemiSphere(SoftObjectStruct);
        
        Eso=SoftObjectStruct.E;
        Vso=SoftObjectStruct.V;
        Fso=SoftObjectStruct.F;
        Fbso=SoftObjectStruct.Fb;
        Cbso=SoftObjectStruct.boundaryMarker;
        
    end
    
    % Cuboid
    if softObject_shape == 2 
        
        boxDim=[SoftObjectLength SoftObjectWidth SoftObjectHeight];
        boxEl=[ceil(SoftObjectLength/SoftObject_PointSpacing) ceil(SoftObjectWidth/SoftObject_PointSpacing) ceil(SoftObjectHeight/SoftObject_PointSpacing)];

        [SoftObjectStruct]=hexMeshBox(boxDim,boxEl);
        
        Eso=SoftObjectStruct.E;
        Vso=SoftObjectStruct.V;
        Fso=SoftObjectStruct.F;
        Fbso=SoftObjectStruct.Fb;
        Cbso=SoftObjectStruct.faceBoundaryMarker;
        
    end
   
    % Shifting the location of the Soft Object
    Vso(:,1)=Vso(:,1)+max(V(:,1))+ObjectCentreDist(1);% X direction shift
    Vso(:,2)=Vso(:,2)+(max(V(:,2))/2)+ObjectCentreDist(2);% Y direction shift - centre of Y
    Vso(:,3)=Vso(:,3)+min(V(:,3))+ObjectCentreDist(3);% Z direction shift
    
    cFigure;
    gpatch(Fb,V,'rw');
    hold on;
    gpatch(Fbso,Vso,Cbso)
    axisGeom;
    icolorbar;
    
    % Taking indices of object boundary conditions
    if softObject_shape ==1 %sphere
        
        bcMidPlane=unique(Fbso(Cbso==2,:))+size(V,1);
        contactObjectSecondary=Fbso(Cbso==1,:)+size(V,1); %outer contact face of hemisphere
    
    elseif softObject_shape == 2 %cuboid
        
        bcMidPlane=unique(Fbso(Cbso==5))+size(V,1);
        contactObjectSecondary=Fbso(Cbso==6 | Cbso==1 | Cbso==2,:)+size(V,1); %outer contact faces of cuboid
        
    end
    
    contactObjectPrimary=Fb(Cb==2 | Cb==5,:); % object contact faces of the Pneunet
    
    % Merging V, E, F matrices
    
    Fbso2=Fbso+size(V,1);
    E3=Eso+size(V,1);
    V=[V;Vso];
    
    cFigure;
    gpatch(Fbso2,V,'bw');
    hold on;
    gpatch(Fb,V,'rw');
    axisGeom;
end

%% .stl file creation
% Pneunet body .stl files
if stl_control==1% if 3D printing .stls are needed
    [F_PneunetBody,~,~]=element2patch(E1,[],'hex8');
    Fb_PneunetBody=F_PneunetBody(tesBoundary(F_PneunetBody),:);%getting boundary of Pneunet(no SLL)
    Fb_PneunetBodyTri=quad2tri(Fb_PneunetBody,V);%changing to tri for .stl
    [Fb_PneunetBodyTri,V_stl_temp]=patchCleanUnused(Fb_PneunetBodyTri,V);%removing unused V to stop warning in command window
    
    TR_PnuenetBody=triangulation(Fb_PneunetBodyTri,V_stl_temp);
    stlwrite(TR_PnuenetBody,STL_Pnuenet_Body);%writing .stl file
    
   if SLL_control==1 % create SLL .stl if required
       [F_PneunetSLL,~,~]=element2patch(E2,[],'hex8');
        Fb_PneunetSLL=F_PneunetBody(tesBoundary(F_PneunetSLL),:);%getting boundary of SLL
        Fb_PneunetSLLTri=quad2tri(Fb_PneunetSLL,V);%changing to tri for .stl
        [Fb_PneunetSLLTri,V_stl_temp]=patchCleanUnused(Fb_PneunetSLLTri,V);%removing unused V to stop warning in command window
    
        TR_PnuenetSLL=triangulation(Fb_PneunetSLLTri,V_stl_temp);
        stlwrite(TR_PnuenetSLL,STL_Pnuenet_SLL);%writing .stl file
       
       
   end
   



% Mould .stl files
elseif stl_control==2
    MouldThickness=3; %thickness outwards of the mould
bottomDimEl=[3 3 2];

[meshStruct]=hexMeshBox(bottomDimEl,bottomDimEl);

E_bar=meshStruct.E;
V_bar=meshStruct.V;
F_bar=meshStruct.F;
Fb_bar=meshStruct.Fb;
Cb_bar=meshStruct.faceBoundaryMarker;

%Getting centre for easier control of elements
VE_STLBase=patchCentre(E_bar,V_bar);

%Shifting data for easier deletion of sink elements
VE_STLBase(:,1)=abs(VE_STLBase(:,1))
VE_STLBase(:,2)=abs(VE_STLBase(:,2))
VE_STLBase(:,3)=VE_STLBase(:,3)-min(VE_STLBase(:,3))+0.5;

%Deleting sink elements for mould pouring
logicDeleteSTLBase= VE_STLBase(:,3)>0 & VE_STLBase(:,2) ==0 & VE_STLBase(:,1) == 0;
E_bottomMould=E_bar(~logicDeleteSTLBase,:);

%scaling the bottom mould
% bottomMouldDimX=[]
% bottomMouldDimY=
% bottomMouldDimZ=

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
material_number=1; %counter for material to reduce loops, easier to read than indexing size

materialName1='Material1';
febio_spec.Material.material{material_number}.ATTR.name=materialName1;
febio_spec.Material.material{material_number}.ATTR.type='neo-Hookean';
febio_spec.Material.material{material_number}.ATTR.id=material_number;
febio_spec.Material.material{material_number}.E=1;
febio_spec.Material.material{material_number}.v=0.49;

if SLL_control == 1
material_number=material_number+1;

materialName2='Material2';
febio_spec.Material.material{material_number}.ATTR.name=materialName2;
febio_spec.Material.material{material_number}.ATTR.type='Ogden';
febio_spec.Material.material{material_number}.ATTR.id=material_number;
febio_spec.Material.material{material_number}.c1=c2;
febio_spec.Material.material{material_number}.m1=m2;
febio_spec.Material.material{material_number}.c2=c2;
febio_spec.Material.material{material_number}.m2=-m2;
febio_spec.Material.material{material_number}.k=k2;

elseif SLL_control >= 2
material_number=material_number+1;
    
materialName2='Material2';
febio_spec.Material.material{material_number}.ATTR.name=materialName2;
febio_spec.Material.material{material_number}.ATTR.type='neo-Hookean';
febio_spec.Material.material{material_number}.ATTR.id=material_number;
febio_spec.Material.material{material_number}.E=E_material2;
febio_spec.Material.material{material_number}.v=Poissons_material2;
end

if  object_control == 2
    material_number=material_number+1;
    rigidID=material_number;
    materialRigidBody='MaterialRigidBody';
    febio_spec.Material.material{material_number}.ATTR.name=materialRigidBody;
    febio_spec.Material.material{material_number}.ATTR.type='rigid body';
    febio_spec.Material.material{material_number}.ATTR.id=material_number;
    febio_spec.Material.material{material_number}.density=1;
    febio_spec.Material.material{material_number}.center_of_mass=center_of_mass_rigid;

    
elseif object_control==3 % soft object material
    material_number=material_number+1;
    
    MaterialSoftObject='Material3';
    febio_spec.Material.material{material_number}=materialBank.SoftMaterialSample1;
    febio_spec.Material.material{material_number}.ATTR.name=MaterialSoftObject;
    febio_spec.Material.material{material_number}.ATTR.id=material_number;
     
end

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates 

% -> Elements
element_number=1; %counter for elements, looks cleaner than more if loops or indexing based on size

partBodyPneunet='Part1';
febio_spec.Mesh.Elements{element_number}.ATTR.name=partBodyPneunet; %Name of this part
febio_spec.Mesh.Elements{element_number}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{element_number}.elem.VAL=E1; %The element matrix

if SLL_control == 1 %Thick SLL -> hexa8
    element_number=element_number+1;
    
partSLL='Part2';
febio_spec.Mesh.Elements{element_number}.ATTR.name=partSLL; %Name of this part
febio_spec.Mesh.Elements{element_number}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{element_number}.elem.VAL=E2; %The element matrix

elseif SLL_control >= 2 %Thin SLL -> quad4
    element_number=element_number+1;
    
partSLL='Part2';
febio_spec.Mesh.Elements{element_number}.ATTR.name=partSLL; %Name of this part
febio_spec.Mesh.Elements{element_number}.ATTR.type='quad4'; %Element type 
febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{element_number}.elem.VAL=E2; %The element matrix

end

if SLL_control==0
    E2=[]; %setting empty E2 to remove need for homogenous Pnuenet loop for object control
end

if  object_control ==2
    element_number=element_number+1;
    
    partRigidBody='Part3';
    febio_spec.Mesh.Elements{element_number}.ATTR.name=partRigidBody; %Name of this part
    febio_spec.Mesh.Elements{element_number}.ATTR.type='tri3'; %Element type 
    febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=size(E1,1)+size(E2,1)+(1:1:size(F_rigid,1))'; %Element id's
    febio_spec.Mesh.Elements{element_number}.elem.VAL=F_rigid; %The element matrix
     
    
elseif object_control ==3 % Assigning part for soft object
    element_number=element_number+1;

    partSoftObject='Part3';
    febio_spec.Mesh.Elements{element_number}.ATTR.name=partSoftObject; %Name of this part
    febio_spec.Mesh.Elements{element_number}.ATTR.type='hex8'; %Element type 
    febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=size(E1,1)+size(E2,1)+(1:1:size(E3,1))'; %Element id's
    febio_spec.Mesh.Elements{element_number}.elem.VAL=E3; %The element matrix
       
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

else %where there is no chamber contact
q=1;%counter for contact surfaces
contactSurfacesValPrimary(q)=0;% dummy value to allow other contact using same input
contactSurfacesValSecondary(q)=1;% dummy value to allow other contact using same input

end

if object_control>=2 %if there is an object/rigid body
contactSurfacesObjectPrimary='contactObjectPrimary';
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)+2}.ATTR.name=contactSurfacesObjectPrimary;
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)+2}.quad4.ATTR.id=(1:1:size(contactObjectPrimary,1))';%size(F_pressure,1)+
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)+2}.quad4.VAL=contactObjectPrimary;

contactSurfacesObjectSecondary='contactObjectSecondary';
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.ATTR.name=contactSurfacesObjectSecondary;

if object_control==3 %if soft object -> quad4
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.quad4.ATTR.id=(1:1:size(contactObjectSecondary,1))';%size(F_pressure,1)+size(F_contactPrimary,1)+
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.quad4.VAL=contactObjectSecondary;
elseif object_control==2 %if soft object -> tri3
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.tri3.ATTR.id=(1:1:size(contactObjectSecondary,1))';%size(F_pressure,1)+size(F_contactPrimary,1)+
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.tri3.VAL=contactObjectSecondary;
end

% -> Surface pairs
if contact_control==1 % shifting to match up q values if no contact
    q=0;
elseif n==1
    q=0;
end

numStringPair=num2str(q+1);%string for number of the name of pair
febio_spec.Mesh.SurfacePair{q+1}.ATTR.name=strcat(leaderStringContact,numStringPair);
febio_spec.Mesh.SurfacePair{q+1}.primary=contactSurfacesObjectPrimary;
febio_spec.Mesh.SurfacePair{q+1}.secondary=contactSurfacesObjectSecondary;
end


% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

nodeSetName2='bcTipList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.node.ATTR.id=indForceNodes(:);

if object_control==3
nodeSetName3='bcSoftObjectList';
febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
febio_spec.Mesh.NodeSet{3}.node.ATTR.id=bcMidPlane(:);
end

%MeshDomains section
solid_number=1;
shell_number=0;

febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.name=partBodyPneunet;
febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.mat=materialName1;

if SLL_control == 1 %SLL domain
    solid_number=solid_number+1;
    
febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.name=partSLL;
febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.mat=materialName2;
elseif SLL_control >= 2
    shell_number=shell_number+1;
    
febio_spec.MeshDomains.ShellDomain{shell_number}.ATTR.name=partSLL;
febio_spec.MeshDomains.ShellDomain{shell_number}.ATTR.mat=materialName2;
end

if object_control==2 %object domain
    shell_number=shell_number+1;
    
febio_spec.MeshDomains.ShellDomain{shell_number}.ATTR.name=partRigidBody;
febio_spec.MeshDomains.ShellDomain{shell_number}.ATTR.mat=materialRigidBody;
    
    
elseif object_control==3
    solid_number=solid_number+1;

    febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.name=partSoftObject;
    febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.mat=MaterialSoftObject;
 
end

%MeshData secion
%-> Element data
if SLL_control >= 2
febio_spec.MeshData.ElementData{1}.ATTR.var='shell thickness';
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partSLL;
febio_spec.MeshData.ElementData{1}.elem.ATTR.lid=(1:1:size(E2,1))';
febio_spec.MeshData.ElementData{1}.elem.VAL=SLLThicknessShell*ones(size(E2,1),size(E2,2));
end

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.type='fix';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dofs='x,y,z';

if object_control==3 %fixing the object nodes
febio_spec.Boundary.bc{2}.ATTR.type='fix';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName3;
febio_spec.Boundary.bc{2}.dofs='x,y,z';
end

if object_control==1
febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.type='prescribe';
febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.node_set=nodeSetName2;
febio_spec.Step.step{2}.Boundary.bc{1}.dof='z';
febio_spec.Step.step{2}.Boundary.bc{1}.scale.ATTR.lc=2;
febio_spec.Step.step{2}.Boundary.bc{1}.scale.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{1}.relative=1;
end

%Loads section
% -> Surface load
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=appliedPressure;
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

if object_control>=2
febio_spec.Contact.contact{q+1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{q+1}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{q+1}.ATTR.name;
febio_spec.Contact.contact{q+1}.two_pass=1;
febio_spec.Contact.contact{q+1}.laugon=laugon;
febio_spec.Contact.contact{q+1}.tolerance=0.2;
febio_spec.Contact.contact{q+1}.gaptol=0;
febio_spec.Contact.contact{q+1}.minaug=minaug;
febio_spec.Contact.contact{q+1}.maxaug=maxaug;
febio_spec.Contact.contact{q+1}.search_tol=0.01;
febio_spec.Contact.contact{q+1}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));
febio_spec.Contact.contact{q+1}.symmetric_stiffness=0;
febio_spec.Contact.contact{q+1}.auto_penalty=1;
febio_spec.Contact.contact{q+1}.penalty=contactPenalty;
febio_spec.Contact.contact{q+1}.fric_coeff=fric_coeff;

end
end

%Rigid section 
% -> Prescribed rigid body boundary conditions
febio_spec.Rigid.rigid_constraint{1}.ATTR.name='RigidFix_1';
febio_spec.Rigid.rigid_constraint{1}.ATTR.type='fix';
febio_spec.Rigid.rigid_constraint{1}.rb=rigidID;
febio_spec.Rigid.rigid_constraint{1}.dofs='Rx,Ry,Rz,Ru,Rv,Rw';

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 0.7; 1.1 0.8; 2 1];

if object_control==1 % Loadcurve for the Z node fix method 
febio_spec.LoadData.load_controller{2}.ATTR.id=2;%loadcurve ID no.
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='STEP'; 
febio_spec.LoadData.load_controller{2}.points.point.VAL=[1 0; 2 1];
end

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
    if object_control==3
febio_spec.Output.logfile.element_data{1}.VAL=1:(size(E{1},1)+size(E{2},1)+size(E3,1));
    else
febio_spec.Output.logfile.element_data{1}.VAL=1:(size(E{1},1)+size(E{2},1));        
    end
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
    plot_number=7;
    if plot_control==1 || ismember(plot_number,plot_vector)==1
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
    end        
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
     elseif SLL_control>=2
    [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,end));
%     CV{1}=faceToVertexMeasure(E1,V,E_stress_mat(:,:,end));
%     CV{2}=faceToVertexMeasure(E2,V,E_stress_mat(:,:,end));
    end
    
    % Create basic view and store graphics handle to initiate animation
    plot_number=8;
    if plot_control==1 || ismember(plot_number,plot_vector)==1
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
        elseif SLL_control==2
        [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,qt));
        end    
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
    end


%% Bending Angle of Each Chamber
Theta=zeros(1,n-2);


if n>4 % results only comparable with 5+ chambers as there are no chambers far enough from an end of the Pneunet
    
    for i=1:1:n-2 
    normLeft=patchNormal(ContactPair.Secondary{i},V_DEF(:,:,size(V_DEF,3)));%normals of left face of chamber wall
    normRight=patchNormal(ContactPair.Primary{i+1},V_DEF(:,:,size(V_DEF,3)));%normals of right face of chamber wall
    
    normLeft=mean(normLeft,1);%average normal for left side
    normRight=mean(normRight,1);%average normal for right side
    
    Theta(i)=rad2deg(acos((dot(normLeft,normRight,2))/((norm(normLeft)*(norm(normRight))))));% angle between chamber wall normals
    
    
    
    end
    
    plot_number=9;
    if plot_control==1 || ismember(plot_number,plot_vector)==1
    cFigure;
    gpatch(Fb,V_DEF(:,:,size(V_DEF,3)),Cb,'k',0.2); hold on;
    title('Contact surface normal vectors')
    for i=1:1:n-2 
        patchNormPlot(ContactPair.Primary{i+1},V_DEF(:,:,size(V_DEF,3)));
        patchNormPlot(ContactPair.Secondary{i},V_DEF(:,:,size(V_DEF,3))); 
    end
    axisGeom;
    disp(Theta)%display chamber angles to command window
    end
    
end


%% Bending Angle using Tip Nodes
BendingAngle=zeros(size(V_DEF,3),3);%empty vector of bending angles over time
BendingAngle(:,3)=timeVec;%assigning time values
YZnormVec=[1 0 0];

BendingAngleRefPoint=mean(V(indForceNodes,:),1);

for i=1:1:size(V_DEF,3)
    
    normEnd=mean(patchNormal(F_end,V_DEF(:,:,i)),1);%avergage normal vector of the last face
    nodeEnd=mean(V_DEF(indForceNodes,:,i),1);%averaged end Node tracking over time
    
    BendingAngle(i,1)=rad2deg(acos((dot(normEnd,YZnormVec,2))/((norm(normEnd)*(norm(YZnormVec))))));%angle of end face normal vector using 2D dot product
    BendingAngle(i,2)=rad2deg(atan((abs(nodeEnd(3)))/(nodeEnd(1))));%saving angle of end node
    
end

plot_number=10;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure;
ha(1)=plot(BendingAngle(:,3),BendingAngle(:,1),'r');
hold on;
ha(2)=plot(BendingAngle(:,3),BendingAngle(:,2),'b');
legend(ha,{'Normal Vector Angle', 'End Point Angle'});%'SLLsurface'
title('Bending Angle');
end
%% Force

% Importing nodal displacements from a log file
dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),1,1);

%Access data
Force_mat=dataStruct.data; %All force values
timeVec=dataStruct.time; %Time

% Extract Forces at Restricted nodes
Force_mean=zeros(1,size(Force_mat,3));
Force_total=zeros(1,size(Force_mat,3));

for i=1:1:size(Force_mat,3)
    Force_nodal=Force_mat(indForceNodes,3,i); %Force only at restricted nodes
    Force_mean(i)=mean(Force_nodal);%mean force at the end nodes due to Z restriction
    Force_total(i)=sum(Force_nodal);% sum of all the end node forces due to Z restriction
end


%Plotting
plot_number=11;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure;
hforce(1)=plot(timeVec,Force_mean,'r');
hold on;
hforce(2)=plot(timeVec,Force_total,'b');
title('Nodal Force ()');
xlabel('Time (s)'); ylabel('Average End-Node Force ()');
end

end % successful run check




