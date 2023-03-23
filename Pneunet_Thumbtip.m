%% Revisions
%% Rev D offshoot Thumbtip
% Notes:
% * 
% * 
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
SLL_control=1;
%==1: Adds a second material layer on the bottom of the pneunet of
%specified thickness
%==2: Adds a thin shell element the base of the structure
%==3: Adds thin shell SLL elements at in between the top and bottom layers
%created in method one. Bottom layer is now the sane material as the top
%and material 2 is the shell.

%Contact
contact_control=2;
%==1: Contact modelling inactive
%==2: Contact modelling active

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

%% FEA control settings
numTimeSteps=50; %Number of time steps desired
opt_iter=25; %Optimum number of iterations
max_refs=opt_iter*4; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
max_retries=10; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps)*4; %Maximum time step size

runMode='external';%'internal';

%% controlContact parameters
contactPenalty=5;%5
laugon=0;
minaug=1;
maxaug=15;
fric_coeff=0;


%% Load Inputs
%Load
appliedPressure=0.1; %0.15

%% Material Properties
%Material parameter set

%EcoFlex 00-50 (Pagoli et al)(Elsayed et al)
c1_mat1=107.9*(10^(-3)); %Shear-modulus-like parameter
m1_mat1=1.55; %Material parameter setting degree of non-linearity
c2_mat1=21.47*(10^(-6));
m2_mat1=7.86;
c3_mat1=-87.1*(10^(-3));
m3_mat1=-1.91;

%k_factor=100; %Bulk modulus factor 
k1=175;%c1*k_factor; %Bulk modulus


if SLL_control==1
c2=c1_mat1*50; %Shear-modulus-like parameter
m2=2; %Material parameter setting degree of non-linearity
k2=c2*100;%k_factor; %Bulk modulus
end

%Paper




%% Geometry Inputs:

n=5; %no. of chambers 
pointSpacing=1.5;%Only active if MeshCreationTool == 2, designated aproximate size of desired cube element length and width if possible (mm) 

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

%Simple mesh based initial geometry
% cFigure;
% gpatch(Fb_bar,V_bar,Cb_bar);
% axisGeom;

VE_bar=patchCentre(E_bar,V_bar);

VE_bar=abs(VE_bar);

%check if n is even or odd
evenTest=mod(n,2); %check if n chambers is even or odd

%Removing unneeded elements to create two chambers, if elements are needed
%to create the divide, they are removed using logic vectors


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

%% Here 09 Feb 2023
VE=patchCentre(E,V);
logicLast=VE(:,1)>(max(VE(:,1))-1);
E_last=E(logicLast,:);

cFigure;
gpatch(E_last,V);
axisGeom;

%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% disc. These labels can be used to define boundary conditions. 

%Define supported node sets
bcSupportList=unique(Fb(Cb==1,:)); %Node set part of selected face

%Get pressure faces
F_pressure=Fb(Cb==0,:); 

%Get end face
F_end=Fb(Cb==2,:);

%Get end face corner nodes
Spine_list=unique(Fb(Cb==5,:));

[LIA,indForceNodes]=ismember(unique(F_end),Spine_list);%unique index when both end conditions are met
indForceNodes=indForceNodes(indForceNodes~=0);
%% Defingin Contact surfaces
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
cFigure;
hold on;
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
logicSLL=VE(:,3)<SLL_thickness_elements;
E1=E(~logicSLL,:);%main body
E2=E(logicSLL,:);%SLL layer
E=[E1;E2];

elseif SLL_control == 2 % Shell SLL material

logic_FSLL=Cb==5;
FSLL=Fb(logic_FSLL,:);

E1=E;
E2=FSLL;
E={};
E{1}=E1;
E{2}=E2;

cFigure;
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; gpatch(E{2},V,'g');
axisGeom;


elseif SLL_control == 3 
FE=patchCentre(F,V); %   
logic_FSLL=FE(:,3)==SLL_thickness_elements;
FSLL=F(logic_FSLL,:);

E1=E;
E2=FSLL;
E={};
E{1}=E1;
E{2}=E2;

cFigure;
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; gpatch(E{2},V,'g');
axisGeom;

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


cFigure; %Displays scaled geometry of two cells
title('Scaled geometry','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; plotV(V,'k.','MarkerSize',markerSize/2);
axisGeom; 
colormap(turbo(250)); icolorbar; 
% camlight headlight; 
gdrawnow; 

%% Thumbtip
% E_last=
Thumb_length=5;
numSteps=2;
layerThickness=Thumb_length/numSteps;

V_basethumb=V(F_end,:);%elements of thumb base
[E_new,VE_new,FP1,FP2]=patchThick(F_end,V,1,layerThickness,numSteps);


cFigure;
gpatch(F,VE_new,1); hold on;
gpatch(E_new,VE_new,2); axisGeom;

cFigure;
gpatch(F,VE_new,1); hold on;
gpatch(FP1,VE_new,'bw'); 
gpatch(FP2,VE_new,'rw');axisGeom;
%%
C_newTemp=1:1:size(E_new,1);
[F_new,C_newTemp,CF_newTemp]=element2patch(E_new,C_newTemp,'hex8');
%%
cFigure;
gpatch(F_new,VE_new); hold on;
axisGeom;
%%
cFigure;
gpatch(F_new,VE_new,CF_newTemp); hold on;
axisGeom;


%% 

CF_newTemp=CF_newTemp+max(Cb);
Cb=[Cb;CF_newTemp];

Fb_unfin=[Fb;F_new];
%%
cFigure;
gpatch(Fb_unfin,VE_new,Cb);
axisGeom; icolorbar;

cFigure;
gpatch(F_new,VE_new,CF_newTemp);
axisGeom; icolorbar;
%%
logicClear=~ismember(Cb,[2; 8]);%removes faces made internal from the thumb tip addition
Fb_unfin2=Fb_unfin(logicClear,:);
%%
Cb_2=Cb(logicClear,:);
%%
Cb_2(Cb_2==9)=2;
Cb_2(Cb_2==10)=5;
Cb_2(Cb_2==11)=6;
Cb_2(Cb_2==12)=4;
Cb_2(Cb_2==13)=3;
%%
cFigure;
gpatch(Fb_unfin2,VE_new,Cb_2)
axisGeom; icolorbar;

%%
cFigure;
gpatch(Fb_unfin2(Cb==9,:),VE_new,9);
hold on;
gpatch(Fb_unfin2(Cb==3,:),VE_new,3);
gpatch(Fb_unfin2(Cb==2,:),VE_new,2);
axisGeom; icolorbar;







% T_maxZ=max(V_basethumb(:,3))/2%max R of thumb in Z
% T_maxY=max(V_basethumb(:,2))/2%max R of thumb in Y
% 
% V_basethumb(:,1)
% 
% Length_thumb=10;
% El_thumb=ceil(Length_thumb/pointSpacing)
% 
% X_co=unique(V_basethumb(:,1))
% Y_co=unique(V_basethumb(:,2))
% Z_co=unique(V_basethumb(:,3))
% 
% numElTY=size(Y_co,1)
% numElTZ=size(Z_co,1)
% 
% boxDim=[numElTZ numElTY El_thumb];
% boxEl=[numElTZ numElTY El_thumb];
% 
% [meshStructT]=hexMeshBox(boxDim,boxEl);
% 
% ET_bar=meshStructT.E;
% VT_bar=meshStructT.V;
% FT_bar=meshStructT.F;
% FbT_bar=meshStructT.Fb;
% CbT_bar=meshStructT.faceBoundaryMarker;
% 
% VT_bar(:,1)=VT_bar(:,1)-min(VT_bar(:,1));
% VT_bar(:,2)=VT_bar(:,2)-min(VT_bar(:,2));
% VT_bar(:,3)=VT_bar(:,3)-min(VT_bar(:,3));
% 
% 
% for i=1:1:size(VT_bar,1)
%    VT_bar(i,1)=X_co+(VT_bar(i,1)*pointSpacing);
%    VT_bar(i,2)=Y_co(VT_bar(i,2)+1)
%    disp('bla')
%    VT_bar(i,3)=Z_co(VT_bar(i,3)+1)
% end
% 
% cFigure;
% gpatch(FbT_bar,VT_bar,CbT_bar)
% axisGeom;
% 


% if mod(Height)==0
% 
%     nStacks=Height/2;
%    
% elseif mod(Height)==1
% 
%     nStacks=(Height-1)/2;
%     
% end
% 
% for i=1:1:nStacks
%     
%     
    
    
    





















% T_tol=pointSpacing/5;%tolerance for the thumb
% 
% %testV=
% 
% xs=0:pointSpacing:ceil(T_maxZ,pointSpacing)*pointSpacing
% xs=repmat(xs,size(V_basethumb,1),1)
% xs=xs(:)
% V_basethumb=repmat(V_basethumb,n,1);
% VS(:,1)=VS(:,1)+xs;
% 
% 
% VS=repmat(V1,n,1);
% VS(:,1)=VS(:,1)+xs;




















