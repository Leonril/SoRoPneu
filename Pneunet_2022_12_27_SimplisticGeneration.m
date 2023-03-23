%% Improving prevuious models mesh sizing based on key dimensions
% this allows:
%
% * More robust analysis
% * Simpler control over the mesh sizings in key areas 
%
% This route was not used previously due to
%
% * Increases complexity of the code
%
%
% Next steps involve:
%
% * Potentially using interpolation to all smoothing of outer mesh to inner
% mesh. This would allow mesh sizing based on major dimensions only and
% inner dimensions would interpolate points to meet the outer - even if not
% matched well in some axes.

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

%% Geometry Inputs:

n=3; %no. of chambers 

% X direction 

Chamber_length=8; % internal chamber length (mm)
Chamber_wt=0.5; % wall thickness expanding (mm)
Gap_length=1; %(mm)

% y direction 

Side_wall_thickness=1;
Chamber_width=8; % internal chamber width (mm)
Channel_width=2; % width of the intenal channel (mm)

% Z direction 

SLL_thickness=1; %Strain Limiting Layer (mm)% for Method 2, this could be made zero and Mat1 increased
Mat1_base_thickness=0.5; %mat 1 basse layer thickness (mm)
Channel_height=1; %height of internal channel between chambers (mm)
Channel_roof_thickness=1;% (mm)
Chamber_roof_thickness=2.5; % (mm)
Chamber_height=10; %internal height of the chamber (mm)

%% Desired Mesh Inputs

% X direction

Chamber_length_elements=2; % no of mesh elements on this length
Chamber_wt_elements=1; % no of mesh elements on this length
Gap_length_elements=1; % no of mesh elements on this length

% y direction 

Side_wall_thickness_elements=1; % no of mesh elements on this length
Chamber_width_elements=3; % no of mesh elements on this length (Inclusive of channel width)
Channel_width_elements=1; % no of mesh elements on this length

% Z direction 

SLL_thickness_elements=1; %Strain Limiting Layer (mm)% for Method 2, this could be made zero and Mat1 increased
Mat1_base_thickness_elements=1; %mat 1 basse layer thickness (mm)
Channel_height_elements=1; %height of internal channel between chambers (mm)
Channel_roof_thickness_elements=1;% (mm)
Chamber_roof_thickness_elements=1; % (mm)
Chamber_height_elements=4; % no of mesh elements on this length (Inclusive of channel parameters)


%% Simplfied geometry creation
% Geometry is defined using mesh elements only

% X direction
Length=(2*((2*Chamber_wt_elements)+Chamber_length_elements))+(Gap_length_elements);


% Y direction
Width=(2*Side_wall_thickness_elements)+Chamber_width_elements;


% Z direction
Height=SLL_thickness_elements+Mat1_base_thickness_elements+Chamber_height_elements+Chamber_roof_thickness_elements;

boxDim=[Length Width Height];
boxE1=[Length Width Height];

[meshStruct]=hexMeshBox(boxDim,boxE1);

E_bar=meshStruct.E;
V_bar=meshStruct.V;
F_bar=meshStruct.F;
Fb_bar=meshStruct.Fb;
Cb_bar=meshStruct.faceBoundaryMarker;

% for i=1:1:3 % shifting everything to positive coordinates 
%    V_bar(:,i)=V_bar(:,i)-min(V_bar(:,i)); 
% end
V_bar(:,3)=V_bar(:,3)-min(V_bar(:,3)); 

%Simple mesh based initial geometry
cFigure;
gpatch(Fb_bar,V_bar,Cb_bar);
axisGeom;

VE_bar=patchCentre(E_bar,V_bar); 
VE_bar=abs(VE_bar);

%Removing unneeded elements to create two chambers, if elements are needed
%to create the diveide, they are removed using logic vectors
logicDeleteOuterX=VE_bar(:,1)<(Gap_length_elements/2);
logicDeleteOuterZ=VE_bar(:,3)>(SLL_thickness_elements+Mat1_base_thickness_elements+Channel_height_elements+Channel_roof_thickness_elements);
logicKeepOuter=~(logicDeleteOuterX.*logicDeleteOuterZ);
logicKeepOuter=logical(logicKeepOuter);


logicDeleteInnerX1=(VE_bar(:,1)>((Gap_length_elements/2)+Chamber_wt_elements) & VE_bar(:,1)<((Gap_length_elements/2)+Chamber_wt_elements+Chamber_length_elements));
logicDeleteInnerX2=(VE_bar(:,1)<((Gap_length_elements/2)+Chamber_wt_elements+Chamber_length_elements));
logicDeleteInnerY1=~(VE_bar(:,2)>(Chamber_width_elements/2));
logicDeleteInnerY2=~(VE_bar(:,2)>(Channel_width_elements/2));
logicDeleteInnerZ1=(VE_bar(:,3)>(SLL_thickness_elements+Mat1_base_thickness_elements) & VE_bar(:,3)<(Height-Chamber_roof_thickness_elements));
logicDeleteInnerZ2=(VE_bar(:,3)>(SLL_thickness_elements+Mat1_base_thickness_elements) & VE_bar(:,3)<(SLL_thickness_elements+Mat1_base_thickness_elements+Channel_height_elements));

logicKeepChamber=~(logicDeleteInnerX1 == 1 & logicDeleteInnerY1 == 1 & logicDeleteInnerZ1 == 1);
logicKeepChannel=~(logicDeleteInnerX2 == 1 & logicDeleteInnerY2 == 1 & logicDeleteInnerZ2 == 1);
logicKeepInner=logicKeepChamber.*logicKeepChannel;

logicKeepInner=logicKeepInner(logicKeepOuter,:);
logicKeepInner=logical(logicKeepInner);

E1=E_bar(logicKeepOuter,:); %removes outer elements
cFigure
gpatch(E1,V_bar)

F1=element2patch(E1);
[indBoundary1]=tesBoundary(F1);

F2=element2patch(E1(logicKeepInner,:));%removes internal elements
[indBoundary2]=tesBoundary(F2);


Fb=F2(indBoundary2,:);
Cb=7*ones(size(Fb,1),1);

for q=1:1:6
    F_Cb1=Fb_bar(Cb_bar==q,:);
    logicNow=all(ismember(Fb,F_Cb1),2);
    Cb(logicNow)=q;
end

Cb(~any(ismember(Fb,F1(indBoundary1,:)),2))=0;

% Remove unused nodes and clean up index matrices
 
[E,V,indFix2]=patchCleanUnused(E1(logicKeepInner,:),V_bar);
Fb=indFix2(Fb);

F=indFix2(F2);

cFigure; 
title('Unscaled simplified geometry','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; plotV(V,'k.','MarkerSize',markerSize/2);
axisGeom; 
colormap(turbo(250)); icolorbar; 
camlight headlight; 
gdrawnow; 

%% Scaling to match desired geometry



















