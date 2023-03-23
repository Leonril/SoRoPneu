%% Revisions
%% Rev A
% Notes:
% * SLL layer nearly finalised - method 3 doesnt converge
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
MeshCreationTool=1;%
% == 1: Uses 'must points' at key variable locations
% == 2: Uses some 'must points' -> uses interpolation to scale internal
% nodes
% == 3: Uses some 'must points' -> uses smoothing to scale internal nodes

%SLL method
SLL_control=1;
%==1: Adds a second material layer on the bottom of the pneunet of
%specified thickness
%==2: Adds a thin shell element the base of the structure
%==3: Adds thin shell SLL elements at in between the top and bottom layers
%created in method one. Bottom layer is now the sane material as the top
%and material 2 is the shell.

%Contact
contact_control=1;
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
numTimeSteps=10; %Number of time steps desired
opt_iter=25; %Optimum number of iterations
max_refs=opt_iter*8; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
max_retries=10; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps)*4; %Maximum time step size

runMode='internal';%'internal';

%% controlContact parameters
contactPenalty=5;%5
laugon=1;
minaug=1;
maxaug=15;
fric_coeff=0;


%% Load Inputs
%Load
appliedPressure=0.15; %0.15

%% Material Properties
%Material parameter set
c1=1; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=100; %Bulk modulus factor 
k1=c1*k_factor; %Bulk modulus

c2=c1*2; %Shear-modulus-like parameter
m2=2; %Material parameter setting degree of non-linearity
k2=c2*k_factor; %Bulk modulus

%% Geometry Inputs:

n=4; %no. of chambers 

% X direction 

Chamber_length=15; % internal chamber length (mm)
Chamber_wt=1.5; % wall thickness expanding (mm)
Gap_length=1.8; %(mm)

% y direction 

Side_wall_thickness=2;
Chamber_width=15; % internal chamber width (mm)
Channel_width=2; % width of the intenal channel (mm)

% Z direction 

SLL_thickness=1.5; %Strain Limiting Layer (mm)% for Method 2, this could be made zero and Mat1 increased
SLLThicknessShell=0.1;

Mat1_base_thickness=3; %mat 1 basse layer thickness (mm)
Channel_height=3; %height of internal channel between chambers (mm)
Channel_roof_thickness=1;% (mm)
Chamber_roof_thickness=4; % (mm)
Chamber_height=25; %internal height of the chamber (mm)

%% Desired Mesh Inputs

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


%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% disc. These labels can be used to define boundary conditions. 

%Define supported node sets
bcSupportList=unique(Fb(Cb==1,:)); %Node set part of selected face

%Get pressure faces
F_pressure=Fb(Cb==0,:); 

% Contact surfaces
logicContactSurf=Cb==7;%chamber wall faces
Normals=patchNormal(Fb,V);%gets normal vector of all facets

logicContactPosX=Normals(:,1)==1;%where norm v is pos x 
logicContactNegX=Normals(:,1)==-1;%where norm v is neg x 

logicContactAllPrimarySets=logical(logicContactSurf.*logicContactPosX);
logicContactAllSecondarySets=logical(logicContactSurf.*logicContactNegX);

% F_contactPrimary=Fb(logicContactPrimary,:);
% F_contactSecondary=Fb(logicContactSecondary,:);


Vm=patchCentre(Fb,V);% centre location of faces, where used here should not affect x axis as contact faces are parrallel, this is used to find the faces on the x ccordinate of interest


if n>1
for q=1:1:n-1
%where the primary and secondary for each contact pair should be found
desiredPrimary=(2*Chamber_wt_elements)+Chamber_length_elements+((q-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
desiredSecondary=(2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements+((q-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));

logicContactCoordinatePrimary=Vm(:,1)==desiredPrimary;
logicContactCoordinateSecondary=Vm(:,1)==desiredSecondary;

logicContactPrimary=logical(logicContactCoordinatePrimary.*logicContactAllPrimarySets);
logicContactSecondary=logical(logicContactCoordinateSecondary.*logicContactAllSecondarySets);

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





%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,Cb,'none',1);%normally set transp to 0.5, 0.05 shows pressure face better

% hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',15);
% hl(2)=gpatch(F_pressure,V,'r','k',1);
% hl(3)=gpatch(F2,V,'c','k',1);

patchNormPlot(F_pressure,V);
%legend(hl,{'BC full support','Pressure surface', 'SLL surface'});

axisGeom(gca,fontSize);
camlight headlight; icolorbar;
gdrawnow; 

cFigure;
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
% hold on; gpatch(E{2},V,'g');
axisGeom;

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='3.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.max_ups=max_ups;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k1;

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.c1=c2;
febio_spec.Material.material{2}.m1=m2;
febio_spec.Material.material{2}.c2=c2;
febio_spec.Material.material{2}.m2=-m2;
febio_spec.Material.material{2}.k=k2;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
if SLL_control == 1 
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E1; %The element matrix
 
partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2; %The element matrix

elseif SLL_control >= 2
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E1; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='quad4'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2; %The element matrix

end



% -> Surfaces
surfaceName1='LoadedSurface';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure;

contactSurfacesValPrimary=2:2:(2*(n-1));% attribute number of contact faces
contactSurfacesValSecondary=3:2:((2*(n-1))+1);% attribute number of contact faces

if n>1
    for q=1:1:n-1
        numStringPrimary=num2str(contactSurfacesValPrimary(q));
        numStringSecondary=num2str(contactSurfacesValSecondary(q));
        numStringPair=num2str(q);
        leaderStringSurf='contactSurface';
        leaderStringContact='Contact';
              
contactSurfacesStringPrimary=strcat(leaderStringSurf,numStringPrimary);%creates string 'ContactSurfaceINT'
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.ATTR.name=contactSurfacesStringPrimary;
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.quad4.ATTR.id=(1:1:size(ContactPair.Primary{q},1))';%size(F_pressure,1)+
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.quad4.VAL=ContactPair.Primary{q};

contactSurfacesStringSecondary=strcat(leaderStringSurf,numStringSecondary);%creates string 'ContactSurfaceINT'
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.ATTR.name=contactSurfacesStringSecondary;
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.quad4.ATTR.id=(1:1:size(ContactPair.Secondary{q},1))';%size(F_pressure,1)+size(F_contactPrimary,1)+
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.quad4.VAL=ContactPair.Secondary{q};

% -> Surface pairs
febio_spec.Mesh.SurfacePair{q}.ATTR.name=strcat(leaderStringSurf,numStringPair);
febio_spec.Mesh.SurfacePair{q}.primary=contactSurfacesStringPrimary;
febio_spec.Mesh.SurfacePair{q}.secondary=contactSurfacesStringSecondary;
    end
end

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

if SLL_control == 1
febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;
elseif SLL_control >= 2
febio_spec.MeshDomains.ShellDomain{1}.ATTR.name=partName2;
febio_spec.MeshDomains.ShellDomain{1}.ATTR.mat=materialName2;
end

%MeshData secion
%-> Element data
if SLL_control >= 2
febio_spec.MeshData.ElementData{1}.ATTR.var='shell thickness';
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName2;
febio_spec.MeshData.ElementData{1}.elem.ATTR.lid=(1:1:size(E2,1))';
febio_spec.MeshData.ElementData{1}.elem.VAL=SLLThicknessShell*ones(size(E2,1),size(E2,2));
end

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.type='fix';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dofs='x,y,z';

%Loads section
% -> Surface load
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=appliedPressure;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;

% %Contact section
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
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 1];

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

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
     elseif SLL_control==2
    [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,end));
%     CV{1}=faceToVertexMeasure(E1,V,E_stress_mat(:,:,end));
%     CV{2}=faceToVertexMeasure(E2,V,E_stress_mat(:,:,end));
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








