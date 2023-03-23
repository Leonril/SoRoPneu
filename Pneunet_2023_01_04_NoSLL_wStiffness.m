%% Due to finding issues in SLL generation, a no SLL version was created
% this allows:
%
% * A cleaner generation of a one material system 
%
% This route was not used previously due to
%
% * Decreases the robustness of the analysis through the investigation of
% less materials
%
%

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

runMode='external';%'internal';

%Contact parameters
contactInitialOffset=0.1;
contactPenalty=5;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0;


%% Load Inputs
%Load
appliedPressure=0.15; 

%Define applied displacement perturbation
prescribedDisplacement_X=2;

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

n=3; %no. of chambers 

% X direction 

Chamber_length=8; % internal chamber length (mm)
Chamber_wt=1.2; % wall thickness expanding (mm)
Gap_length=3; %(mm)

% y direction 

Side_wall_thickness=2;
Chamber_width=8; % internal chamber width (mm)
Channel_width=2; % width of the intenal channel (mm)

% Z direction 

Mat1_base_thickness=1.5; %mat 1 basse layer thickness (mm)
Channel_height=1; %height of internal channel between chambers (mm)
Channel_roof_thickness=1;% (mm)
Chamber_roof_thickness=2.5; % (mm)
Chamber_height=10; %internal height of the chamber (mm)

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

Mat1_base_thickness_elements=2; %mat 1 basse layer thickness (mm)
Channel_height_elements=2; %height of internal channel between chambers (mm)
Channel_roof_thickness_elements=2;% (mm)
Chamber_roof_thickness_elements=3; % (mm)
Chamber_height_elements=12; % no of mesh elements on this length (Inclusive of channel parameters)


%% Simplfied geometry creation
% Simple geometry is defined using mesh elements only

% X direction elements needed
Length=(n*((2*Chamber_wt_elements)+Chamber_length_elements))+((n-1)*Gap_length_elements);

% Y direction elements needed
Width=(2*Side_wall_thickness_elements)+Chamber_width_elements;

% Z direction elements needed
Height=Mat1_base_thickness_elements+Chamber_height_elements+Chamber_roof_thickness_elements;

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
logicDeleteOuterZ=VE_bar(:,3)>(Mat1_base_thickness_elements+Channel_height_elements+Channel_roof_thickness_elements);
logicKeepOuter=~(logicDeleteOuterX.*logicDeleteOuterZ);
logicKeepOuter=logical(logicKeepOuter);

logicDeleteInnerY1=~(VE_bar(:,2)>(Chamber_width_elements/2));
logicDeleteInnerY2=~(VE_bar(:,2)>(Channel_width_elements/2));
logicDeleteInnerZ1=(VE_bar(:,3)>(Mat1_base_thickness_elements) & VE_bar(:,3)<(Height-Chamber_roof_thickness_elements));
logicDeleteInnerZ2=(VE_bar(:,3)>(Mat1_base_thickness_elements) & VE_bar(:,3)<(Mat1_base_thickness_elements+Channel_height_elements));



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
hold on; plotV(V,'k.','MarkerSize',markerSize/2);
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
logicContactSurf=Cb==7;
F_contact=Fb(logicContactSurf,:);

%Get end face for stiffness analysis
bcPrescribeList=unique(Fb(Cb==2,:));%Node set for displacement test


%% Scaling to match desired geometry
if MeshCreationTool == 1 %'Must points' method
%Creating vectors of needed Coordinates
V_height=zeros(1,Height); %creating empty vectors to reduce memory
V_length=zeros(1,Length);
V_width=zeros(1,Width);

%Height scaling vector creation
current_height=0;
for i=1:1:Height
    
    if i <= Mat1_base_thickness_elements
        V_height(i)=current_height+(Mat1_base_thickness*(1/Mat1_base_thickness_elements));
        
    elseif i <= Channel_height_elements+Mat1_base_thickness_elements
        V_height(i)=current_height+(Channel_height*(1/Channel_height_elements));
        
    elseif i <= Channel_roof_thickness_elements+Channel_height_elements+Mat1_base_thickness_elements 
        V_height(i)=current_height+(Channel_roof_thickness*(1/Channel_roof_thickness_elements));
        
    elseif i <= Chamber_height_elements+Mat1_base_thickness_elements
        V_height(i)=current_height+((Chamber_height-Channel_height-Channel_roof_thickness)*(1/(Chamber_height_elements-Channel_height_elements-Channel_roof_thickness_elements)));
        
    elseif  i <= Chamber_roof_thickness_elements+Chamber_height_elements+Mat1_base_thickness_elements
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



%% Interpolation /  Smoothing methods
elseif MeshCreationTool > 1 % 
V_length=zeros(1,Length);
    
% Simplistically scales entire element
V(:,2)=V(:,2)-min(V(:,2));%undoes centre alignment on Y axis
total_chamber_length_elements=((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements);
outer_chamber_length=((2*Chamber_wt)+Chamber_length);
current_length=0;
for i=1:1:Length
    itemp=i;
    total_chamber_length_elements=((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements);
    chamber_count=floor(itemp/total_chamber_length_elements);
    itemp=itemp-(chamber_count*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
   
        if itemp == 0
            V_length(i)=current_length+(Gap_length*(1/Gap_length_elements));
            
        elseif itemp <= (2*Chamber_wt_elements)+Chamber_length_elements
            V_length(i)=current_length+(outer_chamber_length*(1/(total_chamber_length_elements-Gap_length_elements)));

        elseif itemp <= Gap_length_elements+(2*Chamber_wt_elements)+Chamber_length_elements
            V_length(i)=current_length+(Gap_length*(1/Gap_length_elements));

        end
    
    
    current_length=V_length(i);
end
V_length=[0 V_length];



for i=1:1:size(V,1)
    x_current=V(i,1);
    y_current=V(i,2);
    z_current=V(i,3);
    
    v_indexX=V(i,1)+1;%takes the element coordinate as an index to true value vectors
    VXreal=V_length(v_indexX);
    
    
    V(i,1)=VXreal;
    V(i,2)=y_current*(((2*Side_wall_thickness)+Chamber_width)/(Width));
    V(i,3)=z_current*((Mat1_base_thickness+Chamber_height+Chamber_roof_thickness)/Height);

end

%% Create node sets for interpolation
indFs=unique(Fb(Cb==0,:));%index of pressure face
indFb=unique(Fb(Cb~=0,:));%index of outside faces that we do not want to move in warping

indAll=1:1:size(V,1); %All node numbers
indInterp=indAll(~ismember(indAll,indFs) & ~ismember(indAll,indFb));

V_internal=V_element(indFs,:);%nodes of pressure faces that need to be moved to correct position



%% Rudimentary scaling
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

V_internal(:,2)=V_internal(:,2)-min(V_internal(:,2));%moves everything to base 0 for easier control
V_internal(:,3)=V_internal(:,3)-min(V_internal(:,3));%moves everything to base 0 for easier control

for i=1:1:size(V_internal,1)
    x_current=V_internal(i,1);
    y_current=V_internal(i,2);
    z_current=V_internal(i,3);
    
    v_indexX=V_internal(i,1)+1;%takes the element coordinate as an index to true value vectors
    VXreal=V_length(v_indexX);

    V_internal(i,1)=VXreal;
    
    if V_internal(i,2)<= ((Chamber_width_elements-Channel_width_elements)/2)
        V_internal(i,2)=Side_wall_thickness+(y_current*(((Chamber_width-Channel_width)/2)/((Chamber_width_elements-Channel_width_elements)/2)));
        
    elseif V_internal(i,2)<= ((Chamber_width_elements+Channel_width_elements)/2)
        V_internal(i,2)=Side_wall_thickness+((Chamber_width-Channel_width)/2)+((y_current-((Chamber_width_elements-Channel_width_elements)/2))*(Channel_width/Channel_width_elements));
        
    elseif V_internal(i,2)<= (Chamber_width_elements)
        V_internal(i,2)=Side_wall_thickness+((Chamber_width+Channel_width)/2)+((y_current-((Chamber_width_elements+Channel_width_elements)/2))*(((Chamber_width-Channel_width)/2)/((Chamber_width_elements-Channel_width_elements)/2)));
        
    end
    
    if V_internal(i,3)<= Channel_height_elements
        V_internal(i,3)=(Mat1_base_thickness)+(z_current*(Channel_height/Channel_height_elements));
        
    elseif V_internal(i,3) <= Chamber_height_elements
        V_internal(i,3)=(Mat1_base_thickness+Channel_height)+((z_current-Channel_height_elements)*((Chamber_height-Channel_height)/(Chamber_height_elements-Channel_height_elements)));
        
    end
    
    
end
end
if MeshCreationTool == 2
%% Derive desired nodal displacement
Vp=V;%Copy of V for proposed position of internal faces
Vp(indFs,:)=V_internal;%applies correct pressure face positions to rudimentary scaled Vp
Up=Vp-V; %Compute the displacement that occured

%% Distribute displacement effect using interpolation

V_interp=[V(indFs,:); V(indFb,:)]; %Coordinates used in interpolation
Ux_interp=[Up(indFs,1); zeros(numel(indFb),1)]; %"values" at interpolation points
Uy_interp=[Up(indFs,2); zeros(numel(indFb),1)]; %"values" at interpolation points
Uz_interp=[Up(indFs,3); zeros(numel(indFb),1)]; %"values" at interpolation points



%Set-up interpolation function
interpFuncX = scatteredInterpolant(V_interp,Ux_interp,'natural');
interpFuncY = scatteredInterpolant(V_interp,Uy_interp,'natural');
interpFuncZ = scatteredInterpolant(V_interp,Uz_interp,'natural');

%Do interpolation to get a displacement metric at all nodes
Ux_indInterp=interpFuncX(V(indInterp,:));
Uy_indInterp=interpFuncY(V(indInterp,:));
Uz_indInterp=interpFuncZ(V(indInterp,:));

Umag_indInterp=abs(Ux_indInterp)+abs(Uy_indInterp)+abs(Uz_indInterp);

%Apply displacement
V(indInterp,1)=V(indInterp,1)+Ux_indInterp;
V(indInterp,2)=V(indInterp,2)+Uy_indInterp;
V(indInterp,3)=V(indInterp,3)+Uz_indInterp; %Shift non-prescribed
V(indFs,:)=V(indFs,:)+Up(indFs,:); %Shift desired set by desired amount

% Plotting model
cFigure;
title('Box boundaries faces','FontSize',fontSize);
hold on;

gpatch(F,V,'bw','b',0.1,1);
gpatch(F_pressure,V,'rw','k',1);
scatterV(V(indInterp,:),100,Umag_indInterp,'filled')

axisGeom(gca,fontSize); 
camlight headlight; colormap spectral; colorbar;
drawnow; 


cFigure; %Displays scaled geometry of two cells
title('inside face','FontSize',fontSize);
gpatch(F_pressure,V,'','k',0.5);%0.5 transperancy normally hold on
% hold on; plotV(V,'k.','MarkerSize',markerSize/2);
axisGeom; 
colormap(turbo(250)); icolorbar; 
% camlight headlight; 
gdrawnow; 


elseif MeshCreationTool ==3 % Smoothing method

V(indFs,:)=V_internal;%applies correct pressure face positions to rudimentary scaled V
V_presmooth=V;%copys V for comparison with smoothed V
%Smoothing settings
cPar.Method='LAP'; %Laplacian smooth method
cPar.n=500; %Max. number of smoothing iterations
cPar.RigidConstraints=unique([indFs(:);indFb(:)]); %Points not to move during smoothing
cPar.Tolerance=1e-2; %Tolerance (sum of squared coordinate distances based) to stop before max. number of iterations

V=tesSmooth(F,V,[],cPar); %Smooth coordinates

U=V-V_presmooth;%displacement of each node
Umag=abs(U(indInterp,1))+abs(U(indInterp,2))+abs(U(indInterp,3));%magnitude of displacement for display


% Plotting model
cFigure;
title('Box boundaries faces','FontSize',fontSize);
hold on;

gpatch(F,V,'bw','b',0.1,1);
gpatch(F_pressure,V,'rw','k',1);
scatterV(V(indInterp,:),100,Umag,'filled')
axisGeom(gca,fontSize); 
camlight headlight; colormap spectral; colorbar;
drawnow; 

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

hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',15);
hl(2)=gpatch(F_pressure,V,'r','k',1);

patchNormPlot(F_pressure,V);
legend(hl,{'BC full support','Pressure surface'});

axisGeom(gca,fontSize);
camlight headlight; icolorbar;
gdrawnow; 


%% Displacement surface elements
logicTopSurface=Cb==2;
F_top=Fb(logicTopSurface,:);
center_of_mass=mean(V(unique(F_top(:)),:),1);

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
stepStruct.Control.time_stepper.dtmin=dtmin;
stepStruct.Control.time_stepper.dtmax=dtmax; 
stepStruct.Control.time_stepper.max_retries=max_retries;
stepStruct.Control.time_stepper.opt_iter=opt_iter;

%Add template based default settings to proposed control section
[stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

febio_spec.Step.step{1}.Control=stepStruct.Control;
febio_spec.Step.step{1}.ATTR.id=1;
febio_spec.Step.step{2}.Control=stepStruct.Control;
febio_spec.Step.step{2}.ATTR.id=2;
febio_spec.Step.step{3}.Control=stepStruct.Control;
febio_spec.Step.step{3}.ATTR.id=3;
febio_spec.Step.step{4}.Control=stepStruct.Control;
febio_spec.Step.step{4}.ATTR.id=4;

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

materialName3='Material3';
febio_spec.Material.material{3}.ATTR.name=materialName3;
febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1;
febio_spec.Material.material{3}.center_of_mass=center_of_mass;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='quad4'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E,1)+(1:1:size(F_top,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=F_top; %The element matrix


% -> Surfaces
surfaceName1='LoadedSurface1';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure;

surfaceName2='contactSurface1';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=(1:1:size(F_contact,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL=F_contact;

surfaceName3='contactSurface2';
febio_spec.Mesh.Surface{3}.ATTR.name=surfaceName3;
febio_spec.Mesh.Surface{3}.quad4.ATTR.id=(1:1:size(F_contact,1))';
febio_spec.Mesh.Surface{3}.quad4.VAL=F_contact;


% -> Surface pairs
febio_spec.Mesh.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Mesh.SurfacePair{1}.primary=surfaceName2;
febio_spec.Mesh.SurfacePair{1}.secondary=surfaceName3;

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

febio_spec.MeshDomains.ShellDomain.ATTR.name=partName2;
febio_spec.MeshDomains.ShellDomain.ATTR.mat=materialName3;

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

%Contact section
febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{1}.ATTR.name;
febio_spec.Contact.contact{1}.two_pass=0;
febio_spec.Contact.contact{1}.laugon=laugon;
febio_spec.Contact.contact{1}.tolerance=0.2;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.minaug=minaug;
febio_spec.Contact.contact{1}.maxaug=maxaug;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

%Rigid section 
% ->Rigid body fix boundary conditions
febio_spec.Rigid.rigid_constraint{1}.ATTR.name='RigidFix_1';
febio_spec.Rigid.rigid_constraint{1}.ATTR.type='fix';
febio_spec.Rigid.rigid_constraint{1}.rb=3;
febio_spec.Rigid.rigid_constraint{1}.dofs='Ry';

% ->Rigid body prescribe boundary conditions
febio_spec.Step.step{2}.Rigid.rigid_constraint{1}.ATTR.name='RigidPrescribe';
febio_spec.Step.step{2}.Rigid.rigid_constraint{1}.ATTR.type='prescribe';
febio_spec.Step.step{2}.Rigid.rigid_constraint{1}.rb=3;
febio_spec.Step.step{2}.Rigid.rigid_constraint{1}.dof='Rx';
febio_spec.Step.step{2}.Rigid.rigid_constraint{1}.value.ATTR.lc=2;
febio_spec.Step.step{2}.Rigid.rigid_constraint{1}.value.VAL=prescribedDisplacement_X;
febio_spec.Step.step{2}.Rigid.rigid_constraint{1}.relative=1;

febio_spec.Step.step{3}.Rigid.rigid_constraint{1}.ATTR.name='RigidPrescribe';
febio_spec.Step.step{3}.Rigid.rigid_constraint{1}.ATTR.type='prescribe';
febio_spec.Step.step{3}.Rigid.rigid_constraint{1}.rb=3;
febio_spec.Step.step{3}.Rigid.rigid_constraint{1}.dof='Rx';
febio_spec.Step.step{3}.Rigid.rigid_constraint{1}.value.ATTR.lc=3;
febio_spec.Step.step{3}.Rigid.rigid_constraint{1}.value.VAL=-prescribedDisplacement_X;
febio_spec.Step.step{3}.Rigid.rigid_constraint{1}.relative=1;

febio_spec.Step.step{4}.Rigid.rigid_constraint{1}.ATTR.name='RigidPrescribe';
febio_spec.Step.step{4}.Rigid.rigid_constraint{1}.ATTR.type='prescribe';
febio_spec.Step.step{4}.Rigid.rigid_constraint{1}.rb=3;
febio_spec.Step.step{4}.Rigid.rigid_constraint{1}.dof='Rx';
febio_spec.Step.step{4}.Rigid.rigid_constraint{1}.value.ATTR.lc=4;
febio_spec.Step.step{4}.Rigid.rigid_constraint{1}.value.VAL=-prescribedDisplacement_X;
febio_spec.Step.step{4}.Rigid.rigid_constraint{1}.relative=1;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 1; 2 1; 3 1; 4 1];

febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{2}.points.point.VAL=[0 0; 1 0; 2 1;];

febio_spec.LoadData.load_controller{3}.ATTR.id=3;
febio_spec.LoadData.load_controller{3}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{3}.points.point.VAL=[0 0; 1 0; 2 0; 3 1];

febio_spec.LoadData.load_controller{4}.ATTR.id=4;
febio_spec.LoadData.load_controller{4}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{4}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{4}.points.point.VAL=[0 0; 1 0; 2 0; 3 0; 4 1];

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
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

febio_spec.Output.logfile.rigid_body_data{1}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.data='Fx';
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.delim=',';

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
    % Importing nodal displacements from a log file
    dataStructForce=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),1,1);
    timeData=dataStructForce.time(:);
    forceData=dataStructForce.data(:);
    
    logicPush=timeData>=1 & timeData<=2;
    timePush=timeData(logicPush);
    timePush=timePush-1;
    forcePush=forceData(logicPush);
    dispPush=timePush.*prescribedDisplacement_X;
    
    logicPull=timeData>=3 & timeData<=4;
    timePull=timeData(logicPull);
    timePull=timePull-3;
    forcePull=forceData(logicPull);
    dispPull=timePull.*-prescribedDisplacement_X;
    
    Ux=[flipud(dispPull(2:end)); dispPush(2:end)];
    Fx=[flipud(forcePull(2:end)); forcePush(2:end)];
    
    dF=diff(Fx);
    du=diff(Ux);
    S_diff=dF./du;
    u_diff=Ux(1:end-1)+du/2;
            
    ui=linspace(-prescribedDisplacement_X,prescribedDisplacement_X,100);
    S=interp1(u_diff,S_diff,ui,'linear','extrap');

    %%
    
    cFigure; hold on; 
    xlabel('U_x [mm]'); ylabel('F_x [N]');
    hp1=plot(dispPull,forcePull,'b.-','LineWidth',3,'MarkerSize',25);
    hp2=plot(dispPush,forcePush,'r.-','LineWidth',3,'MarkerSize',25);
    legend([hp1 hp2],{'Pull','Push'},'Location','NorthEastOutside')
    set(gca,'FontSize',fontSize);
    axis tight; axis square; grid on; box on;
    drawnow
    
    %%
    
    cFigure; hold on; 
    xlabel('U_x [mm]'); ylabel('S_x [N/mm]');
    hp1=plot(u_diff,S_diff,'k.','MarkerSize',50);
    hp2=plot(ui,S,'k-','LineWidth',3);
    legend([hp1 hp2],{'FEA','Interpolated'},'Location','NorthEastOutside')
    set(gca,'FontSize',fontSize);
    axis tight; axis square; grid on; box on;
    drawnow
    
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
    
    [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,end));
    
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
        
        [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
    
end









