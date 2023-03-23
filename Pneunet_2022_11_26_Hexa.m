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

% SLL control 
SLL_control=1;%1 for thick bottom layer (3D print), 2 for embedded shell (moulding) 

%SLL parameters for method 2 
SLL_ZHeight=0.9;%(mm)
SLLThicknessShell=0.1; %(mm)

% Auto Mesh Size -  if not active, manually specify mesh sizes further below
AutoMeshSize=1; % if 1, will auto choose the biggest size of mesh for XYZ in steps of 1
AutoBaseMaxMesh=1; %(mm) Max size for auto mesh
AutoBaseMinMesh=0.1; %(mm) Min size for auto mesh
AutoMeshIncrement=0.1; %(mm) smallest mesh increment for auto mesh

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


%% Auto Mesh Sizes

if AutoMeshSize ==1 
    XSizes=[Chamber_length Chamber_wt Gap_length];
    YSizes=[Side_wall_thickness Chamber_width (Side_wall_thickness+(Chamber_width/2)-(Channel_width)) (Side_wall_thickness+(Chamber_width/2)+(Channel_width)) Channel_width];
    ZSizes=[SLL_thickness Mat1_base_thickness Channel_height Channel_roof_thickness Chamber_height Chamber_roof_thickness];
    StepX=AutoBaseMaxMesh;
    StepY=AutoBaseMaxMesh;
    StepZ=AutoBaseMaxMesh;
    AutoSizeMeshX=1;
    AutoSizeMeshY=1;
    AutoSizeMeshZ=1;
else
    StepX=0.5;
    StepY=1;
    StepZ=0.5;
end

TestCounter=0;
for steps=AutoBaseMinMesh:AutoMeshIncrement:AutoBaseMaxMesh
    for xx=1:1:size(XSizes,2)
       TestSizeX=mod(XSizes(xx),steps);
       if TestSizeX < (AutoBaseMinMesh/2)
           TestCounter=TestCounter+1;
       end
       if TestCounter == size(XSizes,2)
           StepX=steps;
       end
    end
    TestCounter=0;
    for yy=1:1:size(YSizes,2)
       TestSizeY=mod(YSizes(yy),steps);
       if TestSizeY < (AutoBaseMinMesh/2)
           TestCounter=TestCounter+1;
       end
       if TestCounter == size(YSizes,2)
           StepY=steps;
       end
    end
    TestCounter=0;
    for zz=1:1:size(ZSizes,2)
       TestSizeZ=mod(ZSizes(zz),steps);
       if TestSizeZ < (AutoBaseMinMesh/2)
           TestCounter=TestCounter+1;
       end
       if TestCounter == size(ZSizes,2)
           StepZ=steps;
       end
    end
end


%% Geometry calculations
% X direction
Total_length=((n)*((2*Chamber_wt)+Chamber_length))+((n-1)*(Gap_length));


% Y direction
Total_width=(2*Side_wall_thickness)+Chamber_width;


% Z direction
Total_height=SLL_thickness+Mat1_base_thickness+Chamber_height+Chamber_roof_thickness;

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


%% OUTER GEOMETRY
CX1=VE_bar(:,1);
CX=abs(CX1);
CZCh=VE_bar(:,3);
CZCh=CZCh-min(CZCh)+(StepZ/2);
imax=size(CX,1);
upperlimZCh=SLL_thickness+Mat1_base_thickness+Channel_height+Channel_roof_thickness;

evenTest=mod(n,2); %check if n chambers is even or odd
jmax=n+1;
for i=1:1:imax%checks every element for whether it should exist, 
CXtemp=0;
CZtemp=0;
if evenTest ==0
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
elseif evenTest == 1
        lowerlimX=(Chamber_length/2)+Chamber_wt;
        upperlimX=(Chamber_length/2)+Chamber_wt+Gap_length;
    for j=1:1:jmax
        if CX(i) >= lowerlimX && CX(i) <= upperlimX
            CXtemp=1;
        elseif CZCh(i)<= upperlimZCh
            CZtemp=1;
        end
        lowerlimX=(Chamber_length/2)+Chamber_wt+(j*(Gap_length+(2*Chamber_wt)+Chamber_length));
        upperlimX=(Chamber_length/2)+Chamber_wt+Gap_length+(j*(Gap_length+(2*Chamber_wt)+Chamber_length));
    end
end
    if CXtemp == 1
        CX(i)=0; %deletes unneccessary elements
    end
    if CZtemp==1
        CZCh(i)=0;
    end
end

%% INTERNAL SURFACE
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


for i=1:1:imax
CXtemp=0;
CYtemp=0;
CZtemp=0;
C_Channel_temp=0;
if evenTest == 0
lowerlimX=(Gap_length/2)+Chamber_wt;
upperlimX=(Gap_length/2)+Chamber_wt+Chamber_length;
    for j=1:1:jmax
        if CXI(i) >= lowerlimX && CXI(i) <= upperlimX
            CXtemp=1;
        elseif CXI(i) >= lowerlimX_Chann && CXI(i) <= upperlimX_Chann
            C_Channel_temp=1;
        end
        lowerlimX=(Gap_length/2)+Chamber_wt+(j*(Gap_length+(2*Chamber_wt)+Chamber_length));
        upperlimX=(Gap_length/2)+Chamber_wt+Chamber_length+(j*(Gap_length+(2*Chamber_wt)+Chamber_length));
    end
elseif evenTest == 1
lowerlimX=0;
upperlimX=Chamber_length/2;
    for j=1:1:jmax
        if CXI(i) >= lowerlimX && CXI(i) <= upperlimX
            CXtemp=1;
        elseif CXI(i) >= lowerlimX_Chann && CXI(i) <= upperlimX_Chann
            C_Channel_temp=1;
        end
        lowerlimX=(0.5*Chamber_length)+((j-1)*(Chamber_length))+(j*((2*Chamber_wt)+Gap_length));
        upperlimX=lowerlimX+Chamber_length;
    end
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





%%


logicKeep1=~(CX==0 & CZCh>0); %Outer surface logic
logicKeep3=~(CXI ==0 & CYI==0 & CZI==0  ); %Chamber logic
logicKeep4=~(CXIChannel==0 & CYIChannel==0 & CZIChannel==0); %Channel logic

logicKeep2=logicKeep3.*logicKeep4;
logicKeep2=logicKeep2(logicKeep1,:);
logicKeep2=logical(logicKeep2);


E1=E_bar(logicKeep1,:); %removes outer elements
F1=element2patch(E1);
[indBoundary1]=tesBoundary(F1);

F2=element2patch(E1(logicKeep2,:));%removes internal elements
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

%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% disc. These labels can be used to define boundary conditions. 

%Define supported node sets
bcSupportList=unique(Fb(Cb==1,:)); %Node set part of selected face

%Get pressure faces
F_pressure=Fb(Cb==0,:); 

%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'w','none',0.5);

hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
hl(2)=gpatch(F_pressure,V,'r','k',1);

patchNormPlot(F_pressure,V);
legend(hl,{'BC full support','Pressure surface'});

axisGeom(gca,fontSize);
camlight headlight; 
gdrawnow; 

%% Update this to take all SLL layer


%
if SLL_control == 1 % thick second elestomer bottom layer
    VSLL=patchCentre(E,V);
    CZSLL=VSLL(:,3);
    CZSLL=CZSLL-min(CZSLL)+(StepZ/2);
    logicSLL=CZSLL<SLL_thickness;    
    
    E1=E(~logicSLL,:); %Other elements
    E2=E(logicSLL,:); %SLL element layer
    E=[E1;E2];

    [F1]=element2patch(E1);
    [F2]=element2patch(E2);
    
    cFigure; hold on; 
    gpatch(F1,V,'bw','k',0.5);
    gpatch(F2,V,'rw','k',0.5);
    %plotV(VSLL2(indicesSLL,:),'k.','MarkerSize',markerSize);

    axisGeom; 
    colormap(turbo(250)); icolorbar; 
    camlight headlight; 
    gdrawnow; 

elseif SLL_control == 2 % thin shell element embedded
    logicSLL2=V(:,3)== min(V(:,3));
    VSLL2=V(logicSLL2,:);
    VSLL2(:,3)= min(V(:,3))+SLL_ZHeight;
    logic_FSLL2=Cb==5;
    FSLL2=Fb(logic_FSLL2,:);
    FSLL2=FSLL2+max(Fb(:));
    V=[V; VSLL2];
    
%     FSLL2=[];
%     for i=1:1:size(F,1)
%     m=0;
%     for j=1:1:size(F,2)
%         if V(F(i,j),3) == (min(V(:,3))+SLL_thickness)
%                 m=m+1;   
%         end
%          if m == size(F,2)
%            FSLL2=[FSLL2; F(i,:)];
%         end
%     end
%    
%     end
%     FSLL2=unique(FSLL2, 'rows');

    cFigure; hold on; 
    gpatch(F,V,'bw','k',0);
    gpatch(FSLL2,V,'bw','k',1);
    axisGeom; 
    colormap(turbo(250)); icolorbar; 
    camlight headlight; 
    gdrawnow; 
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

elseif SLL_control == 2
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='quad4'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E,1)+(1:1:size(FSLL2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=FSLL2; %The element matrix
end

% -> Surfaces
surfaceName1='LoadedSurface';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure;

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
elseif SLL_control == 2
febio_spec.MeshDomains.ShellDomain{1}.ATTR.name=partName2;
febio_spec.MeshDomains.ShellDomain{1}.ATTR.mat=materialName2;
end

%MeshData secion
%-> Element data
if SLL_control == 2
febio_spec.MeshData.ElementData{1}.ATTR.var='shell thickness';
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName2;
febio_spec.MeshData.ElementData{1}.elem.ATTR.lid=(1:1:size(FSLL2,1))';
febio_spec.MeshData.ElementData{1}.elem.VAL=SLLThicknessShell*ones(size(FSLL2,1),size(FSLL2,2));
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
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

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



