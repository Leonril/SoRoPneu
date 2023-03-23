%%This is a testing program for the modelling of a soft pneumatic actuator
%%using Gibbon and FeBio. 

% Leon Riley - 18477946

%% Keywords
%
% * febio_spec version 3.0
% * febio, FEBio
% * pressure loading
% * hexahedral elements, hex8
% * pneunet actuator
% * soft robotic

%%


clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=5;
markerSize2=10;
lineWidth=3;
cMap=viridis(20); %colormap 

%% Control parameters
% Load
appliedPressure=-0.1;

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting force


%% Pneunet Geometry Parameters
n=3; % number of chambers

% Height Parameters
h_c=6; %height of the pneunet chambers (mm) from the base
rt=2;%roof thickness
ct=0.5;
ct_r=0.5;%roof of channel thickness
bh=2;%base height (mm)

h_nc=ct+ct_r+bh; %height of the gap between chambers (mm)

% Length Parameters
l_c=7; %length of the pneuent chamber (mm)
wt=1;%wall thickness (mm)

%Strain limiting layer parameters
SLL_thickness=0.1; %thickness of SLL

% Width Parameters
s=2; %spacing between the pneunet chambers (mm)
t=10;%width of the chambers
t_chY=2; %thickness of chambers in the y axis


% Parameter Verification

if h_c <= (rt+h_nc)
    error='Height values are not sensible'
    return
elseif l_c <= (2*wt)
    error='Chamber lenght is too small'
    return
elseif t <= (2*t_chY)
    error='Chamber width is too small'    
end

%% Pneunet Material Parameters
k_factor=50; %Bulk modulus factor 


%material 1 - Pneunet Elastomer
c1=1; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k1=c1*k_factor; %Bulk modulus


%material 2 - Strain Limiting Layer
c2=50*c1; %Shear-modulus-like parameter
m2=2; %Material parameter setting degree of non-linearity
k2=c2*k_factor; %Bulk modulus


%% FEA control settings
numTimeSteps=20; %Number of time steps desired
opt_iter=6; %Optimum number of iterations
max_refs=opt_iter*2; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps)*2; %Maximum time step size

runMode='external';%'internal';

%% Geometry Profile
%% Outer Surface
Vo=[0 h_nc; 0 h_c; l_c h_c; l_c h_nc;]; %geometry of one chamber
w=s+l_c;%x distance of one period

[VS]=copyOffset(Vo,n,w);%extends V for all chambers
VSStart=VS(1,:); VSEnd=VS(end,:);%duplicates VS start and end point 
VSStart(1,2)=0; VSEnd(1,2)=0;%moves start and end point to zero in Z axis
VO=[VSStart; VS; VSEnd];%places start and end points into VS
x=size(VO,1);%taking the number of points per side 

%FO=1:1:x;%creating one large face of one side of the pneunet


%% Pressure Surface
Vi=[wt h_nc-ct_r; wt h_c-rt; l_c-wt h_c-rt; l_c-wt h_nc-ct_r]; %geometry of one chamber

[VI]=copyOffset(Vi,n,w);%extends V for all chambers
VIStart=VI(1,:); VIEnd=VI(end,:);%duplicates VS start and end point 
VIStart(1,2)=bh; VIEnd(1,2)=bh;%moves start and end point to zero in Z axis
VI=[VIStart; VI; VIEnd];%places start and end points into VS
xI=size(VI,1);%taking the number of points per side 

%FI=1:1:xI;%creating one large face of one side of inside the pneunet

%% SLL - Strain Limiting Layer





% Vsll=[SLL_xpos (((bh/2)-(SLL_thickness/2))+zposSLL); SLL_xpos (((bh/2)+(SLL_thickness/2))+zposSLL); ((VO(end,1))-SLL_xpos) (((bh/2)+(SLL_thickness/2))+zposSLL);((VO(end,1))-SLL_xpos) (((bh/2)-(SLL_thickness/2))+zposSLL)];
%geometry of SLL in XY




%% Test draw

cFigure; hold on;
plotV(VS,'b.-','MarkerSize',markerSize,'LineWidth',lineWidth); 
plotV(VI,'b.-','MarkerSize',markerSize,'LineWidth',lineWidth); 
% plotV(Vsll,'y.-','MarkerSize',markerSize,'LineWidth',lineWidth); 

axisGeom; view(2);
drawnow; 


%% Surface Meshing
%% Outer Face
regionCell={VO(:,[1 2])};
resampleCurveOpt=1; %Option to turn on/off resampling of input boundary curves
pointSpacing=0.5;% spacing between mesh points
[FO,VO,boundaryIndO]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,0);
VO(:,3)=0;%adding Y value to V
VOb=VO(boundaryIndO,:);

%
d=t;%thickness
ns=ceil(d./pointSpacing);%number of partitions in the sub
ns=ns+iseven(ns);

cPar.depth=d; 
cPar.numSteps=ns;
cPar.patchType='tri_slash'; 
cPar.dir=1;
cPar.closeLoopOpt=1; 
[FOt,VOt]=polyExtrude(VOb,cPar);

VO2=VO;
VO2(:,3)=d;%adding the Z value to V2
[FO,VO,CO]=joinElementSets({fliplr(FO),FO,fliplr(FOt)},{VO,VO2,VOt});
[FO,VO]=mergeVertices(FO,VO);

%% Inner Surface
regionCell={VI(:,[1 2])};
pointSpacingI=0.5;% spacing between mesh points
[FI,VI,boundaryIndI]=regionTriMesh2D(regionCell,pointSpacingI,resampleCurveOpt,0);
VI(:,3)=wt;%adding Y value to V
VIb=VI(boundaryIndI,:);

%

dI=t-(2*wt);%thickness
nsI=ceil(dI./pointSpacingI);%number of partitions in the sub
nsI=nsI+iseven(nsI);

cPar.depth=dI; 
cPar.numSteps=nsI;
cPar.patchType='tri_slash'; 
cPar.dir=1;
cPar.closeLoopOpt=1; 
[FIt,VIt]=polyExtrude(VIb,cPar);

VI2=VI;
VI2(:,3)=t-wt;%adding the Z value to V2
[FI,VI,CI]=joinElementSets({fliplr(FI),FI,fliplr(FIt)},{VI,VI2,VIt});
[FI,VI]=mergeVertices(FI,VI);



%% SLL
%pointSpacingSLL=0.5;% spacing between mesh points

 

%% Colour ID of Fixed Face
[CO] = XYZColourDesignate(FO,VO,CO,2,0,0); %custom function to set Spine face
[CO] = XYZColourDesignate(FO,VO,CO,1,0,0);%custom function to set BC face colour

%% Colour IDs
% CO=CO+2;%shifting all colour IDs up two
CItemp=CI+max(CO);

cFigure; hold on; 
gpatch(FO,VO,CO);
patchNormPlot(FO,VO);
colormap spectral; icolorbar;
axisGeom; camlight headlight;
drawnow; 


%% SLL
Fsll_top=FO(ismember(CO,2),:);%gets spine face

[Fsll_top, Vsll_top]=patchCleanUnused(Fsll_top,VO);%extracts F and V for spine 


[~,~,dirVec]=patchNormal(Fsll_top,Vsll_top);
dirVec=vecnormalize(dirVec);
[E_SLL,V_SLL]=patchThick(Fsll_top,Vsll_top,dirVec,SLL_thickness,1);
%elements for SLL

CSLL_temp=(1:1:size(E_SLL,1))';
[FSLLSides,CSLLfaces,CSLLsides]=element2patch(E_SLL,CSLL_temp,'penta6');

%Sorting Boundary Colours for SLL
CSLLside2=CSLLsides{1,2};
FSLLside2=FSLLSides{1,2};
% 
% sizeFSSLside2=size(FSLLside2,1);
% % for tempCLLcounter=1:1:sizeFSSLside2
    [CSLLside2] = XYZColourDesignate(FSLLside2,V_SLL,CSLLside2,1,min(VO(:,1)),0);%custom function to set BC face colour
    [CSLLside2] = XYZColourDesignate(FSLLside2,V_SLL,CSLLside2,1,max(VO(:,1)),0);
    [CSLLside2] = XYZColourDesignate(FSLLside2,V_SLL,CSLLside2,3,min(VO(:,3)),0);
    [CSLLside2] = XYZColourDesignate(FSLLside2,V_SLL,CSLLside2,3,max(VO(:,3)),0);

for tempCLLcounter=1:1:size(CSLLside2,1)
   if CSLLside2(tempCLLcounter) > 4 
       CSLLside2(tempCLLcounter)=5;
   end
end
CSLLsides{1,2}=CSLLside2+2;
cFigure;
gpatch(FSLLSides,V_SLL,CSLLsides,'k',0.9); hold on;
axisGeom;icolorbar;

CSLLsides{1,2}=CSLLsides{1,2}+max(CItemp);
CSLLsides{1,1}=CSLLsides{1,1}+max(CItemp);
% Plots to Inv Colours
% cFigure
% gpatch(Fsll_top,Vsll_top); hold on;
% patchNormPlot(Fsll_top,Vsll_top);
% axisGeom;
% 
% cFigure
% gpatch(E_SLL,V_SLL); hold on;
% patchNormPlot(E_SLL,V_SLL);
% axisGeom;


% cFigure;
% subplot(1,2,1); title('Element colors');
% gpatch(FSLLSides,V_SLL,CSLLfaces,'k',0.5);
% axisGeom(gca,fontSize);
% camlight headlight;
% colormap(gca,gjet(250)); colorbar;
% 
% subplot(1,2,2); title('Face side type colors');
% gpatch(FSLLSides,V_SLL,CSLLsides,'k',0.9);
% 
% axisGeom(gca,fontSize);
% camlight headlight;
% colormap(gca,gjet(6)); icolorbar;
% 
% drawnow;
%

%% Defining regions
% Joining surface sets
[F,V,C]=joinElementSets({FO,FI,FSLLSides},{VO,VI,V_SLL},{CO,CItemp,CSLLsides});
tttt

%% could I just manually add the the other sides in {1,1} of the penta?


cFigure; hold on; 
gpatch(F,V,C);
colormap spectral; icolorbar;
axisGeom; camlight headlight;
drawnow; 




ggggg
%Finding interior points
[V_region1]=getInnerPoint({FO,FI},{VO,VI}); 
[V_region2]=getInnerPoint({FSLL},{VSLL}); 

[V_holes]=getInnerPoint(FI,VI); 

V_regions=[V_region1; V_region2];

%Volume parameters

% Volume parameters
[vol1]=tetVolMeanEst(FO,VO);
[vol2]=tetVolMeanEst(FSLL,VSLL);

regionTetVolumes=[vol1 vol2];
stringOpt='-pq1.2AaY'; %Tetgen options

%% Testing the faces of both the inside and outside surfaces
cFigure; hold on; 
gpatch(FO,VO,CO);
patchNormPlot(FO,VO);
colormap spectral; icolorbar;
axisGeom; camlight headlight;
drawnow; 


cFigure; hold on; 
gpatch(FI,VI,CI);
patchNormPlot(FI,VI);
colormap spectral; icolorbar;
axisGeom; camlight headlight;
drawnow; 

cFigure; hold on; 
gpatch(FSLL,VSLL,CSLL,'k');
patchNormPlot(FSLL,VSLL);
colormap spectral; icolorbar;
axisGeom; camlight headlight;
view(2);
drawnow; 


cFigure; hold on; 
gpatch(F,V,C,'k',0.5);
patchNormPlot(F,V);
colormap spectral; icolorbar;
axisGeom; camlight headlight;
drawnow; 

gjhghhg


%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% disc. These labels can be used to define boundary conditions. 

%Define supported node sets
bcSupportList=unique(F(C==1,:)); %Node set part of selected face

%Get pressure faces
F_pressure=F(ismember(C,[5 6 7]),:); 

%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(F,V,C,'none',0.5);%normally set transp to 0.5, 0.05 shows pressure face better

hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',15);
hl(2)=gpatch(F_pressure,V,'r','k',1);

patchNormPlot(F_pressure,V);
legend(hl,{'BC full support','Pressure surface'});

axisGeom(gca,fontSize);
camlight headlight; icolorbar;
gdrawnow; 


%%
% Mesh using TetGen

%Create tetgen input structure
inputStruct.stringOpt=stringOpt; %Tetgen options
inputStruct.Faces=F; %Boundary faces
inputStruct.Nodes=V; %Nodes of boundary
inputStruct.faceBoundaryMarker=C; 
inputStruct.regionPoints=V_regions; %Interior points for regions
inputStruct.holePoints=V_holes; %Interior points for holes
inputStruct.regionA=regionTetVolumes; %Desired tetrahedral volume for each region

% Mesh model using tetrahedral elements using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access mesh output structure

E=meshOutput.elements; %The elements
V=meshOutput.nodes; %The vertices or nodes
CE=meshOutput.elementMaterialID; %Element material or region id
Fb=meshOutput.facesBoundary; %The boundary faces
Cb=meshOutput.boundaryMarker; %The boundary markers

%% Output of mesh results
% Visualization

hf=cFigure; 
subplot(1,2,1); hold on;
title('Input boundaries','FontSize',fontSize);
hp(1)=gpatch(Fb,V,Cb,'k',faceAlpha1);
hp(2)=plotV(V_regions,'r.','MarkerSize',5);
legend(hp,{'Input mesh','Interior point(s)'},'Location','NorthWestOutside');
axisGeom(gca,fontSize); camlight headlight;
colormap(cMap); icolorbar;

hs=subplot(1,2,2); hold on;
title('Tetrahedral mesh','FontSize',fontSize);

% Visualizing using |meshView|
optionStruct.hFig=[hf,hs];
meshView(meshOutput,optionStruct);

axisGeom(gca,fontSize); 
gdrawnow;


%% Splitting element sets

E1=E(meshOutput.elementMaterialID==-3,:);
E2=E(meshOutput.elementMaterialID==-2,:);

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
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E1; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='tet4'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2; %The element matrix

% -> Surfaces
surfaceName1='LoadedSurface';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.tri3.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL=F_pressure;

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;

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
 %febView(febio_spec); %Viewing the febio file|

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
