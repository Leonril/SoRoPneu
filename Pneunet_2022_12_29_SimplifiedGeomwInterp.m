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

%Mesh Tool
MeshCreationTool=1;%
% == 1: Uses 'must points' at key variable locations
% == 2: Uses some 'must points' -> uses interpolation to scale internal
% nodes
% == 3: Uses some 'must points' -> uses smoothing to scale internal nodes


%% Geometry Inputs:

n=2; %no. of chambers 

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

Chamber_length_elements=4; % no of mesh elements on this length
Chamber_wt_elements=1; % no of mesh elements on this length
Gap_length_elements=1; % no of mesh elements on this length

% y direction 

Side_wall_thickness_elements=1; % no of mesh elements on this length
Chamber_width_elements=4; % no of mesh elements on this length (Inclusive of channel width)
Channel_width_elements=2; % no of mesh elements on this length must even if two Chamber width is even, off if chamber width is odd

% Z direction 

SLL_thickness_elements=1; %Strain Limiting Layer (mm)% for Method 2, this could be made zero and Mat1 increased
Mat1_base_thickness_elements=1; %mat 1 basse layer thickness (mm)
Channel_height_elements=1; %height of internal channel between chambers (mm)
Channel_roof_thickness_elements=1;% (mm)
Chamber_roof_thickness_elements=3; % (mm)
Chamber_height_elements=6; % no of mesh elements on this length (Inclusive of channel parameters)


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
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
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




%% Scaling to match desired geometry
if MeshCreationTool == 1 %'Must points' method
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
    V(i,3)=z_current*((SLL_thickness+Mat1_base_thickness+Chamber_height+Chamber_roof_thickness)/Height);

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
        V_internal(i,3)=(SLL_thickness+Mat1_base_thickness)+(z_current*(Channel_height/Channel_height_elements));
        
    elseif V_internal(i,3) <= Chamber_height_elements
        V_internal(i,3)=(SLL_thickness+Mat1_base_thickness+Channel_height)+((z_current-Channel_height_elements)*((Chamber_height-Channel_height)/(Chamber_height_elements-Channel_height_elements)));
        
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

%%



















end




cFigure; %Displays scaled geometry of two cells
title('Scaled geometry','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',0.5);%0.5 transperancy normally hold on
hold on; plotV(V,'k.','MarkerSize',markerSize/2);
axisGeom; 
colormap(turbo(250)); icolorbar; 
% camlight headlight; 
gdrawnow; 
















