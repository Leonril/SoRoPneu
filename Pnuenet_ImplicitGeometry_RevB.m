%% Pneunet Implicit Geometry
%% Revisions
%% Rev B
% Notes:
% * Attempt at pneuent geometry replication
% * 


%%
clear; close all; clc;

%% Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=0.3;
markerSize=40;
lineWidth=3;
cMap=blood(250);

%% Creating shape equation

pointSpacing=0.2;% controlling size of tri mesh

contourLevel=0;% position in eq result where mesh should form

voxelSizeApprox=pointSpacing/3;%approx size of voxels (mm)

nRepeat=4;%number of repeating chambers
x=linspace(-3,3,ceil(6/voxelSizeApprox));%linear vector controlling X space
y=linspace(-3,3,ceil(6/voxelSizeApprox));%linear vector controlling X space
z=linspace(-3,3,ceil(6/voxelSizeApprox));%linear vector controlling X space


[X,Y,Z]=ndgrid(x,y,z);% transforms the x,y,z vectors into a 3D space 


R=sqrt(X.^2+Y.^2); %polar coordinate
V_grid=[X(:) Y(:) Z(:)];% turns X,Y,Z into a grid 

voxelSize=[mean(diff(x)) mean(diff(y)) mean(diff(z))]; %designated Voxel size
theta=atan(Y./X);
M_grid=((R.*sin(theta)-0.5).^4) + (2.*((R.*sin(theta+(pi/2))).^4))+((sin(atan(20*sin(theta-0.5)).^(-1)))+1).*((sin(atan(20*sin(Z)).^(-1)))+1)  %().^4

%sin(nRepeat*Z)-sin(2.8*R)+2.8*R-6
%2.*((R.*sin(theta+(pi/2))).^4))%(5*sin(nRepeat*Z))+%((R.*sin(theta)-0.5).^4) %+ (2.*((R.*sin(theta+(pi/2))).^4))
%((R.*sin(theta)-0.5).^4) + (2.*((R.*sin(theta+(pi/2))).^4))+(sin(tan(20*sin(theta-0.5)).^(-1))).*(sin(tan(20*sin(Z)).^(-1)))  %().^4

%sin(nRepeat*X)-sin(2.8*R)+2.8*R-6
%((R.*sin(theta)-0.5).^4)+(2.*((R.*sin(theta+(pi/2)).^4))+((sin(atan(20*sin(theta-0.5)))).*(sin(atan(20.*sin(Z))))));
%sin(nRepeat*X).*sin(Y)%+sin(Z);
%(0.5.*Y-0.45).^4-((20.*sin(X))./(400.*sin(X).^2+1.^(1./2)+1).*((20.*sin(cos(5.*Z-2.5)))))
%sin(nRepeat*Z)-sin(2.8*R)+2.8*R-6;% Coordinate to attach value to true 

M=reshape(M_grid,size(X)); 

%%

vizStruct.colormap=warmcold(250);
vizStruct.clim=abs(max(M(:)))*[-1 1];
sv3(M,voxelSize,vizStruct); 
camlight headlight; 

%%

imOrigin=min(V_grid,[],1)-voxelSize/2; %image origin

%%

controlPar_isosurface.nSub=[1 1 1];%round(max(v)/2./v);
controlPar_isosurface.capOpt=0; %Option to cap open ended surfaces
controlPar_isosurface.voxelSize=voxelSize;
controlPar_isosurface.contourLevel=contourLevel;
[Fi,Vi]=levelset2isosurface(M,controlPar_isosurface); %Get iso-surface

Fi_sorted=sort(Fi,2);
logicInvalid=any(diff(Fi_sorted,1,2)==0,2);
Fi=Fi(~logicInvalid,:);
[Fi,Vi]=patchCleanUnused(Fi,Vi);
Vi=Vi(:,[2 1 3]);
Vi=Vi+imOrigin;
[Fi,Vi]=triSurfRemoveThreeConnect(Fi,Vi);
[Fi,Vi]=patchCleanUnused(Fi,Vi);

%%
% Remesh using geomgram
optionStructGG.pointSpacing=pointSpacing;
[Fi,Vi]=ggremesh(Fi,Vi,optionStructGG);
Fi=fliplr(Fi); %Invert face orientation

% Eb=patchBoundary(Fi);
% controlPar_smooth.Method='HC';
% controlPar_smooth.Alpha=0.1;
% controlPar_smooth.Beta=0.5;
% controlPar_smooth.n=150;
% controlPar_smooth.RigidConstraints=unique(Eb(:));
% [Vi]=patchSmooth(Fi,Vi,[],controlPar_smooth);

Vi=Vi-imOrigin; 
Vi=Vi(:,[2 1 3]);

%%

gpatch(Fi,Vi,'w','k',1,2);

%%

Eb=patchBoundary(Fi);
groupOptionStruct.outputType='label';
Gb=tesgroup(Eb,groupOptionStruct);

Eb1=Eb(Gb==1,:); %Boundary 1
indB1=edgeListToCurve(Eb1);
indB1=indB1(1:end-1);
Eb2=Eb(Gb==2,:); %Boundary 1
indB2=edgeListToCurve(Eb2);
indB2=indB2(1:end-1);

Z=Vi(:,3);
if mean(Z(indB1))>mean(Z(indB2))
    indTop=indB1;
    indBottom=indB2;
else
    indTop=indB2;
    indBottom=indB1;
end

% patching top edge
[Fb1,Vb1]=regionTriMesh3D({Vi(indTop,:)},pointSpacing,0,'natural');
N1=patchNormal(Fb1,Vb1);
if dot(mean(N1,1),[0 0 1])<0
    Fb1=fliplr(Fb1);
end

% patching bottom edge
[Fb2,Vb2]=regionTriMesh3D({Vi(indBottom,:)},pointSpacing,0,'natural');
N2=patchNormal(Fb2,Vb2);
if dot(mean(N2,1),[0 0 1])>0
    Fb2=fliplr(Fb2);
end

%%

cFigure; hold on;
gpatch(Fi,Vi,'w','k');
patchNormPlot(Fi,Vi);

gpatch(Fb1,Vb1,'rw','k');
patchNormPlot(Fb1,Vb1);

gpatch(Fb2,Vb2,'bw','k');
patchNormPlot(Fb2,Vb2);

plotV(Vi(indTop,:),'r-','LineWidth',3)
plotV(Vi(indBottom,:),'b-','LineWidth',3)

axisGeom; camlight headlight;
drawnow; 

%%

[F,V,C]=joinElementSets({Fi,Fb1,Fb2},{Vi,Vb1,Vb2});
[F,V]=mergeVertices(F,V);

%%

cFigure; hold on;
gpatch(F,V,C,'k');
colormap gjet; icolorbar;
axisGeom; camlight headlight;
drawnow; 
