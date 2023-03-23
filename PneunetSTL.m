%% Stl practice


%% Mould bottom
% general form



bottomDimEl=[4 4 3];

[meshStruct]=hexMeshBox(bottomDimEl,bottomDimEl);

E_bar=meshStruct.E;
V_bar=meshStruct.V;
F_bar=meshStruct.F;
Fb_bar=meshStruct.Fb;
Cb_bar=meshStruct.faceBoundaryMarker;

%Getting centre for easier control of elements
VE_STLBase=patchCentre(E_bar,V_bar);

%Shifting data for easier deletion of sink elements
VE_STLBase(:,1)=abs(VE_STLBase(:,1));
VE_STLBase(:,2)=abs(VE_STLBase(:,2));
VE_STLBase(:,3)=VE_STLBase(:,3)-min(VE_STLBase(:,3))+0.5;

%Deleting sink elements for mould pouring
logicDeleteSTLBase= VE_STLBase(:,3)>2 & VE_STLBase(:,2) < 1 & VE_STLBase(:,1) < 1
E_bottomMould=E_bar(~logicDeleteSTLBase,:);

%scaling the bottom mould
bottomMouldDimX=[]
bottomMouldDimY=
bottomMouldDimZ=