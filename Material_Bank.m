%% Material Bank - developed by Leon Riley for use with GIBBON toolbox
% This code updates file containing common material properties used by
% FEBio as part of thePneunet MEME thesis by Leon Riley for potential
% expansion.

% This allows the user to only define the material name and the ID, which
% allows for simpler analysis where adjusting the material is not a large
% concern, such as geometric analyses.

%Material property templates

% 3rd order Ogden


% febio_spec.Material.material{1}.ATTR.name='SoftmaterialSample1';
% febio_spec.Material.material{1}.ATTR.type='Ogden';
% febio_spec.Material.material{1}.ATTR.id=1;
% febio_spec.Material.material{1}.c1=c1;
% febio_spec.Material.material{1}.m1=m1;
% febio_spec.Material.material{1}.c2=c1;
% febio_spec.Material.material{1}.m2=-m1;
% febio_spec.Material.material{1}.k=k1;



% Neo-Hookean
% febio_spec.Material.material{2}.ATTR.name=materialName2;
% febio_spec.Material.material{2}.ATTR.type='neo-Hookean';
% febio_spec.Material.material{2}.ATTR.id=2;
% febio_spec.Material.material{2}.E=E_material2;
% febio_spec.Material.material{2}.v=Poissons_material2;

% Rigid Body
% materialBank.RigidBody.ATTR.name='RigidBody';
% materialBank.RigidBody.ATTR.type='rigid body';
% materialBank.RigidBody.ATTR.id=1;
% materialBank.RigidBody.density=1;
% materialBank.RigidBody.center_of_mass=[0 0 0]; %will be overwritten

function [outputMaterialBank]=Material_Bank(varargin)


%% Soft Material Sample 1

materialBank.SoftMaterialSample1.ATTR.name='SoftMaterial1';%will be overwritten by user
materialBank.SoftMaterialSample1.ATTR.type='Ogden';
materialBank.SoftMaterialSample1.ATTR.id=1;% will be overwritten by user
materialBank.SoftMaterialSample1.c1=1;
materialBank.SoftMaterialSample1.m1=2;
materialBank.SoftMaterialSample1.c2=1;
materialBank.SoftMaterialSample1.m2=-2;
materialBank.SoftMaterialSample1.k=200;

%% Soft Material Sample 2

materialBank.SoftMaterialSample2.ATTR.name='SoftMaterial2';%will be overwritten by user
materialBank.SoftMaterialSample2.ATTR.type='Ogden';
materialBank.SoftMaterialSample2.ATTR.id=1;% will be overwritten by user
materialBank.SoftMaterialSample2.c1=50;
materialBank.SoftMaterialSample2.m1=2;
materialBank.SoftMaterialSample2.c2=50;
materialBank.SoftMaterialSample2.m2=-2;
materialBank.SoftMaterialSample2.k=5000;

%% Rigid body
materialBank.RigidBody.ATTR.name='RigidBody';
materialBank.RigidBody.ATTR.type='rigid body';
materialBank.RigidBody.ATTR.id=1;
materialBank.RigidBody.density=1;
materialBank.RigidBody.center_of_mass=[0 0 0]; %will be overwritten

%% Ecoflex 00-10

materialBank.Ecoflex0010.ATTR.name='Ecoflex0010';%will be overwritten by user
materialBank.Ecoflex0010.ATTR.type='Ogden';
materialBank.Ecoflex0010.ATTR.id=1;% will be overwritten by user
materialBank.Ecoflex0010.c1=0.039;
materialBank.Ecoflex0010.m1=2.886;
materialBank.Ecoflex0010.c2=0.039;
materialBank.Ecoflex0010.m2=-2.886;
materialBank.Ecoflex0010.k=3.9;%% Soft Material Sample 1

%% Ecoflex 00-30
materialBank.Ecoflex0030.ATTR.name='Ecoflex0030';%will be overwritten by user
materialBank.Ecoflex0030.ATTR.type='Ogden';
materialBank.Ecoflex0030.ATTR.id=1;% will be overwritten by user
materialBank.Ecoflex0030.c1=222;
materialBank.Ecoflex0030.m1=2.686;
materialBank.Ecoflex0030.c2=222;
materialBank.Ecoflex0030.m2=-2.686;
materialBank.Ecoflex0030.k=20000;



%% Outputting all materials

[outputMaterialBank]=materialBank;
