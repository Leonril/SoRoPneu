# SoRoPneu
Open source MATLAB framework employing the use of the GIBBON toolbox and FEA software FEBio3 to allow the rapid optimisation and prototyping of soft robotic pneu-net actuators.
Developed by Leon Riley. Supported by Dr Kevin M. Moerman.


SoRoPneu is only comptaible with GIBBON code, MATLAB and FEBio 3.0. Make sure to install these correctly. 
GIBBON can be downloaded here: 
https://www.gibboncode.org/ or https://github.com/gibbonCode/GIBBON.
FEBio 3.0 can be downoaded here:
https://febio.org/downloads/.

# Contribution Guidelines:
There is no stringent variable name enforcement although try to keep the use consistent thoughout.
Where possible use industry standard variations of F,E,V,C for Face, Element, Vertex and Colour arrays.


# Instructions for use:
Open newest version of main line code in the form 'PneunetFinalRev#'.

Choose your pneu-net design: pick dimensions, materials and strain limiting layer methods.

Choose your FE simulation methods: there are multiple force and contact effects to play with.

Print your actuators: If you have access to a soft material 3D printer or liquid elastomer, use the framework to automatically generate pneu-net .stls or mould .stls.

