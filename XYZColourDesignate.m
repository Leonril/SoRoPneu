%% Original Gibbon Function by Leon Riley

%% Colour Specification tool

%% input C matrix, F matrix and V matrix, and C conditions

%% XYZColourDesignate
%This function checks if all vertices on a face are in plane with an X,  Y
%or Z coordinate, if so it shifts the colours up by one


function [C] = XYZColourDesignate(F,V,C,XYZdirection,XYZvalue,ElseOpt)



for i=1:1:size(F,1)
    m=0;
    for j=1:1:size(F,2)
        if V(F(i,j),XYZdirection) == XYZvalue
                m=m+1;   
        end
         if m == size(F,2)
           C(i)=0; %setting fixed face to ID 0 face  
        end
    end
   
end

if ElseOpt == 1%if this option is active, all other faces are merged to one
   for n=1:1:size(C)
       if C(n) ~= 0
           C(n)=1;
       end
   end
end

C=C+1;


end
