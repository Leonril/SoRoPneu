%clear; close all; clc; 


%% Old version
% function [VS]=copyOffset(V1,n,w)
% % 
% xs=0:w:(n-1)*w;
% xs=repmat(xs,size(V1,1),1);
% xs=xs(:);
% VS=repmat(V1,n,1);
% VS(:,1)=VS(:,1)+xs;

%% With F functionality
function [FS,VS]=copyOffset(F1,V1,n,w)

xs=0:w:(n-1)*w;
xs=repmat(xs,size(V1,1),1);
xs=xs(:);
VS=repmat(V1,n,1);
VS(:,1)=VS(:,1)+xs;

is=0:size(V1,1):(n-1)*size(V1,1);
is=repmat(is,size(F1,1),1);
is=is(:);
FS=repmat(F1,n,1);
FS=FS+repmat(is,1,size(FS,2));

end
