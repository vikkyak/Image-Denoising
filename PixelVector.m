
%% Neighborhood Matrix Generation
function [vector,index]=PixelVector(image,i,j,M)
x=(i-M):(i+M);
y=(j-M):(j+M);
neighborhood_length = (2*M+1)^2;
Window=image(x,y);
vector=reshape(Window,[1,neighborhood_length]);
index=combvec(x,y);
end







