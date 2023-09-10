function [coord_modi,totint_modi] = modified_model(coord,totint)
int = 0;
coord_modi = zeros(totint/2,size(coord,2));
for i =1:totint
    if norm(coord(i,:)) > 0.005
        int = int+1;
        coord_modi(int,:) = coord(i,:);
    end
end
totint_modi = int;