function [nodefam,numfam,pointfam,fail] = neighboornodes(coord,Geome)
delta = Geome.delta;
totnode = size(coord,1);
numfam = int32(zeros(totnode,1));
pointfam = int32(zeros(totnode+1, 1));
fail = ones(totnode,200);
%nodefam = int32(zeros(totnode,totnode));
if Geome.prec % *** The precrack process applies only to 2D problems ***
    TipL = Geome.TipL;
    TipR = Geome.TipR;
    XY1 = [TipL(:,1) TipL(:,2) TipR(:,1) TipR(:,2)];
end
[Idx,~] = rangesearch(coord,coord,delta);
for i = 1:totnode
    Idx{i,1} = sort(Idx{i,1}(2:end));
    %numfam(i,1) = size(Idx{i,1},2);
    numfam(i,1) = size(Idx{i,1},2);
    pointfam(i+1,1) = pointfam(i)+numfam(i,1);
    %nodefam(i,1:numfam(i,1)) = Idx{i,1};
    if Geome.prec && dof == 2
        for k = 1:numfam(i,1)
            j = Idx{i,1}(1,k);
            XY2 = [coord(j, 1) coord(j, 2) coord(i, 1) coord(i, 2)];
            out = lineSegmentIntersect(XY1,XY2);
            Ints = -out.intAdjacencyMatrix; 
            % if the crack line XY1 intersects XY2, Ints = -1, else Ints = 0
            if Ints == -1
                fail(i,j) = 0;
            elseif Ints ~= -1 || Ints ~= 0
                error('ERROR:Check (lineSegmentIntersect function)')
            end
        end
    end
end
nodefam = int32(([Idx{:}])');


