function outputVolume = surface2volume(inputPatches, inputGrid, verboseOutput)

DEFAULT_GRIDSIZE = 1;

if (~exist('verboseOutput','var'))
    outputFID = 0;
elseif (~verboseOutput)
    outputFID = 0;
else
    outputFID = 1;
end

tic; fprintf(outputFID,'Initializing ... ');

nFaces = size(inputPatches.faces,1);

if (size(inputPatches.faces,2) ~= 3)
   error('Matlab:surface2volume','Input faces must be triangles.'); 
end



if (~(exist('gridSize','var')))
    gridSize = DEFAULT_GRIDSIZE;
end

if (~exist('inputGrid','var'))
    inputGrid = [];
end

if (isempty(inputGrid))
    minVPos = min(inputPatches.vertices,[],1);
    maxVPos = max(inputPatches.vertices,[],1);
    [X,Y,Z] = meshgrid( minVPos(1) - gridSize : gridSize : maxVPos(1)+ gridSize, ...
                        minVPos(2) - gridSize : gridSize : maxVPos(2)+ gridSize, ...
                        minVPos(3) - gridSize : gridSize : maxVPos(3)+ gridSize);
    inputGrid = {X Y Z};
else
    % test for gridsize of inputgrid
    grX = diff(inputGrid{1,1}([1,2]));
    grY = diff(inputGrid{1,2}([1,2]));
    grZ = diff(inputGrid{1,3}([1,2]));
    if (~(grX | grY | grZ))
       error('Matlab:surface2volume','The input grid is not equidistant.');
    else
        sizeV = [grX grY grZ];
        pV = find(sizeV);
        gridSize = sizeV(pV(1));
    end
end

startGridPosX = min(inputGrid{1,1}(:));
startGridPosY = min(inputGrid{1,2}(:));
startGridPosZ = min(inputGrid{1,3}(:));

endGridPosX = max(inputGrid{1,1}(:));
endGridPosY = max(inputGrid{1,2}(:));
endGridPosZ = max(inputGrid{1,3}(:));


if (~iscell(inputGrid))
    error('Matlab:surface2volume','Input grid must be a cell array.');
end

initialVSize = size(inputGrid{1,1});

outputVolume = zeros(initialVSize,'uint8');

fprintf(outputFID,'done in %g sec\n',toc);

% raster faces to points

tic; fprintf(outputFID,'Rasterize points in patches to grid points ... ');

for iFace = 1:nFaces
    pointMatrix = inputPatches.vertices([inputPatches.faces(iFace,:)],:);
    % test if triangle is in the search volume 
    % (this cost around 10% performance, if were is nothing to omit, 
    % however it can reduce the amount of cycles, if the search volume 
    % is smaller than the surface object)
    inX = nnz(((pointMatrix(:,1) >= startGridPosX) & (pointMatrix(:,1) <= endGridPosX)));
    inY = nnz(((pointMatrix(:,2) >= startGridPosY) & (pointMatrix(:,2) <= endGridPosY)));
    inZ = nnz(((pointMatrix(:,3) >= startGridPosZ) & (pointMatrix(:,3) <= endGridPosZ)));
    if (inX  & inY & inZ)
        tPoints = pointsontriangle(pointMatrix, gridSize/2);
        % round positions on grid indices
        tPoints(:,1) = round((tPoints(:,1)-startGridPosX+gridSize)./gridSize) + 1;
        tPoints(:,2) = round((tPoints(:,2)-startGridPosY+gridSize)./gridSize) + 1;
        tPoints(:,3) = round((tPoints(:,3)-startGridPosZ+gridSize)./gridSize) + 1;
        % set the indices belonging to the points in the outputvolume to 1
        nTestPoints = size(tPoints,1);
        for iTestPoint = 1:nTestPoints
            outputVolume(tPoints(iTestPoint,1),tPoints(iTestPoint,2),tPoints(iTestPoint,3)) = 1;
        end
    end
end

% crop outputvolume if some point were found outside
if (size(outputVolume) ~= initialVSize);
    outputVolume = outputVolume(1:initialVSize(1),1:initialVSize(2),1:initialVSize(2));
end

fprintf(outputFID,'done in %g sec\n',toc);

% fill the background of the volume using imfill, assuming that the lower
% left corner (beginning of the block is located in background)

tic; fprintf(outputFID,'Fill volume ... ');

volumeShape = outputVolume;

outputVolume = ~(imfill(logical(outputVolume),[1 1 1],6));
if (all(outputVolume) |  ~any(outputVolume))
    warning('Matlab:surface2volume','The output is full flooded, assume that the surface was not close. Retrying slicewise.');
    for i=1:initialVSize(3)
        outputVolume(:,:,i) = ~(imfill(logical(volumeShape(:,:,i)),[1 1],4));
    end
    if (all(outputVolume) |  ~any(outputVolume))
        warning('Matlab:surface2volume','Sorry, failed again.');
    end
end

outputVolume = logical(imdilate(uint8(outputVolume),strel('ball',1,1)));

fprintf(outputFID,'done in %g sec\n',toc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: pointsontriangle                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputPoints = pointsontriangle(triVertices, minElementLength)

DEFAULT_MINELEMENT_LENGTH = 1;

% calculate the 3 vectors of the triangle
edges = triVertices([1,1,2],:)-triVertices([2,3,3],:);

% calculate the lengths of these vectors
el = zeros(3,1);
el(1) = norm(edges(1,:));
el(2) = norm(edges(2,:));
el(3) = norm(edges(3,:));

if (~exist('minElementLength','var'))
    minElementLength = DEFAULT_MINELEMENT_LENGTH;
end

% if none of the edges is larger than the minElementLength return the
% inputVertices
if (~(nnz(el>minElementLength)))
    outputPoints = triVertices;
    return
end

longestEdge = min(el);

A = [(0:minElementLength/longestEdge:1) 1];
B = [(0:minElementLength/longestEdge:1) 1];
[A, B] = combine_vectors(A,B);
C = 1-A-B;
idInTriangle = find(C>=0);
outputPoints = zeros(length(idInTriangle),3);
outputPoints(:,1) = (triVertices(1,1).*A(idInTriangle) + triVertices(2,1).*B(idInTriangle) + triVertices(3,1).*C(idInTriangle))';
outputPoints(:,2) = (triVertices(1,2).*A(idInTriangle) + triVertices(2,2).*B(idInTriangle) + triVertices(3,2).*C(idInTriangle))';
outputPoints(:,3) = (triVertices(1,3).*A(idInTriangle) + triVertices(2,3).*B(idInTriangle) + triVertices(3,3).*C(idInTriangle))';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: combine_vectors                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nv1, nv2] = combine_vectors(v1,v2)

matrix = v1;

matrix = [repmat(matrix,[1 length(v2)]);reshape(repmat(v2,[length(v1), 1]),1,[])];

nv1 = matrix(1,:);
nv2 = matrix(2,:);
