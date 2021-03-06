function [Objects] = ImageAnalysisGolgi3D_PerField(ch1,ch2,InfoTableThis, PreviewPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%vol(ch1)
%vol(ch2, 0, 1000)
    %% Segment nuclei
    NucleiBlurred = imfilter(ch1, fspecial('gaussian',5, 4)); % vol(NucleiBlurred)
    NucleiMask = NucleiBlurred > 50; % vol(NucleiMask)
        if sum(NucleiMask(:))== 0
            Objects = {}; 
            return
        end
    %% Segment Golgi
    GolgiDoG = imfilter(ch2, fspecial('gaussian', 7, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 7, 2), 'symmetric'); % vol(GolgiDoG, 0, 100, 'hot')
    GolgiMaskLocal = GolgiDoG > 10; % vol(GolgiMaskLocal)
    GolgiMask = GolgiMaskLocal; % & GolgiGlobalMask;
    GolgiMask = bwareaopen (GolgiMask, 7);
    if sum(GolgiMask(:))== 0
            Objects = {}; 
            return
        end% vol(GolgiMask)

    %% Morphometrics
    GolgiLM = bwlabeln(GolgiMask);
    GolgiObjects = regionprops('table', GolgiLM, {'Area'});
    %% Extract features

    Objects = table();
    Objects.File = InfoTableThis.files;
    Objects.Barcode = InfoTableThis.Barcode;
    Objects.AreaName = InfoTableThis.AreaName;
    Objects.ROW = InfoTableThis.Row;
    Objects.COL = InfoTableThis.Column;
    Objects.Field = InfoTableThis.field;
    %% Nucleus derived
    Objects.NucArea = sum(NucleiMask(:));
    %% Golgi derived
    Objects.GolgiArea = sum(GolgiMask(:));% ATPB plus Tom20
    Objects.GolgiAreaNorm = sum(GolgiMask(:))/sum(NucleiMask(:));
    voxelSizeX = 0.2152;%Bin2
    voxelSizeY = 0.2152; 
    voxelSizeZ = 1;
    Objects.MinGolgiVol = min(GolgiObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MaxGolgiVol = max(GolgiObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MeanGolgiVol = mean(GolgiObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.StdGolgiVol = std(GolgiObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MedGolgiVol = median(GolgiObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MadGolgiVol = mad(GolgiObjects.Area, 1) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.TotGolgiVolume = sum(GolgiMask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.TotGolgiVolumeNorm = (sum(GolgiMask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ)/sum(NucleiMask(:));
    Objects.CountGolgi = size(GolgiObjects, 1);  
    Objects.CountGolgiNorm = size(GolgiObjects, 1)/sum(NucleiMask(:));
    Objects.AreaVectors = {GolgiObjects.Area};
 
    % Shape
    Conn6Strel = {};
    Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
    Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel = logical(cat(3, Conn6Strel{:}));
    GolgiErodedMask = imerode(GolgiMask, Conn6Strel);
    GolgiPerimMask = (GolgiMask - GolgiErodedMask) > 0;

    % skeleton  
    skel = Skeleton3D(GolgiMask);    
    
    [AdjacencyMatrix, node, link] = Skel2Graph3D(skel,0); % 0 for keeping all branches 

    %% Additional feature images and vectors
    GolgiBodyLabelIm = bwlabeln(GolgiErodedMask, 6);
    Objects.NodeDegreeVector = {sum(AdjacencyMatrix, 1)};

    % Erosion derived
    Objects.GolgiSkelPixels = sum(skel(:));
    Objects.GolgiPerimPixels = sum(GolgiPerimMask(:));
    Objects.GolgiBodyPixels = sum(GolgiErodedMask(:)); 
    Objects.GolgiBodyCount = max(GolgiBodyLabelIm(:)); % Needed for invagination feature
    Objects.GolgiShapeBySurface = Objects.GolgiBodyPixels / Objects.GolgiPerimPixels; % Roundness feature
    Objects.GolgiBodycountByGolgicount = Objects.GolgiBodyCount / Objects.CountGolgi; % Invagination feature

    % Skeleton derived
    Objects.TotalNodeCount = size(node, 2);
    Objects.AverageNodePerGolgi = Objects.TotalNodeCount / Objects.CountGolgi;
    Objects.TotalLinkCount = size(link, 2);
    Objects.NodesPerGolgiMean = Objects.TotalNodeCount / Objects.CountGolgi;
    Objects.LinksPerGolgiMean = Objects.TotalLinkCount / Objects.CountGolgi;
    Objects.AverageNodeDegree = mean(Objects.NodeDegreeVector{:}); %average lenght of branch in pixels
    Objects.MedianNodeDegree = median(Objects.NodeDegreeVector{:});
    Objects.StdNodeDegree = std(Objects.NodeDegreeVector{:});
    Objects.MadNodeDegree = mad(Objects.NodeDegreeVector{:}, 1);
       
    %% 2D previews
    
    % Scalebar
    imSize = size(ch1);
    [BarMask, BarCenter] = f_barMask(20, voxelSizeX, imSize, imSize(1)-100, 100, 10);
    %vol(BarMask)
    RGB = cat(3, imadjust(ch2(:,:,3)),zeros([size(ch1,1),size(ch1,2)], 'uint16'), imadjust(ch1(:,:,3)));
    %RGB = cat(3, imadjust(max(ch2, [], 3), [0; 0.01],[0; 1]),zeros([size(ch1,1),size(ch1,2)], 'uint16'), imadjust(max(ch1, [], 3)));
    RGB = imoverlay(RGB, BarMask, [1 1 1]);
    % imtool(RGB)
  
    GolgiMaskPreview = imoverlay(imadjust(max(ch2, [], 3)), bwperim(max(GolgiMask, [], 3)), [1 0 0]);
    GolgiMaskPreview = imoverlay(GolgiMaskPreview, BarMask, [1 1 1]);
    %imtool(GolgiMaskPreview )
    GolgiRawPreview = imadjust(max(ch2, [], 3));
    GolgiRawPreview = imoverlay(GolgiRawPreview, BarMask, [1 1 1]);
    %imtool(GolgiRawPreview)
    
    SavePathGolgiMaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_Golgi_mask.png'];
    SavePathGolgirawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_Golgi_raw.png'];
    SavePathRGBPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RGB.png'];
   
   
    imwrite(GolgiMaskPreview, SavePathGolgiMaskPreview)
    imwrite(GolgiRawPreview,  SavePathGolgirawPreview)
    imwrite(RGB, SavePathRGBPreview)

end

