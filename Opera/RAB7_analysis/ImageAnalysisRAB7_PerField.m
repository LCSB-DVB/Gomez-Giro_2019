function [Objects] = ImageAnalysisRAB7_PerField(ch1,ch2,InfoTableThis, PreviewPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%vol(ch1)
%vol(ch2, 0, 500)
    %% Segment nuclei
    NucleiBlurred = imfilter(ch1, fspecial('gaussian',5, 4)); % vol(NucleiBlurred)
    NucleiMask = NucleiBlurred > 50; % vol(NucleiMask)
        if sum(NucleiMask(:))== 0
            Objects = {}; 
            return
        end
    %% Segment RAB7somes
    RAB7DoG = imfilter(ch2, fspecial('gaussian', 7, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 7, 2), 'symmetric'); % vol(RAB7DoG, 0, 100, 'hot')
    RAB7MaskLocal = RAB7DoG > 40; % vol(RAB7MaskLocal)
    RAB7Mask = RAB7MaskLocal; % & RAB7GlobalMask;
    RAB7Mask = bwareaopen (RAB7Mask, 7); % vol(RAB7Mask)
    if sum(RAB7Mask(:))== 0
            Objects = {}; 
            return
        end% vol(RAB7Mask)

    %% Morphometrics
    RAB7LM = bwlabeln(RAB7Mask);
    RAB7Objects = regionprops('table', RAB7LM, {'Area'});
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
    %% RAB7 derived
    Objects.RAB7Area = sum(RAB7Mask(:));
    Objects.RAB7AreaNorm = sum(RAB7Mask(:))/sum(NucleiMask(:));
    voxelSizeX = 0.2152;%Bin2
    voxelSizeY = 0.2152; 
    voxelSizeZ = 0.4;
    Objects.MinRAB7Vol = min(RAB7Objects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MaxRAB7Vol = max(RAB7Objects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MeanRAB7Vol = mean(RAB7Objects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.StdRAB7Vol = std(RAB7Objects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MedRAB7Vol = median(RAB7Objects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MadRAB7Vol = mad(RAB7Objects.Area, 1) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.TotRAB7Volume = sum(RAB7Mask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.TotRAB7VolumeNorm = (sum(RAB7Mask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ)/sum(NucleiMask(:));
    Objects.CountRAB7 = size(RAB7Objects, 1);
    Objects.CountRAB7Norm = size(RAB7Objects, 1)/sum(NucleiMask(:));                
    Objects.AreaVectors = {RAB7Objects.Area};

    % Shape
    Conn6Strel = {};
    Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
    Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel = logical(cat(3, Conn6Strel{:}));
    RAB7ErodedMask = imerode(RAB7Mask, Conn6Strel);
    RAB7PerimMask = (RAB7Mask - RAB7ErodedMask) > 0;
    
    %% Additional feature images and vectors
    RAB7BodyLabelIm = bwlabeln(RAB7ErodedMask, 6);
    
    % Erosion derived
    Objects.RAB7PerimPixels = sum(RAB7PerimMask(:));
    Objects.RAB7BodyPixels = sum(RAB7ErodedMask(:)); 
    Objects.RAB7BodyCount = max(RAB7BodyLabelIm(:)); % Needed for invagination feature
    Objects.RAB7ShapeBySurface = Objects.RAB7BodyPixels / Objects.RAB7PerimPixels; % Roundness feature
    Objects.RAB7BodycountByRAB7count = Objects.RAB7BodyCount / Objects.CountRAB7; % Invagination feature

           
    %% 2D previews
    
    % Scalebar
    imSize = size(ch1);
    [BarMask, BarCenter] = f_barMask(20, voxelSizeX, imSize, imSize(1)-100, 100, 10);
    %vol(BarMask)
    %RGB = cat(3, zeros([size(ch1,1),size(ch1,2)], 'uint16'),imadjust(max(ch1, [], 3)), zeros([size(ch1,1),size(ch1,2)], 'uint16'));
    RGB = cat(3, zeros([size(ch1,1),size(ch1,2)], 'uint16'), imadjust(ch2(:,:,2)), imadjust(ch1(:,:,2)));
    RGB = imoverlay(RGB, BarMask, [1 1 1]);
    % imtool(RGB)
  
    RAB7MaskPreview = imoverlay(imadjust(max(ch2, [], 3)), bwperim(max(RAB7Mask, [], 3)), [1 0 0]);
    RAB7MaskPreview = imoverlay(RAB7MaskPreview, BarMask, [1 1 1]);
    %imtool(RAB7MaskPreview)
    RAB7RawPreview = imadjust(max(ch2, [], 3));
    RAB7RawPreview = imoverlay(RAB7RawPreview, BarMask, [1 1 1]);
    %imtool(RAB7RawPreview)
    
    SavePathRAB7MaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RAB7_mask.png'];
    SavePathRAB7rawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RAB7_raw.png'];
    SavePathRGBPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RGB.png'];
   
   
    imwrite(RAB7MaskPreview, SavePathRAB7MaskPreview)
    imwrite(RAB7RawPreview,  SavePathRAB7rawPreview)
    imwrite(RGB, SavePathRGBPreview)

end

