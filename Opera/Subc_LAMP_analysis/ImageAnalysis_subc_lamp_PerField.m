function [Objects] = ImageAnalysis_subc_lamp_PerField(ch1,ch2,ch3,InfoTableThis, PreviewPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%vol(ch1)
%vol(ch2, 0, 500)
%vol(ch3, 0, 500)
    %% Segment nuclei
    NucleiBlurred = imfilter(ch1, fspecial('gaussian',5, 4)); % vol(NucleiBlurred)
    NucleiMask = NucleiBlurred > 50; % vol(NucleiMask)
        if sum(NucleiMask(:))== 0
            Objects = {}; 
            return
        end
    %% Segment LAMPsomes
    LAMPDoG = imfilter(ch2, fspecial('gaussian', 5, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 5, 2), 'symmetric'); % vol(LAMPDoG, 0, 10, 'hot')
    LAMPMaskLocal = LAMPDoG > 10; % vol(LAMPMaskLocal)
    LAMPMask = LAMPMaskLocal; % & LAMPGlobalMask;
    LAMPMask = bwareaopen (LAMPMask, 2); % vol(LAMPMask)
%     if sum(LAMPMask(:))== 0
%             Objects = {}; 
%             return
%         end% vol(GolgiMask)
%     
    %% Segment SUBCsomes
    SUBCDoG = imfilter(ch3, fspecial('gaussian', 5, 1), 'symmetric') - imfilter(ch3, fspecial('gaussian', 5, 2), 'symmetric'); % vol(SUBCDoG, 0, 100, 'hot')
    SUBCMaskLocal = SUBCDoG > 10; % vol(SUBCMaskLocal)
    SUBCMask = SUBCMaskLocal; % & SUBCGlobalMask;
    SUBCMask = bwareaopen (SUBCMask, 2); % vol(SUBCMask)
%     if sum(SUBCMask(:))== 0
%             Objects = {}; 
%             return
%         end% vol(SUBCMask)  
        
    %% Colocalization Mask
    ColoMask = LAMPMask & SUBCMask; %vol(ColoMask)   
    
    %% Morphometrics
    LAMPLM = bwlabeln(LAMPMask);
    LAMPObjects = regionprops('table', LAMPLM, {'Area'});
    SUBCLM = bwlabeln(SUBCMask);
    SUBCObjects = regionprops('table', SUBCLM, {'Area'});
    ColoLM = bwlabeln(ColoMask);
    ColoObjects = regionprops('table', ColoLM, {'Area'});
           
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
    %% LAMP derived
    Objects.LAMPArea = sum(LAMPMask(:));
    Objects.LAMPAreaNorm = sum(LAMPMask(:))/sum(NucleiMask(:));
    voxelSizeX = 0.2152;%Bin2
    voxelSizeY = 0.2152; 
    voxelSizeZ = 0.4;
    Objects.MinLAMPVol = min(LAMPObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MaxLAMPVol = max(LAMPObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MeanLAMPVol = mean(LAMPObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.StdLAMPVol = std(LAMPObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MedLAMPVol = median(LAMPObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MadLAMPVol = mad(LAMPObjects.Area, 1) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.TotLAMPVolume = sum(LAMPMask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.TotLAMPVolumeNorm = (sum(LAMPMask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ)/sum(NucleiMask(:));
    Objects.CountLAMP = size(LAMPObjects, 1);
    Objects.CountLAMPNorm = size(LAMPObjects, 1)/sum(NucleiMask(:));                
    Objects.AreaVectors = {LAMPObjects.Area};

    % Shape
    Conn6Strel = {};
    Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
    Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel = logical(cat(3, Conn6Strel{:}));
    LAMPErodedMask = imerode(LAMPMask, Conn6Strel);
    LAMPPerimMask = (LAMPMask - LAMPErodedMask) > 0;
    
    % Additional feature images and vectors
    LAMPBodyLabelIm = bwlabeln(LAMPErodedMask, 6);
    
    % Erosion derived
    Objects.LAMPPerimPixels = sum(LAMPPerimMask(:));
    Objects.LAMPBodyPixels = sum(LAMPErodedMask(:)); 
    Objects.LAMPBodyCount = max(LAMPBodyLabelIm(:)); % Needed for invagination feature
    Objects.LAMPShapeBySurface = Objects.LAMPBodyPixels / Objects.LAMPPerimPixels; % Roundness feature
    Objects.LAMPBodycountByLAMPcount = Objects.LAMPBodyCount / Objects.CountLAMP; % Invagination feature

    %% SUBC derived
    Objects.SUBCArea = sum(SUBCMask(:));
    Objects.SUBCAreaNorm = sum(SUBCMask(:))/sum(NucleiMask(:));
    voxelSizeX = 0.2152;%Bin2
    voxelSizeY = 0.2152; 
    voxelSizeZ = 0.4;
    Objects.MinSUBCVol = min(SUBCObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MaxSUBCVol = max(SUBCObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MeanSUBCVol = mean(SUBCObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.StdSUBCVol = std(SUBCObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MedSUBCVol = median(SUBCObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.MadSUBCVol = mad(SUBCObjects.Area, 1) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.TotSUBCVolume = sum(SUBCMask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ;
    Objects.TotSUBCVolumeNorm = (sum(SUBCMask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ)/sum(NucleiMask(:));
    Objects.CountSUBC = size(SUBCObjects, 1);
    Objects.CountSUBCNrom = size(SUBCObjects, 1)/sum(NucleiMask(:));                
    Objects.AreaVectors = {SUBCObjects.Area};

    % Shape
    Conn6Strel = {};
    Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
    Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel = logical(cat(3, Conn6Strel{:}));
    SUBCErodedMask = imerode(SUBCMask, Conn6Strel);
    SUBCPerimMask = (SUBCMask - SUBCErodedMask) > 0;
    
    % Additional feature images and vectors
    SUBCBodyLabelIm = bwlabeln(SUBCErodedMask, 6);
    
    % Erosion derived
    Objects.SUBCPerimPixels = sum(SUBCPerimMask(:));
    Objects.SUBCBodyPixels = sum(SUBCErodedMask(:)); 
    Objects.SUBCBodyCount = max(SUBCBodyLabelIm(:)); % Needed for invagination feature
    Objects.SUBCShapeBySurface = Objects.SUBCBodyPixels / Objects.SUBCPerimPixels; % Roundness feature
    Objects.SUBCBodycountBySUBCcount = Objects.SUBCBodyCount / Objects.CountSUBC; % Invagination feature
    
   
    %% 2D previews
    
    % Scalebar
    imSize = size(ch1);
    [BarMask, BarCenter] = f_barMask(20, voxelSizeX, imSize, imSize(1)-50, 75, 10);
    %vol(BarMask)
    %RGB = cat(3, zeros([size(ch1,1),size(ch1,2)], 'uint16'),imadjust(max(ch1, [], 3)), zeros([size(ch1,1),size(ch1,2)], 'uint16'));
    RGB = cat(3, imadjust(ch2(:,:,3)), imadjust(max(ch3, [], 3)), imadjust(ch1(:,:,3)));
    RGB = imoverlay(RGB, BarMask, [1 1 1]);
    % imtool(RGB)
  
    LAMPMaskPreview = imoverlay(imadjust(max(ch2, [], 3)), bwperim(max(LAMPMask, [], 3)), [1 0 0]);
    LAMPMaskPreview = imoverlay(LAMPMaskPreview, BarMask, [1 1 1]);
    %imtool(LAMPMaskPreview )
    LAMPRawPreview = imadjust(max(ch2, [], 3));
    LAMPRawPreview = imoverlay(LAMPRawPreview, BarMask, [1 1 1]);
    %imtool(LAMPRawPreview)
    
    SUBCMaskPreview = imoverlay(imadjust(max(ch3, [], 3)), bwperim(max(SUBCMask, [], 3)), [1 0 0]);
    SUBCMaskPreview = imoverlay(SUBCMaskPreview, BarMask, [1 1 1]);
    %imtool(SUBCMaskPreview )
    SUBCRawPreview = imadjust(max(ch3, [], 3));
    SUBCRawPreview = imoverlay(SUBCRawPreview, BarMask, [1 1 1]);
    %imtool(SUBCRawPreview)
    
   
    SavePathLAMPMaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_LAMP_mask.png'];
    SavePathLAMPrawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_LAMP_raw.png'];
    SavePathSUBCMaskPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_SUBC_mask.png'];
    SavePathSUBCrawPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_SUBC_raw.png'];
    SavePathRGBPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_RGB.png'];
    
   
    imwrite(LAMPMaskPreview, SavePathLAMPMaskPreview)
    imwrite(LAMPRawPreview,  SavePathLAMPrawPreview)
    imwrite(SUBCMaskPreview, SavePathSUBCMaskPreview)
    imwrite(SUBCRawPreview,  SavePathSUBCrawPreview)
    imwrite(RGB, SavePathRGBPreview)
    
    %% Colo derived
    if sum(ColoMask(:))> 0
            Objects.ColoArea = sum(ColoMask(:));
            Objects.ColoAreaNorm = sum(ColoMask(:))/sum(NucleiMask(:));
            voxelSizeX = 0.2152;%Bin2
            voxelSizeY = 0.2152; 
            voxelSizeZ = 0.4;
            Objects.MinColoVol = min(ColoObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
            Objects.MaxColoVol = max(ColoObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
            Objects.MeanColoVol = mean(ColoObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
            Objects.StdColoVol = std(ColoObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
            Objects.MedColoVol = median(ColoObjects.Area) * voxelSizeX * voxelSizeY * voxelSizeZ;
            Objects.MadColoVol = mad(ColoObjects.Area, 1) * voxelSizeX * voxelSizeY * voxelSizeZ;
            Objects.TotColoVolume = sum(ColoMask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ;
            Objects.TotColoVolumeNorm = (sum(ColoMask(:)) * voxelSizeX * voxelSizeY * voxelSizeZ)/sum(NucleiMask(:));
            Objects.CountColo = size(ColoObjects, 1);
            Objects.CountColoNrom = size(ColoObjects, 1)/sum(NucleiMask(:));                
            Objects.AreaVectors = {ColoObjects.Area};

            % Shape
            Conn6Strel = {};
            Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
            Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
            Conn6Strel = logical(cat(3, Conn6Strel{:}));
            ColoErodedMask = imerode(ColoMask, Conn6Strel);
            ColoPerimMask = (ColoMask - ColoErodedMask) > 0;

            % Additional feature images and vectors
            ColoBodyLabelIm = bwlabeln(ColoErodedMask, 6);

            % Erosion derived
            Objects.ColoPerimPixels = sum(ColoPerimMask(:));
            Objects.ColoBodyPixels = sum(ColoErodedMask(:)); 
            Objects.ColoBodyCount = max(ColoBodyLabelIm(:)); % Needed for invagination feature
            Objects.ColoShapeBySurface = Objects.ColoBodyPixels / Objects.ColoPerimPixels; % Roundness feature
            Objects.ColoBodycountByColocount = Objects.ColoBodyCount / Objects.CountColo; % Invagination feature 
            % Previews
            ColoMaskPreview = cat(3, imadjust(max(ch2, [], 3)), imadjust(max(ch3, [], 3)), zeros([size(ch1,1),size(ch1,2)], 'uint16'));
            ColoMaskPreview = imoverlay(ColoMaskPreview, bwperim(max(ColoMask, [], 3)), [0 0 1]);
            %imtool(ColoMaskPreview)
            SavePathColoPreview = [PreviewPath, filesep, Objects.Barcode{:}, '_', num2str(Objects.ROW), '_', num2str(Objects.COL), '_', num2str(Objects.Field), '_Colo.png'];
            imwrite(ColoMaskPreview, SavePathColoPreview)
            else
            Objects = {}; 
            return
    end
end

