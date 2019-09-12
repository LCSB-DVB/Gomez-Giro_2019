function  [ObjectsThisOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch2, ch3, ChannelNames, PreviewPath);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % vol(ch2, 0, 4000) % Alexa 647 >>> GFAP
    % vol(ch3, 0, 5000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))

    %% Initialize variables
    NucleiMask = [];
    GFAPMask = [];
    
    %% Segment nuclei
    %vol(ch3, 0, 10000)
    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 200, 'hot')
    NucleiMask = ch3DoG > 75; %vol(NucleiMask)
    NucleiMask = bwareaopen(NucleiMask, 20);%vol(NucleiMask)
    ch3LP = imfilter(ch3, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch3LP, 0, 4000, 'hot')
    NucMaskHigh =  (ch3LP > 1500) .* NucleiMask; %vol(NucMaskHigh, 0, 1) % tried 3000 (includes lot of death nuclei)
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
    
    %% GFAP (ch2)

    %vol(ch2, 0, 2000)
    ch2MedFilt = []; 
    SizeZ = size(ch2, 3);
    parfor p = 1:SizeZ
        ch2MedFilt(:,:,p) = medfilt2(ch2(:,:,p));
    end
    %vol(ch2MedFilt, 0, 4000, 'hot')
    GFAPMask = ch2MedFilt > 1500; %previously 2000
    GFAPDoG = imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(GFAPDoG, 0, 300, 'hot')
    GFAPDoGMask = GFAPDoG > 300;
    %vol(GFAPDoGMask, 0, 1)
    GFAPMask = GFAPMask & GFAPDoGMask;
    %vol(GFAPMask, 0, 1)
    GFAPMask = bwareaopen(GFAPMask, 300);

    %% Perinuclear Volume (to detect amount of cells positive for GFAP)
    
    %vol(NucleiMask)
    NucleiMaskSingleCells = f_RemoveBigObjects (NucleiMask, 10000); 
    NucDil = imdilate(imdilate(NucleiMaskSingleCells, strel('disk', 4)), strel('sphere',1));
    NucPerim = logical(NucDil) & ~logical(NucleiMaskSingleCells);
    %vol(NucPerim)
    GFAPMaskinNucPerim = GFAPMask & NucPerim;% vol(GFAPMaskinNucPerim)
    
    %% Percent GFAP pos
    %split perinuc
    D = bwdist(NucleiMaskSingleCells);
    %vol(D, 0, 20, 'hot')
    %it(D(:,:,1))
    disp('start watershed')
    tic
    W = watershed(D);
    toc
    disp('watershed done')
    %vol(W)
    NucPerimStencil = uint16(W) .* uint16(imreconstruct(logical(imdilate(NucPerim, strel('disk', 1))), logical(NucleiMaskSingleCells))); % This line was causing the error 20171207 % Function imreconstruct expected MARKER and MASK to have the same class.
    %vol(NucPerimStencil)
    %vol(NucPerim)
    %vol(NucleiMaskSingleCells)
    %toto = imreconstruct(logical(NucPerim), logical(NucleiMaskSingleCells));
       
    PeriNucMask = logical(NucPerimStencil);
    PeriNucMask = bwareaopen(PeriNucMask, 500);
%     vol(PeriNucMask)
    PerinucLM = bwlabeln(PeriNucMask);
    
    % GFAP 
    PeriNucObjectsGFAP = regionprops('table', PerinucLM, double(GFAPMask), 'PixelValues');
    GFAPproportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjectsGFAP, 'InputVariables', 'PixelValues');
    GFAPPos = array2table(table2array(GFAPproportions) > 0.01);
    GFAPPos.Properties.VariableNames(end) = {'GFAPpos'};
    PeriNucObjectsGFAP = [PeriNucObjectsGFAP, GFAPproportions, GFAPPos];
    PeriNucObjectsGFAP.Properties.VariableNames(end-1) = {'GFAPproportion'};
    PeriNucObjectsCompactGFAP = PeriNucObjectsGFAP(:, {'GFAPproportion','GFAPpos'});
    GFAPPercent = (sum(PeriNucObjectsCompactGFAP.GFAPpos)/height(PeriNucObjectsCompactGFAP))*100;

    %% Previews 
    
    % Scalebar
    imSize = [size(ch2, 1), size(ch2, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)

    PreviewGFAP = imoverlay2(imadjust(max(ch2,[],3),[0 0.07]), bwperim(max(GFAPMask,[],3)), [1 0 0]);
    PreviewGFAP = imoverlay2(PreviewGFAP, BarMask, [1 1 1]);
    %imtool(PreviewGFAP)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)
    
    PreviewNucMaskAlive = imoverlay2(imadjust(max(ch3, [], 3), [0 0.07]), bwperim(max(NucMaskAlive,[],3)), [0 0 1]);
    PreviewNucMaskAlive = imoverlay2(PreviewNucMaskAlive, BarMask, [1 1 1]);
    %imtool(PreviewNucMaskAlive)
    
    PreviewNucMaskHigh = imoverlay2(imadjust(max(ch3, [], 3), [0 0.07]), bwperim(max(NucMaskHigh,[],3)), [0 0 1]);
    PreviewNucMaskHigh = imoverlay2(PreviewNucMaskHigh, BarMask, [1 1 1]);
    %imtool(PreviewNucMaskHigh)
    
    
    
    %imwrite(PreviewTuj1, [PreviewPath, filesep, CellLine, '_', 'Tuj1', '_', num2str(s), '.png'])
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewGFAP, [PreviewPath, filesep, IdentityString, '_', 'GFAP', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewNucMaskAlive, [PreviewPath, filesep, IdentityString, '_', 'PreviewNucMaskAlive', '.png'])
    imwrite(PreviewNucMaskHigh, [PreviewPath, filesep, IdentityString, '_', 'PreviewNucMaskHigh', '.png'])
   
    
    %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.GFAPMaskSum = sum(GFAPMask(:));
    ObjectsThisOrganoid.GFAPMaskByNuc = sum(GFAPMask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskHigh = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.GFAPPercent = GFAPPercent;

    

end

