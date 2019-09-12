function  [ObjectsThisOrganoid] = f_ImageAnalysisPerOperetta_cell_count(Label, ch1, ch2, ChannelNames, PreviewPath);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

   
    % vol(ch1, 0, 10000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch2, 0, 2000) % TRITC >>> Tunel

    %% Initialize variables
    NucleiMask = [];
    TunelMask = [];
    TunelMaskinNucMaskHigh = [];
    NucMaskAlive = [];
    %% Segment nuclei
    %vol(ch1, 0, 3000)
    ch1BlurSmall = imfilter(double(ch1), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch1BlurSmall)
    ch1BlurBig = imfilter(double(ch1), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch1BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch1DoG = ch1BlurSmall - ch1BlurBig; %vol(ch1DoG, 0, 200, 'hot')
    NucleiMask = ch1DoG > 75; %vol(NucleiMask)
    NucleiMask = bwareaopen(NucleiMask, 40);%vol(NucleiMask)
    ch1LP = imfilter(ch1, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch1LP, 0, 20000, 'hot')
    NucMaskHigh =  (ch1LP > 8000) .* NucleiMask; %vol(NucMaskHigh, 0, 1)
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
    
     %% Tunel (ch2)
    
    %vol(ch2, 0, 2000)
    ch2MedFilt = [];
    parfor p = 1:size(ch2, 3)
        ch2MedFilt(:,:,p) = medfilt2(ch2(:,:,p));
    end
    %vol(ch2MedFilt, 0, 1500, 'hot')    
    TunelMask = ch2MedFilt > 600; %vol(TunelMask, 0, 1)
%     TunelMask = TunelMask & NucleiMask;
    TunelMask = bwareaopen(TunelMask, 200);    
    TunelMaskinNucMaskHigh = TunelMask & NucMaskHigh; %vol(TunelMaskinNucMaskHigh)
    TunelMaskinNucMaskAlive = TunelMask & NucMaskAlive;
    %% Perinuclear Volume (to detect amount of cells positive for Tunel)
    NucDil = imdilate(imdilate(NucleiMask, strel('disk', 4)), strel('sphere',1));
    NucPerim = NucDil & ~NucleiMask;
    %vol(NucleiMask)
    NucleiMaskSingleCells = f_RemoveBigObjects (NucleiMask, 10000); 
  
    %% Percent Tunel pos
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
    
    PerinucLM = bwlabeln(PeriNucMask);%vol(PerinucLM); vol(uint16(PeriNucMask) .* uint16(THMask), 0,1); vol(THMask +2*PeriNucMask)
    PeriNucObjects = regionprops('table', PerinucLM, double(TunelMask), 'PixelValues');
    Tunelproportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjects, 'InputVariables', 'PixelValues');
    TunelPos = array2table(table2array(Tunelproportions) > 0.01);
    TunelPos.Properties.VariableNames(end) = {'TunelPos'};
    PeriNucObjects = [PeriNucObjects, Tunelproportions, TunelPos];
    PeriNucObjects.Properties.VariableNames(end-1) = {'Tunelproportion'};
    PeriNucObjectsCompact = PeriNucObjects(:, {'Tunelproportion','TunelPos'});
    TunelPercent = (sum(PeriNucObjectsCompact.TunelPos)/height(PeriNucObjectsCompact))*100;
    
    [~, NucCountSingleCells] = bwlabeln(NucleiMaskSingleCells, 6); % clumps removed
    [~, NucCountAll] = bwlabeln(NucleiMask, 6);% clumps included
%     ObjectsThisOrganoid.NucCountSingleCells = NucCountSingleCells;
  
    %% Previews 
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)

    PreviewTunelMask = imoverlay2(imadjust(max(ch2,[],3),[0 0.01]), bwperim(max(TunelMask,[],3)), [1 0 0]);
    PreviewTunelMask = imoverlay2(PreviewTunelMask, BarMask, [1 1 1]);
    %imtool(PreviewTunelMask)
    PreviewHoechst = imoverlay2(imadjust(max(ch1,[],3),[0 0.06]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)

   
    %imwrite(PreviewTuj1, [PreviewPath, filesep, CellLine, '_', 'Tuj1', '_', num2str(s), '.png'])
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewTunelMask, [PreviewPath, filesep, IdentityString, '_', 'Tunel', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    
    %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.TunelMaskSum = sum(TunelMask(:));
    ObjectsThisOrganoid.TunelMaskByNuc = sum(TunelMask(:)) / sum(NucleiMask(:)); 
    ObjectsThisOrganoid.TunelPercent = TunelPercent;
    ObjectsThisOrganoid.NucCountAll = NucCountAll;
    ObjectsThisOrganoid.NucCountSingleCells = NucCountSingleCells;

end

