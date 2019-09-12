%% User inputs
clear
clc

SetupMode = 1; % 1 for creating numeric organoid labels OR 0 for linking the final analysis to human labels
SlideLayout = 'GG_20180417_WT_CO55days_Hoechst_Tunnel568reimaging.txt';

SavePath = 'S:\HCS_Platform\Data\GemmaGomez\Tunel\GG_20180417_WT_CO55days_Hoechst_Tunnel568reimaging';
mkdir(SavePath)
%PreviewSavePath = 'S:\HCS_Platform\Data\SilviaBolognin\OperettaBeginnings\Previews';
PreviewSavePath = [SavePath, filesep, 'Previews'];
mkdir(PreviewSavePath)

%% Parallel pool control
delete(gcp('nocreate'))
myCluster = parcluster;
Workers = myCluster.NumWorkers;
parpool(28) %for HURRICANE
% parpool(Workers) % for MEGATRON

%% Run mode control
if SetupMode
    RunMode = 0;
else
    RunMode = 1;
end

%% Common part

    channelID = 1; % channel to use for overview
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\GG_20180416_Mut_CO55days_Hoechst_Tunnel568TRITC\94bc4588-bbaf-4683-b010-c1cc24c16c7c\metadata.csv');
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\GG_20180410_Mut_CO55days_Hoechst_Tunnel488\fbed441f-022b-4907-b238-5252009739e6\metadata.csv');
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\GG_20180417_Mut_CO55days_Hoechst_Tunnel568TRITC_rerun\aceab072-3850-419c-998b-694e6a74b344\metadata.csv');
InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\GG_20180417_WT_CO55days_Hoechst_Tunnel568reimaging\9f46ff35-de8f-4599-8d6f-5f907dcbcc10\metadata.csv');
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\GG_20180416_Wt_CO55days_Hoechst_Tunnel568TRITC\8ecd8988-f624-4997-b89b-b2bf0bba963e\metadata.csv');

    ChannelNames = unique(InfoTable.Channel);
    Channels = length(unique(InfoTable.Channel));
    Planes = unique(InfoTable.Plane)';
    Timepoints = unique(InfoTable.Timepoint)' + 1;
    [GroupsTable, GroupsIm5DCellArray] = FindGroups(InfoTable); % it(GroupsTable)
    
%% Setup mode

if SetupMode == 1
    
    Preview = CreateLabelHelpPreview(GroupsTable, PreviewSavePath);
    imwrite(Preview, [PreviewSavePath, filesep, 'layout.png'])
    Message = ['The plate Layout has been saved at ', [PreviewSavePath, filesep, 'layout.png'], '. Please save a text file without header, using tab separation, where the first column is the index number as shown in the preview and the second is the area name. Save the text file as SlideLayout_Date.txt in your working directory and set the variable SetupMode to 0 >>> Run'];
    h = msgbox(Message);
    
else    
%% Analysis mode
    
    % Load annotations
    Layout = readtable(SlideLayout)
    Layout.Properties.VariableNames = {'Idx', 'AreaName'};
    
    % Load images and organize in an XYC array
    Groups = unique(GroupsTable(GroupsTable > 0))';
    GroupPreviews = {};
    ObjectsAll = {};

    for g = Groups % Number of organoids
    %for g = 1:2 % Number of organoids
    %for g = 7:11 % Number of organoids
        XYMosaicCells = {};
        GroupZone = GroupsTable == g;
        [GroupIdxRowVec, GroupIdxColVec] = find(GroupZone); % linear indexes
        Elements = sum(GroupZone(:));
        InfoTablesThisGroup = {};
        for e = 1:Elements % Fields of a given organoid
            for c = 1:Channels
            InfoTableThisField = GroupsIm5DCellArray{GroupIdxRowVec(e), GroupIdxColVec(e)};
            InfoTablesThisGroup{e} = InfoTableThisField;
            InfoTableThisChannel = InfoTableThisField(strcmp(InfoTableThisField.Channel, ChannelNames{c}), :);
                clear Im4D
                for t = Timepoints
                    for p = Planes
                        InfoTableThisChannelThisPlane = InfoTableThisChannel(InfoTableThisChannel.Plane == p, :);
                        ImPathThisPlane = InfoTableThisChannelThisPlane.Path{:};   
                        Im4D(:,:,t,p) = imread(ImPathThisPlane); % it(Im4D(:,:,t,p))
                    end
                end
               XYMosaicCells{c}{GroupIdxRowVec(e), GroupIdxColVec(e)} = Im4D; % Operetta counterpart of XYmosaicCells for Opera
            end
        end

        InfoTableThisGroup = vertcat(InfoTablesThisGroup{:});

        %% Remove empty cells
        XYMosaicCells = cellfun(@(x) GroupClipper(x),  XYMosaicCells, 'UniformOutput', false);

        %% Stitch
        XYmosaicContourCell = cellfun(@(x) stdfilt(x, ones(3)), XYMosaicCells{1}, 'UniformOutput', false);
        XPositions = unique(InfoTableThisGroup.PositionX); % m
        YPositions = unique(InfoTableThisGroup.PositionY); % m
        ResolutionXY = 675 / 1360; % um per pixel
        MaxPixelDrift = 30;
        PreviewChannel = 1;
        ShowProgress = 0;
        [CroppedMosaic, StitchedIm] = f_stitching_operetta(XYMosaicCells, XYmosaicContourCell, XPositions, YPositions, ResolutionXY, MaxPixelDrift, PreviewChannel, ShowProgress);
        GroupPreviews{g} = max(CroppedMosaic{channelID},[],3); %it(GroupPreviews{g})
        
        %% Image analysis
        Label = Layout(g,:);
        %try
            ObjectsThisOrganoid = f_ImageAnalysisPerOperetta_cell_count(Label, CroppedMosaic{1}, CroppedMosaic{2}, ChannelNames, PreviewSavePath);
            %ObjectsThisOrganoid = Copy_of_f_ImageAnalysisPerOperettaOrganoid_20171129(Label, CroppedMosaic{1}(1:500, 1:500, 10:30), CroppedMosaic{2}(1:500, 1:500, 10:30), CroppedMosaic{3}(1:500, 1:500, 10:30), ChannelNames, PreviewSavePath);
%         catch
%             Errors{g} = 'Image analysis failed';
%             continue % next group g
%         end
            ObjectsAll{g} = ObjectsThisOrganoid;

    end
    
    Objects = vertcat(ObjectsAll{:});
    save([SavePath, filesep, 'Objects.mat'], 'Objects');
    writetable(Objects, [SavePath, filesep, 'Objects.csv'])
    writetable(Objects, [SavePath, filesep, 'Objects.xlsx'])

    %% Preview of the whole slide

    SizeSingleIm = size(XYMosaicCells{1}{1,1});
    SizeSingleIm = SizeSingleIm(1:2);
    RescanGridSize = size(GroupsTable);
    GreatPreview = zeros(SizeSingleIm(1)*RescanGridSize(1), SizeSingleIm(2)*RescanGridSize(2), 'uint16');
    ImHeight = SizeSingleIm(1);
    ImWidth = SizeSingleIm(2);
    StartRCell = {};
    StartCCell = {};

    for g = Groups
        g
        StitchedGroupSize = size(GroupPreviews{g});
        ZoneNow = GroupsTable == g;
        [R,C] = find(ZoneNow)
        StartR = min(R);
        StartC = min(C);
        StartRPixel = ((StartR-1) * ImHeight) + 1;
        %EndRPixel = StartRPixel + (3 * ImHeight) - 1;
        EndRPixel = StartRPixel + StitchedGroupSize(1) - 1;
        StartCPixel = ((StartC-1) * ImWidth) + 1;
        %EndCPixel = StartCPixel + (3 * ImWidth) - 1;
        EndCPixel = StartCPixel + StitchedGroupSize(2) - 1;
        GreatPreview(StartRPixel:EndRPixel, StartCPixel:EndCPixel) = GroupPreviews{g};
        StartRCell{g} = StartRPixel;
        StartCCell{g} = StartCPixel;
    end

    Zoomfactor = 50;
    %GreatPreviewResized = imresize(imadjust(GreatPreview), 1/Zoomfactor);
    GreatPreviewResized = imresize(imadjust(GreatPreview, [0 0.02], [0 1]), 1/Zoomfactor);

    for g = Groups
        GreatPreviewResized = insertText(GreatPreviewResized, [round(StartCCell{g}/Zoomfactor), round(StartRCell{g}/Zoomfactor)], num2str(g), 'FontSize', 12, 'BoxColor', 'red', 'TextColor', 'white');
    end
    
    %imtool(GreatPreview)
    %imtool(GreatPreviewResized)
    imwrite(GreatPreviewResized, [SavePath, filesep, 'GreatPreview.png'])
    
    % it(GreatPreviewResized)
    % save([SavePath, filesep, 'WorkspaceIncludingObjects.mat'])
    
end




