%% Mitochondrial Morphometrics Opera 3D 
clear
clc

%% User inputs

Barcodes = {'S:\OperaQEHS\OperaDB\Gemma Giro\GG_20180411_60X_Endo_SubC_LAMP_2\GG_20180411_60X_Endo_SubC_LAMP_2',...
    'S:\OperaQEHS\OperaDB\Gemma Giro\GG_20180411_60X_Endo_SubC_LAMP_3\GG_20180411_60X_Endo_SubC_LAMP_3'};


%% Run documentation

FolderThisAnalysis = ['S:\HCS_Platform\Data\JavierJarazo\For_Gemma\Epithelial_SubC_LAMP_', datestr(clock,'yyyymmdd_HH_MM_SS')];
mkdir(FolderThisAnalysis)
FileNameShort = mfilename;
newbackup = sprintf('%s_log.m',[FolderThisAnalysis, '\', FileNameShort]);
FileNameAndLocation=[mfilename('fullpath')];

currentfile = strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);
PreviewPath = [FolderThisAnalysis, '\Previews'];
mkdir(PreviewPath)
VersionMatlab = version;
save([FolderThisAnalysis, '\', 'MatlabVersion.mat'], 'VersionMatlab');

f_LogDependencies(FileNameShort, FolderThisAnalysis); 
ObjectsAllFields = {};

% delete(gcp('nocreate'))
% parpool()

for n = 1:size(Barcodes,2)


        Barcode = Barcodes{n};
        InfoTable = f_InfoTable(Barcodes{n});
        ObjectsAllFields = {};
        parfor l = 1:height(InfoTable)
        %for l = 1:2
            cube = readflexcube(InfoTable.files{l}, 'PlaneCount', 6); % Read 4-D image cube
            ch1 = cube.data(:, :, :, 1);%Hoechst
            ch2 = cube.data(:, :, :, 2);%LAMP
            ch3 = cube.data(:, :, :, 3);%SUBC
            %vol(ch1)
            %vol(ch2)
            %vol(ch3)
               
            InfoTableThis = InfoTable(l,:);
            ObjectsAllFields{n,l} = ImageAnalysis_subc_lamp_PerField(ch1,ch2,ch3,InfoTableThis, PreviewPath);
            
        end
        
ObjectsAll = vertcat(ObjectsAllFields{:});
ObjectsName = InfoTable.Barcode{1};
save([FolderThisAnalysis, filesep, ObjectsName], 'ObjectsAll')

AreaNamesToKeep = ObjectsAll.Properties.VariableNames(~ismember(ObjectsAll.Properties.VariableNames, {'AreaVectors', 'NodeDegreeVector'}));
Summary = ObjectsAll(:, AreaNamesToKeep); % remove vectors
writetable(Summary, [FolderThisAnalysis, filesep, ObjectsName, '.csv'], 'WriteVariableNames', true); % Saving as comma separated file     
        
end

