classdef STARMapDataset
    % STARMAPDATASET Summary of this class goes here
    % Detailed explanation goes here
    
    properties
        basePath = "";
        prefix;
        
        % barcode info
        seqToName
        nameToSeq;
        
        barcodeCS; % matrix of barcodes in color space
        barcodeNames; % cell array of barcode names in same order as barcodeCS
        barcodeSeqs; % cell array of barcode sequences in colorspace (as string)
        
        allPoints; % all poitns found by algorithm
        ccStats; % all information of cc found by algorithm
        
        pixelQual; % pixel quality for each round 
        baseCall; % all refined base calling location 
        
        Nround;
        imgs; % unregistered images
        registeredImages; % registered images
        
        bases; % raw bases for every point ie: '121345'
        baseCS; % matrix of colorseqs for each point ie: [1 2 1 3 4 5]
        qualScores; % base scores for every point

        goodReads; % indices of good points
        goodPoints; % just good points
        goodBases; % bases that are in codebook
        goodScores; % scores of good bases
        
        avgReg;
        outPath;     
        
        nissl;
        dapi;
        cellLocs;
    end
    
    methods
        
        % Pipeline Obj
        function obj = STARMapDataset(basePath, Nround, prefix)
            obj.basePath = basePath;
            obj.outPath = fullfile(obj.basePath, 'output');
            obj.Nround = Nround;
            obj.prefix = prefix;
            fprintf('Pipeline Obj is generated...\n');
            
            if ~exist(obj.outPath, 'dir')
                mkdir(obj.outPath)
            end
        end
        
        % Load Raw Images
        function obj = LoadRawImages(obj, prefix, doNormalize, Nchannel, zrange)
            
            if nargin < 3
                doNormalize = true;
            end
            
            if nargin < 4
                Nchannel = 4;
            end
            
            if nargin < 5
                zrange = [];
            end
            
            % load tiff image stacks
            obj.imgs = LoadImageStacks(obj.basePath, prefix, obj.Nround, Nchannel, zrange);   
            
            % equalize image intensities for each round (fancier normalization)
            if doNormalize
                obj.imgs = EqualizeHist3D(obj.imgs);
            end
            
        end
        
        % Arbitrary Thresholding
        function obj = Thresholding(obj, varargin)
            
            % Input parser
            p = inputParser;
            % Defaults
            defaultfactor = 0.2;
            defaultDoSM = true;
            defaultSigma = 1;
            defaultFilter = [5 5 3];

            addOptional(p,'factor',defaultfactor);
            addOptional(p,'DoSmoothing',defaultDoSM);
            addParameter(p,'Sigma',defaultSigma);
            addParameter(p,'FilterSize',defaultFilter);

            parse(p, varargin{:});
            
            [obj.imgs] = ArbitraryThresholding(obj.imgs, p.Results.factor, 'DoSmoothing', p.Results.DoSmoothing, 'Sigma', p.Results.Sigma, 'FilterSize', p.Results.FilterSize);
        end
        
        % Morphological Reconstruction
        function obj = MorphoRecon(obj, varargin)
            
            % Input parser
            p = inputParser;
            % Defaults
            defaultRadius = 6;
            defaultHeight = 3;


            %addRequired(p,'input_img');
            addOptional(p,'radius',defaultRadius);
            addOptional(p,'height',defaultHeight);


            parse(p, varargin{:});
            
            [obj.imgs] = MorphologicalReconstruction(obj.imgs, p.Results.radius, p.Results.height);
        end
        
        % Global Image Registration 
        function obj = RegisterImages(obj, Ntile)
            fprintf('Registering images...\n');
            [obj.registeredImages,obj.avgReg, ~] = JointRegister3D(obj.imgs, Ntile, false);
            %for i=1:2
                %[registeredImages,avgReg, params] = JointRegister3D(registeredImages,1,true);
            %end

            %[registeredImages,avgReg, params] = JointRegister3D(registeredImages,3,true);
        end
        
        % Local Image Registration 
        function obj = LocalRegisterImages(obj, varargin)
            fprintf('Registering images locally...\n');
            
            % Input parser
            p = inputParser;
            % Defaults
            defaultRef = 1;
            defaultMethod = "max";
            defaultIter = 60;
            defaultAFS = 1;

            %addRequired(p,'object');
            addOptional(p,'ref_round',defaultRef);
            addOptional(p,'method',defaultMethod);
            addParameter(p,'Iterations',defaultIter);
            addParameter(p,'AccumulatedFieldSmoothing',defaultAFS);

            parse(p, varargin{:});
            
            [obj.registeredImages] = LocalRegistration(obj.registeredImages, p.Results.ref_round, 'method', p.Results.method, 'Iterations', p.Results.Iterations, 'AccumulatedFieldSmoothing', p.Results.AccumulatedFieldSmoothing);

        end
        
        % Min-Max nromalization
        function obj = MinMaxNormalize(obj)
            
            obj.imgs = MinMaxNorm(obj.imgs);
            
        end    
        
        % Load Registered Images (Optional)
        function obj = LoadRegisteredImages(obj, reverseChs)
            if nargin < 2
                reverseChs = false;
            end
            obj.registeredImages = LoadRegistered(fullfile(obj.outPath, 'reg3d'), obj.Nround);
            
            if reverseChs                
                for i=1:numel(obj.registeredImages)
                    currRound = obj.registeredImages{i};
                    temp2 = currRound(:,:,:,2);
                    temp3 = currRound(:,:,:,3);
                    currRound(:,:,:,3) = temp2;
                    currRound(:,:,:,2) = temp3;
                    obj.registeredImages{i} = currRound;
                end
            end
        end
        
        % Swap Channel 2 & 3
        function obj = SwapChannels2And3(obj)
            % swap channels 2 and 3 for sequential scane
            for i=1:numel(obj.registeredImages)
    
                temp = obj.registeredImages{i};
                ch2 = temp(:,:,:,2);
                ch3 = temp(:,:,:,3);
                temp(:,:,:,3) = ch2;
                temp(:,:,:,2) = ch3;
                obj.registeredImages{i} = temp;
            end
            fprintf('Channel 2 & 3 have been swapped...\n');
        end
        
        % Load Code Book
        function obj = LoadCodebook(obj, doReverse)
            nargin
            if nargin < 2
                doReverse = true;
            end
            % load hash tables of name -> seq and seq -> name
            % where 'seq' is the string representation of the barcode in colorspace
            [obj.nameToSeq, obj.seqToName] = LoadCodebook(obj.basePath, doReverse);            
            seqStrs = obj.seqToName.keys;
            seqCS = []; % color sequences in matrix form for computing hamming distances
            for i=1:numel(seqStrs)
                seqCS = [seqCS; Str2Colorseq(seqStrs{i})];
            end
            obj.barcodeCS = seqCS; % matrix of color seqs
            obj.barcodeNames = obj.seqToName.values; % cell array of seq names
            obj.barcodeSeqs = seqStrs; % str color seqs
        end
        
        % Save Registered Images (Optional)
        function obj = SaveRegisteredImages(obj)
            reg3DPath = fullfile(obj.outPath, 'reg3d');
            if ~exist(reg3DPath, 'dir')
                mkdir(reg3DPath);
            end
            
            SaveColorMax(reg3DPath, obj.registeredImages);
            SaveRegistered(reg3DPath, obj.registeredImages);            
        end
        
        % Save Thres Images (Optional)
        function obj = SaveThresImages(obj)
            thres3DPath = fullfile(obj.outPath, 'thres');
            if ~exist(thres3DPath, 'dir')
                mkdir(thres3DPath);
            end
            
            SaveThres(thres3DPath, obj.registeredImages);            
        end
        
        % Load CMTK Registered Images (null)
        function obj = LoadCMTKRegisteredImages(obj)
        end
        
        % Extract Reads by Color
        function obj = ExtractReadsByColor(obj, thresholdRel, ksize, kdims, scoreThresh, regionSize, Nchannel)
            if nargin < 6
                regionSize = [3 3 0];
            end
            
            if nargin < 7
                Nchannel = 4;
            end
            
            fprintf('Extracting reads...\n');
            % first find the spots
            % using round 1
            points = FindAllSpots3D({obj.registeredImages{1}},'thresholdRel', thresholdRel, 'ksize', ksize, 'kdims', kdims);
            % using all rounds 
            %points = FindAllSpots3D(obj.registeredImages,'thresholdRel', thresholdRel, 'ksize', ksize, 'kdims', kdims);
            %points = FindAllSpots3D({obj.avgReg},'thresholdRel', thresholdRel, 'ksize', ksize, 'kdims', kdims);            
            obj.allPoints = points.centroids;
            obj.ccStats = points.ccStats;
            %points = FindAllSpots3D(regMax, 'thresholdRel', 0.2, 'ksize', 3.5, 'kdims',[11 11]);
            %obj.allPoints = FindSpots3D(obj.avgReg, 'thresholdRel', thresholdRel, 'ksize', ksize, 'kdims', kdims);
            
            % Determine Base Calling location 
            [obj.pixelQual, obj.baseCall] = RefineBaseCall(obj.registeredImages, obj.ccStats);
            
            % extract base colors from each
            cs = ExtractColorSeq(obj.registeredImages, obj.allPoints, regionSize, Nchannel);
            
            % debug
            %size(cs)
            %cs(1,1,:)

            % convert colors into bases and scores
            [obj.bases, obj.baseCS, obj.qualScores] = GetBaseSeq(cs);
            
            % remove reads with any infinite values
            finiteScores = ~any(isinf(obj.qualScores),2);            
            
            % remove reads with bad quality scores
            totalScores = mean(obj.qualScores(finiteScores,:),2);
            figure(1);
            hist(mean(obj.qualScores,2),100)
            xlabel('Average scores'); ylabel('Count');
            
            % export_fig is an add-on
            export_fig(fullfile(obj.outPath, 'average_scores.png'));

            %
            belowScoreThresh = mean(obj.qualScores,2) < scoreThresh; 
            s = sprintf('%f [%d / %d] percent of reads are below score thresh\n',...
                sum(belowScoreThresh)/numel(belowScoreThresh),...
                sum(belowScoreThresh), ...
                numel(belowScoreThresh));
            fid = fopen(fullfile(obj.outPath, 'stats.txt'), 'w');
            fprintf(s); fprintf(fid, s);
            fclose(fid);

            toKeep = belowScoreThresh & finiteScores; % points to keep;
            obj.bases = obj.bases(toKeep);
            obj.qualScores = obj.qualScores(toKeep,:);
            obj.allPoints = obj.allPoints(toKeep,:);
            obj.baseCS = obj.baseCS(toKeep,:);
            figure(1);
            errorbar(mean(obj.qualScores), std(obj.qualScores),'ko-'); 
            xlim([0 obj.Nround+1]); 
            xlabel('Round'); ylabel('Average qual score');
            export_fig(fullfile(obj.outPath, 'average_qual_score_belowThresh.png'));
            
            bases = obj.bases;
            qualScores = obj.qualScores;
            allPoints = obj.allPoints;
            save(fullfile(obj.outPath, 'points.mat'), 'bases', 'qualScores', 'allPoints');
        end
        
        % Save all gene images
        function obj = SaveAllGeneImages(obj, xsize, ysize, genesToShow)
            
            if nargin < 4                
                genesToShow = obj.nameToSeq.keys;
            end
            if nargin > 2
                figure('position', [100 100 xsize ysize]);           
            else                
                figure(1);
                clf;
            end
            imagesPath = fullfile(obj.outPath, 'images');
            if ~exist(imagesPath, 'dir')
                mkdir(imagesPath);
            end
            for i=1:numel(genesToShow)            
                gene = genesToShow{i};
                clf;
                fprintf('%s\n', gene);
                obj.PlotGene(gene); 
                export_fig(fullfile(imagesPath, [gene '_gene.png']));
            end
            
        end
        
        %
        % functions to actually assign reads to cells
        %
        
        % Load Nissl DAPI
        function obj = LoadNisslDapi(obj,zrange)
            if nargin < 2
                zrange = [];
            end
            nisslImagePath = fullfile(obj.basePath, 'nissl', [obj.prefix 'nissl_decon_ch01.tif']);
            dapiImagePath = fullfile(obj.basePath, 'nissl', [obj.prefix 'nissl_decon_ch00.tif']);
            if isempty(zrange)
                obj.nissl = max(loadtiff(nisslImagePath),[],3);
                obj.dapi = max(loadtiff(dapiImagePath),[],3);
            else
                nissl = loadtiff(nisslImagePath);
                dapi = loadtiff(dapiImagePath);
                obj.nissl = max(nissl(:,:,zrange),[],3);
                obj.dapi = max(dapi(:,:,zrange),[],3);
            end
            imwrite(obj.nissl, fullfile(obj.basePath, 'nissl', 'nissl_maxproj.tif'));
            imwrite(obj.dapi, fullfile(obj.basePath, 'nissl', 'dapi_maxproj.tif'));
        end
        
        % Auto-select DAPI cells 
        function obj = AutomaticallySelectDapiCells(obj)
            t = graythresh(obj.dapi);
            bw = imbinarize(imgaussfilt(obj.dapi,2), t);
            bw = bwareaopen(bw, 100);
            for i=1:3
                bw = imerode(bw, strel('disk', 8));
            end
            rp = regionprops(bw);
            figure, imagesc(obj.dapi); colormap gray; hold on;
            cellLocs = [];
            for i=1:numel(rp)
                plot(rp(i).Centroid(1), rp(i).Centroid(2),'r.'); hold on;
                cellLocs = [cellLocs; rp(i).Centroid];
            end
            obj.cellLocs = cellLocs;
            save(fullfile(obj.basePath, 'output', 'cellLocs.mat'), 'cellLocs');
        end
        
        % Manual-select DAPI cells
        function obj = ManuallySelectDapiCells(obj)
            avgImg = obj.dapi;
            imagesc(avgImg); colormap gray;
            cellLocs = [];
            while true
                [x,y] = ginput(1);
                cellLocs = [cellLocs; x y];
                if x < size(avgImg,2) && y < size(avgImg,1)
                    imagesc(avgImg); colormap gray; hold on;        
                    plot(cellLocs(:,1),cellLocs(:,2),'r.');        
                else
                    break;
                end    
            end            
            obj.cellLocs = cellLocs;
            save(fullfile(obj.basePath, 'output', 'cellLocs.mat'), 'cellLocs');
        end
        
        % Find Cells' Nissl
        function obj = FindCellsNissl(obj)
            t = graythresh(obj.nissl);
            bw = imbinarize(imgaussfilt(obj.nissl,2), t);
            bw = bwareaopen(bw, 100);
            bw = imdilate(bw, strel('disk', 8));
            D = bwdist(bw);
            
            markers = zeros(size(obj.nissl));
            for i=1:size(obj.cellLocs,1)
                markers(round(obj.cellLocs(i,1)), round(obj.cellLocs(i,2))) = 255;
            end
            D = imimposemin(D, markers);
            L = watershed(D);
            imshow(label2rgb(L, 'jet', 'w','shuffle'));
            
        end
        
        % Assign Reads (null)
        function obj = AssignReadsNissl(obj)
        end
        
        % Load points 
        function obj = LoadPoints(obj)
            
            S = load(fullfile(obj.outPath, 'points.mat'));
            obj.bases = S.bases;
            obj.qualScores = S.qualScores;
            obj.allPoints = S.allPoints;

            S = load(fullfile(obj.outPath, 'goodPoints.mat'));
            obj.goodReads = S.goodReads;            
            obj.goodPoints = S.goodPoints; 
            obj.goodBases = S.goodBases;
            obj.goodScores = S.goodScores;            
        end

        % Load numpy points
        function obj = LoadNumpyPoints(obj,regionSize, scoreThresh)
            if nargin < 2
                regionSize = [2 2 1];
            end
            if nargin < 3
                scoreThresh = 0.2;
            end
            ch1 = readNPY(fullfile(obj.outPath, 'ch1_points.npy'));
            ch2 = readNPY(fullfile(obj.outPath, 'ch2_points.npy'));
            ch3 = readNPY(fullfile(obj.outPath, 'ch3_points.npy'));
            ch4 = readNPY(fullfile(obj.outPath, 'ch4_points.npy'));
            temp = [ch1; ch2; ch3; ch4];
            allPoints = zeros(size(temp,1),3);
            allPoints(:,1) = temp(:,3);
            allPoints(:,2) = temp(:,2);
            allPoints(:,3) = temp(:,1);
            obj.allPoints = allPoints; 
            
            cs = ExtractColorSeq(obj.registeredImages, obj.allPoints, regionSize);

            % convert colors into bases and scores
            [obj.bases, obj.baseCS, obj.qualScores] = GetBaseSeq(cs);
            
            % remove reads with any infinite values
            finiteScores = ~any(isinf(obj.qualScores),2);            
            
            % remove reads with bad quality scores
            totalScores = mean(obj.qualScores(finiteScores,:),2);
            figure(1);
            hist(mean(obj.qualScores,2),100)
            xlabel('Average scores'); ylabel('Count');
            export_fig(fullfile(obj.outPath, 'average_scores.png'));

            %
            belowScoreThresh = mean(obj.qualScores,2) < scoreThresh; 
            s = sprintf('%f [%d / %d] percent of reads are below score thresh\n',...
                sum(belowScoreThresh)/numel(belowScoreThresh),...
                sum(belowScoreThresh), ...
                numel(belowScoreThresh));
            fid = fopen(fullfile(obj.outPath, 'stats.txt'), 'w');
            fprintf(s); fprintf(fid, s);
            fclose(fid);

            toKeep = belowScoreThresh & finiteScores; % points to keep;
            obj.bases = obj.bases(toKeep);
            obj.qualScores = obj.qualScores(toKeep,:);
            obj.allPoints = obj.allPoints(toKeep,:);
            obj.baseCS = obj.baseCS(toKeep,:);
            figure(1);
            errorbar(mean(obj.qualScores), std(obj.qualScores),'ko-'); 
            xlim([0 obj.Nround+1]); 
            xlabel('Round'); ylabel('Average qual score');
            export_fig(fullfile(obj.outPath, 'average_qual_score_belowThresh.png'));
            
            bases = obj.bases;
            qualScores = obj.qualScores;
            allPoints = obj.allPoints;
            save(fullfile(obj.outPath, 'points.mat'), 'bases', 'qualScores', 'allPoints');            
        end
        
        % Extract Reads by Position
        function obj = ExtractReadsByPosition(obj, thresholdRel, ksize, kdims)
            % try to extract reads by matching dots across rounds
            % doesn't work as well as intensity based matching
            points = FindAllSpots3D(obj.registeredImages,...
                'thresholdRel', thresholdRel,...
                'ksize', ksize,...
                'kdims',kdims);
            Nround = max(points.seqRound);
            Npoint = sum(points.seqRound==1);
            seqVals = zeros(Npoint, Nround);
            seqVals(:,1) = points.colorCh(points.seqRound==1); % get color vals for first sequencing round
            distThresh = 4;
            firstRound = points.centroids(points.seqRound==1,:);
            badIdx= zeros(Npoint,1);
            for i=2:Nround
                currIdx = find(points.seqRound==i);
                currRound = points.centroids(currIdx,:);
                [matchIdx,dist] = knnsearch(currRound, firstRound, 'k', 1); % find matching index in currRound and distance for each point in firstRound
                disp(['Round ' num2str(i) ' ' num2str(sum(dist < distThresh)/numel(dist)) '% match']);
                for j=1:numel(matchIdx)
                    if dist(j) < distThresh
                        seqVals(j,i) = points.colorCh(currIdx(matchIdx(j)));                    
                    else
                        badIdx(j) = 1;
                    end
                end
            end

            badIdx = boolean(badIdx);
            fprintf('%f total spots match from first round\n', 1-sum(badIdx)/length(badIdx));            
            obj.allPoints = firstRound(badIdx == 0,:);
        end
        
        % Read Filtration 
        function obj = FilterReads(obj, endBases, correctErrors)
            if nargin < 2
                endBases = ['C'; 'C'];
            end            
            if nargin < 4
                correctErrors = false;
            end
            fid = fopen(fullfile(obj.outPath, 'stats.txt'),'a');

            fprintf('Filtering reads...\n');
            % Filter reads by whether they are in thecodebook
            % this is just as a sanity check; reads are actually filtered
            % by whether they are in the codebook
            filtBases = obj.baseCS;
            correctSeqs = zeros(size(filtBases,1),1); % count sequences that are of correct form
            
            
            % older code for filtering reads in sequence space
            bases = DecodeBases(filtBases, endBases(1));
            for i=1:numel(bases)
                currSeq = bases{i};
                currSeq;
                if currSeq(1) == endBases(1) && currSeq(end) == endBases(2)
                    correctSeqs(i) = 1;
                end
            end
            s = sprintf('%f [%d / %d] percent of good reads are %sNNNNN%s\n',...
                sum(correctSeqs)/size(filtBases,1),...
                sum(correctSeqs),...
                size(filtBases,1),...
                endBases(1),...
                endBases(2));
            fprintf(s); fprintf(fid, s);
            
            
            % filter reads based on codebook
            Nbases = numel(obj.bases);
            inCodebook = zeros(Nbases,1);
            codebookSeqs = obj.seqToName.keys;
            for s=1:Nbases
                str = obj.bases{s};
                if ismember(str, codebookSeqs)         
                    inCodebook(s) = 1;
                else
                    inCodebook(s) = 0;
                end
            end
            
            
            % correct things in colospace
            if correctErrors
                Ncorrected = 0;            
                upd = textprogressbar(sum(inCodebook==0));
                correctionCount = 0;

                qualScores = obj.qualScores;
                baseCS = obj.baseCS;
                barcodeCS = obj.barcodeCS;
                Nround = obj.Nround;

                correctedIdx = zeros(Nbases,1);
                for s=1:Nbases                   
                    if inCodebook(s) == 0
                        correctionCount = correctionCount + 1;
                        if mod(s,100) == 0
                            upd(correctionCount);
                        end

                        currBase = baseCS(s,:);
                        
                        dists = pdist2(currBase, barcodeCS, 'hamming')*Nround; % calculate hamming distance between currBase and each barcode
                        minDist = min(dists);
                        if minDist == 1
                            potentialMatches = find(dists == minDist);
                            if numel(potentialMatches) > 1 % if multiple matches, choose the one with the lowest quality score for the mismatched position
                                %{
                                currQuals = qualScores(s,:);
                                mismatchedPositionScores = zeros(numel(potentialMatches),1);       
                                for matchIdx=1:numel(potentialMatches)
                                    currBarcode = barcodeCS(potentialMatches(matchIdx),:);
                                    for baseIdx=1:numel(currBarcode)
                                        if currBarcode(baseIdx) ~= currBase(baseIdx)
                                            mismatchedPositionScores(matchIdx) = currQuals(baseIdx);
                                        end
                                    end
                                    % compare 
                                end
                                correctedIdx(s) = potentialMatches(mismatchedPositionScores == min(mismatchedPositionScores));
                                inCodebook(s) = 1;
                                %}
                            else
                                correctedIdx(s) = potentialMatches; % look up correct seq                                                                        
                                inCodebook(s) = 1;
                            end                
                        end
                    end
                end
                Ncorrected = sum(correctedIdx > 0);
                disp(Ncorrected);
            
                % actually reassign bases
                reassigned = find(inCodebook(s) == 0 & correctedIdx(s) > 0);
                for s=reassigned                    
                    obj.bases{s} = obj.barcodeSeqs(correctedIdx(s));                    
                end
            end
            % only save reads that are below the score thresh and are in
            % the codebook
            goodReads = inCodebook==1;
            obj.goodReads = goodReads;
            obj.goodPoints = obj.allPoints(goodReads,:);
            obj.goodBases = obj.bases(goodReads);
            obj.goodScores = obj.qualScores(goodReads,:);

            goodPoints = obj.goodPoints; 
            goodBases = obj.goodBases;
            goodScores = obj.goodScores;
            save(fullfile(obj.outPath, 'goodPoints.mat'), 'goodReads', 'goodPoints', 'goodBases', 'goodScores');

            figure(1);
            errorbar(mean(obj.goodScores), std(obj.goodScores),'ko-'); 
            xlim([0 obj.Nround+1]); 
            xlabel('Round'); ylabel('Average qual score');
            export_fig(fullfile(obj.outPath, 'average_qual_score_inCodebook.png'));


            s = sprintf('%f [%d / %d] percent of good reads are in codebook\n',...
                sum(goodReads)/Nbases,...
                sum(goodReads),...
                Nbases);
            fprintf(s); fprintf(fid, s);

            if correctErrors
                s = sprintf('%f [%d / %d] reads in codebook were corrected\n', ...
                    Ncorrected/sum(inCodebook), ...
                    Ncorrected, ...
                    sum(inCodebook));
                fprintf(s); fprintf(fid, s);
            end


                s = sprintf('%f [%d / %d] percent of %sNNNNN%s reads are in codebook\n',...
                    sum(goodReads)/sum(correctSeqs),...
                    sum(goodReads), ...
                    sum(correctSeqs),...
                    endBases(1),...
                    endBases(2));            

            fprintf(s); fprintf(fid, s);
            fclose(fid);

             
        end
        
        % Get points for gene
        function p = GetPointsForGene(obj,geneName)            
            matchIdx = GetSeqIdx(obj.goodBases,obj.nameToSeq(geneName));            
            p = obj.goodPoints(matchIdx==1,:);
        end
        
        % Plot barcodes
        function PlotBarcodes(obj)
            seqs = obj.seqToName.keys;
            N  = ceil(sqrt(numel(seqs)));
            for i=1:numel(seqs)
                matchIdx = GetSeqIdx(obj.goodBases,seqs{i});
                %imagesc(max(avgReg,[],3));colormap gray; hold on;
                subplot(N,N,i);    
                %imagesc(max(avgReg,[],3)); colormap gray; hold on;
                plot(obj.goodPoints(matchIdx==1,1), obj.goodPoints(matchIdx==1,2),'k.','markersize',0.1); 
                title(obj.seqToName(seqs{i}));
                axis off;
            end        
            set(gcf,'color','w');
        end
        
        % Plot gene
        function PlotGene(obj, geneName, color, showBackground)
            if nargin < 3
                color = 'k.';
            end
            if nargin < 4
                showBackground = true;
            end
            matchIdx = GetSeqIdx(obj.goodBases,obj.nameToSeq(geneName));            
            if showBackground
                imagesc(obj.nissl/5,[0 255]); colormap(flipud(gray)); hold on;
            else
            %imagesc(max(obj.avgReg,[],3)); colormap(flipud(gray)); hold on;
                imagesc(max(max(obj.registeredImages{1},[],4),[],3)); colormap(flipud(gray)); hold on;
            end
            plot(obj.goodPoints(matchIdx==1,1), obj.goodPoints(matchIdx==1,2),'k.','markersize',15); 
            %plot(obj.goodPoints(matchIdx==1,1), obj.goodPoints(matchIdx==1,2),'r.','markersize',10); 
            %title(geneName);
            axis off;
        end

        % Save top genes
        function counts = SaveTopGenes(obj)
            geneNames = obj.nameToSeq.keys;
            counts = containers.Map(geneNames, zeros(numel(geneNames),1));
            for i=1:numel(obj.goodBases)               
                counts(obj.seqToName(obj.goodBases{i})) = counts(obj.seqToName(obj.goodBases{i})) + 1;
            end
            fid = fopen(fullfile(obj.outPath, 'gene_counts.csv'),'w');
            for j=counts.keys
                fwrite(fid, sprintf('%s,%s,%d\n', j{1}, obj.nameToSeq(j{1}), counts(j{1})));
            end
            fclose(fid);
        end
        
        % Make average image
        function obj = MakeAverageImage(obj)
            [X,Y,Z,~] = size(obj.registeredImages{1});
            obj.avgReg = zeros(X,Y,Z);
            for i=1:numel(obj.registeredImages)
                obj.avgReg = obj.avgReg + single(squeeze(max(obj.registeredImages{i},[],4)));
            end
            obj.avgReg = obj.avgReg / numel(obj.registeredImages);    
        end
        
        % Save Points
        function SavePoints(obj)
        end
        
        % Change all.Points
        function obj = changePoints(obj, input_val, regionSize, scoreThresh)
            obj.allPoints = input_val;
                        
            % extract base colors from each
            cs = ExtractColorSeq(obj.registeredImages, obj.allPoints, regionSize);
            
            % debug
            %size(cs)
            %cs(1,1,:)

            % convert colors into bases and scores
            [obj.bases, obj.baseCS, obj.qualScores] = GetBaseSeq(cs);
            
            % remove reads with any infinite values
            finiteScores = ~any(isinf(obj.qualScores),2);            
            
            % remove reads with bad quality scores
            totalScores = mean(obj.qualScores(finiteScores,:),2);
            figure(1);
            hist(mean(obj.qualScores,2),100)
            xlabel('Average scores'); ylabel('Count');
            
            % export_fig is an add-on
            export_fig(fullfile(obj.outPath, 'average_scores.png'));

            %
            belowScoreThresh = mean(obj.qualScores,2) < scoreThresh; 
            s = sprintf('%f [%d / %d] percent of reads are below score thresh\n',...
                sum(belowScoreThresh)/numel(belowScoreThresh),...
                sum(belowScoreThresh), ...
                numel(belowScoreThresh));
            fid = fopen(fullfile(obj.outPath, 'stats.txt'), 'w');
            fprintf(s); fprintf(fid, s);
            fclose(fid);

            toKeep = belowScoreThresh & finiteScores; % points to keep;
            obj.bases = obj.bases(toKeep);
            obj.qualScores = obj.qualScores(toKeep,:);
            obj.allPoints = obj.allPoints(toKeep,:);
            obj.baseCS = obj.baseCS(toKeep,:);
            figure(1);
            errorbar(mean(obj.qualScores), std(obj.qualScores),'ko-'); 
            xlim([0 obj.Nround+1]); 
            xlabel('Round'); ylabel('Average qual score');
            export_fig(fullfile(obj.outPath, 'average_qual_score_belowThresh.png'));
            
            bases = obj.bases;
            qualScores = obj.qualScores;
            allPoints = obj.allPoints;
            save(fullfile(obj.outPath, 'points.mat'), 'bases', 'qualScores', 'allPoints');
        end
        
    end
end
