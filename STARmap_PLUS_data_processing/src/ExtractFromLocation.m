function obj = ExtractFromLocation( obj, voxelSize )
%ExtractFromLocation

    % get dims
    [x, y, z, Nchannel, Nround] = size(obj.registeredImages);
    Npoint = size(obj.allSpots,1);
    colorSeq = zeros(Npoint, Nround, Nchannel); % "color value" of each dot in each channel of each sequencing round 

    fprintf('Geting color sequence for each voxel...\n');

    for i=1:Npoint
        
        % Get voxel for each dot
        curr_p = obj.allSpots(i,:);
        extentsX = GetExtents(curr_p(2), voxelSize(1), x);
        extentsY = GetExtents(curr_p(1), voxelSize(2), y);                    
        extentsZ = GetExtents(curr_p(3), voxelSize(3), z);    
        
        
        for r=1:Nround
            currVol = obj.registeredImages(extentsX, extentsY, extentsZ, :, r); % 4-D array
            
            colorSeq(i,r,:) = single(squeeze(sum(currVol, [1 2 3]))); % sum along row,col,z
            colorSeq(i,r,:) = colorSeq(i,r,:) ./ (sqrt(sum(squeeze(colorSeq(i,r,:)).^2)) + 1E-6); % +1E-6 avoids denominator equaling to 0
        end
        

    end
    
    
    fprintf('\n');
    
    [obj.allReads, obj.basecsMat, obj.allScores] = new_GetBaseSeq(colorSeq);

    % remove reads with any infinite values
    finiteScores = ~any(isinf(obj.allScores), 2);    
    Ndropped = sum(any(isinf(obj.allScores), 2));
    
    obj.allReads = obj.allReads(finiteScores);
    obj.allScores = obj.allScores(finiteScores,:);
    obj.allSpots = obj.allSpots(finiteScores,:);
    obj.basecsMat = obj.basecsMat(finiteScores,:);
    
    s = sprintf('Number of spots were dropped because of multi-max color: %d\n', Ndropped);
    fprintf(s);
    fprintf(obj.log, s);
    
% 
%     % remove reads with bad quality scores
%     if showPlots
%         figure(1);
%         histogram(mean(allScores, 2), 100)
%         xlabel('Average scores'); ylabel('Count');
%     end
%     
%     % export_fig is an add-on
%     % export_fig(fullfile(obj.outPath, 'average_scores.png'));
% 
% 
%     belowScoreThresh = mean(allScores, 2) < q_score_thers; % 0.5; 
%     s = sprintf('%f [%d / %d] percent of reads are below score thresh %d\n',...
%         sum(belowScoreThresh)/numel(belowScoreThresh),...
%         sum(belowScoreThresh), ...
%         numel(belowScoreThresh), ...
%         q_score_thers);
%     fprintf(s);
%     
%     if ~isempty(obj_log)
%         fprintf(obj_log, s);
%     end
%     
%     
% %     fid = fopen(fullfile(obj.outPath, 'stats.txt'), 'w');
% %     fprintf(s); fprintf(fid, s);
% %     fclose(fid);
%     
% %     if ~isempty(dividerBase)
% %         k = {'AT','CT','GT','TT',...
% %             'AG', 'CG', 'GG', 'TG',...
% %             'AC', 'CC', 'GC', 'TC',...
% %             'AA', 'CA', 'GA', 'TA'};
% %         v = {4,3,2,1,3,4,1,2,2,1,4,3,1,2,3,4};
% %         codeMap = containers.Map(k, v);  
% %         
% %         dividerColor = codeMap(dividerBase);
% %         divideLoc = Nround / 2;
% %         allReads = cellfun(@(x) insertAfter(x, divideLoc, int2str(dividerColor)), allReads, 'UniformOutput', false);
% %         divideMat = repmat(dividerColor, size(csMat, 1), 1);
% %         csMat = [csMat(:, 1:divideLoc) divideMat csMat(:, divideLoc+1:end)];
% %     end
%     
%    
%     if showPlots
%         figure(2);
%         errorbar(mean(allScores), std(allScores),'ko-'); 
%         xlim([0 Nround+1]); 
%         xlabel('Round'); ylabel('Average qual score');
%     end
%     
% %     export_fig(fullfile(obj.outPath, 'average_qual_score_belowThresh.png'));
% %     save(fullfile(obj.outPath, 'points.mat'), 'allReads', 'qualScores', 'allPoints');

end

function e = GetExtents(pos, voxelSize, lim)

if pos-voxelSize < 1 
    e1 = 1;
else
    e1 = pos-voxelSize;
end

if pos+voxelSize > lim
    e2 = lim;
else
    e2 = pos+voxelSize;
end

e = e1:e2;

end

