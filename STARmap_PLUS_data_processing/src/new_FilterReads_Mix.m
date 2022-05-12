function obj = new_FilterReads_Mix( obj, endBases, endBases_mix, showPlots )
% FilterReads
  
%     if nargin < 4
%         correctErrors = false;
%     end

    % Filter reads by whether they are in the codebook
    % This is just as a sanity check; reads are actually filtered
    % by whether they are in the codebook
    filtBases = obj.basecsMat;
    correctSeqs = zeros(size(filtBases, 1), 1); % count sequences that are of correct form
    correctSeqs_mix = zeros(size(filtBases, 1), 1);
    correctSeqs_T = zeros(size(filtBases, 1), 1);
    correctSeqs_G = zeros(size(filtBases, 1), 1);
    
    % older code for filtering reads in sequence space
    bases = new_DecodeCS(filtBases, endBases(1));
    
    for i=1:numel(bases)
        currSeq = bases{i};
        if currSeq(1) == endBases(1) && currSeq(end) == endBases(2)
            correctSeqs(i) = 1;
        elseif currSeq(1) == endBases_mix(1) && currSeq(end) == endBases_mix(2)
            correctSeqs_mix(i) = 1;
        elseif currSeq(1) == endBases_mix(1) && currSeq(end) == 'T'
            correctSeqs_T(i) = 1;
        elseif currSeq(1) == endBases_mix(1) && currSeq(end) == 'G'
            correctSeqs_G(i) = 1;
        end
    end
    
    % form 1
    fprintf('Filtration Statistics:\n');
    fprintf(obj.log, 'Filtration Statistics:\n');
    score_1 = sum(correctSeqs)/size(filtBases,1);
    s = sprintf('%f [%d / %d] percent of good reads are %sNNNNN%s\n',...
        sum(correctSeqs)/size(filtBases,1),...
        sum(correctSeqs),...
        size(filtBases,1),...
        endBases(1),...
        endBases(2));
    fprintf(s);
    fprintf(obj.log, s);

    % form 2
    % score_2 = sum(correctSeqs_mix)/size(filtBases,1);
    s = sprintf('%f [%d / %d] percent of good reads are %sNNNNN%s\n',...
        sum(correctSeqs_mix)/size(filtBases,1),...
        sum(correctSeqs_mix),...
        size(filtBases,1),...
        endBases_mix(1),...
        endBases_mix(2));
    fprintf(s);
    fprintf(obj.log, s);
    
    % form 3
    % score_2 = sum(correctSeqs_mix)/size(filtBases,1);
    s = sprintf('%f [%d / %d] percent of good reads are %sNNNNNT\n',...
        sum(correctSeqs_T)/size(filtBases,1),...
        sum(correctSeqs_T),...
        size(filtBases,1),...
        endBases_mix(1));
    fprintf(s);
    fprintf(obj.log, s);
    
    % form 4
    % score_2 = sum(correctSeqs_mix)/size(filtBases,1);
    s = sprintf('%f [%d / %d] percent of good reads are %sNNNNNG\n',...
        sum(correctSeqs_G)/size(filtBases,1),...
        sum(correctSeqs_G),...
        size(filtBases,1),...
        endBases_mix(1));
    fprintf(s);
    fprintf(obj.log, s);
    
    % filter reads based on codebook
    Nreads = numel(obj.allReads);
    inCodebook = zeros(Nreads,1);
    codebookSeqs = obj.barcodeSeqs;
    for s=1:Nreads
        str = obj.allReads{s};
        if ismember(str, codebookSeqs)         
            inCodebook(s) = 1;
        else
            inCodebook(s) = 0;
        end
    end

    % only save reads that are below the score thresh and are in
    % the codebook
    readsToKeep = inCodebook==1;
    obj.goodSpots = obj.allSpots(readsToKeep,:);
    obj.goodReads = obj.allReads(readsToKeep);
    obj.goodScores = obj.allScores(readsToKeep,:);

    if showPlots
        figure(1);
        errorbar(mean(obj.goodScores), std(obj.goodScores),'ko-'); 
        xlim([0 obj.Nround+1]); 
        xlabel('Round'); ylabel('Average qual score');
    end

    score_2 = sum(readsToKeep)/Nreads;
    s = sprintf('%f [%d / %d] percent of good reads are in codebook\n',...
        sum(readsToKeep)/Nreads,...
        sum(readsToKeep),...
        Nreads);
    fprintf(s);
    fprintf(obj.log, s);

%     if correctErrors
%         s = sprintf('%f [%d / %d] reads in codebook were corrected\n', ...
%             Ncorrected/sum(inCodebook), ...
%             Ncorrected, ...
%             sum(inCodebook));
%         fprintf(s); fprintf(fid, s);
%     end

    score_3 = sum(readsToKeep)/sum(correctSeqs);
    form1_goodReads = inCodebook == 1 & correctSeqs == 1;
    s = sprintf('%f [%d / %d] percent of %sNNNNN%s reads are in codebook\n',...
        sum(form1_goodReads)/sum(correctSeqs),...
        sum(form1_goodReads), ...
        sum(correctSeqs),...
        endBases(1),...
        endBases(2));            
    fprintf(s);
    fprintf(obj.log, s);
    
    % score_3 = sum(readsToKeep)/sum(correctSeqs);
    form2_goodReads = inCodebook == 1 & correctSeqs_mix == 1;
    s = sprintf('%f [%d / %d] percent of %sNNNNN%s reads are in codebook\n',...
        sum(form2_goodReads)/sum(correctSeqs_mix),...
        sum(form2_goodReads), ...
        sum(correctSeqs_mix),...
        endBases_mix(1),...
        endBases_mix(2));            
    fprintf(s);
    fprintf(obj.log, s);
    
    obj.FilterScores = [score_1 score_2 score_3]; 
    
end