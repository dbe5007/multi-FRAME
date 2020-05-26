% Calculate accuracy by behavioral trial type (e.g., Hit acc, Miss acc)
            switch taskInfo.accuracyFlag
                case 'Yes'
                    
                    % Pull trial information from Mask by target condition
                    Cond1 = ~cellfun(@isempty, strfind(currDataset.sa.labels,...
                        taskInfo.Conditions{1,1}));
                    Cond1idx = find(Cond1 == 1);
                    Cond2 = ~cellfun(@isempty, strfind(currDataset.sa.labels,...
                        taskInfo.Conditions{1,2}));
                    Cond2idx = find(Cond2 == 1);
                    trialInfo.Cond1 = currDataset.sa.labels(Cond1idx);
                    trialInfo.Cond2 = currDataset.sa.labels(Cond2idx);
                    
                    switch bootstrap.flag
                        case 'No'
                            
                            % Pull accuracy information by target condition
                            acc.Cond1 = correct.(regionName)(Cond1idx');
                            acc.Cond2 = correct.(regionName)(Cond2idx');
                            acc.Cond1 = acc.Cond1';
                            acc.Cond2 = acc.Cond2';
                            
                            % Find behavioral trials for each condition
                            acc.rec1 = ~cellfun(@isempty, ...
                                strfind(trialInfo.Cond1,'Hit'));
                            acc.fam1 = ~cellfun(@isempty, ...
                                strfind(trialInfo.Cond1,'Familiar'));
                            acc.miss1 = ~cellfun(@isempty, ...
                                strfind(trialInfo.Cond1,'Miss'));
                            
                            acc.rec2 = ~cellfun(@isempty, strfind(trialInfo.Cond2,'Hit'));
                            acc.fam2 = ~cellfun(@isempty, strfind(trialInfo.Cond2,'Familiar'));
                            acc.miss2 = ~cellfun(@isempty, strfind(trialInfo.Cond2,'Miss'));
                            
                            acc.recFinal1 = acc.Cond1(acc.rec1==1);
                            acc.famFinal1 = acc.Cond1(acc.fam1==1);
                            acc.missFinal1 = acc.Cond1(acc.miss1==1);
                            
                            acc.recFinal2 = acc.Cond2(acc.rec2==1);
                            acc.famFinal2 = acc.Cond2(acc.fam2==1);
                            acc.missFinal2 = acc.Cond2(acc.miss2==1);
                            
                            if iteration==1
                                finalTableTrialTypeAcc=cell...
                                    (length(subjects)+1,length(rois)+2);
                                finalTableTrialTypeAcc{1,1}='subjectid';
                                finalTableTrialTypeAcc{1,2}='Trial Type';
                                tmpCnt=3;
                                
                                for headerTrialType=1:length(rois)
                                    finalTableTrialTypeAcc{1,tmpCnt}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_trainedRecAccuracy');
                                    finalTableTrialTypeAcc{1,tmpCnt+1}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_novelRecAccuracy');
                                    
                                    finalTableTrialTypeAcc{1,tmpCnt+2}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_trainedRecNumTrials');
                                    finalTableTrialTypeAcc{1,tmpCnt+3}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_novelRecNumTrials');
                                    
                                    finalTableTrialTypeAcc{1,tmpCnt+4}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_trainedFamAccuracy');
                                    finalTableTrialTypeAcc{1,tmpCnt+5}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_novelFamAccuracy');
                                    
                                    finalTableTrialTypeAcc{1,tmpCnt+6}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_trainedFamNumTrials');
                                    finalTableTrialTypeAcc{1,tmpCnt+7}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_novelFamNumTrials');
                                    
                                    finalTableTrialTypeAcc{1,tmpCnt+8}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_trainedMissAccuracy');
                                    finalTableTrialTypeAcc{1,tmpCnt+9}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_novelMissAccuracy');
                                    
                                    finalTableTrialTypeAcc{1,tmpCnt+10}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_trainedMissNumTrials');
                                    finalTableTrialTypeAcc{1,tmpCnt+11}=strcat...
                                        (rois{1,headerTrialType}(1:end-4),'_novelMissNumTrials');
                                    
                                    tmpCnt=tmpCnt+12;
                                end
                                
                                row=2;
                                headerTrialType=3;
                                clear tempcount;
                                TrialTypeCombo = {strcat(taskInfo.Conditions{1},...
                                    '_v_',taskInfo.Conditions{2})};
                            end
                            
                            % Add subject, trial type, and accuracy to table
                            finalTableTrialTypeAcc{row,1}=subject;
                            finalTableTrialTypeAcc{row,2}=TrialTypeCombo{1,1};
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (nanmean(acc.recFinal1));
                            headerTrialType=headerTrialType+1;
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (nanmean(acc.recFinal2));
                            headerTrialType=headerTrialType+1;
                            
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (length(acc.recFinal1) - sum(isnan(acc.recFinal1)));
                            headerTrialType=headerTrialType+1;
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (length(acc.recFinal2) - sum(isnan(acc.recFinal2)));
                            headerTrialType=headerTrialType+1;
                            
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (nanmean(acc.famFinal1));
                            headerTrialType=headerTrialType+1;
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (nanmean(acc.famFinal2));
                            headerTrialType=headerTrialType+1;
                            
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (length(acc.famFinal1) - sum(isnan(acc.famFinal1)));
                            headerTrialType=headerTrialType+1;
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (length(acc.famFinal2) - sum(isnan(acc.famFinal2)));
                            headerTrialType=headerTrialType+1;
                            
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (nanmean(acc.missFinal1));
                            headerTrialType=headerTrialType+1;
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (nanmean(acc.missFinal2));
                            headerTrialType=headerTrialType+1;
                            
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (length(acc.missFinal1) - sum(isnan(acc.missFinal1)));
                            headerTrialType=headerTrialType+1;
                            finalTableTrialTypeAcc{row,headerTrialType}=num2str...
                                (length(acc.missFinal2) - sum(isnan(acc.missFinal2)));
                            headerTrialType=headerTrialType+1;
                            
                    end
                    
            end
            
            %         switch bootstrap.flag
%             case 'No'
%                 save([fileparts(outputPath) filesep 'finalTableTrialType.mat'],...
%                     'finalTableTrialTypeAcc')
%         end

%                 % Save classifier accuracy by trial type/behavior to CSV
%                 file = fopen([out_path '_MVPA' filesep 'allAccuraciesByTrialType.csv'], 'w');
%                 
%                 for a=1:size(finalTableTrialTypeAcc,1)
%                     for b=1:size(finalTableTrialTypeAcc,2)
%                         var = eval('finalTableTrialTypeAcc{a,b}');
%                         try
%                             fprintf(file, '%s', var);
%                         end
%                         fprintf(file, ',');
%                     end
%                     fprintf(file, '\n');
%                 end
%                 fclose(file);