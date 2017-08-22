%Create a task regressor file for auditory oddball that is happening in
%some of the visual attention streams (needs to be same number of TPs and
%one column
clear



condnames = {'PREPOSTFIX','AUD_LEFT','AUD_RIGHT','VIS_LEFT','VIS_RIGHT','PASSIVE','FIXATION','CUE'};

SUBJ = '130217CC';
%mkdir(['./Paradigms/' SUBJ]);


for PN = 1:4

    clearvars -except PN SUBJ condnames

    %PAR_FILE_NAME = './Paradigms/modlocsw-244670.par';
    %PAR_OUT =       './Paradigms/modlocsw-paradigm1.par';

    %fn = strcat('Results/OUTPUT_', SUBJ, '_',num2str(PN), '*.mat');

    [r,r]= system(strcat('ls Results/MODLOC_OUTPUT_', SUBJ, '_',num2str(PN), '_', '*.mat'));
    r=strread(r,'%s');
    r=sort(r);
    rr = [r{2} ' ' r{1}];
    load(rr);

    PAR_OUT =       ['./Paradigms/' SUBJ '/bmlsaudoddball-par' num2str(PN) '.par'];
    
    audoddball = zeros(NUM_TP,1);
    
    for b = 1:length(COND_ALL)
        if COND_ALL(b) == 3 || COND_ALL(b) == 4
            al = TRI_CHAR_IND(b,1,:);
            ar = TRI_CHAR_IND(b,2,:);
            if sum(al==ar)>0              
                fprintf('block is %i',b)
                bad = find(al==ar)
                audoddball(ceil(stimprestimes(b,bad)/TR))=1;
            end

        end
    end
    
    fprintf('total bad timepoints is %i\n',sum(audoddball))

    fid = fopen(PAR_OUT, 'w');



    fprintf(fid, '%i\n',audoddball);
    

    fclose(fid);



end

