%Create paradigm file for each stream 
clear



condnames = {'PREPOSTFIX','AUD_LEFT','AUD_RIGHT','VIS_LEFT','VIS_RIGHT','PASSIVE','FIXATION','CUE'};

SUBJ = '130306WW';
%mkdir(['./Paradigms/' SUBJ]);


for PN = 5:5

    clearvars -except PN SUBJ condnames

    %PAR_FILE_NAME = './Paradigms/modlocsw-244670.par';
    %PAR_OUT =       './Paradigms/modlocsw-paradigm1.par';

    %fn = strcat('Results/OUTPUT_', SUBJ, '_',num2str(PN), '*.mat');

    [r,r]= system(strcat('ls Results/MODLOC_OUTPUT_', SUBJ, '_',num2str(PN), '_', '*.mat'));
    r=strread(r,'%s');
    r=sort(r);
    rr = [r{2} ' ' r{1}];
    load(rr);

    PAR_OUT =       ['./Paradigms/' SUBJ '/blockmls-par' num2str(PN) '.par'];


    % TR and NUM_TP from loaded file
    TASK_TIME = NUM_TP * TR;


    %First go through and convert the paradigm to be every 2.6s
    newparconds = 7*ones(1,length(COND_ALL)*2);
    newparconds(2:2:end)=COND_ALL;
    newparconds = [0 newparconds 0];
    duration = []; %Duration of each stimulus (column 3
    onset = 0;


    for i = 1:length(newparconds)
        condtext{i} = condnames{1+newparconds(i)};

        if newparconds(i) == 0
            %as many times as FIX_TP for fixation
            duration(i) = FIXATION_TIME;

        elseif newparconds(i) == 7
            %Once for cue
            duration(i) = CUE_TIME;
        else
            %the rest are task blocks, so use BLOCK_TP
            duration(i) = BLOCK_TIME;
        end

        if i>1
            onset(i) = onset(i-1)+duration(i-1);
        end
    end

    fid = fopen(PAR_OUT, 'w');

    for i=1:length(newparconds)

        fprintf(fid, '%3.1f\t%i\t%3.1f\t%1.1f\t%s\n',onset(i),newparconds(i),duration(i),1.0,condtext{i});
    end

    fclose(fid);



end

