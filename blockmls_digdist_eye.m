function [ ] = blockmls_digdist_eye(varargin )
%MODLOCsw Modality and location attention switching study
% BLOCKED VERSION WITH DISTRACTORS THAT ARE ALL DIGITS
%   Input: subjectID, run number, paradigm number (1 or 2, default = 1), initial stream
%   (ex: 'AL' for auditory left, defaul is 'VL')
%   TP  = 138 (130 with 4 TRs of fix on both sides)
%   TR = 2.6 

%Warn people that there are no Q's or I's

%NOTES: 
% Moved aud streams farther apart, may even do this a little more
% Could increase number of targets, or make it more flexible
% There is currently an error in this code that results in two streams
% presenting the same digit on occasion. 



%%%%%%%%%%%%%%%%%%%%%
%Timing:
% Paradigm file is 140 TR    

% Task: Press a button to match the number in the attended stream
% 
% 6 blocks at 10 TRs each with 1 TR cue
% 12 blocks total
% switching in all directions
% button presses for attended digit
% target every 3-4 seconds 
% 
% need to record which block you are in on a mistake


% Record results
% Record auditory cues
% Make new paradigm files for the analysis

%Results:
%  Correct/Incorrect by block, condition, overall
%  FP, CR, H, M, and k score 
%  Reaction times by block, condition, and overall
% Record all stimulus and response times




% Options to be turned into function arguments later
extravis  = 1; %1 for extra visual distractor streams
overlap =0;


options.location = 'laptop'; %'Harvard';
PAR_NUM = 0;
subjID = 'TEST';
FIX_TP = 4; %4; %Initial fixation, not including block of fixation
CUE_TP = 1;
BLOCK_TP = 10;
NUM_TP = 140; 
TR = 2.6;

options.eyetracker  = false;
options.eyecalibration = false; %Run calibration on eyetracker (every other run)
blocksovertime = 0;


for arg = 1:2:nargin    
  switch varargin{arg}
      % Override simulation parameters
      case {'p','par','parnum','PAR_NUM'} 
          PAR_NUM   = varargin{arg+1};
      case {'s' ,'subj','subject'}
          subjID    = varargin{arg+1};
      case 'FIX_TP'
          FIX_TP = varargin{arg+1};
      case {'Location','location','L','l'}
          options.location = varargin{arg+1};
      case {'eye', 'e', 'eyetracker','Eyelink'}
          options.eyetracker = varargin{arg+1};
      case {'eyecal','calibration','eyecalibration'}
          options.eyecalibration = varargin{arg+1};
  end
end

if options.eyetracker == true
    edfFile = input('Input unique name for saving Eyelink data, only 8 or less letters and num allowed:', 's');
    edfFile = strcat(edfFile,'.edf');
    options.driftcorrect = false; % Drift correct
    p.v_dist 		= 88;	% 88cm for 12 and 32 channel,  viewing distance behavioral: 60(cm)
    p.mon_width  	= 40;	% 40 scanner, 28 mac air, 43 tammy office, 51 my office horizontal dimension of viewable Screen (cm)
end

FIXATION_TIME = FIX_TP * TR;
CUE_TIME = CUE_TP * TR;
BLOCK_TIME = BLOCK_TP * TR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % Data collection initalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kbnumber = initializeKeyboard();  

%These will be the defaults, but overwritten by load 'MGH_responsesettings.mat'

ans1 = KbName('1!'); 
ans2 = KbName('2@');   
ans3 = KbName('3#');  
ans4 = KbName('4$'); 

kbnumber = 1; %initializeKeyboard(); 

if strcmp(options.location,'laptop')
    %load exproom_responsesettings.mat
    load LAPTOP_responsesettings.mat
elseif strcmp(options.location,'desktop')
    load SamKey_responsesettings.mat
elseif strcmp(options.location,'MGH')
    load MGH_responsesettings.mat
elseif strcmp(options.location,'Harvard')
    load('Harvard_responsesettings.mat');
elseif strcmp(options.location,'DeskLaptop')
    load('DeskLaptop_responsesettings.mat');
end


stopkey = KbName('q');
trigger = KbName('=+');
trigger2 = KbName('='); %Backup trigger
buttons = [ans1 ans2 ans3 ans4]; %1 and 2 for targ, 3 for cue


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          %File locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fiximgroot = './Stimuli/';
fiximgfilenames = 'fixaw.png';
AUD_FILE_ROOT = './Stimuli/stim300/'; %Aud files should be char.wav


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       % Paradigms and Timing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0 = FIX
% 1 = AL, 2 = AR, 3 = VL, 4= VR, 5=passive, 6=fix

% PARADIGMS
PAR_ALL_OPTIONS = [1 5 4 3 6 2 3 1 6 4 5 2 ;
                   2 3 6 1 5 4 3 5 2 4 1 6 ;
                   3 2 5 1 4 6 4 1 6 5 2 3 ;
                   4 1 3 6 2 5 2 6 3 5 1 4 ;
                   5 4 2 6 3 1 6 3 2 1 4 5];
               
                            
if PAR_NUM > size(PAR_ALL_OPTIONS,1)
    error('Invalid Paradigm Number, try a smaller value');
else  
    if PAR_NUM == 0
        COND_ALL = [1 4]; %Run two blocks for testing
    else
        COND_ALL = PAR_ALL_OPTIONS(PAR_NUM,:); %All conditions
    end 
end

NUM_BLOCKS = length(COND_ALL);
if NUM_TP ~= (2*FIX_TP + NUM_BLOCKS*(CUE_TP + BLOCK_TP))
    input('NUM_TP does not match paradigm. Continue?');
end


NUM_INPUTS = 4; %Two aud Two vid
STREAM_AL =1; STREAM_AR=2; STREAM_VL=3; STREAM_VR=4;

NUM_TRIALS_IN_TR = 4; 
STIM_PRES_TIME = .3;  %This should equal the length of the audio file / FS
VIS_WAIT = (TR/NUM_TRIALS_IN_TR)-STIM_PRES_TIME;
PRESTIM_TIME = .25; %VIS_WAIT/2;

TRI_PER_BLOCK = NUM_TRIALS_IN_TR * BLOCK_TP;
NUM_TRIALS = TRI_PER_BLOCK * NUM_BLOCKS;

ALLOWED_RESP_DELAY = 1.5; %s response buffer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Stimulus Type Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define target, cue, and distractor stimuli
%All auditory files should be named same as stimuli .wav

CUE_CHAR = {'Listen Left', 'Listen Right', 'Watch Left', 'Watch Right', 'Passive', 'Fixation'};
CUE_CHAR_NS = {'ListenLeft', 'ListenRight', 'WatchLeft', 'WatchRight', 'Passive', 'Fixation'};

load(strcat(fiximgroot','AUD_CUES_blockmls.mat')); %File is AUD_CUES(num cues, length of signal 40001, 2)

TARG_CHAR =     {'1','2','3','4'}; %This never comes up in first paradigm
TARG_CORR_ANS   = [ 1 2 3 4]; %Correct answers corresponding to TARG_CHAR
DIST_CHAR = {'1','2','3','4', '5', '6', '8', '9'};
TARG_STREAM_DIST_CHAR = {'A','F','G','H','J','K','L','M','N','P','R','X','Y'}; %Only needs to be as long as DIST_CHAR
%DIST_CHAR = {'A','C','F','G','H','J','K','L','M','N','P','R','T','V','X','Y'};
%EXT_DIST_CHAR = {'&', '$', '@', '{'};
EXT_DIST_CHAR = {'1','2','3','4', '5', '6', '8', '9', '1','2','3','4', '5', '6', '8', '9'}; % Doubled to see if code below will work, but this made repeats in aud
%EXT_DIST_CHAR = {'1','2','3','4', '5', '6', '8', '9', '0'}; % Added 0 to see if enough
%EXT_DIST_CHAR = {'A','C','F','G','H','J','K','L','M','N','P','R','T','V','X','Y'};
%EXT_DIST_CHAR = DIST_CHAR;

STIM_CHAR = [TARG_CHAR DIST_CHAR TARG_STREAM_DIST_CHAR]; %List of all characters in set order
STIM_CORR_ANS = TARG_CORR_ANS;


if extravis == 1  %If visual distractor streams 
    NUM_EXTRA_VIS = 6;
    EXTRA_VIS_IND = zeros(NUM_BLOCKS, NUM_EXTRA_VIS, TRI_PER_BLOCK);   
    for b = 1:NUM_BLOCKS
        [EXTRA_VIS_IND(b,:,:)] = createDistListInd(NUM_EXTRA_VIS, TRI_PER_BLOCK, length(EXT_DIST_CHAR));  %SWM: might need to change this, could be impossible
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     % Target Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep track of if there is a target presented at all


TARG_PER_BLOCK = 3;  %SWM: may want to fix how we determine how many targets


%Create list of non-repeating, non-matching distractors
%   Some will be unused, but assures no repeats without recalculating
DIST_LIST_IND = zeros(NUM_BLOCKS, NUM_INPUTS, TRI_PER_BLOCK); %blocks, stream, tri per block
 for b = 1:NUM_BLOCKS
    [DIST_LIST_IND(b,:,:)] = createDistListInd(NUM_INPUTS,TRI_PER_BLOCK, length(DIST_CHAR));
    % Write over the attended stream with new random values that will come
    % from TARG_STREAM_DIST_CHAR if not passive
    if COND_ALL(b) < 5
        [DIST_LIST_IND(b,COND_ALL(b),:)] = createDistListInd(1,TRI_PER_BLOCK, length(TARG_STREAM_DIST_CHAR));
        DIST_LIST_IND(b,COND_ALL(b),:) = DIST_LIST_IND(b,COND_ALL(b),:) + length(DIST_CHAR); %Offset by DIST_CHAR because of order of STIM_CHAR
    end
 end
 
% Create basis for each stimulus presentation using STIM_CHAR (targ then
% distractors)
TRI_CHAR_IND = DIST_LIST_IND + length(TARG_CHAR); %Shift values to match STIM_CHAR

% Generate the Target index List: index in TARG_CHAR for # of targs
% presented
TARG_LIST = randi(length(TARG_CHAR),[NUM_BLOCKS,TARG_PER_BLOCK]); 


%Determine where targets will fall for each block (should not be first
%trial or last three trials 
TARG_FALT_IND = zeros(NUM_BLOCKS, TRI_PER_BLOCK); %blocks, stream, tri per block
for b = 1:NUM_BLOCKS
    if COND_ALL(b) < 5
    %Determine targets
        rand_base = 2:4:(TRI_PER_BLOCK-3); %This means that the targets will always be a multiple of 4 apart
        rand_tp = rand_base(randperm(length(rand_base)));   
        rand_t = rand_tp(1:TARG_PER_BLOCK);
        rand_off = randi(TARG_PER_BLOCK); % random one in the list of offset by 1
        rand_t(rand_off) = rand_t(rand_off) + 1; %Selecting one to move off by 1, so not all multiples of 4 
        TARG_FALT_IND(b,rand_t) = 1; %Targets are 1
        
        % Set the trial character index to include tarets and false targets
        TRI_CHAR_IND(b,COND_ALL(b),rand_t) = TARG_LIST(b,:);
        
        % Make sure that the targets don't have distractor digits that are
        % the same
        for r = 1:length(rand_t)
            temp_target_value = TRI_CHAR_IND(b,COND_ALL(b),rand_t(r));
            %Detect if one of the targets and distractor digits is the same
            targdistdigsame = find(TRI_CHAR_IND(b,:,rand_t(r)) == (temp_target_value + length(TARG_CHAR)));
            if length(targdistdigsame) == 1 %If one
                %If so, replace with another distractor that does not match
                %any of the distractors in this list
                newdistval = 1:length(DIST_CHAR); 
                newdistval(newdistval==temp_target_value) = []; %Remove the target (this relies on the target values overlapping with the distractor digits
                newdistval = newdistval + length(TARG_CHAR); % increase this by target char value so stim_char is a distractor
                otherdistval = TRI_CHAR_IND(b,~COND_ALL(b),rand_t(r)); %list of other distractor values (these are + length(TARG_CHAR)
                for o = 1:length(otherdistval)
                    newdistval(newdistval==otherdistval(o)) = [];  %Remove the other distractor values
                end
                TRI_CHAR_IND(b,targdistdigsame,rand_t(r)) = newdistval(randi(length(newdistval))); %Replace this value with a random new distractor
            elseif  length(targdistdigsame) > 1
                error('Revisit this code, problem with distractors');
            end
        end
        
 
% DON'T do anything if passive, just leave these all as number distractors
%        
     end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                           %Auditory Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrchannels =2;
FS = 44100;
freq = FS; %48000;  %SWM: WHY IS THIS DIFFERENT THAN FS??
AUD_STIM_LEN = STIM_PRES_TIME * FS;  %Number of points in the aud stim file (may be issue 
%AUD_HEMI_DELAY_LEN = round(.000660 * FS); %Number of points to offset for a 500 us delay  
%SWM: CHANGED THIS TO SEE IF WE COULD GET BETTER STREAM SEPARATION
AUD_HEMI_DELAY_LEN = round(.000800 * FS); %Number of points to offset for a 500 us delay  

%Load auditory stimuli; Output format: [:,2, x], where the 2nd is
%delayed by AUD_HEMI_DELAY (so stereo will be left located) and x is the
%number of TARG, CUE, or DIST
[AUD_STIM] = loadAUD_STIM(AUD_FILE_ROOT, AUD_STIM_LEN, AUD_HEMI_DELAY_LEN, STIM_CHAR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                           %Visual Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digit stimuli settings
if overlap ==1
    hdist = 10; %dist from stim to center plus other adjustments
    fontsize = 100;
    digwidth = 56/2; %SWM - Need proper way to determine this
    voffset = 50;
else
    hdist = 150; %dist from stim to center plus other adjustments
    fontsize = 50;
    digwidth = 35; %50; %SWM - Need proper way to determine this
    voffset = 30;
end
EV_Toffset = 50; %Top offset
EV_Boffset = 50; %Bottom offset
EV_Foffset = hdist - 50; %Far offset   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         % Setup Responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

responses = zeros(2,NUM_TRIALS); %Store response button and time of response
stimprestimes = zeros(NUM_BLOCKS, TRI_PER_BLOCK); %Store time of each stimulus presentatio
stimprestimesGet = zeros(NUM_BLOCKS, TRI_PER_BLOCK); %Store time of each stimulus presentatio

for c = 1:length(COND_ALL)
    COND_ALL_NAMES{c} = num2str(COND_ALL(c));
end

R = struct('Condition',COND_ALL_NAMES,'StimTime',[],'RespTime',[],'ReactTime', [], 'StimChar',[], 'CorrResp',[],...
    'SubjResp',[], 'Correct',[],'Incorrect',[],'Miss',[], 'FalsePos', [], 'Valid', []);
 

try 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % Set Up Screens and Colors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    screens=Screen('Screens');
    fprintf('trying from the start')
    %screenNumber = 1;
    screenNumber=max(screens);   
	oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    
    % Color settings

    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    %gray=GrayIndex(screenNumber);
    
    cuecolor = [0 200 255];
    bluecolor = white - 50;  % Make stimuli white/gray instead of blue
    
    
    %ListenChar(2); %Suppress keyboard 
 
    HideCursor;
    [w, wRect]=Screen('OpenWindow',screenNumber, black);
    % Set up alpha-blending for smooth (anti-aliased) drawing of dots:
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    priorityLevel=MaxPriority(w);
    Priority(priorityLevel);
           
    cx = wRect(3)/2;
    cy = wRect(4)/2;
    Screen('TextFont', w, 'Helvetica')
    Screen('TextStyle', w, 1); %bold
    Screen('TextSize', w, fontsize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             %Load fixation cue   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fiximgfile= [fiximgroot fiximgfilenames];
    [fiximgdata fiximgmap fiximgalpha] = imread(fiximgfile);
    fiximgdata(:,:,4) = fiximgalpha(:,:) - 50;        
    fixtex = Screen('MakeTexture',w, fiximgdata);
    clear fiximgdata fiximgfile fiximgmap fiximgalpha

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Eye tracker settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.eyetracker == true 
        %set up pixesl per degree
        pix_per_deg = pi * wRect(3) / atan(p.mon_width/p.v_dist/2) / 360;	% pixels per degree
    
        %convert stim size to pixels
        box=2*pix_per_deg; % size of 2 degree box in pixels
       
        EYE.WINDOW = 60;
        EYE.OUT_WINDOW = 100;
        EYE.DISTANCE  = 200;
        xo=round(wRect(3)/2);      
        yo=round(wRect(4)/2);
        %Rects
        fixcenterx1 = xo-EYE.DISTANCE /2;
        fixcentery1 = yo;
        fixcenterx2 = xo+EYE.DISTANCE /2;
        fixcentery2 = yo;
       
        el=EyelinkInitDefaults(w); 

        if ~EyelinkInit(0)
            fprintf('Eyelink Init aborted.\n');
            % Shutdown Eyelink:
            Eyelink('Shutdown');
            return;
        end

        connected = Eyelink('IsConnected')%#ok<NASGU> %%Just to verify connection
        [vEyelink vsEyelink]=Eyelink('GetTrackerVersion');
        fprintf('Running experiment on a ''%s'' tracker.\n', vsEyelink );

        % open file to record data to
        tempeye = Eyelink('Openfile', edfFile); %%doesn't return
        if tempeye~=0
            error('Cannot create EDF file ''%s'' ', edfFile);
            Eyelink( 'Shutdown');
            return;
        end
        
        Eyelink('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox tempASTM experiment''');

        %   SET UP TRACKER CONFIGURATION
        if options.eyecalibration
            Eyelink('command', 'calibration_type = HV5'); %Using HV5 instead of 9, since just fixation
        end
        %	set parser (conservative saccade thresholds)
        Eyelink('command', 'recording_parse_type = GAZE');
        Eyelink('command', 'saccade_acceleration_threshold = 8000');
        Eyelink('command', 'saccade_velocity_threshold = 30');
        Eyelink('command', 'saccade_motion_threshold = 0.0');
        Eyelink('command', 'saccade_pursuit_fixup = 60');
        Eyelink('command', 'fixation_update_interval = 0');

        %	set EDF file contents
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,GAZERES,HREF,AREA');
        Eyelink('command', 'file_event_data  = GAZE,GAZERES,AREA,VELOCITY,HREF');

        %	set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,HREF,AREA');
        Eyelink('command', 'link_event_data  = GAZE,GAZERES,AREA,VELOCITY,HREF');
        
        
        % make sure we're still connected.
        if Eyelink('IsConnected')~=1
            Eyelink( 'Shutdown');
            display('Shutting down Eyelink');
            return;
        end;
        
        % Calibrate the eye tracker
        if options.eyecalibration
            EyelinkDoTrackerSetup(el);
        end
        eye_used = Eyelink('EyeAvailable');
               
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         %Do drift correction again ? then wait??? 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.driftcorrect
            black = BlackIndex(screenNumber);
            Screen('FillRect', w, black); 
            Screen('DrawTexture',w,fixtex);
            Screen('Flip',w);
            EyelinkDoDriftCorrect(el, xo, yo, 0, 1);  %SWM EYE: May need to present fixation cross here
            input('Done drift correction. Press button to continue \n')
        end

        % setup the proper calibration foreground and background colors
        el.backgroundcolour = black;
        el.foregroundcolour = bluecolor;
        
        % Write some info to the eyetracker about conditions 
        Eyelink('Message', 'TRIAL_VAR_LABELS ListenLeft ListenRight WatchLeft WatchRight Passive Fixation');
        Eyelink('Message', 'V_TRIAL_GROUPING ListenLeft ListenRight WatchLeft WatchRight Passive Fixation');
        
        % Before recording, we place reference graphics on the host display
        % Must be offline to draw to EyeLink screen
        Eyelink('Command', 'set_idle_mode'); %toggle to off
        % clear tracker display and draw box at fix point
        Eyelink('Command', 'clear_screen 0');
        Eyelink('command', 'draw_box %d %d %d %d 15', xo-box, yo-box, xo+box, yo+box); %box in pixels

        Eyelink('Command', 'set_idle_mode'); %toggle back on
        WaitSecs(0.01);   
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          %Present fixation cross  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    black = BlackIndex(screenNumber);
    Screen('FillRect', w, black); 
    Screen('DrawTexture',w,fixtex);
    Screen('Flip',w);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %check for forced stop or trigger  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while 1
    [keyIsDown, secs, keyCode, deltaSecs] =KbCheck(-1);
        if keyCode(stopkey)>0
            Screen('CloseAll');
            PsychPortAudio('Close');
            disp('stop key pressed');
            break;
        elseif keyCode(trigger)>0 || keyCode(trigger2)>0
            disp('trigger pressed');
            break;
        end
    end
    
    startSecs = tic;  
    startSecsGet = GetSecs;  
    
    % Initilize data collection
    KbQueueCreate(kbnumber);
    KbQueueStart;
    
    % Perform basic initialization of the sound driver:
    InitializePsychSound;
    pahandle = PsychPortAudio('Open', [], [], 1, freq, nrchannels);  %needed?    
    
    
    Screen('DrawTexture',w,fixtex);
    Screen('Flip',w);
       
    
    %%%%%%%%%%%%%%%%%%
    %Initialize stim variables
    last_STIM_time = GetSecs - startSecsGet; %Time relative to start of last relevant stimulus presentation
    corr_STIM_button = 0; %Index in buttons for correct resp

    
    while (toc(startSecs) < FIXATION_TIME); end;
        
        
    %Main loop
    for curBlock = 1:NUM_BLOCKS
        blockSecs = tic;
        active_STIM = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    Turn on eyetracker
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.eyetracker == true
            % Send trial id message to eyelink file
            % Sending a 'TRIALID' message to mark the start of a trial in Data 
            % Viewer.  This is different than the start of recording message 
            % START that is logged when the trial recording begins. The viewer
            % will not parse any messages, events, or samples, that exist in 
            % the data file prior to this message. 
            Eyelink('Message', 'TRIALID %d', COND_ALL(curBlock));

            % This supplies the title at the bottom of the eyetracker display
            Eyelink('command', 'record_status_message "TRIAL %d"', COND_ALL(curBlock)); 
            
            
            if options.driftcorrect  %SWM EYE: not sure if this will work
                EyelinkDoDriftCorrect(el, xo, yo, 0, 1);  
 
            end

            %%%actually start recording
            Eyelink('StartRecording');
            
            Screen('DrawTexture',w,fixtex); %Not sure if you need to flip the screen again for timing, but doing it here just in case
            Screen('Flip',w);
            Eyelink('Message', 'SYNCTIME');
        end
        
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Cue block type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        %Present Auditory and Visual cue
        audStimC = squeeze(AUD_CUES(COND_ALL(curBlock),:,:));   
        pahandleC = PsychPortAudio('Open', [], [], 1, freq, size(audStimC',1));         
        PsychPortAudio('FillBuffer', pahandleC, audStimC');
        %Draw visual cue
        Screen('DrawTexture',w,fixtex); 

        
        while (toc(blockSecs) < 0.1); end; %Wait very briefly to make sure buffer gets filled
        
        DrawFormattedText(w,CUE_CHAR{COND_ALL(curBlock)},'center','center',cuecolor); 
        t3 = PsychPortAudio('Start', pahandleC, 1, 0, 1);
        Screen('Flip',w); %Display
        
        Screen('DrawTexture',w,fixtex); 
        while (toc(blockSecs) < (CUE_TIME-.600)); end
        Screen('Flip',w); %Display
        PsychPortAudio('Close',pahandleC);
        
        KbQueueFlush(); %Reset last button across blocks
        
        
        % Load first stimulus
        if COND_ALL(curBlock) ~= 6 % If not fixation next
            % Load auditory stimuli
            audLeft = AUD_STIM(:,:,TRI_CHAR_IND(curBlock, STREAM_AL,1));
            audRight = [AUD_STIM(:,2,TRI_CHAR_IND(curBlock, STREAM_AR,1))  AUD_STIM(:,1,TRI_CHAR_IND(curBlock, STREAM_AR,1))];
            audStimB = audLeft + audRight;  %SWM: Should I be normalizing?  
            pahandleB = PsychPortAudio('Open', [], [], 1, freq, size(audStimB',1));         
            PsychPortAudio('FillBuffer', pahandleB, audStimB');
            % Load visual stimuli
            Screen('DrawTexture',w,fixtex); %select cue for this trial per paradigm
            DrawFormattedText(w,STIM_CHAR{TRI_CHAR_IND(curBlock, STREAM_VL,1)},cx-hdist-digwidth, 'center',bluecolor); %Draw VL
            DrawFormattedText(w,STIM_CHAR{TRI_CHAR_IND(curBlock, STREAM_VR,1)},cx+hdist,         'center',bluecolor); %Draw VR 

            if extravis == 1  %If visual distractor streams 
                DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 1,1)},cx-EV_Foffset-digwidth, 'center',bluecolor); %Draw L Far
                DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 2,1)},cx+ EV_Foffset,   'center' ,bluecolor); %Draw R Far 
                DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 3,1)},cx-hdist-digwidth, cy-voffset-EV_Toffset,bluecolor); %Draw L Top
                DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 4,1)},cx+hdist,          cy-voffset-EV_Toffset,bluecolor); %Draw R Top
                DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 5,1)},cx-hdist-digwidth, cy-voffset+EV_Boffset,bluecolor); %Draw L Bottom
                DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 6,1)},cx+hdist,          cy-voffset+EV_Boffset,bluecolor); %Draw R Bottom 
            end
              
        end
        
        while (toc(blockSecs) < CUE_TIME); end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Task block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If FIXATION
        if COND_ALL(curBlock) == 6 
            while (toc(blockSecs) < (CUE_TIME + BLOCK_TIME)); end
        else
            
        % If STIMULI BLOCKS
            for curIndex = 1:TRI_PER_BLOCK
                trialSecs = tic; 

                % Wait for prestim buffer time and present stimuli
                while (toc(trialSecs) < PRESTIM_TIME); end;          
                t3 = PsychPortAudio('Start', pahandleB, 1, 0, 1);
                Screen('Flip',w); %Display 
                stimprestimes(curBlock,curIndex) = toc(startSecs);
                stimprestimesGet(curBlock,curIndex) = GetSecs-startSecsGet;
                
                %If a relevant stimulus is presented
                if TARG_FALT_IND(curBlock, curIndex) == 1
                    active_STIM = 1; %
                    last_STIM_time = stimprestimesGet(curBlock,curIndex); %toc(startSecs); %Time relative to start of last relevant stimulus presentation
                    last_CHAR_IND = TRI_CHAR_IND(curBlock,COND_ALL(curBlock), curIndex); %Block, stream, index (only works b/c streams match par block num
                    corr_STIM_button = STIM_CORR_ANS(last_CHAR_IND); %Index in buttons for correct resp
                end

                % Turn off stimulus after appropriate amount of time
                while (toc(trialSecs) < PRESTIM_TIME + STIM_PRES_TIME); end; 
                Screen('DrawTexture',w,fixtex); 
                Screen('Flip',w);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             %Collect responses and record data    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                [responses lastbutton lbtime dp] = collectResponses_check(responses, stopkey, buttons,startSecsGet);

                %If correct resonse in proper time window
                
                if lastbutton > 0 && active_STIM == 1 && (lbtime - last_STIM_time) > 0 %Added to make sure no neg reaction times
                    %Check to see if they match and update the appropriate
                    %condition
                    R(1,curBlock).StimTime(end+1)  = last_STIM_time;
                    R(1,curBlock).RespTime(end+1)  = lbtime;
                    R(1,curBlock).ReactTime(end+1) = lbtime - last_STIM_time;
                    R(1,curBlock).StimChar(end+1)  = last_CHAR_IND;
                    R(1,curBlock).CorrResp(end+1)  = corr_STIM_button;
                    R(1,curBlock).SubjResp(end+1)  = lastbutton;
                    R(1,curBlock).Correct(end+1)   = (corr_STIM_button ==lastbutton);
                    R(1,curBlock).Incorrect(end+1) = (corr_STIM_button ~=lastbutton);
                    R(1,curBlock).Miss(end+1)      = 0;
                    R(1,curBlock).FalsePos(end+1)  = 0;
                    R(1,curBlock).Valid(end+1)      = active_STIM;
                    active_STIM = 0; %Reset active_STIM

                %If no response yet in proper time window therefore a miss
                elseif active_STIM == 1 && lastbutton == 0 && ((GetSecs - startSecsGet  - last_STIM_time) > ALLOWED_RESP_DELAY)
                    %Make a note of a miss
                    R(1,curBlock).StimTime(end+1)  = last_STIM_time;
                    R(1,curBlock).RespTime(end+1)  = 0; %Did not response
                    R(1,curBlock).ReactTime(end+1) = nan; 
                    R(1,curBlock).StimChar(end+1)  = last_CHAR_IND;
                    R(1,curBlock).CorrResp(end+1)  = corr_STIM_button;
                    R(1,curBlock).SubjResp(end+1)  = 0;
                    R(1,curBlock).Correct(end+1)   = 0;
                    R(1,curBlock).Incorrect(end+1) = 0;
                    R(1,curBlock).Miss(end+1)      = 1;
                    R(1,curBlock).FalsePos(end+1)  = 0;
                    R(1,curBlock).Valid(end+1)     = active_STIM;

                    active_STIM = 0; %Reset active_STIM
                    
                % If false positive    
                elseif lastbutton > 0 
                    %Make note of a false positive... you can do this in the
                    %current distractor condition, since we made room for these
                    %anway
                    R(1,curBlock).StimTime(end+1)  = last_STIM_time;
                    R(1,curBlock).RespTime(end+1)  = lbtime;
                    R(1,curBlock).ReactTime(end+1) = lbtime - last_STIM_time;  %How late from last actual stimulus... not useful but might record late misses
                    R(1,curBlock).StimChar(end+1)  = nan; %If there was no FALT
                    R(1,curBlock).CorrResp(end+1)  = nan;

                    R(1,curBlock).Valid(end+1)      = active_STIM;
                    R(1,curBlock).SubjResp(end+1)  = lastbutton;
                    R(1,curBlock).Correct(end+1)   = 0;
                    R(1,curBlock).Incorrect(end+1) = 0;
                    R(1,curBlock).Miss(end+1)      = 0;
                    R(1,curBlock).FalsePos(end+1)  = 1;

                end           

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    % End Trial
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                PsychPortAudio('Close',pahandleB);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           % Load next stimulus if not last trial
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if curIndex < TRI_PER_BLOCK 
                    % Load auditory stimuli
                    audLeft = AUD_STIM(:,:,TRI_CHAR_IND(curBlock, STREAM_AL,curIndex + 1));
                    audRight = [AUD_STIM(:,2,TRI_CHAR_IND(curBlock, STREAM_AR,curIndex + 1))  AUD_STIM(:,1,TRI_CHAR_IND(curBlock, STREAM_AR,curIndex + 1))];
                    audStimB = audLeft + audRight;  %SWM: Should I be normalizing?  
                    pahandleB = PsychPortAudio('Open', [], [], 1, freq, size(audStimB',1));         
                    PsychPortAudio('FillBuffer', pahandleB, audStimB');
                    % Load visual stimuli
                    Screen('DrawTexture',w,fixtex); %select cue for this trial per paradigm
                    DrawFormattedText(w,STIM_CHAR{TRI_CHAR_IND(curBlock, STREAM_VL,curIndex + 1)},cx-hdist-digwidth, 'center',bluecolor); %Draw VL
                    DrawFormattedText(w,STIM_CHAR{TRI_CHAR_IND(curBlock, STREAM_VR,curIndex + 1)},cx+hdist,         'center',bluecolor); %Draw VR 

                    if extravis == 1  %If visual distractor streams 
                        DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 1,curIndex + 1)},cx-EV_Foffset-digwidth, 'center',bluecolor); %Draw L Far
                        DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 2,curIndex + 1)},cx+ EV_Foffset,   'center' ,bluecolor); %Draw R Far 
                        DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 3,curIndex + 1)},cx-hdist-digwidth, cy-voffset-EV_Toffset,bluecolor); %Draw L Top
                        DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 4,curIndex + 1)},cx+hdist,          cy-voffset-EV_Toffset,bluecolor); %Draw R Top
                        DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 5,curIndex + 1)},cx-hdist-digwidth, cy-voffset+EV_Boffset,bluecolor); %Draw L Bottom
                        DrawFormattedText(w,EXT_DIST_CHAR{EXTRA_VIS_IND(curBlock, 6,curIndex + 1)},cx+hdist,          cy-voffset+EV_Boffset,bluecolor); %Draw R Bottom 
                    end
                end
                              

                while (toc(trialSecs) < VIS_WAIT + STIM_PRES_TIME); end; 

            end 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     % Stop and record eyetracking   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.eyetracker == true
            % adds 100 msec of data to catch final events
            %WaitSecs(0.1);
            % stop the recording of eye-movements for the current block
            Eyelink('StopRecording');
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, xo-box, yo-box, xo+box, yo+box, CUE_CHAR_NS{COND_ALL(curBlock)});

            Eyelink('Message', '!V TRIAL_VAR index %d', curBlock); %Eventually replace with index for each trial and move up into trial section
            Eyelink('Message', '!V TRIAL_VAR type %d', COND_ALL(curBlock));

            %%give it info about what group belongs to
            Eyelink('Message', '!V TRIAL_VAR_DATA %s', CUE_CHAR_NS{COND_ALL(curBlock)});

            % Sending a 'TRIAL_RESULT' message to mark the end of a trial in 
            % Data Viewer. This is different than the end of recording message 
            % END that is logged when the trial recording ends. The viewer will
            % not parse any messages, events, or samples that exist in the data 
            % file after this message.
            Eyelink('Message', 'TRIAL_RESULT 0')
        end
        if (toc(blockSecs) >= (CUE_TIME + BLOCK_TIME))
            blocksovertime = blocksovertime + 1;
        else
           while (toc(blockSecs) < (CUE_TIME + BLOCK_TIME)); end 
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              % Final fixation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Screen('DrawTexture',w,fixtex);
    Screen('Flip',w);
    Totaltime = toc(startSecs) + FIXATION_TIME
    
    WaitSecs(FIXATION_TIME);
    [responses lastbutton lbtime dp] = collectResponses_check( responses, stopkey, buttons,startSecsGet);
    
    
    if lastbutton > 0 && active_STIM == 1 && (lbtime - last_STIM_time) > 0 %Added to make sure no neg reaction times
        %Check to see if they match and update the appropriate
        %condition
        R(1,curBlock).StimTime(end+1)  = last_STIM_time;
        R(1,curBlock).RespTime(end+1)  = lbtime;
        R(1,curBlock).ReactTime(end+1) = lbtime - last_STIM_time;
        R(1,curBlock).StimChar(end+1)  = last_CHAR_IND;
        R(1,curBlock).CorrResp(end+1)  = corr_STIM_button;
        R(1,curBlock).SubjResp(end+1)  = lastbutton;
        R(1,curBlock).Correct(end+1)   = (corr_STIM_button ==lastbutton);
        R(1,curBlock).Incorrect(end+1) = (corr_STIM_button ~=lastbutton);
        R(1,curBlock).Miss(end+1)      = 0;
        R(1,curBlock).FalsePos(end+1)  = 0;
        R(1,curBlock).Valid(end+1)      = active_STIM;

    %If no response yet in proper time window therefore a miss
    elseif active_STIM == 1 && lastbutton == 0 && ((GetSecs - startSecsGet  - last_STIM_time) > ALLOWED_RESP_DELAY)
        %Make a note of a miss
        R(1,curBlock).StimTime(end+1)  = last_STIM_time;
        R(1,curBlock).RespTime(end+1)  = 0; %Did not response
        R(1,curBlock).ReactTime(end+1) = nan; 
        R(1,curBlock).StimChar(end+1)  = last_CHAR_IND;
        R(1,curBlock).CorrResp(end+1)  = corr_STIM_button;
        R(1,curBlock).SubjResp(end+1)  = 0;
        R(1,curBlock).Correct(end+1)   = 0;
        R(1,curBlock).Incorrect(end+1) = 0;
        R(1,curBlock).Miss(end+1)      = 1;
        R(1,curBlock).FalsePos(end+1)  = 0;
        R(1,curBlock).Valid(end+1)     = active_STIM;

    % If false positive    
    elseif lastbutton > 0 
        %Make note of a false positive... you can do this in the
        %current distractor condition, since we made room for these
        %anway
        R(1,curBlock).StimTime(end+1)  = last_STIM_time;
        R(1,curBlock).RespTime(end+1)  = lbtime;
        R(1,curBlock).ReactTime(end+1) = lbtime - last_STIM_time;  %How late from last actual stimulus... not useful but might record late misses
        R(1,curBlock).StimChar(end+1)  = nan; %If there was no FALT
        R(1,curBlock).CorrResp(end+1)  = nan;

        R(1,curBlock).Valid(end+1)      = active_STIM;
        R(1,curBlock).SubjResp(end+1)  = lastbutton;
        R(1,curBlock).Correct(end+1)   = 0;
        R(1,curBlock).Incorrect(end+1) = 0;
        R(1,curBlock).Miss(end+1)      = 0;
        R(1,curBlock).FalsePos(end+1)  = 1;
    end            
                
    Screen('CloseAll');
    PsychPortAudio('Close');
    
    
    KbQueueRelease;
    % Revive the mouse cursor.
    ShowCursor;
    ListenChar(0); %Unsuppress keyboard
       
    %Analyze the data
    Correct = zeros(1,max(COND_ALL));
    Incorrect = zeros(1,max(COND_ALL));
    Miss = zeros(1,max(COND_ALL));
    FalsePos = zeros(1,max(COND_ALL));

    MeanBlockRT = zeros(1,max(COND_ALL));
    

%     %By Condition
    for i = 1:(NUM_BLOCKS)

        Correct(str2num(R(1,i).Condition))       = sum(R(1,i).Correct) + Correct(str2num(R(1,i).Condition));
        Incorrect(str2num(R(1,i).Condition))     = sum(R(1,i).Incorrect) + Incorrect(str2num(R(1,i).Condition));
        Miss(str2num(R(1,i).Condition))          = sum(R(1,i).Miss) + Miss(str2num(R(1,i).Condition));
        FalsePos(str2num(R(1,i).Condition))      = sum(R(1,i).FalsePos) + FalsePos(str2num(R(1,i).Condition));
        MeanBlockRT(str2num(R(1,i).Condition))   = sum(R(1,i).ReactTime(~R(1,i).Miss)) + MeanBlockRT(str2num(R(1,i).Condition));   
    end
    NUM_DIV = TARG_PER_BLOCK * 2;
%     NUM_DIV_FALT = FALT_PER_BLOCK * 2;
     Correct./NUM_DIV
     Incorrect./NUM_DIV
     Miss./NUM_DIV
     FalsePos

    
    save(strcat('./Results/MODLOC_OUTPUT_',subjID,'_',datestr(now)));
       
        
    % Restore preferences
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     Shut down eyetracker
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.eyetracker == true 
        Eyelink('Command', 'set_idle_mode'); % Set back to idle mode to close trial
        WaitSecs(0.5);
        Eyelink('CloseFile');
        % download data file
        try
            fprintf('Receiving data file ''%s''\n', edfFile );
            status=Eyelink('ReceiveFile'); %SWM: fix this to write to proper place
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            end
        catch ETE
            fprintf('Problem receiving data file ''%s''\n', edfFile );
            disp(ETE.message)
        end

        %%%%%%%%%%%%%%%%shut it down
        Eyelink('ShutDown');
    end
    
    
    
    
    
catch SE
    ListenChar(0); %Unsuppress keyboard
    Screen('CloseAll');
    PsychPortAudio('Close');
    disp('caught some errors: See SE');
    disp(SE.message);
    SE.stack.line
    %psychrethrow(psychlasterror);
    KbQueueRelease;
    save(strcat('./Results/MODLOC_OUTPUT_',subjID,'ERROR_',datestr(now)));
    %save(strcat('OUTPUT_',subjID,num2str(runNum)),'responses', 'paradigm', 'tripar', 'triparoffset','triresults','LOGstim','LOGcurStream','LOGfutStream', 'LOG ');
    if options.eyetracker == true
        Eyelink('Command', 'set_idle_mode'); % Set back to idle mode to close trial
        WaitSecs(0.5);
        Eyelink('CloseFile');
        % download data file
        try
            fprintf('Receiving data file ''%s''\n', edfFile );
            status=Eyelink('ReceiveFile'); %SWM: fix this to write to proper place
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            end
        catch ETE
            fprintf('Problem receiving data file ''%s''\n', edfFile );
            disp(ETE.message)
        end

        %%%%%%%%%%%%%%%%shut it down
        Eyelink('ShutDown');
    end

end


end


function [AUD_STIM] = loadAUD_STIM(AUD_FILE_ROOT, AUD_STIM_LEN, AUD_HEMI_DELAY_LEN, STIM_CHAR)
%loadAUD_STIM:  Load the auditory stimuli
%   Output format: [:,2], where the 2nd is delayed by AUD_HEMI_DELAY (so stereo will be left located)
%   Note: The delayed version will be cropped by AUD_HEMI_DELAY/FS, make
%   sure there is no clipping because of this.  If so, will need to extend
%   the AUD_STIM_LEN
%   AUD_STIM is indexed in the third dimension by the STIM_CHAR index

%Initialize outputs
AUD_STIM = zeros(AUD_STIM_LEN,2,length(STIM_CHAR));

try
    for i = 1:length(STIM_CHAR)
        audStimFile = strcat(AUD_FILE_ROOT, STIM_CHAR{i}, '.wav');
        AUD_STIM(:,1,i) = audioread(audStimFile);
        AUD_STIM((AUD_HEMI_DELAY_LEN+1):end,2,i)= AUD_STIM(1:(end-AUD_HEMI_DELAY_LEN),1,i);
    end
catch ERROR_loadAudStim
    disp('ERROR: See ERROR_loadAudStim below: \n');
    disp(ERROR_loadAudStim.message)
end
    
end

function [ distListInd ] = createDistListInd(numInputs,numTrials,lenDistStim)
%createDistList: Create list of non-repeating, non-matching distractors
%   Some will be unused, but assures no repeats without recalculating

distListInd = zeros(numInputs,numTrials);
for nT = 1:numTrials
    sV = 1:lenDistStim;
    for nI = 1:numInputs
        %Remove the matching value prior to prevent repeat
        if nT>1; toRem = distListInd(nI,nT-1); sV(sV==toRem)=[]; end;
        ri = randi(length(sV));
        distListInd(nI,nT)=sV(ri);
        sV(ri)=[];
    end
end

end

function [ responses lastbutton lbtime dp] = collectResponses_check( responses, stopkey, buttons,startSecsGet)
%collectResponses: Record responses and times
%   Note: lastbutton defaults to the higher index if two buttons are
%   pressed simultaneously.

    %Check for button press
    [pressed, firstPress, junk, lastPress]=KbQueueCheck();       
    if pressed == 1
        fpt = firstPress(buttons); 
        lpt = lastPress(buttons);
        lbtime = max(lpt); %Time of last button press
        lastbutton = find(lpt==lbtime, 1, 'last' );  %Index in buttons of last button press
        lp = find(lpt>fpt);  %Find the last press indices in buttons for repeats only
        if ~isempty(lp) 
            lpsecs = lpt(lp)-startSecsGet; 
        else
            lpsecs = [];
        end
        
        ind1 = find(responses(1,:),1,'last')+1; %Index for start of response recording
        if isempty(ind1); ind1=1;end;
        dp = length(find(fpt));  %Find the button indicies for pressed
        indf = ind1 + length(find(fpt)) + length(lp) -1;
    

        responses(1,ind1:indf) = [find(fpt) lp];         % Button indicies        
        responses(2,ind1:indf) = [(fpt(fpt>0)-startSecsGet) lpsecs];   % Press times

    else
        lastbutton = 0;
        lbtime = 0;
        dp = 0;
    end
    lbtime = lbtime - startSecsGet;
    KbQueueFlush();
    
    [keyIsDown, secs, keyCode, deltaSecs] =KbCheck(-3);
    if keyCode(stopkey)>0
            Screen('CloseAll');
            PsychPortAudio('Close');
            error('STOP key pressed');          
    end
    

end



