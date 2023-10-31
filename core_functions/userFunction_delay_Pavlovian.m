function [logger , logger_valsROIs] =...
    userFunction_BMIv11_withZ(source, event, varargin)
%% Variable stuff
tic
global logger loggerNames...
    pe shifter...
    params data rois sm...
    frameNum...
    baselineStuff trialStuff...
    
% persistent  

%% IMPORT DATA
% TODO: Create NI-DAQ output to control optogenetic trigger
% TODO: Create NI-DAQ input (?) to read teensy output
% TODO: Catch trial: uncued reward session and cued omission session

frameNum = source.hSI.hStackManager.framesDone;

if frameNum == 1
    params = struct();
    data = struct();
    rois = struct();
    sm = struct();
end

data.currentImage = int32(source.hSI.hDisplay.lastFrame{1});
% disp(class(data.currentImage))
data.currentImage_gpu = gpuArray(data.currentImage);
data.hash_image = simple_image_hash(data.currentImage);  %% slower on gpu
data.MC.current_position_z = source.hSI.hFastZ.currentFastZs{1}.targetPosition;

%% == USER SETTINGS ==
if frameNum == 1
    % SETTINGS: General
    params.paths.directory = 'D:\RH_local\data\cage_0914\mouse_0914\20231030\analysis_data';
%     params.paths.expSettings = 'D:\RH_local\data\cage_0415\mouse_0415R\20230816\analysis_data\day0_analysis\expSettings.mat';
    params.paths.expSettings = false;
    
    % SETTINGS: TIMING
    params.timing.frameRate          = 30;
    params.timing.duration_plotting  = 30 * params.timing.frameRate; % ADJUSTABLE: change number value (in seconds). Duration of x axis in plots
    params.timing.duration_session   = source.hSI.hStackManager.framesPerSlice;
    
    if params.timing.duration_session == inf % when you click FOCUS...
        params.timing.duration_session = params.timing.frameRate*60*60; % Just a random number that is long enough
    end
    
    disp(['duration_session :  ', num2str(params.timing.duration_session)])
    
    % 20230509 Log Frames / File check
    disp(['Log Frames per File :  ', num2str(source.hSI.hScan2D.logFramesPerFile)])
    if source.hSI.hScan2D.logFramesPerFile ~= 1000
        error('Frames per File is not 1000')
    end
    
    % SETTINGS: Motion correction
    params.MC.numFrames_avgWin_zCorr      = 30*2;
    params.MC.intervalFrames_zCorr        = 5;
    params.MC.min_interval_z_correction   = 20*params.timing.frameRate;
    params.MC.max_delta_z_correction      = 0.5;
    params.MC.bandpass_freqs              = [1/64, 1/4];
    params.MC.bandpass_orderButter        = 3;
    params.MC.device                      = 'cuda';
    params.MC.frame_shape_yx              = int64([512,512]);
    
    % SETTINGS: Cursor: Only works for params.trial.block_trial == false
    params.cursor.factor_to_use = 2;
    params.cursor.angle_power = 2.0;
 
    params.cursor.threshold_reward     = 1.5;
    params.cursor.thresh_quiescence_cursorDecoder = 0.15;
    params.cursor.thresh_quiescence_cursorMag = 0;
    params.cursor.win_smooth_cursor    = 1; % smoothing window (in frames)
    params.cursor.bounds_cursor        = [-params.cursor.threshold_reward , params.cursor.threshold_reward *1.5];
    params.cursor.range_freqOutput     = [1000 18000]; % this is set in the teensy code (only here for logging purposes)
    params.cursor.voltage_at_threshold = 3.1; % this will be the maximum output voltage ([0:voltage_at_threshold])
    
    % SETTINGS: Mode
    params.mode = 'sm_only';
%     params.mode = 'imaging';
    
    % SETTINGS: Trials
    % below in unit seconds
    params.trial.reward_duration = 52; % 01/14/2023 in ms calibrated to 2.5 uL/reward
    params.trial.reward_delay = 0; % in ms
    % params.trial.LED_duration = 0.2; % in s
    % params.trial.LED_ramp_duration = 0.1; % in s
    params.trial.duration_trial          = 20;
    params.trial.duration_timeout        = 4;
    
%     params.trial.duration_threshold      = 0.066; % in seconds
    params.trial.duration_threshold      = 3; % in frames
    
    params.trial.duration_ITI            = 3; % in seconds
    params.trial.duration_rewardDelivery = 1.00; % before it was .2 10/10/22: Controls how long the reward tone is
    
%     params.trial.duration_quiescenceHold = 0.5;  % in seconds
    params.trial.duration_quiescenceHold      = 3; % in frames
    
    params.trial.duration_buildingUpStats    = round(params.timing.frameRate * 60 * 1);
    
    if isscalar(params.paths.expSettings)==false
        if isfile(params.paths.expSettings)
            clear expSettings
%             load(params.paths.expSettings);
            load(params.paths.expSettings, 'expSettings'); % To avoid adding variable to a static workspace
%             disp(num2str(exist('expSettings') > 0))
            assert(exist('expSettings') > 0, 'RH ERROR: Imported file from params.paths.expSettings, but variable name is not expSettings');
            
            % overwrite all parameters contained within expSettings
            disp('Inherit expSettings parameters')
            params = overwrite_struct_fields(params, expSettings);
        else
            error('RH ERROR: params.paths.expSettings is not false, but does not point to a valid file')
        end
    end
end

%% INITIALIZE EXPERIMENT
if frameNum == 1 
    % 20230724 Experiment Length Sanity Check
    if source.hSI.hStackManager.framesPerSlice ~= inf % Does not run this part when you hit FOCUS
        if params.timing.duration_session ~= (baselineStuff.block_sequence.expected_duration_session * params.timing.frameRate)
            error(['Expected Session Length does not match: Scanimage ',...
            num2str(params.timing.duration_session),...
            ', baselineStuff ',...
            num2str(baselineStuff.block_sequence.expected_duration_session)])
        end
    end
    
    path_trialStuff = [params.paths.directory , '\trialStuff.mat'];
    load(path_trialStuff);
    disp(['LOADED trialStuff from:  ' , path_trialStuff])
    
    if strcmp(params.mode, 'imaging')
        % type_stack = 'stack_warped';
        type_stack = 'stack_sparse';
        zstack = load([params.paths.directory , '\', type_stack, '.mat']);
        
        %%% Motion correction python code prep
        try
            pe = pyenv('Version', 'C:\ProgramData\Miniconda\envs\matlab_env\python');  %% prepare python environment
        catch
            disp('failed to initalize Python environment. The environment may already by loaded')
        end
        py.importlib.import_module('bph.motion_correction');
        py.importlib.import_module('rp.rolling_percentile');
        
        s_y = floor((size(data.currentImage,1)-params.MC.frame_shape_yx(1))/2) + 1;
        s_x = floor((size(data.currentImage,2)-params.MC.frame_shape_yx(2))/2) + 1;
        data.MC.idx_im_MC_crop_y = s_y:s_y+params.MC.frame_shape_yx(1)-1;
        data.MC.idx_im_MC_crop_x = s_x:s_x+params.MC.frame_shape_yx(2)-1;
        
        data.im_zstack = eval(['zstack.', type_stack, '.stack_avg']);
        data.im_zstack = single(data.im_zstack(:, data.MC.idx_im_MC_crop_y, data.MC.idx_im_MC_crop_x));
        data.MC.stepSize_zstack = eval(['zstack.', type_stack, '.step_size_um']);
        data.MC.n_slices_zstack = size(data.im_zstack, 1);
        data.MC.idx_middle_frame = ceil(data.MC.n_slices_zstack/2);
    
        im = squeeze(data.im_zstack(data.MC.idx_middle_frame, :,:));
        data.MC.im_refIm_MC_2D = gpuArray(single(im));
        
        % Initialize the shifter class
        shifter = py.bph.motion_correction.Shifter_rigid(params.MC.device);
        shifter.make_mask(py.tuple(params.MC.frame_shape_yx), py.tuple(params.MC.bandpass_freqs), params.MC.bandpass_orderButter);
        shifter.preprocess_template_images(gather(single(cat(1, permute(data.MC.im_refIm_MC_2D, [3,1,2]), data.im_zstack))), py.int(0));
        
        data.MC.im_buffer_rolling_z = gpuArray(zeros([size(data.MC.im_refIm_MC_2D) , params.MC.numFrames_avgWin_zCorr], 'int32'));
    %     data.MC.im_buffer_rolling_z = gpuArray(zeros([size(data.MC.im_refIm_MC_2D) , params.MC.numFrames_avgWin_zCorr], 'int16'));
        data.MC.counter_buffer_rolling_z = 0;
    end
    
end

%% == Session Starting & counting ==

% == Start Session ==
if frameNum == 1
    disp('hi. NEW SESSION STARTED')
    startSession
    disp('frameNum = 1')
end
% ======== COMMENT THIS IN/OUT TO START SESSION ===========================
% startSession
% =========================================================================

%% == MOTION CORRECTION ==
% % FASTER ON CPU THAN GPU
% % Make a cropped version of the current image for use in motion correction
if strcmp(params.mode, 'imaging')
    data.MC.img_MC_2d_moving = data.currentImage_gpu(data.MC.idx_im_MC_crop_y(1):data.MC.idx_im_MC_crop_y(end), data.MC.idx_im_MC_crop_x(1):data.MC.idx_im_MC_crop_x(end));
    
    % Track frames for slow Z-axis motion correction
    if (frameNum >= 0) && (mod(frameNum, params.MC.intervalFrames_zCorr) == 0)
        data.MC.counter_buffer_rolling_z = data.MC.counter_buffer_rolling_z + 1;
        data.MC.im_buffer_rolling_z_mean = rolling_z_mean_obj.update_mean(rolling_z_mean_obj.idx_new, data.MC.img_MC_2d_moving, data.MC.im_buffer_rolling_z(:,:,mod(data.MC.counter_buffer_rolling_z ,params.MC.numFrames_avgWin_zCorr)+1), rolling_z_mean_obj.win_size, rolling_z_mean_obj.mean_old);
        
        rolling_z_mean_obj.mean_old = data.MC.im_buffer_rolling_z_mean;
        rolling_z_mean_obj.idx_new = rolling_z_mean_obj.idx_new + 1;
        
        data.MC.im_buffer_rolling_z(:,:,mod(data.MC.counter_buffer_rolling_z ,params.MC.numFrames_avgWin_zCorr)+1) = data.MC.img_MC_2d_moving;
    elseif frameNum < params.MC.intervalFrames_zCorr
        data.MC.im_buffer_rolling_z_mean = data.MC.im_buffer_rolling_z(:,:,1);
    end
    
    out = shifter.find_translation_shifts(gather(data.MC.img_MC_2d_moving), py.int(0));  %% 0-indexed
    shifts_yx          = single(int32(out{1}.numpy()));
    % shifts_yx          = single(int16(out{1}.numpy()));
    data.MC.yShift     = shifts_yx(1);
    data.MC.xShift     = shifts_yx(2);
    data.MC.maxCorr_2d = single(out{2}.numpy());
    
    out = shifter.find_translation_shifts(gather(data.MC.im_buffer_rolling_z_mean), py.list(int64([1:data.MC.n_slices_zstack])));  %% 0-indexed
    data.MC.maxCorr_z = single(out{2}.numpy());
    [maxVal, maxArg] = max(data.MC.maxCorr_z);
    data.MC.delta_z = (data.MC.idx_middle_frame-maxArg) * data.MC.stepSize_zstack;
    
    % data.MC.xShift = 0;
    % data.MC.yShift = 0;
    % data.MC.maxCorr_2d = 0;
    % data.MC.maxCorr_z = zeros(1,5);
    % data.MC.delta_z = 0;
    
    if abs(data.MC.xShift) >60
        data.MC.xShift = 0;
    end
    if abs(data.MC.yShift) >60
        data.MC.yShift = 0;
    end
end

%% == ROLLING STATS ==
sm.fakeFeedback_inUse = NaN;
if sm.CE_experimentRunning
    % 03/27/2023 block-trial preparation
    if params.blocks.block_trial
        if ~sm.CE_trial && ((frameNum - sm.blocks.frameNum_blockStart >= sm.blocks.block_timecap * params.timing.frameRate) || (sm.NumOfRewardsAcquired - sm.blocks.rewardNum_blockStart >= sm.blocks.block_rewardcap))
            disp(['Starting new block.'])
            temp_blockNum = sm.blocks.blockNum + 1;
            if temp_blockNum > params.blocks.blockNum_cap
                sm.blocks.blockNum = temp_blockNum;
                disp(['All Blocks Completed!'])
            else
                sm.blocks = baselineStuff.block_sequence.blocks(temp_blockNum);
                sm.cursor = baselineStuff.cursors(sm.blocks.cursor_to_use);
                sm.blocks.blockNum = temp_blockNum;
                disp(['Current Block: ', num2str(sm.blocks.blockNum), ' / ', num2str(params.blocks.blockNum_cap)])
            end
            sm.blocks.frameNum_blockStart = frameNum;
            sm.blocks.rewardNum_blockStart = sm.NumOfRewardsAcquired;
            disp(sm.blocks)
            disp(sm.cursor)
        end
    else
        % Set to defaults
        sm.cursor = params.cursor;
%         sm.blocks = params.blocks;
    end    
    %% Trial prep 
    sm.trialType_CSOn = trialStuff.condTrialBool(sm.trialNum,1);
    sm.trialType_optostim = trialStuff.condTrialBool(sm.trialNum,2);
    sm.trialType_rewardOn = trialStuff.condTrialBool(sm.trialNum,3);
    
    %%  ===== TRIAL STRUCTURE =====
    % CE = current epoch
    % ET = epoch transition signal
    % CS = current state
    
    % START BUILDING UP STATS
    if frameNum == 1
        sm.CE_buildingUpStats = 1;
        %     sm.soundVolume = 0;
        %     setSoundVolumeTeensy(sm.soundVolume);
    end
    % BUILDING UP STATS
    if sm.CE_buildingUpStats
        source.hSI.task_cursorAmplitude.writeDigitalData(0);
        source.hSI.task_goalAmplitude.writeDigitalData(0);
    end
    % END BUILDING UP STATS
    if sm.CE_buildingUpStats && frameNum > params.trial.duration_buildingUpStats
        sm.CE_buildingUpStats = 0;
        sm.ET_waitForBaseline = 1;
    end
        
    % START WAIT FOR BASELINE
    if sm.ET_waitForBaseline
        sm.ET_waitForBaseline = 0;
        sm.CE_waitForBaseline = 1;
        sm.soundVolume = 0;
        %     setSoundVolumeTeensy(sm.soundVolume);
        sm.counter_baseline_ENL = 0;
    end
    % WAIT FOR BASELINE: ENL required
    if sm.CE_waitForBaseline
        sm.counter_baseline_ENL = sm.counter_baseline_ENL + 1;
    end
    % END WAIT FOR BASELINE (QUIESCENCE ACHIEVED)
    if sm.counter_baseline_ENL > params.trial.duration_baseline_ENL
        % if readDigitalPin(lick_teensy, 'Pin')
        sm.baseline_ENL_violation = source.hSI.lick_detection.readDigitalData();
        if ~sm.baseline_ENL_violation
            sm.CE_waitForBaseline = 0;
            sm.ET_trialStart = 1;
            sm.baseline_ENL_violation = 1;
        end
    end
    
    % START TRIAL
    if sm.ET_trialStart
        sm.ET_trialStart = 0;
        sm.CE_trial = 1;
        sm.counter_trialIdx = 0;
        
        updateLoggerTrials_START
        if sm.trialType_CSOn
            sm.soundVolume = 1;
        else
            sm.soundVolume = 0;
        end

        if sm.trialType_optostim
            sm.counter_optostim = 0;
        else
            sm.counter_optostim = nan;
        end

        %     setSoundVolumeTeensy(sm.soundVolume);
        source.hSI.task_cursorAmplitude.writeDigitalData(sm.soundVolume);
        
        sm.frequencyOverride = 0;
    end
    % COUNT TRIAL DURATION
    if sm.CE_trial
        sm.counter_trialIdx = sm.counter_trialIdx + 1;
    end

    % Optogenetic stimulation
    if sm.counter_trialIdx >= round(params.timing.frameRate * params.trial.delay_optostim)
        giveoptostim()
    end
    
    % ENL violation check before reward delivery
    if sm.CE_trial && sm.counter_trialIdx >= round(params.timing.frameRate * params.trial.duration_postCS_ENL)
        sm.postCS_ENL_violation = source.hSI.lick_detection.readDigitalData();

        % END TRIAL: FAILURE
        if sm.postCS_ENL_violation
            sm.CE_trial = 0;
            sm.ET_timeout = 1;
            sm.counter_trialIdx = NaN;
            updateLoggerTrials_END(0)
            sm.trialNum = sm.trialNum+1;
            sm.postCS_ENL_violation = 1;  
    
        % END TRIAL: Success
        else
            updateLoggerTrials_END(1)
            sm.CE_trial = 0;
            %     ET_rewardToneHold = 1;
            sm.ET_rewardDelivery = 1;
            sm.trialNum = sm.trialNum+1;
            sm.counter_trialIdx = NaN;
        end
    end

    % ENL violation -> START TIMEOUT
    if sm.ET_timeout
        sm.ET_timeout = 0;
        sm.CE_timeout = 1;
        sm.counter_timeout = 0;
        sm.soundVolume = 0;
        sm.NumOfTimeouts = sm.NumOfTimeouts + 1;
        source.hSI.task_cursorAmplitude.writeDigitalData(sm.soundVolume);
    end
    % COUNT TIMEOUT DURATION
    if sm.CE_timeout
        sm.counter_timeout = sm.counter_timeout + 1;
    end
    % END TIMEOUT
    if sm.CE_timeout && sm.counter_timeout >= round(params.timing.frameRate * params.trial.duration_timeout)
        sm.CE_timeout = 0;
        sm.ET_ITI_withZ = 1;
    end

    
    % ENL Pass -> START DELIVER REWARD
    if sm.ET_rewardDelivery
        sm.ET_rewardDelivery = 0;
        sm.CE_rewardDelivery = 1;
        sm.counter_rewardDelivery = 0;
        sm.frequencyOverride = 1;
        %     sm.soundVolume = 0;
        %     setSoundVolumeTeensy(sm.soundVolume);
        if sm.trialType_rewardOn
            %     giveReward2(1, 1, reward_duration, reward_delay, LED_duration, 1, LED_ramp_duration); % in ms. This is in the START section so that it only delivers once
            giveReward3(source, 1, 0, params.trial.reward_duration, params.trial.reward_delay, params.trial.LED_duration, 1, params.trial.LED_ramp_duration); % in ms. This is in the START section so that it only delivers once
            %         giveReward3(source, 1, 0, 500, reward_delay, LED_duration, 1, LED_ramp_duration); % in ms. This is in the START section so that it only delivers once
        end
        sm.NumOfRewardsAcquired = sm.NumOfRewardsAcquired + 1;
        
        %         save([params.paths.directory , '\logger.mat'], 'logger')
        %         saveParams(params.paths.directory)
        %         disp(['Logger & Params Saved: frameCounter = ' num2str(frameNum)]);
    end
    % COUNT DELIVER REWARD
    if sm.CE_rewardDelivery
        sm.counter_rewardDelivery = sm.counter_rewardDelivery + 1;
    end
    % END DELIVER REWARD
    if sm.CE_rewardDelivery && sm.counter_rewardDelivery >= round(params.timing.frameRate * params.trial.duration_rewardDelivery)
        sm.CE_rewardDelivery = 0;
        sm.ET_ITI_withZ = 1;
        sm.frequencyOverride = 0;
    end
    
    sm.delta_moved = 0; % place holder to potentially be overwritten by 'moveFastZ' function below
    % START INTER-TRIAL-INTERVAL (POST-REWARD): WITH Z-CORRECTION
    if sm.ET_ITI_withZ
        sm.ET_ITI_withZ = 0;
        sm.CE_ITI_withZ = 1;
        sm.counter_ITI_withZ = 0;
        sm.soundVolume = 0;
        source.hSI.task_cursorAmplitude.writeDigitalData(sm.soundVolume);
        source.hSI.task_goalAmplitude.writeDigitalData(0);
        
        if strcmp(params.mode, 'imaging') && (frameNum > (sm.counter_last_z_correction + params.MC.min_interval_z_correction))
            if data.MC.delta_z ~=0
                clampedDelta = sign(data.MC.delta_z) * min(abs(data.MC.delta_z), params.MC.max_delta_z_correction);
                [~] = moveFastZ(source, [], clampedDelta, [], [20,380]);
                disp(['moving fast Z by one step: ', num2str(clampedDelta), ', new position: ', num2str(data.MC.current_position_z + clampedDelta)]) %num2str(data.MC.delta_z)])
                sm.delta_moved = clampedDelta;
                sm.counter_last_z_correction = frameNum;
            elseif abs(data.MC.maxCorr_z(data.MC.idx_middle_frame + 1) - data.MC.maxCorr_z(data.MC.idx_middle_frame - 1)) > abs(max([data.MC.maxCorr_z(data.MC.idx_middle_frame - 1), data.MC.maxCorr_z(data.MC.idx_middle_frame + 1)] - data.MC.maxCorr_z(data.MC.idx_middle_frame)))
                clampedDelta = sign(data.MC.maxCorr_z(data.MC.idx_middle_frame - 1) - data.MC.maxCorr_z(data.MC.idx_middle_frame + 1)) * params.MC.max_delta_z_correction;
                [~] = moveFastZ(source, [], clampedDelta, [], [20,380]);
                disp(['moving fast Z by one step: ', num2str(clampedDelta), ', new position: ', num2str(data.MC.current_position_z + clampedDelta)]) %num2str(data.MC.delta_z)])
                sm.delta_moved = clampedDelta;
                sm.counter_last_z_correction = frameNum;
            end
        end
    end
    % COUNT INTER-TRIAL-INTERVAL (POST-REWARD)
    if sm.CE_ITI_withZ
        sm.counter_ITI_withZ = sm.counter_ITI_withZ + 1;
    end
    % END INTER-TRIAL-INTERVAL
    if sm.CE_ITI_withZ && sm.counter_ITI_withZ >= round(params.timing.frameRate * params.trial.duration_ITI)
        sm.counter_ITI_withZ = NaN;
        sm.CE_ITI_withZ = 0;
        sm.ET_waitForBaseline = 1;
    end
    
end
%% Teensy Output calculations: 2.8 kHz cue

if sm.CE_experimentRunning
    % data.voltage_cursorCurrentPos = convert_cursor_to_voltage(sm.cursor.threshold_reward , sm.cursor.bounds_cursor , sm.cursor.voltage_at_threshold);
    data.voltage_cursorCurrentPos = 1.2;
    source.hSI.task_cursorCurrentPos.writeAnalogData(double(data.voltage_cursorCurrentPos));
    % data.freqToOutput = convert_voltage_to_frequency(data.voltage_cursorCurrentPos , 3.3 , sm.cursor.range_freqOutput); % for logging purposes only. function should mimic (exactly) the voltage to frequency transformation on teensy
end

%% Plotting

if frameNum>1
    
    plotUpdatedOutput2([sm.CE_waitForBaseline*0.1, sm.CE_trial*0.2,...
        sm.CE_rewardDelivery*0.3, sm.CE_timeout*0.4, sm.CE_buildingUpStats*0.5, sm.fakeFeedback_inUse*0.6],...
        params.timing.duration_plotting, params.timing.frameRate, 'Rewards', 10, 22, ['# Rewards: ' , num2str(sm.NumOfRewardsAcquired) , ' ; # Timeouts: ' , num2str(sm.NumOfTimeouts)])
    plotUpdatedOutput3([data.MC.xShift' data.MC.yShift'], params.timing.duration_plotting, params.timing.frameRate, 'Motion Correction Shifts', 10, 11)
    
    if frameNum > 25
        plotUpdatedOutput4(nanmean(logger.motionCorrection(frameNum-15:frameNum,3),1),...
            params.timing.duration_plotting, params.timing.frameRate, 'Motion Correction Correlation Rolling', 10, 12)
    end
    
    plotUpdatedOutput6([data.ROIs.cursor_output, data.ROIs.cursor_brain],...
        params.timing.duration_plotting, params.timing.frameRate, ['cursor_output, ', 'cursor_brain'] , 10,3)
    
    
    if mod(frameNum,30) == 0 && frameNum > 300
        plotUpdatedOutput5([nanmean(logger.motionCorrection(frameNum-300:frameNum-1,3),1)],...
            params.timing.duration_session, params.timing.frameRate, 'Motion Correction Correlation All', 10, 1)
    end
    
    plotUpdatedOutput7(data.MC.maxCorr_z,...
        params.timing.duration_plotting, params.timing.frameRate, 'Z Frame Correlations', 10, 10)
    
    % % Below are too slow for online calculation
    %     plotUpdatedOutput(data.ROIs.dFoF(:,1:10) + [1:10]*2, params.timing.duration_plotting, params.timing.frameRate, 'neurons', 1, 1),
    %     plotUpdatedOutput8(data.ROIs.decoder_magnitudes + [1:length(data.ROIs.decoder_magnitudes)]*100, params.timing.duration_plotting, params.timing.frameRate, 'magnitudes', 5, 1),
    %     plotUpdatedOutput9(data.ROIs.decoder_angles + [1:length(data.ROIs.decoder_angles)]*1.5, params.timing.duration_plotting, params.timing.frameRate, 'angles', 5, 1),
    
    %     drawnow
end

%% DATA LOGGING
if ~isnan(frameNum)
    logger.timeSeries(frameNum,1) = frameNum;
    logger.timeSeries(frameNum,2) = sm.CS_quiescence;
    logger.timeSeries(frameNum,3) = sm.ET_trialStart;
    logger.timeSeries(frameNum,4) = sm.CE_trial;
    logger.timeSeries(frameNum,5) = sm.soundVolume;
    logger.timeSeries(frameNum,6) = sm.counter_trialIdx;
    logger.timeSeries(frameNum,7) = sm.CS_threshold;
    logger.timeSeries(frameNum,8) = sm.ET_rewardToneHold; % reward signals
    logger.timeSeries(frameNum,9) = sm.CE_rewardToneHold;
    logger.timeSeries(frameNum,10) = sm.counter_rewardToneHold;
    logger.timeSeries(frameNum,11) = sm.frequencyOverride;
    logger.timeSeries(frameNum,12) = sm.ET_rewardDelivery;
    logger.timeSeries(frameNum,13) = sm.CE_rewardDelivery;
    logger.timeSeries(frameNum,14) = sm.counter_rewardDelivery;
    logger.timeSeries(frameNum,15) = sm.ET_ITI_withZ;
    logger.timeSeries(frameNum,16) = sm.CE_ITI_withZ;
    logger.timeSeries(frameNum,17) = sm.counter_ITI_withZ;
    logger.timeSeries(frameNum,18) = sm.ET_waitForBaseline;
    logger.timeSeries(frameNum,19) = sm.CE_waitForBaseline;
    logger.timeSeries(frameNum,20) = sm.ET_timeout;
    logger.timeSeries(frameNum,21) = sm.CE_timeout;
    logger.timeSeries(frameNum,22) = sm.counter_timeout;
    logger.timeSeries(frameNum,23) = sm.CE_buildingUpStats;
    logger.timeSeries(frameNum,24) = sm.CE_experimentRunning;
    logger.timeSeries(frameNum,25) = sm.NumOfRewardsAcquired;
    logger.timeSeries(frameNum,26) = sm.NumOfTimeouts;
    logger.timeSeries(frameNum,27) = data.hash_image;
    logger.timeSeries(frameNum,28) = sm.trialNum;
    logger.timeSeries(frameNum,29) = sm.trialNum * sm.CE_trial;
    logger.timeSeries(frameNum,30) = sm.fakeFeedback_inUse;
    logger.timeSeries(frameNum,31) = sm.trialType_cursorOn;
    logger.timeSeries(frameNum,32) = sm.trialType_feedbackLinked;
    logger.timeSeries(frameNum,33) = sm.trialType_rewardOn;
    logger.timeSeries(frameNum,34) = sm.counter_last_z_correction;
    logger.timeSeries(frameNum,35) = sm.delta_moved;
    logger.timeSeries(frameNum,36) = sm.blocks.blockNum;
    
    logger.timers(frameNum,1) = now;
    logger.timers(frameNum,2) = toc;
    
    if strcmp(params.mode, 'imaging')
        logger.motionCorrection(frameNum,1) = gather(data.MC.xShift);
        logger.motionCorrection(frameNum,2) = gather(data.MC.yShift);
        logger.motionCorrection(frameNum,3) = gather(data.MC.maxCorr_2d(1));
        logger.motionCorrection(frameNum,4) = data.MC.current_position_z;
        logger.motionCorrection(frameNum,5) = sm.delta_moved;
        logger.motionCorrection(frameNum,6:10) = data.MC.maxCorr_z(1:end);
    end
end

%% End Session
% if  sm.CE_experimentRunning && (frameNum == round(params.timing.duration_session * 0.90))
%     %     source.hSI.task_cursorAmplitude.writeDigitalData(0);
%     %     source.hSI.task_goalAmplitude.writeDigitalData(0);
%     saveLogger(params.paths.directory);
%     saveParams(params.paths.directory);
%     %     source.hSI.task_cursorAmplitude.writeDigitalData(1);
%     %     source.hSI.task_goalAmplitude.writeDigitalData(1);
% end


if sm.CE_experimentRunning && ((frameNum == round(params.timing.duration_session * 0.98)) || (sm.blocks.blockNum > params.blocks.blockNum_cap))
    sm.CE_experimentRunning = 0;
    endSession
end

%% FUNCTIONS
    function updateLoggerTrials_START % calls at beginning of a trial
        logger.trials(sm.trialNum,1) = sm.trialNum;
        logger.trials(sm.trialNum,2) = now;
        logger.trials(sm.trialNum,3) = frameNum;
        logger.trials(sm.trialNum,4) = sm.trialType_CSOn;
        logger.trials(sm.trialNum,5) = sm.trialType_optostim;
        logger.trials(sm.trialNum,6) = sm.trialType_rewardOn;
        logger.trials(sm.trialNum,7) = sm.blocks.blockNum;
    end
    function updateLoggerTrials_END(success_outcome) % calls at end of a trial
        logger.trials(sm.trialNum,9) = sm.trialNum;
        logger.trials(sm.trialNum,10) = now;
        logger.trials(sm.trialNum,11) = frameNum;
        logger.trials(sm.trialNum,12) = success_outcome;
    end
    function startSession
        % If block-trial, INITIALIZE PARAMETERS
        if strcmp(params.mode, 'BMI') && params.blocks.block_trial
            params.blocks.blockNum_cap       = baselineStuff.block_sequence.blockNum_cap; % After this number of block trials, finish the exp.

            sm.blocks = baselineStuff.block_sequence.blocks(1);
            sm.cursor = baselineStuff.cursors(sm.blocks.cursor_to_use);
            sm.blocks.blockNum = 1; % 1-INDEXED!!!
            sm.blocks.frameNum_blockStart =  params.trial.duration_buildingUpStats;
            sm.blocks.rewardNum_blockStart = 0;
            
            disp(sm.blocks)
            disp(sm.cursor)
            disp("dFoF parameters")
            disp(params.dFoF)
        else
            sm.blocks  = struct();
            sm.blocks.blockNum = 0;
            params.blocks.blockNum_cap = inf;
            sm.cursor = params.cursor();
            disp("No Block Trial: Load Default Parameter")
            disp(sm.cursor)
        end
        
        % INITIALIZE VARIABLES
        sm.CE_waitForBaseline = 0;
        sm.CS_quiescence = 0;
        sm.ET_trialStart = 0;
        sm.CE_trial = 0;
        sm.soundVolume = 0;
        sm.counter_trialIdx = 0;
        sm.CS_threshold = 0;
        sm.ET_rewardToneHold = 0; % reward signals
        sm.CE_rewardToneHold = 0;
        sm.counter_rewardToneHold = 0;
        sm.frequencyOverride = 0;
        sm.ET_rewardDelivery = 0;
        sm.CE_rewardDelivery = 0;
        sm.counter_rewardDelivery = 0;
        sm.ET_ITI_withZ = 0;
        sm.CE_ITI_withZ = 0;
        sm.counter_ITI_withZ = 0;
        sm.ET_waitForBaseline = 0;
        sm.CE_waitForBaseline = 0;
        sm.ET_timeout = 0;
        sm.CE_timeout = 0;
        sm.counter_timeout = 0;
        sm.counter_last_z_correction = 0;
        sm.CE_buildingUpStats = 1;
        sm.CE_experimentRunning = 1;


        sm.baseline_ENL_violation = 1;
        sm.postCS_ENL_violation = 1;

        
        sm.NumOfRewardsAcquired = 0;
        sm.NumOfTimeouts = 0;
        sm.trialNum = 1;
        
        loggerNames.timeSeries{1} = 'frameNum';
        loggerNames.timeSeries{2} = 'CS_quiescence';
        loggerNames.timeSeries{3} = 'ET_trialStart';
        loggerNames.timeSeries{4} = 'CE_trial';
        loggerNames.timeSeries{5} = 'soundVolume';
        loggerNames.timeSeries{6} = 'counter_trialIdx';
        loggerNames.timeSeries{7} = 'CS_threshold';
        loggerNames.timeSeries{8} = 'ET_rewardToneHold'; % reward signals
        loggerNames.timeSeries{9} = 'CE_rewardToneHold';
        loggerNames.timeSeries{10} = 'counter_rewardToneHold';
        loggerNames.timeSeries{11} = 'frequencyOverride';
        loggerNames.timeSeries{12} = 'ET_rewardDelivery';
        loggerNames.timeSeries{13} = 'CE_rewardDelivery';
        loggerNames.timeSeries{14} = 'counter_rewardDelivery';
        loggerNames.timeSeries{15} = 'ET_ITI_withZ';
        loggerNames.timeSeries{16} = 'CE_ITI_withZ';
        loggerNames.timeSeries{17} = 'counter_ITI_withZ';
        loggerNames.timeSeries{18} = 'ET_waitForBaseline';
        loggerNames.timeSeries{19} = 'CE_waitForBaseline';
        loggerNames.timeSeries{20} = 'ET_timeout';
        loggerNames.timeSeries{21} = 'CE_timeout';
        loggerNames.timeSeries{22} = 'counter_timeout';
        loggerNames.timeSeries{23} = 'CE_buildingUpStats';
        loggerNames.timeSeries{24} = 'CE_experimentRunning';
        loggerNames.timeSeries{25} = 'NumOfRewardsAcquired';
        loggerNames.timeSeries{26} = 'NumOfTimeouts';
        loggerNames.timeSeries{27} = 'image_hash';
        loggerNames.timeSeries{28} = 'trialNum';
        loggerNames.timeSeries{29} = 'trialNum*CE_trial';
        loggerNames.timeSeries{30} = 'fakeFeedback_inUse';
        loggerNames.timeSeries{31} = 'trialType_cursorOn';
        loggerNames.timeSeries{32} = 'trialType_feedbackLinked';
        loggerNames.timeSeries{33} = 'trialType_rewardOn';
        loggerNames.timeSeries{34} = 'counter_last_z_correction';
        loggerNames.timeSeries{35} = 'delta_moved';
        loggerNames.timeSeries{36} = 'blockNum';
        
        loggerNames.timers{1} = 'time_now';
        loggerNames.timers{2} = 'tic_toc';
        
        loggerNames.motionCorrection{1} = 'xShift';
        loggerNames.motionCorrection{2} = 'yShift';
        loggerNames.motionCorrection{3} = 'MC_correlation';
        loggerNames.motionCorrection{4} = 'current_position_z';
        loggerNames.motionCorrection{5} = 'deltaMoved';
        loggerNames.motionCorrection{6} = 'z_correlation_1';
        loggerNames.motionCorrection{7} = 'z_correlation_2';
        loggerNames.motionCorrection{8} = 'z_correlation_3';
        loggerNames.motionCorrection{9} = 'z_correlation_4';
        loggerNames.motionCorrection{10} = 'z_correlation_5';
        
        loggerNames.trials{1} = 'trialNum_trialStart';
        loggerNames.trials{2} = 'time_now_trialStart';
        loggerNames.trials{3} = 'frameNum_trialStart';
        loggerNames.trials{4} = 'trialType_cursorOn';
        loggerNames.trials{5} = 'trialType_feedbackLinked';
        loggerNames.trials{6} = 'trialType_rewardOn';
        loggerNames.trials{7} = 'blockNum';
        loggerNames.trials{8} = 'factor_to_use';
        loggerNames.trials{9} = 'trialNum_trialEnd';
        loggerNames.trials{10} = 'time_now_trialEnd';
        loggerNames.trials{11} = 'frameNum_trialEnd';
        loggerNames.trials{12} = 'success_outcome';
        
        %         clear logger
        logger.timeSeries = NaN(params.timing.duration_session, length(loggerNames.timeSeries));
        logger.timers = NaN(params.timing.duration_session, length(loggerNames.timers));
        logger.motionCorrection = NaN(params.timing.duration_session,  length(loggerNames.motionCorrection));
        logger.trials = NaN(size(trialStuff.condTrials,1),  length(loggerNames.trials));

        saveParams(params.paths.directory)
    end

    function endSession
        disp('SESSION OVER')
        frameNum = NaN;
        sm.CE_experimentRunning = 0;
        
        saveLogger(params.paths.directory)
        saveParams(params.paths.directory)
        disp('=== Loggers and expParams saved ===')
        
        
        sm.CE_waitForBaseline = 0;
        sm.CS_quiescence = 0;
        sm.ET_trialStart = 0;
        sm.CE_trial = 0;
        sm.soundVolume = 0;
        sm.counter_trialIdx = 0;
        sm.CS_threshold = 0;
        sm.ET_rewardToneHold = 0; % reward signals
        sm.CE_rewardToneHold = 0;
        sm.counter_rewardToneHold = 0;
        sm.frequencyOverride = 0;
        sm.ET_rewardDelivery = 0;
        sm.CE_rewardDelivery = 0;
        sm.counter_rewardDelivery = 0;
        sm.ET_ITI_withZ = 0;
        sm.CE_ITI_withZ = 0;
        sm.counter_ITI_withZ = 0;
        sm.ET_waitForBaseline = 0;
        sm.CE_waitForBaseline = 0;
        sm.ET_timeout = 0;
        sm.CE_timeout = 0;
        sm.counter_timeout = 0;
        
        %         frameNum = 0;
        sm.CE_buildingUpStats = 0;
        %         sm.CE_experimentRunning = 0;
        data.ROIs.cursor_output = 0;
        
        %         setSoundVolumeTeensy(0);
        source.hSI.task_cursorAmplitude.writeDigitalData(0);
        source.hSI.task_goalAmplitude.writeDigitalData(0);
        
    end

    function saveLogger(directory)
        disp("Saving logger.mat...")
        save([directory , '\logger.mat'], 'logger');
    end

    function saveParams(directory)
        expParams.params = params;
        expParams.directory = directory;
        
        expParams.image_hash_function = 'hash = sum(sum(image,1).^2)';
        
        expParams.loggerNames = loggerNames;
        
        expParams.baselineStuff = baselineStuff; % Too Big!
        
%         save([directory , '\expParams.mat'], 'expParams')
%         save([directory , '\expParams.mat'], 'expParams', '-v7')
        save([directory , '\expParams.mat'], 'expParams', '-v7.3')
        %         save([directory , '\motionCorrectionRefImages.mat'], 'motionCorrectionRefImages')
    end

%     function [refIm_crop_conjFFT_shift, refIm_crop, indRange_y_crop, indRange_x_crop] = make_fft_for_MC(refIm)
%         refIm = single(refIm);
%         % crop_factor = 5;
%         crop_size = 256; % MAKE A POWER OF 2! eg 32,64,128,256,512
%
%         length_x = size(refIm,2);
%         length_y = size(refIm,1);
%         middle_x = size(refIm,2)/2;
%         middle_y = size(refIm,1)/2;
%
%         % indRange_y_crop = [round(middle_y - length_y/crop_factor) , round(middle_y + length_y/crop_factor) ];
%         % indRange_x_crop = [round(middle_x - length_y/crop_factor) , round(middle_x + length_y/crop_factor) ];
%
%         indRange_y_crop = [round(middle_y - (crop_size/2-1)) , round(middle_y + (crop_size/2)) ];
%         indRange_x_crop = [round(middle_x - (crop_size/2-1)) , round(middle_x + (crop_size/2)) ];
%
%         refIm_crop = refIm(indRange_y_crop(1) : indRange_y_crop(2) , indRange_x_crop(1) : indRange_x_crop(2)) ;
%
%         refIm_crop_conjFFT = conj(fft2(refIm_crop));
%         refIm_crop_conjFFT_shift = fftshift(refIm_crop_conjFFT);
%
%         % size(refIm_crop_conjFFT_shift,1);
%         % if mod(size(refIm_crop_conjFFT_shift,1) , 2) == 0
%         %     disp('RH WARNING: y length of refIm_crop_conjFFT_shift is even. Something is very wrong')
%         % end
%         % if mod(size(refIm_crop_conjFFT_shift,2) , 2) == 0
%         %     disp('RH WARNING: x length of refIm_crop_conjFFT_shift is even. Something is very wrong')
%         % end
%
% %         refIm_crop_conjFFT_shift_centerIdx = ceil(size(refIm_crop_conjFFT_shift)/2);
%     end

%     function [delta,frame_corrs, xShifts, yShifts] = calculate_z_position(img_MC_moving_rolling_z, registrationImage, refIm_crop_conjFFT_shift, referenceDiffs, maskPref, borderOuter, borderInner)
%         image_toUse = mean(img_MC_moving_rolling_z, 3);
%         [delta, frame_corrs, xShifts, yShifts] = zCorrection(image_toUse, registrationImage, ...
%             refIm_crop_conjFFT_shift, referenceDiffs, maskPref, borderOuter, borderInner);
%
%     end

    function currentPosition = moveFastZ(source, evt, delta, position, range_position)
        
        if ~exist('range_position')
            range_position = [0, 200];
        end
        
        fastZDevice = source.hSI.hFastZ.currentFastZs{1};
        %Select the FastZ device (you likely have just one, so index at 1)
        
        currentPosition = fastZDevice.targetPosition;
        
        if ~exist('position')
            newPosition = currentPosition + delta;
        elseif isempty(position)
            newPosition = currentPosition + delta;
        else
            newPosition = position;
        end
        % scalar finite number indicating depth within the lower and upper travel bounds set by the user.
        
        if range_position(1) > newPosition | range_position(2) < newPosition
            error(['RH ERROR: newPosition if out of range. Range: ', range_position, ' Attempted position: ', newPosition])
        end
        
                force = true;
        % do the move even if it is a grab or loop acquisition. Don't try using this with a stack acquisition.
        
        %     source.hSI.hFastZ.move(fastZDevice, position,  force);
        
        source.hSI.hFastZ.move(fastZDevice, newPosition);
        
    end
%     toc
end