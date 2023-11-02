function optostim_train(source, train_delay_in_ms, train_length_in_ms, train_frequency, pulse_length_in_ms)

traceDuration = ceil(((train_length_in_ms*1.1 + train_delay_in_ms)/1000) * source.hSI.sesh_reward.Rate);

train_t = 0 : 1/source.hSI.sesh_reward.Rate : round(train_length_in_ms/1000 * 1.1);
train_repetition = 0 + pulse_length_in_ms/1000 : 1/train_frequency : round(train_length_in_ms/1000 * 1.1);
pulse_train = pulstran(train_t, train_repetition, @rectpuls, pulse_length_in_ms/1000);

output_pulse_train = zeros(traceDuration,1);
delayFrames = round((train_delay_in_ms/1000) * source.hSI.sesh_reward.Rate);
output_pulse_train(delayFrames+1 : delayFrames + numel(pulse_train)) = pulse_train;

% queueOutputData(source.hSI.sesh_reward, [ outputTrace_solenoid     outputTrace_5VForLickDetection   outputTrace_LED ]); % [ (reward solenoid) (lick detection voltage)]
queueOutputData(source.hSI.sesh_optostim, [ output_pulse_train ]); % [ (reward solenoid) (lick detection voltage)]

listener_in = addlistener(source.hSI.sesh_optostim,'DataAvailable',@(sesh,event)nullFunction); % check my old function liveTF for how to use this type of function. Here it is necessary to use because we have an analog channel (being used purely as a clock for the digital output channel). This requires a listener for some reason.

startBackground(source.hSI.sesh_optostim);

end