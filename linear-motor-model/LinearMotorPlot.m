% Run this to connect to the model
mphopen Linearmotor
accData = model.result.table("tbl5").getReal();
oldFrequency = 0;
for samples = 1:length(accData)
    frequency = accData(samples);
    if oldFrequency ~= frequency
        samples = samples - 1;
        break
    end
end

clf;
hold on;
for i = 1:(length(accData)/samples)
    startSample = (i-1)*samples+1;
    stopSample = startSample + samples-1;
    sample = accData(startSample:stopSample,:);
    frequency = sample(:,1);
    force = sample(:,3);
    secondHalf = floor(length(force)*1/2):length(force);
    plot(frequency(secondHalf), force(secondHalf));
    plot(frequency(secondHalf), mean(force(secondHalf)), "o");
end
xlabel("Frequency")
ylabel("Force")
title("50cm stator")