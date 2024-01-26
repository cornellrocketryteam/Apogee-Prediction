% RUN THIS SCRIPT AFTER RUNNING CDR_flight_modeling








function delta = ChangeRate(prediction_resids)
% gives the differences between adjacent entries of a vector

delta = []; 

for i = 2:length(prediction_resids)
    delta(i-1) = prediction_resids(i) - prediction_resids(i-1);
end

end

function SMA_vector = SMA(vector,L)
% gives the SMA of a vector
% vector is the vector you want to obtain the SMA for
% L is the length of the SMA

SMA_vector = [];

for i = L:length(vector)
    SMA_vector(i-(L-1)) = mean(vector(i-(L-1):L));
end


end
