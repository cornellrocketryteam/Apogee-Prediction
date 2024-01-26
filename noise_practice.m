close all
sigma = 6.448; % standard deviation
mean = 30;

points = zeros(1,10000);

for i = 1:length(points)
    r = rand();
    z = -sqrt(2) * erfcinv(r*2);
    points(i) = mean + z*sigma;
end

histogram(points)