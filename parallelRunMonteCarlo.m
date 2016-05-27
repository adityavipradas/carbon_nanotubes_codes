%Run Monte Carlo simulations
function [] = parallelRunMonteCarlo(film_length, CNT_length, num_trials, CNT_increment)
tic;
CNT_numList = [];
numTrials = [];
shortDist = [];
matlabpool(2);
parfor i = 1:1:num_trials
    disp(i);
    mindist = 0;
    CNT_num = 0;
    while mindist == 0
        CNT_num = CNT_num + CNT_increment;
        mindist = parallel_cnt_shortest_path(CNT_num, film_length, CNT_length);
    CNT_numList(i) = CNT_num;
    numTrials(i) = i;
    shortDist(i) = mindist;
    end
end
disp(CNT_numList);
disp(shortDist);

mkdir('monteCarloResults');
cd('monteCarloResults');
vid = 'montecarlo';

h1 = figure();
hold on;
grid on;
scatter(numTrials, CNT_numList, 'filled')
%set(gca,'XTick',1:1:num_trials);
title(sprintf('Number of %d um long CNTs in %d x %d domain for %d increment (mean = %f, mode = %f)',...
    CNT_length, film_length, film_length, CNT_increment, mean(CNT_numList), mode(CNT_numList)))
xlabel('Trial number')
ylabel('Number of CNTs giving connectivity')
vname = [vid 'Num' 'jpg'];
saveas(h1, vname, 'jpg');
set(clf,'visible','off');

h2 = figure();
hold on;
grid on;
scatter(numTrials, shortDist, 'filled')
%set(gca,'XTick',1:1:num_trials);
title(sprintf('Connectivity distance for %d um long CNTs in %d x %d domain for %d CNT increment (mean = %f um)',...
    CNT_length, film_length, film_length, CNT_increment, mean(shortDist)))
xlabel('Trial number')
ylabel('Connectivity distance between top and bottom edges(um)')
vname = [vid 'Dist' 'jpg'];
saveas(h2, vname, 'jpg');
set(clf,'visible','off');
toc;






