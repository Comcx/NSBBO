%% The problem that the EA is going to solve

function [init, cost, feasible] = problem

init = @problemInit;
feasible = @problemFeasible;
cost = @problemCost;
end



%% Init Problem
function Population = problemInit(Config)

global r;
global k;     k = 0.005;
global y;

% Initialize population
for popIndex = 1 : Config.popSize
   for chromIndex = 1 : Config.chromLength
      Population(popIndex).chrom(chromIndex) = rand();
   end
end

end


function Population = problemCost(Population, Config)

global k;
global r;
global y;
global a;     a = 0.05;
global c;

for popIndex = 1 : Config.popSize
   vector = Population(popIndex).chrom' * r;
   punish = max( 0 , abs( sum(Population(popIndex).chrom) -1+c ) - 10^(-4) ) + ...
            max( 0 , sum( Population(popIndex).chrom-y ) ) + ...
            max( 0 , -( cvar(vector) + k*sum(vector) ) );
   Population(popIndex).object(1) = 1 / sum(vector);
   Population(popIndex).object(2) = 1 / skewness(vector);
   Population(popIndex).object(3) = kurtosis(vector);
   Population(popIndex).object(4) = length( find(vector<k) );
   Population(popIndex).object(5) = punish;
end



end


%% Make each population feasible
function Population = problemFeasible(Population, Config)

for i = 1 : length(Population)
    for k = 1 : length(Population(1).chrom)
       if Population(i).chrom(k) < 0 || Population(i).chrom(k) > 1
          Population(i).chrom(k) = rand();
       end
    end
    
end
end

