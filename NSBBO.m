
%% NSBBO ¡ª¡ª Non-domination Sort EA

function  ans  = NSBBO( Problem )

% Non-domination Sort Biogeography-based optimization (NSBBO) software for
% minimizing a few continuous functions

% INPUTS: ProblemFunction = the handle of the function that returns the handles of the initialization, cost, and feasibility functions.
% DisplayFlag = true or false, whether or not to display and plot results.
% RandSeed = random number seed
% OUTPUTS: MinCost = array of best solution, one element for each generation



%% Initialization start

Config.maxGen = 50;               % generation count limit 
Config.pModify = 0.90;             % habitat modification probability (between 0 and 1)0.95
Config.pMutate = 0.35;             % mutation probability0.25
Config.keep = 30;                  % elitism parameter: how many of the best habitats to keep from one generation to the next
Config.dyList = [ 0.4, 0.95, 0.4, 0.05 ];
                                   % dynamic list to adjust
Config.popSize = 100;              % population size
Config.chromLength = m;           % number of variables in each population chrom (i.e., problem dimension)
Config.objectNum = 4;              % number of objects

[ cost, feasible, Population ] = NSBBOInit( Problem, Config );
Island = Population;
record = zeros(1, Config.maxGen);


%% Starting to evolve

% A few pre-steps
% Sort from best to worst
Population = cost(Population, Config);
Population = popSort(Population, Config);

for genIndex = 1 : Config.maxGen
    
    % To save the elite ppulation
    Prior = keepPrior(Population, Config.keep);
    
    %% Modify
    [ mu , lambda ] = getRatio( Population );
    
    for popIndex = 1 : Config.popSize
           if rand() > Config.pModify
                  continue;
           end
           for positionIndex = 1 : Config.chromLength
                   if rand() < lambda(popIndex)
                          randomNum = rand * sum(mu);
                          Select = mu(1);
                          selectIndex = 1;
                          while (randomNum > Select) && (selectIndex < Config.popSize)
                                 selectIndex = selectIndex + 1;
                                 Select = Select + mu(selectIndex);
                          end
                          Island(popIndex).chrom(positionIndex) = (genIndex/Config.maxGen) * Population(popIndex).chrom(positionIndex) + ...
                                                               ( 1-genIndex/Config.maxGen) * Population(selectIndex).chrom(positionIndex);
                   end
           end
    end
    
    
    %% Mutate
    F = 0.5;
    for popIndex = 1 : Config.popSize
        for positionIndex = 1 : Config.chromLength
            if rand() < Config.pMutate
                  Island(popIndex).chrom(positionIndex) =  Population(popIndex).chrom(positionIndex) + ...
                                                             F*(Population(1).chrom(positionIndex) - Population(popIndex).chrom(positionIndex))+ ...
                                                             F*(Population(ceil(rand()*Config.popSize)).chrom(positionIndex) - ...
                                                                Population(ceil(rand()*Config.popSize)).chrom(positionIndex));
                       
            end
        end
    end
    
    
    %% After-steps
    
    Population = Island;
    
    % Make sure the population does not have duplicates. 
    Population = clearDups(Population);
    Population = feasible(Population, Config);
    Population = cost(Population, Config);
    
    % Sort from best to worst
    Population = popSort(Population, Config);
    
    Population = mix(Population, Prior);
    Population = popSort(Population, Config);
    Population = clearDups(Population);
    
    % Dynamic adjustment
    avgCost = 0;
    maxCost = sum(Population(1).object) / 2;
    for popIndex = 1 : Config.popSize
              avgCost = avgCost + sum(Population(popIndex).object) / 2;
              if avgCost > maxCost
                     maxCost = avgCost;
              end
    end
    avgCost = avgCost / Config.popSize;
    record(genIndex) = avgCost;
    
%     selectIndex = ceil( Config.popSize * rand );
%     selectCost = sum(Population(selectIndex).object) / 2;
    
    if genIndex > 1
              if avgCost >= record(genIndex-1)
                     Config.pModify = Config.dyList(1);
                     Config.pMutate = Config.dyList(3);
              
              else
                     Config.pModify = Config.dyList(2);
                     Config.pMutate = Config.dyList(4);
              end
    end
    
    disp([ '>> Generation # ', num2str(genIndex) ]);
    disp([ '-> avgCost ', num2str(avgCost) ]);
    
    
end

Population = cost(Population, Config);
Population = popSort(Population, Config);
ans = conclude(Population, Config);
stop = 1;

return



%% To init NSBBO

function [ cost, feasible, Population ] = NSBBOInit( Problem, Config )

% Initialize population-based optimization software.
if ~exist('RandSeed', 'var')
    RandSeed = round(sum(100*clock));
end
rand('state', RandSeed); % initialize random number generator

disp(['random # seed = ', num2str(RandSeed)]);

% Get the addresses of the initialization, cost, and feasibility functions.
[init, cost, feasible] = Problem();
% Initialize the population.
Population = init(Config);
Population = feasible( Population, Config );
% Make sure the population does not have duplicates. 
Population = clearDups(Population);
Population = feasible( Population, Config );


disp(['Inited Generation # 0 ']);

return



%% To clear duplications
function [Population] = clearDups(Population)

% Make sure there are no duplicate individuals in the population.
% This logic does not make 100% sure that no duplicates exist, but any duplicates that are found are
% randomly mutated, so there should be a good chance that there are no duplicates after this procedure.
for i = 1 : length(Population)
    Chrom1 = sort(Population(i).chrom);
    
    for j = i+1 : length(Population)
        Chrom2 = sort(Population(j).chrom);
        
        if isequal(Chrom1, Chrom2)
            parnum = ceil( ( length(Population(j).chrom) )* rand );
            Population(j).chrom(:,parnum) = rand;
        end
    end
    
end

return



%% To get the ratio
function [lambda, mu] = getRatio(Population)

% Compute immigration rate and extinction rate for each species count.
% lambda(i) is the immigration rate for individual i.
% mu(i) is the extinction rate for individual i.
% This routine assumes the population is sorted from most fit to least fit.

popSize = length(Population);

for popIndex = 1 : popSize
%        if Population.pop(popIndex).distance == Inf
%               Distance = 1;
%        else
%               Distance = Population.pop(popIndex).distance;
%        end
    mu(popIndex) = sin(  (popSize - popIndex) / popSize  );
    lambda(popIndex) = 1 - mu(popIndex);
end

return





%% The strategy to keep elite

function Prior = keepPrior(Population , size)
% to keep Prior

for popIndex = 1 : size
       Prior(popIndex) = Population(popIndex);
       
end

return


function Population = mix(Population, Prior)
% to put elite
popSize = length(Population);

for keepIndex = 1 : length(Prior)
      Population(popSize-keepIndex+1) = Prior(keepIndex);
      
end

return



%%

function conclusion = conclude( Population, Config )


return





