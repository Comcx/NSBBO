
%% Population Sort

function population = popSort( Population, Config )
%POPSORT

popSize = length( Population );

population = nonDominationSort(Population);
population = crowdingDistance(population);

sortList = zeros(popSize, 2);
for popIndex = 1 : popSize
       sortList(popIndex, :) = [population(popIndex).rank, population(popIndex).distance];
end

[none, indices] = sortrows(sortList, [1, -2]);
tempPopulation = population;
for popIndex = 1 : popSize
       population(popIndex) = tempPopulation( indices(popIndex) );
end


end


%% Whether Dominated or not

function res = dominated(pop1, pop2)

flag = 1;
less = 0;

for objectIndex = 1 : length(pop1.object)
       if pop1.object(objectIndex) < pop2.object(objectIndex)
             flag = 0;
       elseif pop1.object(objectIndex) > pop2.object(objectIndex)
             less = 1;
       end
end
res = flag && less;

end



%% Non-Domination Sort 

% This function sort the current popultion based on non-domination. All the 
% individuals in the first front are given a rank of 1, the second front 
% individuals are assigned rank 2 and so on. After assigning the rank the 
% crowding in each front is calculated. 

function Population = nonDominationSort(Population)

popSize = length(Population);
leftCounter = popSize;
rankCounter = 1;
visited = zeros(1, popSize);
N_p = zeros(1, popSize);
S_p = cell(1, popSize);

for i = 1 : popSize
       %visited(i) = 1;
       counter = 1;
       S_p{i}.set(counter) = 0;
       for j = 1 : popSize
              if i~=j && dominated( Population(i), Population(j) )
                     N_p(i) = N_p(i) + 1;
              elseif i~=j && dominated( Population(j), Population(i) )
                     S_p{i}.set(counter) = j;
                     counter = counter + 1;
              end
       end
       counter = 0;
end

% To rank the Population
targetIndex = 0;            % find the target to sort
while( leftCounter > 0 )    % Continue until there's no pop left
       for popIndex = 1 : popSize
             if  N_p(popIndex) == 0 && ~visited(popIndex)
                    visited(popIndex) = 1;
                    targetIndex = popIndex;
                    leftCounter = leftCounter - 1;
                    Population(popIndex).rank = rankCounter;
                    
                    for setIndex = 1 : length( S_p{targetIndex}.set )
                           if S_p{targetIndex}.set(setIndex) ~= 0
                              N_p(S_p{targetIndex}.set(setIndex)) = N_p(S_p{targetIndex}.set(setIndex)) - 1;    
                           end
                    end
             end
       end

       rankCounter = rankCounter + 1;
       
end


end % end nonDominationSort function



%% Crowding-distance Assignment

function population = crowdingDistance(Population)

popSize = length(Population);

for popIndex = 1 : popSize
       Population(popIndex).distance = 0;
end

for objectIndex = 1 : length(Population(1).object)
       population = getSorted(Population, objectIndex);
       f_min = population(1).object(objectIndex);
       f_max = population(popSize).object(objectIndex);
       population(1).distance = Inf;
       population(popSize).distance = Inf;
       
       for popIndex = 2 : popSize-1
              population(popIndex).distance = population(popIndex).distance + ...
                                              abs( population(popIndex-1).object(objectIndex) - population(popIndex+1).object(objectIndex) ) / ...
                                              (f_max - f_min + 1);
       end
end

end


% Para-function to get sorted population by the object [objectIndex]
function population = getSorted(Population, objectIndex)

objectList = zeros( length(Population), length(Population(1).object) );

for popIndex = 1 : length(Population)
       objectList(popIndex, :) = Population(popIndex).object;
end

[none, indices] = sortrows(objectList, objectIndex);
population = Population;

for popIndex = 1 : length(Population)
       population(popIndex) = Population(indices(popIndex));
end

end




