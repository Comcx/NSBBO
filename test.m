function res = test(Problem)
%TEST to test a small unit
%


for time = 1 : 3
       
       disp(['->> Time # ', num2str(time)]);

global skelePath;
global freeSet;
skelePath = PRM(2);
[m, n] = size(skelePath);

skeleDist = 0;
for pathIndex = 2 : m
       skeleDist = skeleDist + ...
              sqrt( (skelePath(pathIndex,1)-skelePath(pathIndex-1,1)).^2 + ...
                    (skelePath(pathIndex,2)-skelePath(pathIndex-1,2)).^2 );
end

skelePath = skelePath(2:m-1, :);
freeSet = getFreeSpace(skelePath);
skeleDist

pathSet(time) = NSBBO(Problem);

end


end

