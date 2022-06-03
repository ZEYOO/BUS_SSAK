% Define the step function here to devide the region
% into two pieces 
function step = step(x)
     x(x>0) = 1;
     x(x<=0) = 0;
     step = x;
end