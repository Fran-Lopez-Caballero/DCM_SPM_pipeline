% This function generates the model space of reduced models 

function [output] = gen_model_space_NMDA()

% Models for repetition effects (none, extr, intr, both)
%--------------------------------------------------------------------------
F =   [ 0 0 0 0 0 0; ...
     	0 0 0 0 0 0; ...
     	1 0 0 0 0 0; ...
        0 1 0 0 0 0; ...
      	0 0 1 0 0 0; ...
      	0 0 0 1 0 0];  
    
B =   [ 0 0 1 0 0 0; ...
     	0 0 0 1 0 0; ...
     	0 0 0 0 1 0; ...
        0 0 0 0 0 1; ...
      	0 0 0 0 0 0; ...
      	0 0 0 0 0 0];   
    
output{1}.name = '0';
output{1}.matrix = zeros(6);

output{2}.name = 'F';
output{2}.matrix = F;

output{3}.name = 'B';
output{3}.matrix = B;

output{4}.name = 'FB';
output{4}.matrix = F + B;
