function [ auxi_var ] = auxiliary_var_selection( vq )
% This function is to solve the selection of auxiliary variable
% Define the maximum allowed rate and the control parameter V
[nSF, nBS] = size(vq); % number of BSs and subflows
auxi_var = zeros(nSF, nBS);
V = 5;
max_rate = log(1 + 10); % maximum rate of each link
min_rate = 0;
for sub_f = 1:nSF
    for bs = 1:nBS
        if vq(sub_f,bs) ~= 0
            auxi_var(sub_f,bs) = V./vq(sub_f,bs);
        else
            auxi_var(sub_f,bs) =  max_rate;
        end
        % Select 
        auxi_var(sub_f,bs) = max(auxi_var(sub_f,bs),min_rate);
        auxi_var(sub_f,bs) = min(auxi_var(sub_f,bs),max_rate);
    end     
end

end

