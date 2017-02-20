function [ nq ] = network_queue_update( selectedaction, pre_nq, outgoing_data, incoming_data )
global alpha1 % fraction of data, user 1 divides for subflow1

[nSF, nBS] = size(pre_nq); % number of BSs and subflows
nq = zeros(nSF, nBS);

indicator_BSs = routingtable (selectedaction);
    % Only update the arrival rate for MBS b=1
%     bs = 1;
%         nq(1, bs) = max(pre_nq(1, bs) - indicator_BSs(bs,1) * ( outgoing_data(1, bs) ), 0) ...
%             +  indicator_BSs(bs,1)* ( incoming_data(bs, 1) );
%         nq(2, bs) = max(pre_nq(2, bs) - indicator_BSs(bs,2) * ( outgoing_data(2, bs) ), 0) + alpha1;
%         nq(3, bs) = max(pre_nq(3, bs) - indicator_BSs(bs,3) * ( outgoing_data(3, bs) ), 0) + alpha1;
%         nq(4, bs) = max(pre_nq(4, bs) - indicator_BSs(bs,4) * ( outgoing_data(4, bs) ), 0) + alpha1; 
    for bs = 1:nBS
        nq(1, bs) = max(pre_nq(1, bs) - indicator_BSs(bs,1)* ( outgoing_data(1, bs) ), 0) ...
            +  indicator_BSs(bs,1)* ( incoming_data(bs, 1) );
        nq(2, bs) = max(pre_nq(2, bs) - indicator_BSs(bs,2)* ( outgoing_data(2, bs) ), 0) ...
            + indicator_BSs(bs,2)* ( incoming_data(bs, 2) );
        nq(3, bs) = max(pre_nq(3, bs) - indicator_BSs(bs,3)* ( outgoing_data(3, bs) ), 0) ...
            + indicator_BSs(bs,3)* ( incoming_data(bs, 3) );
        nq(4, bs) = max(pre_nq(4, bs) - indicator_BSs(bs,4)* ( outgoing_data(4, bs) ), 0) ...
            + indicator_BSs(bs,4)* ( incoming_data(bs, 4) );
    end

end

