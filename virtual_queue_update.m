function [ vq ] = virtual_queue_update( pre_vq,ax_var, rate )
    % Update the value of virtual queue
    vq = max(pre_vq + ax_var - rate, 0);
end

