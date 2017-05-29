function [ isthatsame ] = check_action( a,b )
% This function is used to check if two players play the same action.
isthatsame = (a == 1 && b == 6) || (a == 6 && b == 1) ...
                || (a == 2 && b == 5)|| (a == 5 && b == 2) ...
                            || (a == 3 && b == 4)|| (a == 4 && b == 3);

end

