function [ ] = loadsolvers( )
% This function is used to setup the yalmip tool for optimization model,
% and select the solver 
fprintf('\nAdding solvers\n');
% % For the desktop
addpath K:/yalmip;
addpath K:/yalmip/extras;
addpath K:/yalmip/demos;
addpath K:/yalmip/solvers;
addpath K:/yalmip/modules;
addpath K:/yalmip/modules/parametric;
addpath K:/yalmip/modules/moment;
addpath K:/yalmip/modules/global;
addpath K:/yalmip/modules/sos;
addpath K:/yalmip/operators;
addpath K:/yalmip/solvers/sedumi;
install_sedumi;
addpath K:/yalmip/solvers/sdpt3;
install_sdpt3;

addpath 'C:/Program Files/Mosek/7';
addpath 'C:/Program Files/Mosek/7/toolbox/r2013a';


% % % Note that for Server
% 
%  addpath /home/vkien/cvx/structures;
%  addpath /home/vkien/cvx/lib;
%  addpath /home/vkien/cvx/functions;
%  addpath /home/vkien/cvx/functions/vec;
%  addpath /home/vkien/cvx/commands;
%  addpath /home/vkien/cvx/builtins; 
%  addpath /home/vkien/cvx;
%  cvx_setup;
% 
%  addpath /home/vkien/yalmip;
%  addpath /home/vkien/yalmip/extras;
%  addpath /home/vkien/yalmip/demos;
%  addpath /home/vkien/yalmip/solvers;
%  addpath /home/vkien/yalmip/modules;
%  addpath /home/vkien/yalmip/modules/parametric;
%  addpath /home/vkien/yalmip/modules/moment;
%  addpath /home/vkien/yalmip/modules/global;
%  addpath /home/vkien/yalmip/modules/sos;
%  addpath /home/vkien/yalmip/operators;
%  addpath /home/vkien/yalmip/solvers/sedumi;
%  install_sedumi;
%  addpath /home/vkien/yalmip/solvers/sdpt3;
%  install_sdpt3;
%  % for Mosek solver
%  
%  addpath /home/vkien/mosek;
%  addpath /home/vkien/mosek/7;
%  addpath /home/vkien/mosek/7/toolbox/r2013a;

end

