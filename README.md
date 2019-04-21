# The Matlab code for the concurrent topology optimization

%======================================================================================================================%
% Function ConTop2D:                                                                                                   %
% A compact and efficient MATLAB code for Concurrent topology optimization of multiscale composite structures          %
% in Matlab.                                                                                                           %
%                                                                                                                      %
% Developed by: Jie Gao, Zhen Luo, Liang Xia and Liang Gao*                                                            %
% Email: gaoliang@mail.hust.edu.cn (GabrielJie_Tian@163.com)                                                           %
%                                                                                                                      %
% Main references:                                                                                                     %
%                                                                                                                      %
% (1) Jie Gao, Zhen Luo, Liang Xia, Liang Gao. Topology optimization for concurrent design of material and structure   %
% in Matlab. submitted to the journal "Structural and multidisciplinary optimization". 2018                            %
%                                                                                                                      %
% (2) Xia L, Breitkopf P. Design of materials using topology optimization and energy-based homogenization approach in  %
% Matlab. % Structural and multidisciplinary optimization, 2015, 52(6): 1229-1241.                                     %
%                                                                                                                      %
% *********************************************   Disclaimer   ******************************************************* %
% The authors reserve all rights for the programs. The programs may be distributed and used for academic and           %
% educational purposes. The authors do not guarantee that the code is free from errors,and they shall not be liable    %
% in any event caused by the use of the program.                                                                       %
%======================================================================================================================%

including:
Contop2D.m is the main function for two dimensional structures.
Contop3D.m is the mian function for three dimensional structures.
EBHM2D.m is the sub function to evaluate the effective macroscopic properties in 2D.
EBHM3D.m is the sub function to evaluate the effective macroscopic properties in 3D.
elementMatVec2D.m is the sub function to compute the 2D element stiffness matrix.
elementMatVec3D.m is the sub function to compute the 3D element stiffness matrix.
