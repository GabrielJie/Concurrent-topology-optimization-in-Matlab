# The Matlab code for the concurrent topology optimization


This paper presents the compact and efficient Matlab codes for the concurrent topology optimization of multiscale composite structures not only in 2D scenario, but also considering 3D cases. 
A modified SIMP approach (Sigmund 2007) is employed to implement the concurrent topological design, with an energy-based homogenization method (EBHM) to evaluate the macroscopic effective properties of the microstructure. 
The 2D and 3D Matlab codes in the paper are developed, using the 88-line 2D SIMP code (Struct Multidisc Optim 43(1): 1-16, 2011) and the 169-line 3D topology optimization code (Struct Multidisc Optim 50(6): 1175-1196, 2014), respectively. 
This paper mainly contributes to the following four aspects: (1) the code architecture for the topology optimization of cellular composite structures (ConTop2D.m and ConTop3D.m); (2) the code to compute the 3D isoparametric element stiffness matrix (elementMatVec3D.m); (3) the EBHM to predict the macroscopic effective properties of 2D and 3D material microstructures (EBHM2D.m and EBHM3D.m); and (4) the code to calculate the sensitivities of the objective function with respect to the design variables at two scales.


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
% in Matlab. submitted to the journal "Structural and multidisciplinary optimization". 2019                            %
% https://doi.org/10.1007/s00158-019-02323-6                                                                           %
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
filtering2d.m is the 2D density filtering function for topology optimization
filtering3d.m is the 3D density filtering function for topology optimization
OC.m is the Optimality Criteria function to evolve design variables.
display_3D.m is the function to present the optimized topologies in 3D.




