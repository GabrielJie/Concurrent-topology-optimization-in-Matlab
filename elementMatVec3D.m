function Ke = elementMatVec3D(a, b, c, DH)
GN_x=[-1/sqrt(3),1/sqrt(3)]; GN_y=GN_x; GN_z=GN_x; GaussWeigh=[1,1];
Ke = zeros(24,24); L = zeros(6,9);
L(1,1) = 1; L(2,5) = 1; L(3,9) = 1;
L(4,2) = 1; L(4,4) = 1; L(5,6) = 1;
L(5,8) = 1; L(6,3) = 1; L(6,7) = 1;
for i=1:length(GN_x)
    for j=1:length(GN_y)
        for k=1:length(GN_z)
            x = GN_x(i);y = GN_y(j);z = GN_z(k);
            dNx = 1/8*[-(1-y)*(1-z)  (1-y)*(1-z)  (1+y)*(1-z) -(1+y)*(1-z) -(1-y)*(1+z)  (1-y)*(1+z)  (1+y)*(1+z) -(1+y)*(1+z)];
            dNy = 1/8*[-(1-x)*(1-z) -(1+x)*(1-z)  (1+x)*(1-z)  (1-x)*(1-z) -(1-x)*(1+z) -(1+x)*(1+z)  (1+x)*(1+z)  (1-x)*(1+z)];
            dNz = 1/8*[-(1-x)*(1-y) -(1+x)*(1-y) -(1+x)*(1+y) -(1-x)*(1+y)  (1-x)*(1-y)  (1+x)*(1-y)  (1+x)*(1+y)  (1-x)*(1+y)];
            J = [dNx;dNy;dNz]*[ -a  a  a  -a  -a  a  a  -a ;  -b  -b  b  b  -b  -b  b  b; -c -c -c -c  c  c  c  c]';
            G = [inv(J) zeros(3) zeros(3);zeros(3) inv(J) zeros(3);zeros(3) zeros(3) inv(J)];
            dN(1,1:3:24) = dNx; dN(2,1:3:24) = dNy; dN(3,1:3:24) = dNz;
            dN(4,2:3:24) = dNx; dN(5,2:3:24) = dNy; dN(6,2:3:24) = dNz;
            dN(7,3:3:24) = dNx; dN(8,3:3:24) = dNy; dN(9,3:3:24) = dNz;
            Be = L*G*dN;
            Ke = Ke + GaussWeigh(i)*GaussWeigh(j)*GaussWeigh(k)*det(J)*(Be'*DH*Be);
        end
    end
end
end
%======================================================================================================================%
% Subfunction elementMatVec3D:                                                                                         %
% A compact and efficient MATLAB code to evaluate the 3D isoparametric element stiffness matrix                        %
%                                                                                                                      %
% Developed by: Jie Gao, Zhen Luo, Liang Xia and Liang Gao*                                                            %
% Email: gaoliang@mail.hust.edu.cn (GabrielJie_Tian@163.com)                                                           %
%                                                                                                                      %
% Main references:                                                                                                     %
%                                                                                                                      %
% (1) Jie Gao, Zhen Luo, Liang Xia, Liang Gao. Concurrent topology optimization of multiscale composite structures     %
% in Matlab. Accepted in Structural and multidisciplinary optimization.                                                %
%                                                                                                                      %
% (2) Xia L, Breitkopf P. Design of materials using topology optimization and energy-based homogenization approach in  %
% Matlab. % Structural and multidisciplinary optimization, 2015, 52(6): 1229-1241.                                     %
%                                                                                                                      %
% *********************************************   Disclaimer   ******************************************************* %
% The authors reserve all rights for the programs. The programs may be distributed and used for academic and           %
% educational purposes. The authors do not guarantee that the code is free from errors,and they shall not be liable    %
% in any event caused by the use of the program.                                                                       %
%======================================================================================================================%