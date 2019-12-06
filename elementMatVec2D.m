function Ke = elementMatVec2D(a, b, DH)
GaussNodes = [-1/sqrt(3); 1/sqrt(3)]; GaussWeigh = [1 1];
L = [1 0 0 0; 0 0 0 1; 0 1 1 0];
Ke = zeros(8,8);
for i = 1:2
    for j = 1:2
        GN_x = GaussNodes(i); GN_y = GaussNodes(j);
        dN_x = 1/4*[-(1-GN_x)  (1-GN_x) (1+GN_x) -(1+GN_x)];
        dN_y = 1/4*[-(1-GN_y) -(1+GN_y) (1+GN_y)  (1-GN_y)];
        J = [dN_x; dN_y]*[ -a  a  a  -a;  -b  -b  b  b]';
        G = [inv(J) zeros(size(J)); zeros(size(J)) inv(J)];
        dN(1,1:2:8) = dN_x; dN(2,1:2:8) = dN_y;
        dN(3,2:2:8) = dN_x; dN(4,2:2:8) = dN_y;
        Be = L*G*dN;
        Ke = Ke + GaussWeigh(i)*GaussWeigh(j)*det(J)*Be'*DH*Be;
    end
end
end
%======================================================================================================================%
% Subfunction elementMatVec2D:                                                                                         %
% A compact and efficient MATLAB code to evaluate the 2D isoparametric element stiffness matrix                        %
%                                                                                                                      %
% Developed by: Erik Andreassen and Casper Schousboe Andreasen                                                         %
% Email: erand@mek.dtu.dk (E.Andreassen).                                                                              %
%                                                                                                                      %
% Main references:                                                                                                     %
%                                                                                                                      %
% (1) Andreassen E, Andreasen C S. How to determine composite material properties using numerical homogenization.      %
% Computational Materials Science, 2014, 83: 488-495.                                                                  %
%                                                                                                                      %
% *********************************************   Disclaimer   ******************************************************* %
% The authors reserve all rights for the programs. The programs may be distributed and used for academic and           %
% educational purposes. The authors do not guarantee that the code is free from errors,and they shall not be liable    %
% in any event caused by the use of the program.                                                                       %
%======================================================================================================================%