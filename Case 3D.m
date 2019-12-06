%% load 1
Macro_struct = [6, 2, 0.4, 60, 20, 4, 0.25];
Micro_struct = [0.1, 0.1, 0.1, 20, 20, 20, 0.3];
Load_CASE = 1; penal = 3; rmin = 1.5; ft = 3;
ConTop3D(Macro_struct, Micro_struct, Load_CASE, penal, rmin, ft)
%% load 2
Macro_struct = [2, 1.6, 2, 20, 16, 20, 0.2];
Micro_struct = [0.1, 0.1, 0.1, 20, 20, 20, 0.3];
Load_CASE = 2; penal = 3; rmin = 2.0; ft = 3;
ConTop3D(Macro_struct, Micro_struct, Load_CASE, penal, rmin, ft)
Macro_struct = [2, 1.6, 2, 20, 16, 20, 0.2];
Micro_struct = [0.1, 0.1, 0.1, 20, 20, 20, 0.3];
Load_CASE = 2; penal = 4; rmin = 1.5; ft = 3;
ConTop3D(Macro_struct, Micro_struct, Load_CASE, penal, rmin, ft)
Macro_struct = [2, 1.6, 2, 20, 16, 20, 0.2];
Micro_struct = [0.1, 0.1, 0.1, 20, 20, 20, 0.3];
Load_CASE = 2; penal = 4; rmin = 2.0; ft = 3;
ConTop3D(Macro_struct, Micro_struct, Load_CASE, penal, rmin, ft)