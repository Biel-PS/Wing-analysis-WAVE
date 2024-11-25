%% BEAM STRUCTURE: WINGBOX
% BIEL PUJADAS SURIOL / JUDITH BAILEN LOZANO

clear
close all

%% 1) PREPROCESS

% 1.1 Input data (define your input parameters here)
data.ni = 2;  % Degrees of freedom per node

% Material properties matrix
m = [% Each column corresponds to a material property (thickness, young modulus, shear modulus...)
    22e-3   210e9   80e9;
    15e-3   210e9   80e9;
    3.5e-3  210e9   80e9;
];

% 1.2 Build geometry (mesh)
data.c = 2; % chord length

h1 = 0.25*data.c; h2 = 0.15*data.c;
d = 0.3*data.c;

Section_area = 2*m(3,1)*sqrt(d^2+((h1-h2)/2)^2) + m(1,1)*h1 + m(2,1)*h2;

nnodes = 1000;
open = 0;


% Nodal coordinates matrix, Nodal connectivities matrix, Material connectivities matrix
[x, Tn, Tm] = node_pos(h1,h2,d,nnodes,open);

data.nnod = size(x,1); % Number of nodes 
data.nd = size(x,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom
data.nel = size(Tn,1); % Number of elements 
data.nne = size(Tn,2); % Number of nodes in a bar

% Create degrees of freedom connectivities matrix
Td = connectDOF(data,Tn);


% 1.3 Input boundary conditions
% Fixed nodes matrix
p = [% Each row is a prescribed degree of freedom | column_1 = node, column_2 = direction, column_3 = value of prescribed displacement
    1   1   0;
    1   2   0;
    1   3   0
];

%% 2) CROSS-SECTION ANALYSIS
% Section Properties
[Xo,Yo,Xs,Ys,Atot,Ixx,Iyy,Ixy,J,Ain] = SectionProperties(data,x,Tn,m,Tm,open);

I = [Ixx Ixy; Ixy Iyy];
% Normal stress distribution
Mx_p = -1;
My_p = 0;
[sigma,s] = normalStress(data,x,Tn,Xo,Yo,Ixx,Iyy,Ixy,Mx_p,My_p);

figure(1)
plot(s(1,:),sigma(1,:)); grid on; xlabel('Longitud de arco (m)'); ylabel('Tensión (Pa)');
title('Distribución de la tensión normal a lo largo de la sección');

% Tangential stress distribution (Shear)
Sx_p = 0; Sy_p = 1;
[tau_s,s] = TangentialShear(data,x,Tn,m,Tm,Xo,Yo,Xs,Ys,Ain,Ixx,Iyy,Ixy,Sx_p,Sy_p,open);

figure(2)
plot(s(1,:),tau_s(1,:)); grid on; xlabel('Longitud de arco (m)'); ylabel('Tensión (Pa)');
title('Distribución de la tensión por cortante a lo largo de la sección');

% Tangential stress distribution (Torsion)
Mz_p = 1;
[tau_t,s] = TangentialTorsion(data,x,Tn,m,Tm,Mz_p,J,Ain,open);

figure(3)
plot(s(1,:),tau_t(1,:)); grid on; xlabel('Longitud de arco (m)'); ylabel('Tensión (Pa)');
title('Distribución de la tensión por torsion a lo largo de la sección');

%% 3) BEAM ANALYSIS

b = 16; 
be = 0.25*b; 
zm = 0.48*data.c; 
za = 0.25*data.c; 
ze = 0.3*data.c; 
v_inf = 750*10/36; 
rho = 1.225; 
Cl = 0.1; 
W =  140*b*9.81; 
We = 2100*9.81;
nnodes = 513; 
xi_S = d + 0.3*data.c - Xs; % Distància del LE al SC

[fe,me,x_nodes_biga] = Element_function (be,b,ze,za,zm,data,v_inf,rho,Cl,nnodes,W,We,xi_S);

[Tn_b,Tm_b] = Tn_biga(nnodes); % Càlcul Tn biga

data_b.ni = 3;
data_b.nnod = size(x_nodes_biga,1); % Number of nodes 
data_b.nd = size(x,2);   % Problem dimension
data_b.ndof = data_b.nnod*data_b.ni;  % Total number of degrees of freedom
data_b.nel = size(Tn_b,1); % Number of elements 
data_b.nne = size(Tn_b,2); % Number of nodes in a bar

Td_b = connectDOF(data_b,Tn_b);

m_b = [210e9 Section_area 80e9];

% 2.1.1 Compute element stiffness matrices
Kel = StiffnessFunction(data_b,x_nodes_biga,Tn_b,m_b,Tm_b,J,Ixx);

% 2.1.2 Compute element force vectors
fel = ForceFunction(data_b,x_nodes_biga,Tn_b,fe,me);

% 2.2 Assemble global stiffness matrix
[K,f] = assemblyFunction(data_b,Td_b,Kel,fel);

% 2.3.1 Apply prescribed DOFs
[up,vp] = applyBC(data_b,p);

% 2.4 Solve system
[u,r] = solveSystem(data_b,K,f,up,vp);

% 2.5 Compute stress
sig = stressFunction(data_b,x_nodes_biga,Tn_b,m_b,Tm_b,Td_b,u);

% 2.6 Compute internal forces
[xel,Sel,Mbel,Mtel] = InternalForces(data_b,x_nodes_biga,Tn_b,Td_b,Kel,u);

% Figures
% figure(4)
% plot(x_nodes_biga(1:end-1,1),Sel);
% title('Fuerza de cortante a lo largo de la envergadura'); grid on;
% xlabel('Envergadura (m)'); ylabel('Fuerza de Cortante (N)');
% 
% figure(5)
% plot(x_nodes_biga(1:end-1,1),Mbel);
% title('Momento flector a lo largo de la envergadura'); grid on;
% xlabel('Envergadura (m)'); ylabel('Momento Flector (Nm)');
% 
% figure(6)
% plot(x_nodes_biga(1:end-1,1),Mtel);
% title('Momento torson a lo largo de la envergadura'); grid on;
% xlabel('Envergadura (m)'); ylabel('Fuerza Torsor (Nm)');
% 
% figure(7)
% plot(x_nodes_biga(:,1),u(1:3:end));
% title('Deflexión vertical a lo largo de la envergadura'); grid on;
% xlabel('Envergadura (m)'); ylabel('Deflexión (m)');
% 
% figure(8)
% plot(x_nodes_biga(:,1),u(2:3:end));
% title('Rotación por flexión de la sección'); grid on;
% xlabel('Envergadura (m)'); ylabel('Rotación (m)');
% 
% figure(9)
% plot(x_nodes_biga(:,1),u(3:3:end));
% title('Rotación por torsión de la sección'); grid on;
% xlabel('Envergadura (m)'); ylabel('Rotación (m)');

%% 3.2) CONVERGENCE
%     lim = 10;
%     nnodes_used = zeros(1,lim-1);
%     utip = zeros(3,lim-1);
%     uref = u;
% 
% for i = 2:lim
% 
%     nnodes = 2^i+1;
%     nnodes_used(1,i-1) = nnodes;
%     [fe,me,x_nodes_biga] = Element_function (be,b,ze,za,zm,data,v_inf,rho,Cl,nnodes,W,We,xi_S);
% 
%     [Tn_b,Tm_b] = Tn_biga(nnodes);
% 
%     data_b.ni = 3;
%     data_b.nnod = size(x_nodes_biga,1); % Number of nodes 
%     data_b.nd = size(x_nodes_biga,2);   % Problem dimension
%     data_b.ndof = data_b.nnod*data_b.ni;  % Total number of degrees of freedom
%     data_b.nel = size(Tn_b,1); % Number of elements 
%     data_b.nne = size(Tn_b,2); % Number of nodes in a bar
% 
%     Td_b = connectDOF(data_b,Tn_b);
% 
%     Kel = StiffnessFunction(data_b,x_nodes_biga,Tn_b,m_b,Tm_b,J,Ixx);
% 
%     fel = ForceFunction(data_b,x_nodes_biga,Tn_b,fe,me);
% 
%     [K,f] = assemblyFunction(data_b,Td_b,Kel,fel);
% 
%     [up,vp] = applyBC(data_b,p);
% 
%     [u,r] = solveSystem(data_b,K,f,up,vp);
%     utip (:,i-1) = u(size(u)-2:end,1);
% 
% end
% 
% figure(10)
% plot(nnodes_used,utip(1,:));
% title('Convergencia');
% 
% figure(11)
% plot(x_nodes_biga(1:end-1,1),fe(2,:));
% title('Fuerza de cada elemento');
% 
% figure(12)
% plot(x_nodes_biga(1:end-1,1),me(1,:));
% title('Momento de cada elemento');

%% 4) VON MISES
nnodes = 1000;

% Open
open = 1;
[sigVM_o,posmax_o] = VonMises(open,h1,h2,d,nnodes,data,m,Sel,Mbel,Mtel);

%A = [sigVM_o]

% Closed
open = 0;
[sigVM_c,posmax_c] = VonMises(open,h1,h2,d,nnodes,data,m,Sel,Mbel,Mtel);


