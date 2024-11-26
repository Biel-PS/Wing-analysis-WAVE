%% BEAM STRUCTURE: WINGBOX
% BIEL PUJADAS SURIOL / JUDITH BAILEN LOZANO

clear
close all

%% 1) PREPROCESS

% 1.1 Input data (define your input parameters here)
data.ni = 3;  % Degrees of freedom per node

% Material properties matrix
m = [% Each column corresponds to a material property (thickness, young modulus, shear modulus...)
    22e-3   210e9   80e9;
    15e-3   210e9   80e9;
    3.5e-3  210e9   80e9;
];

% 1.2 Build geometry (mesh)
data.c = 2.08; % chord length

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

[Xo,Yo,Xs,Ys,Atot,Ixx,Iyy,Ixy,J,Ain] = SectionProperties(data,x,Tn,m,Tm,open);

I = [Ixx Ixy; Ixy Iyy];
% Normal stress distribution
Mx_p = -1;
My_p = 0;
[sigma,s] = normalStress(data,x,Tn,Xo,Yo,Ixx,Iyy,Ixy,Mx_p,My_p);

% figure(1)
% plot(s(1,:),sigma(1,:)); grid on; xlabel('Longitud de arco (m)'); ylabel('Tensión (Pa)');
% title('Distribución de la tensión normal a lo largo de la sección');

% Tangential stress distribution (Shear)
Sx_p = 0; Sy_p = 1;
[tau_s,s] = TangentialShear(data,x,Tn,m,Tm,Xo,Yo,Xs,Ys,Ain,Ixx,Iyy,Ixy,Sx_p,Sy_p,open);

% figure(2)
% plot(s(1,:),tau_s(1,:)); grid on; xlabel('Longitud de arco (m)'); ylabel('Tensión (Pa)');
% title('Distribución de la tensión por cortante a lo largo de la sección');

% Tangential stress distribution (Torsion)
Mz_p = 1;
[tau_t,s] = TangentialTorsion(data,x,Tn,m,Tm,Mz_p,J,Ain,open);

% figure(3)
% plot(s(1,:),tau_t(1,:)); grid on; xlabel('Longitud de arco (m)'); ylabel('Tensión (Pa)');
% title('Distribución de la tensión por torsion a lo largo de la sección');

%% 3) BEAM ANALYSIS

b = 9.115; 
be = 0.25*b; 
zm = 0.48*data.c; 
za = 0.25*data.c; 
ze = 0.3*data.c; 
v_inf = 152.83*10/36; 
rho = 1.20; 
Cl = 2.263; 
W =  (9500/2)*9.81; 
We = 0*9.81;
nnodes = 513; 
nnode_mot = round(nnodes/b * be);
xi_S = d + 0.3*data.c - Xs; % Distància del LE al SC

F = [nnode_mot 1 -We;
    nnode_mot 3 -We*(xi_S-ze)];

[fe,me,x_nodes_biga] = Element_function (be,b,ze,za,zm,data,v_inf,rho,Cl,nnodes,W,We,xi_S);

[Tn_b,Tm_b] = Tn_biga(nnodes); % Càlcul Tn biga

data_b.ni = 3;
data_b.nnod = size(x_nodes_biga,1); % Number of nodes 
data_b.nd = size(x,2);   % Problem dimension
data_b.ndof = data_b.nnod*data_b.ni;  % Total number of degrees of freedom
data_b.nel = size(Tn_b,1); % Number of elements 
data_b.nne = size(Tn_b,2); % Number of nodes in a bar

Td_b = connectDOF(data_b,Tn_b);

m_b = [210e9 Atot 80e9];

% 2.1.1 Compute element stiffness matrices
Kel = StiffnessFunction(data_b,x_nodes_biga,Tn_b,m_b,Tm_b,J,Ixx);

% 2.1.2 Compute element force vectors
fel = ForceFunction(data_b,x_nodes_biga,Tn_b,fe,me);

% 2.2 Assemble global stiffness matrix
[K,f] = assemblyFunction(data_b,Td_b,Kel,fel);

% 2.3.1 Apply prescribed DOFs
[up,vp] = applyBC(data_b,p);

f = pointLoads(data_b,Td_b,f,F);
% 2.4 Solve system
[u,r] = solveSystem(data_b,K,f,up,vp);

u_tip_ref = u(end-2:end,1);

% 2.5 Compute internal forces
[xel,Sel,Mbel,Mtel] = InternalForces(data_b,x_nodes_biga,Tn_b,Td_b,Kel,u);

%load open.mat % generar el codi per secció oberta i guardar les dades de u, sel, mbel, mtel en aquest arxiu per posteriorment, comparar
% % Figures
% 
% Mtelopen = Mtel;
% Mbelopen = Mbel;
% Selopen = Sel;
% uopen = u;
% 
% save open.mat Mtelopen Selopen Mbelopen uopen;

%load open.mat

figure(4)
plot(x_nodes_biga(1:end-1,1),Sel);
title('Fuerza de cortante a lo largo de la envergadura'); grid on;
xlabel('Envergadura (m)'); ylabel('Fuerza de Cortante (N)');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

figure(5)
plot(x_nodes_biga(1:end-1,1),Mbel);
title('Momento flector a lo largo de la envergadura'); grid on;
xlabel('Envergadura (m)'); ylabel('Momento Flector (Nm)');


figure(6)
plot(x_nodes_biga(1:end-1,1),Mtel);
title('Momento torsor a lo largo de la envergadura'); grid on;
xlabel('Envergadura (m)'); ylabel('Momento Torsor (Nm)');

figure(7)
plot(x_nodes_biga(:,1),u(1:3:end));
title('Deflexión vertical a lo largo de la envergadura'); grid on;
xlabel('Envergadura (m)'); ylabel('Deflexión (m)');


figure(8)
plot(x_nodes_biga(:,1),u(2:3:end));
title('Rotación por flexión de la sección'); grid on;
xlabel('Envergadura (m)'); ylabel('Rotación (rad)');


figure(9)
% plot(x_nodes_biga(:,1),u(3:3:end));
% title('Rotación por torsión de la sección'); grid on;
xlabel('Envergadura (m)'); ylabel('Rotación (rad)');
% legend ("Sección cerrada", "Sección abierta");

% Top plot
nexttile;
plot(x_nodes_biga(:,1),u(3:3:end));
title('Rotación por torsión perfil cerrado');
grid('on');




%% 3.2) CONVERGENCE
    lim = 6;
    nnodes_used = zeros(1,lim-1);
    utip = zeros(3,lim-1);
    uref = u;

for i = 2:lim

    nnodes = 2^i+1;
    nnodes_used(1,i-1) = nnodes;
    nnode_mot = round(nnodes/b * be);

    F = [nnode_mot 1 -We;
         nnode_mot 3 -We*(xi_S-ze)];

    [fe,me,x_nodes_biga] = Element_function (be,b,ze,za,zm,data,v_inf,rho,Cl,nnodes,W,We,xi_S);

    [Tn_b,Tm_b] = Tn_biga(nnodes);
    
    data_b.ni = 3;
    data_b.nnod = size(x_nodes_biga,1); % Number of nodes 
    data_b.nd = size(x_nodes_biga,2);   % Problem dimension
    data_b.ndof = data_b.nnod*data_b.ni;  % Total number of degrees of freedom
    data_b.nel = size(Tn_b,1); % Number of elements 
    data_b.nne = size(Tn_b,2); % Number of nodes in a bar
    
    Td_b = connectDOF(data_b,Tn_b);
    
    Kel = StiffnessFunction(data_b,x_nodes_biga,Tn_b,m_b,Tm_b,J,Ixx);
    
    fel = ForceFunction(data_b,x_nodes_biga,Tn_b,fe,me);

    [K,f] = assemblyFunction(data_b,Td_b,Kel,fel);
    
    [up,vp] = applyBC(data_b,p);

    f = pointLoads(data_b,Td_b,f,F);

    [u,r] = solveSystem(data_b,K,f,up,vp);
    utip (:,i-1) = u(size(u)-2:end,1);

    [xel,Sel,Mbel,Mtel] = InternalForces(data_b,x_nodes_biga,Tn_b,Td_b,Kel,u);

end

figure(10)
plot(nnodes_used,abs((utip(1,:)-u_tip_ref(1))/u_tip_ref(1)),nnodes_used,abs((utip(2,:)-u_tip_ref(2))/u_tip_ref(2)),nnodes_used,abs((utip(3,:)-u_tip_ref(3))/u_tip_ref(3)));
title('Convergencia del método con 512 elementos como valor de referencia en ejes logarítmicos');
grid;
xlabel ("Número d'elementos");
ylabel ("Error relativo");
legend ("Desplazamineto vertical", "Rotación de flexión", "Rotación de torsión");
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

figure(13)
plot(nnodes_used,abs((utip(1,:)-u_tip_ref(1))/u_tip_ref(1)),nnodes_used,abs((utip(2,:)-u_tip_ref(2))/u_tip_ref(2)),nnodes_used,abs((utip(3,:)-u_tip_ref(3))/u_tip_ref(3)));
title('Convergencia del método con 512 elementos como valor de referencia');
grid;
xlabel ("Número d'elementos");
ylabel ("Error relativo");
legend ("Desplazamineto vertical", "Rotación de flexión", "Rotación de torsión");

% figure(11)
% plot(x_nodes_biga(1:end-1,1),fe(2,:));
% title('Fuerza de cada elemento');
% 
% figure(12)
% plot(x_nodes_biga(1:end-1,1),me(1,:));
% title('Momento de cada elemento');

%% Von Mises

nnodes = 1000;

% Open
open = 1;
[sigVM_o,pos_o] = VonMises(open,h1,h2,d,nnodes,data,m,Sel,Mbel,Mtel);

[sigmax_o, index_o] = max(sigVM_o);
posmax_o = pos_o(:,index_o);

% Closed
open = 0;
[sigVM_c,pos_c] = VonMises(open,h1,h2,d,nnodes,data,m,Sel,Mbel,Mtel);

[sigmax_c, index_c] = max(sigVM_c);
posmax_c = pos_c(:,index_c);