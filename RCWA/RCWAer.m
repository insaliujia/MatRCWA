function [R,T] = RCWAer(ObjRCWA,source,device,wavenumber)
% This function caculate the reflection and transmission 
% Inside this function there is two sub-functions 
% one named Dev_build construct the convolution matrices
% one named Global_marix caculate the device matrices and connect it with 
% ouside region to construct global matrices
%
% INPUT ARGUMENTS
% Ref contains the parameters in the reflective region 
% Trn contains the parameters in the tranmission regin 
% Source contains the parameters of the source
% Dev contains all the parameters of the device
% PQR is number of spatial harmonics along x and y (HAS TO BE ODD NUM)
% Nx is the number of point along x in real-space grid Ny will be caculated
% by using Nx
%
% OUTPUT ARGUMENTS
%
% R is the reflection 
% T is the transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERIFY INPUT/OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VERIFY NUMBER OF INPUT ARGUMENTS
error(nargchk(4,4,nargin));

% VERIFY NUMBER OF OUTPUT ARGUMENTS
error(nargchk(2,2,nargout));

%% EXTRACT PARAMETERS
% Reflective and reflective region parameters
Ref.ur1 = ObjRCWA.referur(2);                         %permeability in reflection region 
Ref.er1 = ObjRCWA.referur(1);                         %permittivity in reflection region 
Trn.ur2 = ObjRCWA.trnerur(2);                         %permeability in transmission region 
Trn.er2 = ObjRCWA.trnerur(1);                         %permittivity in transmission region
URC = device.URC;                                     %convolution matrix of ur
ERC = device.ERC;                                     %convolution matrix of er

PQR = device.PQR;

CondRecordField = ObjRCWA.RecordField;                %condiction of recording the field
%% CACULATE DEVICE MATRICES 
% EXTRACT PARAMETERS
NH  = PQR(1)*PQR(2);                   %total number of spatial harmonics
L = device.GetLengthLayer;                 %lenght of each layer of the devices
k0 = source.sk(wavenumber);                     %caculate source vector which will be caculated 


%Compute wave vector expansion

I0 = eye(NH);
Z0 = zeros(NH);
Ref.n = sqrt(Ref.er1*Ref.ur1);                      % reflective index in the reflection region
Trn.n = sqrt(Trn.er2*Trn.ur2);                      % reflective index in the transmission region

% Compute Source Vector
kinc = source.skinc(Ref.n); 

%Compute the spatial harmonics to compute wave vector components
m = [-floor(PQR(1)/2):floor(PQR(1)/2)]';
n = [-floor(PQR(2)/2):floor(PQR(2)/2)]';

%Compute wave vector components 
kx = kinc(1) - 2*pi*m/(k0*device.xydimension(1));
ky = kinc(2) - 2*pi*n/(k0*device.xydimension(2));


%Defensive code for the case when lam equals period
kx(kx==0)=0.0000001;
ky(ky==0)=0.0000001;



[ky,kx] = meshgrid(ky,kx);
% Kx = diag(kx(:));
% Ky = diag(ky(:));
Kx = diag(sparse(kx(:)));
Ky = diag(sparse(ky(:)));

%Compute the kz componets
KzR = -conj(sqrt((Ref.er1*Ref.ur1)*I0 - Kx.^2 - Ky.^2));
KzT = conj(sqrt((Trn.er2*Trn.ur2)*I0 - Kx.^2 - Ky.^2));

%Caculate the P and Q in free space
Kz0 = conj(sqrt(I0 - Kx*Kx - Ky*Ky));
Q0 = [Kx*Ky, I0-Kx*Kx; Ky*Ky-I0, -Kx*Ky];                 %eigen-modes for the electric fields in free space
% P0 = Q0;
% OMEGA0 = P0*Q0;
% [W0,LAM0] = eig(OMEGA0);                              % compute eigen-modes 
% LAM0 = sqrt(LAM0);
% V0 = Q0*W0/LAM0;                                        % compute the V concerning magnetic  part 
W0 = [I0 Z0;Z0 I0];
LAM0 = [1i*Kz0 Z0; Z0 1i*Kz0];
V0 = Q0/LAM0;                                   %eigen-modes for the magnetic fields in free space

%defaut W and V in the reflective and transmission region will be in the
%free space
W_ref = W0;
W_trn = W0;
Kz_ref = -Kz0;
Kz_trn = Kz0;
V_ref = V0;
V_trn = V0;

% Initialize Global Scattering Matix
I = eye(NH*2);
Z = zeros(NH*2);
S_G_11 = Z;
S_G_12 = I;
S_G_22 = Z;
S_G_21 = I;


% Main loop for each layer 
for n = 1: sum(device.ilayer) 
    % Caculate Parameters for each Layer 
    P = [Kx/ERC(:,:,n)*Ky, URC(:,:,n)-Kx/ERC(:,:,n)*Kx ;...
            Ky/ERC(:,:,n)*Ky-URC(:,:,n), -Ky/ERC(:,:,n)*Kx];     % matrix concerning magnetic  part 
    Q = [Kx/URC(:,:,n)*Ky, ERC(:,:,n)-Kx/URC(:,:,n)*Kx ;...
            Ky/URC(:,:,n)*Ky-ERC(:,:,n), -Ky/URC(:,:,n)*Kx];     % matrix concerning electrical part 
    OMEGA2 = P*Q;
    [W,LAM] = eig(OMEGA2);                              % compute eigen-modes 
    LAM = sqrt(LAM);
    V = Q*W/LAM;                                        % compute the V concerning magnetic  part 
    
    % Record the parameters for reconstructe E and H field
    if CondRecordField == 1
        ObjRCWA.field.W_V{1,n} = [W,W;-V,V];
        ObjRCWA.field.LAM{1,n} = sparse(LAM);
    end
    
    % Calculate Scattering Matrix for each layer
    A = W\W0 + V\V0;
    B = W\W0 - V\V0;
%     X = expm(-LAM*k0*L(n));
    X = diag(exp(-diag(LAM)*k0*L(n)));
    S11 = (A-X*B/A*X*B)\(X*B/A*X*A-B); % caculate reflection matrix
    S12 = (A-X*B/A*X*B)\X*(A-B/A*B);      % caculate transmission matrix 
    S22 = S11;
    S21 = S12;
    
    % Update Device Scattering Matrix 
    D = S_G_12/(I-S11*S_G_22);
    F = S21/(I-S_G_22*S11);
    S_G_11 = S_G_11 + D*S11*S_G_21;
    S_G_12 = D*S12;
    S_G_21 = F*S_G_21;
    S_G_22 = S22 + F*S_G_22*S12;
end

%% Connect the device matrix to external regions (This part is not necessary if ref and trn parts is in free space)
if sum(ObjRCWA.referur ~= [1,1]) + sum(ObjRCWA.trnerur ~= [1,1]) ~= 0
    % Caculate reflection-side scattering matrix
    % by defaut the external regions will be free space
    Q_ref = 1/Ref.ur1*[Kx*Ky, Ref.ur1*Ref.er1*I0-Kx*Kx;...
        Ky*Ky-Ref.ur1*Ref.er1*I0, -Kx*Ky];     % matrix concerning electrical part
    Kz_ref = -conj(sqrt(Ref.ur1*Ref.er1*I0 - Kx*Kx - Ky*Ky));
    W_ref = [I0 Z0;Z0 I0];
    LAM_ref = [-1i*Kz_ref Z0; Z0 -1i*Kz_ref];
    V_ref = Q_ref/LAM_ref;                                   %eigen-modes for the magnetic fields in free space
    % Calculate Scattering Matrix for reflection region
    A_ref = W0\W_ref + V0\V_ref;
    B_ref = W0\W_ref - V0\V_ref;
    S11_ref = -A_ref\B_ref;                                         % caculate reflection matrix in the reflection region
    S12_ref = 2*eye(2*NH)/A_ref;                                              % caculate transmission matrix in the reflection region
    S21_ref = 1/2*(A_ref-B_ref/A_ref*B_ref);
    S22_ref = B_ref/A_ref;
    
    % Caculate transmission-side scattering matrix
    Q_trn = 1/Trn.ur2*[Kx*Ky, Trn.er2*Trn.ur2*I0-Kx*Kx;...
        Ky*Ky-Trn.er2*Trn.ur2*I0, -Kx*Ky];     % matrix concerning electrical part
    Kz_trn = conj(sqrt(Trn.ur2*Trn.er2*I0 - Kx*Kx-Ky*Ky));
    W_trn = [I0 Z0;Z0 I0];
    LAM_trn = [1i*Kz_trn Z0; Z0 1i*Kz_trn];
    V_trn = Q_trn/LAM_trn;
    
    % Calculate Scattering Matrix for transmission region
    A_trn = W0\W_trn + V0\V_trn;
    B_trn = W0\W_trn - V0\V_trn;
    S11_trn = B_trn/A_trn;                                         % caculate reflection matrix in the reflection region
    S12_trn = 1/2*(A_trn-B_trn/A_trn*B_trn);                                              % caculate transmission matrix in the reflection region
    S21_trn = 2*inv(A_trn);
    S22_trn = -A_trn\B_trn;
    
    % Connect the reflection region scattering matrix and transmission region
    % scattering matrix with the device scattering matrix to build the globle
    % matrix
    % For the reflection region
    D_ref = S12_ref/(I-S_G_11*S22_ref);
    F_ref = S_G_21/(I-S22_ref*S_G_11);
    S_G_22 = S_G_22 + F_ref*S22_ref*S_G_12;
    S_G_21 = F_ref*S21_ref;
    S_G_12 = D_ref*S_G_12;
    S_G_11 = S11_ref + D_ref*S_G_11*S21_ref;
    
    % For the tranmission region
    D_trn = S_G_12/(I-S11_trn*S_G_22);
    F_trn = S21_trn/(I-S_G_22*S11_trn);
    S_G_11 = S_G_11 + D_trn*S11_trn*S_G_21;
    S_G_12 = D_trn*S12_trn;
    S_G_21 = F_trn*S_G_21;
    S_G_22 = S22_trn + F_trn*S_G_22*S12_trn;
end

%% Compute Reflected and Transmitted Fields
% Construct Delta Vector
delta = zeros(NH,1);
delta(ceil(NH/2),1) = 1;

% Compute source feild
sP = source.sP(wavenumber,Ref.n);


% % caculate polarization vector
% Sur_nor = [0; 0; -1];                     % define surface normal
% 
% % caculate vector along polarizations
% if Source.theta == 0
%     a_te = [0;1;0];
% else
%     a_te = cross((k0*kinc),Sur_nor)/norm(cross((k0*kinc),Sur_nor));
% end
% a_tm = cross(a_te,(k0*kinc))/norm(cross(a_te,(k0*kinc)));
% 
% % Composite polarization vector 
% Source.P = Source.TETM(1)*a_te + Source.TETM(2)*a_tm;
% Source.P = Source.P/norm(Source.P);

% Compute source feild
E_src = [sP(1)*delta;sP(2)*delta];

% Compute source modal coefficients
C_src = W_ref\E_src;

% Compute Transmission and Reflection Modal Coefficients
C_ref = S_G_11*C_src;
C_trn = S_G_21*C_src;

% Caculate transmitted and reflected fields
E_ref = W_ref*C_ref;
E_trn = W_trn*C_trn;

% Compute x and y components in the reflection and tranmission part
X_ref = E_ref(1:NH);
Y_ref = E_ref(NH+1:2*NH);
X_trn = E_trn(1:NH);
Y_trn = E_trn(NH+1:2*NH);

% Caculate longitudinal field components
Z_ref = -Kz_ref\(Kx*X_ref + Ky*Y_ref);
Z_trn = -Kz_trn\(Kx*X_trn + Ky*Y_trn);

% Caculate Reflected and Tranmitted power
R = real(-Kz_ref/kinc(3))*(abs(X_ref).^2+abs(Y_ref).^2+abs(Z_ref).^2);
T = real((Ref.ur1/Trn.ur2)*(Kz_trn/kinc(3)))*(abs(X_trn).^2+abs(Y_trn).^2+abs(Z_trn).^2);

%% Record parameters for caculating field in the device

% Record the parameters for reconstructe E and H field
if CondRecordField == 1
    ObjRCWA.field.wavenumber = wavenumber;
    ObjRCWA.field.L = L;
    ObjRCWA.field.xydimension = device.xydimension;
    ObjRCWA.field.Kx = Kx;
    ObjRCWA.field.Ky = Ky;
    ObjRCWA.field.k0 = k0;
    ObjRCWA.field.URC = URC;
    ObjRCWA.field.ERC = ERC;
    ObjRCWA.field.W_V_ref = sparse([W_ref,W_ref;-V_ref,V_ref]);
    ObjRCWA.field.W_V_trn = sparse([W_trn,W_trn;-V_trn,V_trn]);
    ObjRCWA.field.C_input = [C_src;C_ref];    % input for the first layer
end
    
