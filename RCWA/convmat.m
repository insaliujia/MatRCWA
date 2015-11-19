function C = convmat(A,PQR)
% CONVMAT Rectangular Convolution Matrix
% This code can construct 3D problems
% 
% numel(PQR) = 1; for 1D problems
% numel(PQR) = 2; for 2D problems
% numel(PQR) = 3; for 3D problems
%
% This MATLAB function constructs convolution matrices
% from a real space grid.
%% HANDLE INPUT AND OUTPUT ARGUMENTS
% DETERMINE SIZE OF A
[Nx,Ny,Nz] = size(A);

% HANDLE NUMBER OF HARMONICS FOR ALL DIMENSIONS
if numel(PQR) == 1
    P = PQR(1);
    Q = 1;
    R = 1;
elseif numel(PQR) == 2
    P = PQR(1);
    Q = PQR(2);
    R = 1;
elseif numel(PQR) == 3
    P = PQR(1);
    Q = PQR(2);
    R = PQR(3);
end
% COMPUTE INDICES OF SPATIAL HARMONICS
NH = P*Q*R; %total number
p = [-floor(P/2):+floor(P/2)]; %indices along x
q = [-floor(Q/2):+floor(Q/2)]; %indices along y
r = [-floor(R/2):+floor(R/2)]; %indices along z
% COMPUTE FOURIER COEFFICIENTS OF A
A = fftshift(fftn(A)) / (Nx*Ny*Nz);
% COMPUTE ARRAY INDICES OF CENTER HARMONIC
p0 = 1 + floor(Nx/2);
q0 = 1 + floor(Ny/2);
r0 = 1 + floor(Nz/2);
C = ones(NH,NH);
for rrow = 1 : R
    for qrow = 1 : Q
        for prow = 1 : P
            row = (rrow-1)*Q*P + (qrow-1)*P + prow;
            for rcol = 1 : R
                for qcol = 1 : Q
                    for pcol = 1 : P
                        col = (rcol-1)*Q*P + (qcol-1)*P + pcol;
                        pfft = p(prow) - p(pcol);
                        qfft = q(qrow) - q(qcol);
                        rfft = r(rrow) - r(rcol);
                        C(row,col) = A(p0+pfft,q0+qfft,r0+rfft);
                    end
                end
            end
        end
    end
end