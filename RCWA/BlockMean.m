function Y = BlockMean(X, V, W)
% 2D block mean over 1st and 2nd dim [MEX]
% The mean of V*W elements along the 1st and 2nd dimensions is calculated.
% Y = BlockMean(X, V, W)
% INPUT:
%   X: array of any size.
%   V, W: Scalar numerical with integer value. Each element of the output is
%      the mean over V*W neighbouring elements of the input.
%      V and W are limited to 256 to limit memory usage.
%      A square V*V block is used, if W is omitted.
% OUTPUT:
%   Y: UINT8 or DOUBLE array, the 1st and 2nd dimensions are V and W times
%      shorter: [FLOOR(X / V) x FLOOR(Y / W) x (further dims...)]
%      If the size of the 1st or 2nd dimension is not a multiple of V and W,
%      the remaining elements at the end are skipped.
%      The empty matrix is replied for empty inputs or if the 1st or 2nd
%      dimension is shorter than V or W.
%

% Get size of X and calculate the reduced sizes for the 1st and 2nd dimension:
if nargin < 3
   W = V;
end
S = size(X);
M = S(1) - mod(S(1), V);
N = S(2) - mod(S(2), W);
if M * N == 0
   Y = X([]);  % Copy type of X
   return;
end
MV = M / V;
NW = N / W;

% Cut and reshape input such that the 1st and 3rd dimension have the lengths V
% and W:
XM = reshape(X(1:M, 1:N, :), V, MV, W, NW, []);

% Different methods depending on the type of the input:
Y = sum(sum(XM, 1), 3) .* (1.0 / (V * W));


% Remove singleton dimensions:
S(1) = MV;
S(2) = NW;
Y    = reshape(Y, S);

return;