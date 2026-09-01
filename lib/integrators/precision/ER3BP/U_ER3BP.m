function U = U_ER3BP(q,mu)

if isrow(q) || iscolumn(q)
        q = reshape(q, 1, []); % Ensure q is a 2D matrix
end

% Extract components
x = q(:, 1);
y = q(:, 2);
z = q(:, 3);

% Mass parameters
mu1 = 1 - mu; % Mass of larger primary
mu2 = mu;     % Mass of smaller primary

% Distances (squared)
d = (x + mu2).^2 + y.^2 + z.^2; % Distance to larger primary
r = (x - mu1).^2 + y.^2 + z.^2; % Distance to smaller primary

% Compute U
U = -mu1 ./ sqrt(d) - mu ./ sqrt(r);