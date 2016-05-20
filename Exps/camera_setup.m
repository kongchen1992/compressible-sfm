% This script is used to generate synthetic orthonormal camera matrix.
function Ri = camera_setup(F, type)

if nargin < 2
    type = 'continuous';
end

Ri = cell(F, 1);

switch type
    case 'random'
        for i = 1:F
            R = orth(randn(3));
            Ri{i} = R(1:2, :);
        end
    case 'continuous'
        Ux = @(u) [0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];
        theta = randn(1); % create an angle randomly
        u = randn(3, 1); u = u/norm(u);% create a unit axis randomly
        for i = 1:F
            theta = theta + 1;
            % create a rotation matrix from angle and axis above.
            orthR = cos(theta)*eye(3) + sin(theta)*Ux(u) + (1-cos(theta))*kron(u,u');
            Ri{i} = orthR(1:2, :);
        end
    case 'vertical'
        Ux = @(u) [0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];
        theta = randn(1);
        u = [0;1;0];
        for i = 1:F
            theta = theta + 0.1;
            orthR = cos(theta)*eye(3) + sin(theta)*Ux(u) + (1-cos(theta))*kron(u,u');
            Ri{i} = orthR(1:2, :);
        end
    otherwise
        error('%s is a wrong type!', type);
end
