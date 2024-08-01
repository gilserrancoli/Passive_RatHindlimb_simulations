function [N_matrix, Np_matrix, Npp_matrix] = bspline_creator(varargin)

% BSPLINE_CREATOR creates the matrices necessary for a B-spline curve and
%    its first and second time derivative.
%
%    [N, NP, NPP] = BSPLINE_CREATOR computes and saves to disk the B-spline
%    N coefficent matrix and its first and second time derivatives (NP and
%    NPP) given the default B-spline order (5th), number of nodes (20), and
%    number of output time frames (100).
%
%    [N, NP, NPP] = BSPLINE_CREATOR('ORDER', 'NODES', 'FRAMES','SAMPLEINGRATE') implements
%    the user defined B-spline 'ORDER', number of 'NODES', and number of
%    output time 'FRAMES'.
% 
%    Example
%        ORDER, NODES, and FRAMES can be specified as:
%        [N, Np, Npp] = bspline_creator('5', '20', '100');
%
%    [N, NP, NPP] = BSPLINE_CREATOR('MAZZA') uses 15 predefined B-spline
%    nodes published by Mazza and Cappozzo (2004) as both an example and
%    formulation verification.
%
%    Example
%        The user may see this example using the following command line:
%        bspline_creator('Mazza');
%
%    [N, NP, NPP] = BSPLINE_CREATOR('SINE') uses 15 predefined B-spline
%    nodes to recreate a Sine curve and its first and second time
%    derivatives as both an example and error quantification.
%
%    Example
%        The user may see this example using the following command line:
%        bspline_creator('Sine');

% -------------------------------------------------------------------------

% Check to see if input arguments are present; otherwise, use defaults.
if nargin >= 1
    Order = varargin{1};
else
    Order = 5;
end
if nargin >= 2
    Nodes = varargin{2};
else
    Nodes = 20;
end
if nargin >= 3
    TimeFrames = varargin{3};
else
    TimeFrames = 100;
end
if nargin >= 4
    SamplingRate = varargin{4};
else
    SamplingRate = 1/120;
end

% -------------------------------------------------------------------------

% Define time vector for chosen number of time frames.
XTrue = linspace(0, (SamplingRate)*(TimeFrames-1), TimeFrames)'; % instead of:  XTrue = linspace(0, (2*pi), TimeFrames)';

% -------------------------------------------------------------------------

% Example inputs for 5th order B-spline using 15 nodes from:
%
% Mazza, C and Cappozzo, A, 2004.  An optimization algorithm for human 
% joint angle time-history generation using external force data.  Journal 
% of Biomedical Engineering, Volume 32, Number 5, Pages 764-772.

if nargin == 1 && strcmpi(varargin{1}, 'mazza'),
    
    % Redefine input arguments.
    Order = 5;
    Nodes = 15;
    TimeFrames = 50;

    % Define B-spline node values.
    XIn = linspace(0, (2*pi), 15)'; 
    YIn = [90; 89; 89; 90; 95; 105; 120; 98; 65; 15; 8; 3; 2; 1; 0];

    % Redefine time vector for chosen number of time frames.
    XTrue = linspace(0, (2*pi), TimeFrames)';
    
end

% -------------------------------------------------------------------------

if nargin == 1 && strcmpi(varargin{1}, 'sine'),
    
    % Example input and output comparison for 5th order B-spline using 15 nodes
    % to parameterize y = sin(x) and its first and second derivatives.

    % Redefine input arguments.
    Order = 5;
    Nodes = 15;
    TimeFrames = 50;

    % Define B-spline node values.
    XIn = linspace(0, (2*pi), 15)';
    YIn = [-0.00219805895929
        0.12259977610370
        0.37749238177445
        0.72864778856466
        1.05292045242948
        1.05213726141776
        0.65075564876474
        -0.00000000005385
        -0.65075564887636
        -1.05213726167670
        -1.05292045296986
        -0.72864778927913
        -0.37749238282162
        -0.12259977736280
        0.00219805760280];

    % Redefine time vector for chosen number of time frames.
    XTrue = linspace(0, (2*pi), TimeFrames)';

    % Define true values.
    YTrue = sin(XTrue);
    YpTrue = cos(XTrue);
    YppTrue = -sin(XTrue);
    
end

% -------------------------------------------------------------------------

% Check if input control points are present; otherwise, use zeros.
if exist('XIn') == 0,
    XIn = linspace(0, (2*pi), Nodes)';
    YIn = zeros(Nodes, 1);
end

% Compile the input control points (including row of zeros for later matrix
% multiplication).
Points = [XIn'; YIn'];
Points(3, :) = 0;

% Define the number of control points.
numberPoints = size(Points, 2);

% Check if the degree of the curve is greater than number of control points
% minus one.
if Order > (numberPoints-1),
    error(sprintf('Error: the maximum order of the B-spline cannot exceed %d', numberPoints-1));
end

% Define the parameters n (allowing indices to start at zero for n+1
% control points) and d (controlling degree of polynomial in u and curve
% continuity).
n = numberPoints - 1;
d = Order + 1;

% Define parametric knots (or knot values) Uj for open curve B-spline,
% which relate the parametric variable u to the Pj control points.
for tempJ = 1:(n + d + 1),
    j = tempJ - 1;
    if j < d
        Uj(tempJ) = 0;
    elseif j >= d & j <= n
        Uj(tempJ) = j - d + 1;
    elseif j > n
        Uj(tempJ) = n - d + 2;
    end
end

Ujmax = max(Uj);
tf = 2*pi;
Uj = Uj * (tf/Ujmax);

% Define the range of the parametric variable u.
umax = n - d + 2;
u_range = 0:(umax/(TimeFrames-1)):umax;

u_range = u_range * (tf/umax);
umax = tf;

% Allocate memory for the output points matrix Pout and B-spline basis
% blending functions N (a weighting function.
Pout = zeros(3, (size(u_range, 2)));
N = zeros(1, (n + d));
Np = zeros(1, (n + d));

% Instantiate the N matrix for the shape functions.
N_matrix = [];

% Loop through the entire range of the parametric variable u.
for u = u_range,

    % Loop through each normalized polynomial degree in u.
    for k = 1:d

        % Check if the normalized polynomial degree in u is equal to one.
        if k == 1,

            % Loop through each parametric knot (or knot values).
            for i = 1:numberPoints,

                % Calculate the first order basis functions.
                if (u >= Uj(i)) & (u < Uj(i+1)),
                    N(i) = 1;
                elseif (abs(u - umax) < 1e-8) && (abs(u - Uj(i) - (tf/Ujmax)) < 1e-8) && (u - Uj(i+1) <= 1e-8),
                    N(i) = 1;
                else,
                    N(i) = 0;
                end

            end	% End loop through each parametric knot (or knot values).

        % Otherwise if the normalized polynomial degree in u is not equal 
        % to one.
        else

            % Loop through each parametric knot (or knot values).
            for i=1:numberPoints

                % Calculate the higher order basis functions.
                if ((Uj(i+k-1) - Uj(i)) == 0) && ((Uj(i+k) - Uj(i+1)) == 0),
                    N(i) = 0;
                elseif ((Uj(i+k-1) - Uj(i)) == 0) && ((Uj(i+k) - Uj(i+1)) ~= 0),
                    N(i) = ((Uj(i+k) - u)*N(i+1)) / (Uj(k+i) - Uj(i+1));
                elseif ((Uj(i+k-1) - Uj(i)) ~= 0) && ((Uj(i+k) - Uj(i+1)) == 0),
                    N(i) = ((u - Uj(i))*N(i)) / (Uj(i+k-1) - Uj(i));
                else,
                    N(i) = ((u - Uj(i))*N(i) / (Uj(i+k-1) - Uj(i))) + ((Uj(i+k) - u)*N(i+1) / (Uj(k+i) - Uj(i+1)));
                end

            end % End loop through each parametric knot (or knot values).

        end % End check if the normalized polynomial degree in u is equal to one or not.

    end % End of for loop through each polynomial degree in u.

    % Compile all time frames of the shape functions.
    N_matrix = [N_matrix; N(1, 1:numberPoints)];

end % End of for loop through the entire range of the parametric variable u.

% Compute the output points.
YOut = N_matrix*YIn;

% Compute first derivative matrix and output points.
[Fx, Np_matrix] = gradient(N_matrix, XTrue, XTrue);
YpOut = Np_matrix*YIn;

% Compute second derivative matrix and output points.
[Fxx, Npp_matrix] = gradient(Np_matrix, XTrue, XTrue);
YppOut = Npp_matrix*YIn;

% % Normalize time for following plots.
% XIn = (XIn/(2*pi)) * 100;
% XTrue = (XTrue/(2*pi)) * 100;
% 
% % Plot the control points and resulting B-spline curve.
% if exist('YTrue'),
%     
%     figure
%     clf
% 
%     subplot(3, 1, 1)
%     plot(XTrue, YTrue, 'g.-', 'LineWidth', 0.5);
%     hold on
%     plot(XIn, YIn, 'ko', 'LineWidth', 2);
%     plot(XTrue, YOut, 'bo-', 'LineWidth', 0.5, 'MarkerSize', 5)
%     legend('True', 'Node', 'B-Spline', -1)
%     rms = sqrt(sum((YOut - YTrue).^2)/length(YOut));
%     title(['Generated B-Spline, Q(t) = N(t) * q_i, RMS Error = ' num2str(rms)], 'FontWeight', 'bold')
%     set(gca, 'XTick', [])
%     ylabel('Control')
%     axis tight
% 
%     subplot(3, 1, 2)
%     plot(XTrue, YpTrue, 'g.-', 'LineWidth', 0.5);
%     hold on
%     plot(XTrue, YpOut, 'bo-', 'LineWidth', 0.5, 'MarkerSize', 5)
%     legend('True', 'B-Spline', -1)
%     rms = sqrt(sum((YpOut - YpTrue).^2)/length(YpOut));
%     title(['First Derivative of B-Spline, Q''(t) = N''(t) * q_i, RMS Error = ' num2str(rms)], 'FontWeight', 'bold')
%     set(gca, 'XTick', [])
%     ylabel('Control')
%     axis tight
% 
%     subplot(3, 1, 3)
%     plot(XTrue, YppTrue, 'g.-', 'LineWidth', 0.5);
%     hold on
%     plot(XTrue, YppOut, 'bo-', 'LineWidth', 0.5, 'MarkerSize', 5)
%     legend('True', 'B-Spline', -1)
%     rms = sqrt(sum((YppOut - YppTrue).^2)/length(YppOut));
%     title(['Second Derivative of B-Spline, Q''''(t) = N''''(t) * q_i, RMS Error = ' num2str(rms)], 'FontWeight', 'bold')
%     xlabel('Normalized Time (%)')
%     ylabel('Control')
%     axis tight
%     
% else
% 
%     figure 
%     clf 
%     hold on 
%     plot(XIn, YIn, 'ko', 'LineWidth', 2); 
%     plot(XTrue, YOut, 'bo-', 'LineWidth', 0.5, 'MarkerSize', 5) 
%     legend('Node', 'B-Spline', -1) 
%     title('Generated B-Spline, Q(t) = N(t) * q_i', 'FontWeight', 'bold')
%     xlabel('Normalized Time (%)') 
%     ylabel('Control')
%     axis tight
%     
% end

% Save N matrix and its first and second derivatives for future
% applications.
% save(['N_order' num2str(Order) '_nodes' num2str(Nodes) '_points' num2str(TimeFrames) '_matrix.txt'], 'N_matrix', '-ascii', '-double')
% save(['Np_order' num2str(Order) '_nodes' num2str(Nodes) '_points' num2str(TimeFrames) '_matrix.txt'], 'Np_matrix', '-ascii', '-double')
% save(['Npp_order' num2str(Order) '_nodes' num2str(Nodes) '_points' num2str(TimeFrames) '_matrix.txt'], 'Npp_matrix', '-ascii', '-double')

