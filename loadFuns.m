function message = loadFuns
  assignin('base','VelQuadraticForces',@VelQuadraticForces);
  assignin('base','JacobianBody',@JacobianBody);
  assignin('base','VecTose3',@VecTose3);
  assignin('base','FKinSpace',@FKinSpace);
  assignin('base','Adjoint',@Adjoint);
  assignin('base','ProjectToSE3',@ProjectToSE3);
  assignin('base','JacobianSpace',@JacobianSpace);
  assignin('base','FKinBody',@FKinBody);
  assignin('base','RpToTrans',@RpToTrans);
  assignin('base','ForwardDynamicsTrajectory',@ForwardDynamicsTrajectory);
  assignin('base','CartesianTrajectory',@CartesianTrajectory);
  assignin('base','JointTrajectory',@JointTrajectory);
  assignin('base','ScrewToAxis',@ScrewToAxis);
  assignin('base','VecToso3',@VecToso3);
  assignin('base','AxisAng3',@AxisAng3);
  assignin('base','MatrixLog6',@MatrixLog6);
  assignin('base','MatrixExp3',@MatrixExp3);
  assignin('base','MassMatrix',@MassMatrix);
  assignin('base','ad',@ad);
  assignin('base','InverseDynamics',@InverseDynamics);
  assignin('base','ProjectToSO3',@ProjectToSO3);
  assignin('base','IKinSpace',@IKinSpace);
  assignin('base','AxisAng6',@AxisAng6);
  assignin('base','TestIfSO3',@TestIfSO3);
  assignin('base','DistanceToSO3',@DistanceToSO3);
  assignin('base','ComputedTorque',@ComputedTorque);
  assignin('base','InverseDynamicsTrajectory',@InverseDynamicsTrajectory);
  assignin('base','NearZero',@NearZero);
  assignin('base','so3ToVec',@so3ToVec);
  assignin('base','TransToRp',@TransToRp);
  assignin('base','SimulateControl',@SimulateControl);
  assignin('base','MatrixLog3',@MatrixLog3);
  assignin('base','MatrixExp6',@MatrixExp6);
  assignin('base','Normalize',@Normalize);
  assignin('base','RotInv',@RotInv);
  assignin('base','EndEffectorForces',@EndEffectorForces);
  assignin('base','TestIfSE3',@TestIfSE3);
  assignin('base','TransInv',@TransInv);
  assignin('base','CubicTimeScaling',@CubicTimeScaling);
  assignin('base','GravityForces',@GravityForces);
  assignin('base','EulerStep',@EulerStep);
  assignin('base','DistanceToSE3',@DistanceToSE3);
  assignin('base','se3ToVec',@se3ToVec);
  assignin('base','IKinBody',@IKinBody);
  assignin('base','QuinticTimeScaling',@QuinticTimeScaling);
  assignin('base','ForwardDynamics',@ForwardDynamics);
  assignin('base','ScrewTrajectory',@ScrewTrajectory);
  message='Done importing functions to workspace';
end
function c = VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: A list of joint variables,
%       dthetalist: A list of joint rates,
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
% Returns c: The vector c(thetalist,dthetalist) of Coriolis and centripetal
%            terms for a given thetalist and dthetalist.
% This function calls InverseDynamics with g = 0, Ftip = 0, and 
% ddthetalist = 0.
% Example Input (3 Link Robot):
% 
% clear; clc;
% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3];
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% c = VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)
% 
% Output:
% c =
%    0.2645
%   -0.0551
%   -0.0069

c = InverseDynamics(thetalist, dthetalist, ...
                    zeros(size(thetalist, 1), 1), [0; 0; 0], ...
                    [0; 0; 0; 0; 0; 0], Mlist, Glist, Slist);
end
function Jb = JacobianBody(Blist, thetalist)
% *** CHAPTER 5: VELOCITY KINEMATICS AND STATICS ***
% Takes Blist: The joint screw axes in the end-effector frame when the
%              manipulator is at the home position, in the format of a 
%              matrix with the screw axes as the columns,
%       thetalist: A list of joint coordinates.
% Returns the corresponding body Jacobian (6xn real numbers).
% Example Input:
% 
% clear; clc;
% Blist = [[0; 0; 1;   0; 0.2; 0.2], ...
%        [1; 0; 0;   2;   0;   3], ...
%        [0; 1; 0;   0;   2;   1], ...
%        [1; 0; 0; 0.2; 0.3; 0.4]];
% thetalist = [0.2; 1.1; 0.1; 1.2];
% Jb = JacobianBody(Blist, thetalist)
% 
% Output:
% Jb =
%   -0.0453    0.9950         0    1.0000
%    0.7436    0.0930    0.3624         0
%   -0.6671    0.0362   -0.9320         0
%    2.3259    1.6681    0.5641    0.2000
%   -1.4432    2.9456    1.4331    0.3000
%   -2.0664    1.8288   -1.5887    0.4000

Jb = Blist;
T = eye(4);
for i = length(thetalist) - 1: -1: 1   
    T = T * MatrixExp6(VecTose3(-1 * Blist(:, i + 1) * thetalist(i + 1)));
	Jb(:, i) = Adjoint(T) * Blist(:, i);
end
end
function se3mat = VecTose3(V)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a 6-vector (representing a spatial velocity).
% Returns the corresponding 4x4 se(3) matrix.
% Example Input:
% 
% clear; clc;
% V = [1; 2; 3; 4; 5; 6];
% se3mat = VecTose3(V)
% 
% Output:
% se3mat =
%     0    -3     2     4
%     3     0    -1     5
%    -2     1     0     6
%     0     0     0     0 

se3mat = [VecToso3(V(1: 3)), V(4: 6); 0, 0, 0, 0];
end
function T = FKinSpace(M, Slist, thetalist)
% *** CHAPTER 4: FORWARD KINEMATICS ***
% Takes M: the home configuration (position and orientation) of the 
%          end-effector,
%       Slist: The joint screw axes in the space frame when the manipulator
%              is at the home position,
%       thetalist: A list of joint coordinates.
% Returns T in SE(3) representing the end-effector frame, when the joints 
% are at the specified coordinates (i.t.o Space Frame).
% Example Inputs:
% 
% clear; clc;
% M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
% Slist = [[0; 0;  1;  4; 0;    0], ...
%        [0; 0;  0;  0; 1;    0], ...
%        [0; 0; -1; -6; 0; -0.1]];
% thetalist =[pi / 2; 3; pi];
% T = FKinSpace(M, Slist, thetalist)
% 
% Output:
% T =
%   -0.0000    1.0000         0   -5.0000
%    1.0000    0.0000         0    4.0000
%         0         0   -1.0000    1.6858
%         0         0         0    1.0000

T = M;
for i = size(thetalist): -1: 1
    T = MatrixExp6(VecTose3(Slist(:, i) * thetalist(i))) * T;
end
end
function AdT = Adjoint(T)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes T a transformation matrix SE3. 
% Returns the corresponding 6x6 adjoint representation [AdT].
% Example Input:
% 
% clear; clc;
% T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
% AdT = Adjoint(T)
% 
% Output:
% AdT =
%     1     0     0     0     0     0
%     0     0    -1     0     0     0
%     0     1     0     0     0     0
%     0     0     3     1     0     0
%     3     0     0     0     0    -1
%     0     0     0     0     1     0

[R, p] = TransToRp(T);
AdT = [R, zeros(3); VecToso3(p) * R, R];
end
function T = ProjectToSE3(mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes mat: A matrix near SE(3) to project to SE(3).
% Returns T representing the closest rotation matrix that is in SE(3).
% This function uses singular-value decomposition (see
% http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review)
% and is only appropriate for matrices close to SE(3).
% Example Inputs:
% 
% clear; clc;
% mat = [ 0.675, 0.150,  0.720, 1.2;
%         0.370, 0.771, -0.511, 5.4;
%        -0.630, 0.619,  0.472, 3.6;
%         0.003, 0.002,  0.010, 0.9];
% T = ProjectToSE3(mat)
% 
% Output:
% T =
%     0.6790    0.1489    0.7189    1.2000
%     0.3732    0.7732   -0.5127    5.4000
%    -0.6322    0.6164    0.4694    3.6000
%          0         0         0    1.0000

T = RpToTrans(ProjectToSO3(mat(1: 3, 1: 3)), mat(1: 3, 4));
end
function Js = JacobianSpace(Slist, thetalist)
% *** CHAPTER 5: VELOCITY KINEMATICS AND STATICS ***
% Takes Slist: The joint screw axes in the space frame when the manipulator
%              is at the home position, in the format of a matrix with the
%              screw axes as the columns,
%       thetalist: A list of joint coordinates. 
% Returns the corresponding space Jacobian (6xn real numbers).
% Example Input:
% 
% clear; clc;
% Slist = [[0; 0; 1;   0; 0.2; 0.2], ...
%        [1; 0; 0;   2;   0;   3], ...
%        [0; 1; 0;   0;   2;   1], ...
%        [1; 0; 0; 0.2; 0.3; 0.4]];
% thetalist = [0.2; 1.1; 0.1; 1.2];
% Js = JacobianSpace(Slist, thetalist)
% 
% Output:
% Js =
%         0    0.9801   -0.0901    0.9575
%         0    0.1987    0.4446    0.2849
%    1.0000         0    0.8912   -0.0453
%         0    1.9522   -2.2164   -0.5116
%    0.2000    0.4365   -2.4371    2.7754
%    0.2000    2.9603    3.2357    2.2251

Js = Slist;
T = eye(4);
for i = 2: length(thetalist)
    T = T * MatrixExp6(VecTose3(Slist(:, i - 1) * thetalist(i - 1)));
	Js(:, i) = Adjoint(T) * Slist(:, i);
end
end
function T = FKinBody(M, Blist, thetalist)
% *** CHAPTER 4: FORWARD KINEMATICS ***
% Takes M: the home configuration (position and orientation) of the
%          end-effector,
%       Blist: The joint screw axes in the end-effector frame when the 
%              manipulator is at the home position,
%       thetalist: A list of joint coordinates.
% Returns T in SE(3) representing the end-effector frame when the joints 
% are at the specified coordinates (i.t.o Body Frame).
% Example Inputs:
% 
% clear; clc;
% M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
% Blist = [[0; 0; -1; 2; 0; 0], [0; 0; 0; 0; 1; 0], [0; 0; 1; 0; 0; 0.1]];
% thetalist = [pi / 2; 3; pi];
% T = FKinBody(M, Blist, thetalist)
% 
% Output:
% T =
%   -0.0000    1.0000         0   -5.0000
%    1.0000    0.0000         0    4.0000
%         0         0   -1.0000    1.6858
%         0         0         0    1.0000

T = M;
for i = 1: size(thetalist)
    T = T * MatrixExp6(VecTose3(Blist(:, i) * thetalist(i)));
end
end
function T = RpToTrans(R, p)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes rotation matrix R and position p.
% Returns the corresponding homogeneous transformation matrix T in SE(3).
% Example Input:
% 
% clear; clc;
% R = [[1, 0, 0]; [0, 0, -1]; [0, 1, 0]];
% p = [1; 2; 5];
% T = RpToTrans(R, p)
% 
% Output:
% T =
%     1     0     0     1
%     0     0    -1     2
%     0     1     0     5
%     0     0     0     1

T = [R, p; 0, 0, 0, 1];
end
function [thetamat, dthetamat] ...
         = ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, ...
                                     Ftipmat, Mlist, Glist, Slist, dt, ...
                                     intRes)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: n-vector of initial joint variables,
%       dthetalist: n-vector of initial joint rates,
%       taumat: An N x n matrix of joint forces/torques, where each row is 
%               the joint effort at any time step,
%       g: Gravity vector g,
%       Ftipmat: An N x 6 matrix of spatial forces applied by the 
%                end-effector (If there are no tip forces, the user should 
%                input a zero and a zero matrix will be used),
%       Mlist: List of link frames {i} relative to {i-1} at the home
%              position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
%       dt: The timestep between consecutive joint forces/torques,
%       intRes: Integration resolution is the number of times integration
%               (Euler) takes places between each time step. Must be an 
%               integer value greater than or equal to 1.
% Returns thetamat: The N x n matrix of robot joint angles resulting from 
%                   the specified joint forces/torques,
%         dthetamat: The N x n matrix of robot joint velocities.
% This function simulates the motion of a serial chain given an open-loop 
% history of joint forces/torques. It calls a numerical integration 
% procedure that uses ForwardDynamics.
% Example Inputs (3 Link Robot):
% 
% clc; clear;
% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3];
% taumat = [[3.63, -6.58, -5.57]; [3.74, -5.55, -5.5]; ...
%         [4.31, -0.68, -5.19]; [5.18, 5.63, -4.31]; ...
%         [5.85, 8.17, -2.59]; [5.78, 2.79, -1.7]; ...
%         [4.99, -5.3, -1.19]; [4.08, -9.41, 0.07]; ...
%         [3.56, -10.1, 0.97]; [3.49, -9.41, 1.23]];
% %Initialise robot description (Example with 3 links)
% g = [0; 0; -9.8];
% Ftipmat = ones(size(taumat, 1), 6);
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34);  
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% dt = 0.1;
% intRes = 8;
% [thetamat, dthetamat] ...
% = ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, ...
%                           Ftipmat, Mlist, Glist, Slist, dt, intRes);
% %Output using matplotlib to plot the joint forces/torques
% Tf = size(taumat, 1);
% time=0: (Tf / size(thetamat, 1)): (Tf - (Tf / size(thetamat, 1)));
% plot(time,thetamat(:, 1),'b')
% hold on
% plot(time,thetamat(:, 2), 'g')
% plot(time,thetamat(:, 3), 'r')
% plot(time,dthetamat(:, 1), 'c')
% plot(time,dthetamat(:, 2), 'm')
% plot(time,dthetamat(:, 3), 'y')
% title('Plot of Joint Angles and Joint Velocities')
% xlabel('Time')
% ylabel('Joint Angles/Velocities')
% legend('Theta1', 'Theta2', 'Theta3', 'DTheta1', 'DTheta2', 'DTheta3')
%

taumat = taumat';
Ftipmat = Ftipmat';
thetamat = taumat;
thetamat(:, 1) = thetalist;
dthetamat = taumat;
dthetamat(:, 1) = dthetalist;
for i = 1: size(taumat, 2) - 1
    for j = 1: intRes
       ddthetalist ...
       = ForwardDynamics(thetalist, dthetalist, taumat(:,i), g, ...
                         Ftipmat(:, i), Mlist, Glist, Slist);     
       [thetalist, dthetalist] = EulerStep(thetalist, dthetalist, ...
                                           ddthetalist, dt / intRes);
    end
    thetamat(:, i + 1) = thetalist;
    dthetamat(:, i + 1) = dthetalist;
end
thetamat = thetamat';
dthetamat = dthetamat';
end
function traj = CartesianTrajectory(Xstart, Xend, Tf, N, method)
% *** CHAPTER 9: TRAJECTORY GENERATION ***
% Takes Xstart: The initial end-effector configuration,
%       Xend: The final end-effector configuration,
%       Tf: Total time of the motion in seconds from rest to rest,
%       N: The number of points N > 1 (Start and stop) in the discrete 
%          representation of the trajectory,
%       method: The time-scaling method, where 3 indicates cubic 
%               (third-order polynomial) time scaling and 5 indicates 
%               quintic (fifth-order polynomial) time scaling.
% Returns traj: The discretized trajectory as a list of N matrices in SE(3)
%               separated in time by Tf/(N-1). The first in the list is 
%               Xstart and the Nth is Xend .
% This function is similar to ScrewTrajectory, except the origin of the 
% end-effector frame follows a straight line, decoupled from the rotational
% motion.
% Example Input:
% 
% clear; clc;
% Xstart = [[1, 0, 0, 1]; [0, 1, 0, 0]; [0, 0, 1, 1]; [0, 0, 0, 1]];
% Xend = [[0, 0, 1, 0.1]; [1, 0, 0, 0]; [0, 1, 0, 4.1]; [0, 0, 0, 1]];
% Tf = 5;
% N = 4;
% method = 5;
% traj = CartesianTrajectory(Xstart, Xend, Tf, N, method)
% 
% Output:
% traj =
%    1.0000         0         0    1.0000
%         0    1.0000         0         0
%         0         0    1.0000    1.0000
%         0         0         0    1.0000
%
%    0.9366   -0.2140    0.2774    0.8111
%    0.2774    0.9366   -0.2140         0
%   -0.2140    0.2774    0.9366    1.6506
%         0         0         0    1.0000
%
%    0.2774   -0.2140    0.9366    0.2889
%    0.9366    0.2774   -0.2140         0
%   -0.2140    0.9366    0.2774    3.4494
%         0         0         0    1.0000
%
%   -0.0000    0.0000    1.0000    0.1000
%    1.0000   -0.0000    0.0000         0
%    0.0000    1.0000   -0.0000    4.1000
%         0         0         0    1.0000

timegap = Tf / (N - 1);
traj = cell(1, N);
[Rstart, pstart] = TransToRp(Xstart);
[Rend, pend] = TransToRp(Xend);
for i = 1: N
    if method == 3
        s = CubicTimeScaling(Tf,timegap * (i - 1));
    else
        s = QuinticTimeScaling(Tf,timegap * (i - 1));
    end
    traj{i} ...
    = [Rstart * MatrixExp3(MatrixLog3(Rstart' * Rend) * s), ...
       pstart + s * (pend - pstart); 0, 0, 0, 1];
end
end
function traj = JointTrajectory(thetastart, thetaend, Tf, N, method)
% *** CHAPTER 9: TRAJECTORY GENERATION ***
% Takes thetastart: The initial joint variables,
%       thetaend: The final joint variables,
%       Tf: Total time of the motion in seconds from rest to rest,
%       N: The number of points N > 1 (Start and stop) in the discrete 
%          representation of the trajectory,
%       method: The time-scaling method, where 3 indicates cubic 
%               (third-order polynomial) time scaling and 5 indicates 
%               quintic (fifth-order polynomial) time scaling.
% Returns traj: A trajectory as an N x n matrix, where each row is an
%               n-vector of joint variables at an instant in time. The 
%               first row is thetastart and the Nth row is thetaend . The 
%               elapsed time between each row is Tf/(N - 1).
% The returned trajectory is a straight-line motion in joint space.
% Example Input:
% 
% clear; clc;
% thetastart = [1; 0; 0; 1; 1; 0.2; 0; 1];
% thetaend = [1.2; 0.5; 0.6; 1.1; 2;2; 0.9; 1];
% Tf = 4;
% N = 6;
% method = 3;
% traj = JointTrajectory(thetastart, thetaend, Tf, N, method)
% 
% Output:
% traj =
%   1.0000        0        0   1.0000   1.0000   0.2000        0   1.0000
%   1.0208   0.0520   0.0624   1.0104   1.1040   0.3872   0.0936   1.0000
%   1.0704   0.1760   0.2112   1.0352   1.3520   0.8336   0.3168   1.0000
%   1.1296   0.3240   0.3888   1.0648   1.6480   1.3664   0.5832   1.0000
%   1.1792   0.4480   0.5376   1.0896   1.8960   1.8128   0.8064   1.0000
%   1.2000   0.5000   0.6000   1.1000   2.0000   2.0000   0.9000   1.0000

timegap = Tf / (N - 1);
traj = zeros(size(thetastart, 1), N);
for i = 1: N
    if method == 3
        s = CubicTimeScaling(Tf, timegap * (i - 1));
    else
        s = QuinticTimeScaling(Tf, timegap * (i - 1));
    end
    traj(:, i) = thetastart + s * (thetaend - thetastart);
end
traj = traj';
end
function S = ScrewToAxis(q, s, h)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes q: a point lying on the screw axis,
%       s: a unit vector in the direction of the screw axis, 
%       h: the pitch of the screw axis.
% Returns the corresponding normalized screw axis.
% Example Input:
% 
% clear; clc;
% q = [3; 0; 0];
% s = [0; 0; 1];
% h = 2;
% S = ScrewToAxis(q, s, h)
% 
% Output:
% S =
%     0
%     0
%     1
%     0
%    -3
%     2

S = [s; cross(q, s) + h * s];
end
function so3mat = VecToso3(omg)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a 3-vector (angular velocity).
% Returns the skew symmetric matrix in so(3).
% Example Input:
% 
% clear; clc;
% omg = [1; 2; 3];
% so3mat = VecToso3(omg)
% 
% Output:
% so3mat =
%     0    -3     2
%     3     0    -1
%    -2     1     0

so3mat = [0, -omg(3), omg(2); omg(3), 0, -omg(1); -omg(2), omg(1), 0];
end
function [omghat, theta] = AxisAng3(expc3)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes A 3-vector of exponential coordinates for rotation.
% Returns the unit rotation axis omghat and the corresponding rotation 
% angle theta.
% Example Input:
% 
% clear; clc;
% expc3 = [1; 2; 3];
% [omghat, theta] = AxisAng3(expc3)  
% 
% Output:
% omghat =
%    0.2673
%    0.5345
%    0.8018
% theta =
%    3.7417

theta = norm(expc3);
omghat = expc3 / theta;
end
function expmat = MatrixLog6(T)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a transformation matrix T in SE(3).
% Returns the corresponding se(3) representation of exponential 
% coordinates.
% Example Input:
% 
% clear; clc;
% T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
% expmat = MatrixLog6(T)
% 
% Output:
% expc6 =
%         0         0         0         0
%         0         0   -1.5708    2.3562
%         0    1.5708         0    2.3562
%         0         0         0         0

[R, p] = TransToRp(T);
omgmat = MatrixLog3(R);
if isequal(omgmat, zeros(3))
    expmat = [zeros(3), T(1: 3, 4); 0, 0, 0, 0];
else
    theta = acos((trace(R) - 1) / 2);
    expmat = [omgmat, (eye(3) - omgmat / 2 ...
                      + (1 / theta - cot(theta / 2) / 2) ...
                        * omgmat * omgmat / theta) * p;
              0, 0, 0, 0];    
end
end
function  R = MatrixExp3(so3mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a 3x3 so(3) representation of exponential coordinates.
% Returns R in SO(3) that is achieved by rotating about omghat by theta 
% from an initial orientation R = I.
% Example Input:
% 
% clear; clc;
% so3mat = [[0, -3, 2]; [3, 0, -1]; [-2, 1, 0]];
% R = MatrixExp3(so3mat)  
% 
% Output:
% R =
%   -0.6949    0.7135    0.0893
%   -0.1920   -0.3038    0.9332
%    0.6930    0.6313    0.3481

omgtheta = so3ToVec(so3mat);
if NearZero(norm(omgtheta))
    R = eye(3);
else
    [omghat, theta] = AxisAng3(omgtheta);
    omgmat = so3mat / theta;
    R = eye(3) + sin(theta) * omgmat + (1 - cos(theta)) * omgmat * omgmat;
end
end
function M = MassMatrix(thetalist, Mlist, Glist, Slist)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: A list of joint variables,
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns M: The numerical inertia matrix M(thetalist) of an n-joint serial
%            chain at the given configuration thetalist.
% This function calls InverseDynamics n times, each time passing a 
% ddthetalist vector with a single element equal to one and all other 
% inputs set to zero. Each call of InverseDynamics generates a single 
% column, and these columns are assembled to create the inertia matrix.
% Example Input (3 Link Robot):
% 
% clear; clc;
% thetalist = [0.1; 0.1; 0.1];
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34);  
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% M = MassMatrix(thetalist, Mlist, Glist, Slist)
% 
% Output:
% M =
%   22.5433   -0.3071   -0.0072
%   -0.3071    1.9685    0.4322
%   -0.0072    0.4322    0.1916

n = size(thetalist, 1);
M = zeros(n);
for i = 1: n
   ddthetalist = zeros(n, 1);
   ddthetalist(i) = 1;
   M(:, i) = InverseDynamics(thetalist, zeros(n, 1), ddthetalist, ...
                             [0; 0; 0], [0; 0; 0; 0; 0; 0],Mlist, ...
                             Glist, Slist);
end
end
function adV = ad(V)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes V: 6-vector spatial velocity.
% Returns adV: The corresponding 6x6 matrix.
% Used to calculate the Lie bracket [V1, V2] = [adV1]V2
% Example Input:
%  
% clear; clc;
% V = [1; 2; 3; 4; 5; 6];
% adV = ad(V)
% 
% Output:
% adV =
%     0    -3     2     0     0     0
%     3     0    -1     0     0     0
%    -2     1     0     0     0     0
%     0    -6     5     0    -3     2
%     6     0    -4     3     0    -1
%    -5     4     0    -2     1     0

omgmat = VecToso3(V(1: 3));
adV = [omgmat, zeros(3); VecToso3(V(4: 6)), omgmat];
end
function taulist = InverseDynamics(thetalist, dthetalist, ddthetalist, ...
                                   g, Ftip, Mlist, Glist, Slist)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: n-vector of joint variables,
%       dthetalist: n-vector of joint rates,
%       ddthetalist: n-vector of joint accelerations,
%       g: Gravity vector g,
%       Ftip: Spatial force applied by the end-effector expressed in frame 
%             {n+1},
%       Mlist: List of link frames {i} relative to {i-1} at the home 
%              position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns taulist: The n-vector of required joint forces/torques.
% This function uses forward-backward Newton-Euler iterations to solve the 
% equation:
% taulist = Mlist(thetalist) * ddthetalist + c(thetalist, dthetalist) ...
%           + g(thetalist) + Jtr(thetalist) * Ftip
% Example Input (3 Link Robot):
% 
% clear; clc;
% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3];
% ddthetalist = [2; 1.5; 1];
% g = [0; 0; -9.8];
% Ftip = [1; 1; 1; 1; 1; 1];
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% taulist = InverseDynamics(thetalist, dthetalist, ddthetalist, g, ...
%                         Ftip, Mlist, Glist, Slist)
% 
% Output:
% taulist =
%   74.6962
%  -33.0677
%   -3.2306

n = size(thetalist, 1);
Mi = eye(4);
Ai = zeros(6, n);
AdTi = zeros(6, 6, n + 1);
Vi = zeros(6, n + 1);
Vdi = zeros(6, n + 1);
Vdi(4: 6, 1) = -g;
AdTi(:, :, n + 1) = Adjoint(TransInv(Mlist(:, :, n + 1)));
Fi = Ftip;
taulist = zeros(n, 1);
for i=1: n    
    Mi = Mi * Mlist(:, :, i);
    Ai(:, i) = Adjoint(TransInv(Mi)) * Slist(:, i);    
    AdTi(:, :, i) = Adjoint(MatrixExp6(VecTose3(Ai(:, i) ...
                    * -thetalist(i))) * TransInv(Mlist(:, :, i)));    
    Vi(:, i + 1) = AdTi(:, :, i) * Vi(:, i) + Ai(:, i) * dthetalist(i);
    Vdi(:, i + 1) = AdTi(:, :, i) * Vdi(:, i) ...
                    + Ai(:, i) * ddthetalist(i) ...
                    + ad(Vi(:, i + 1)) * Ai(:, i) * dthetalist(i);    
end
for i = n: -1: 1
    Fi = AdTi(:, :, i + 1)' * Fi + Glist(:, :, i) * Vdi(:, i + 1) ...
         - ad(Vi(:, i + 1))' * (Glist(:, :, i) * Vi(:, i + 1));
    taulist(i) = Fi' * Ai(:, i);
end
end
function R = ProjectToSO3(mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes mat: A matrix near SO(3) to project to SO(3).
% Returns R representing the closest rotation matrix that is in SO(3).
% This function uses singular-value decomposition (see
% http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review)
% and is only appropriate for matrices close to SO(3).
% Example Inputs:
% 
% clear; clc;
% mat = [ 0.675, 0.150,  0.720;
%         0.370, 0.771, -0.511;
%        -0.630, 0.619,  0.472];
% R = ProjectToSO3(mat)
% 
% Output:
% R =
%    0.6790    0.1489    0.7189
%    0.3732    0.7732   -0.5127
%   -0.6322    0.6164    0.4694

[U, S, V] = svd(mat);
R = U * V';
if det(R) < 0
    % In this case the result may be far from mat.
    R = [R(:, 1: 2); -R(:, 3)];
end
end
function [thetalist, success] ...
         = IKinSpace(Slist, M, T, thetalist0, eomg, ev)
% *** CHAPTER 6: INVERSE KINEMATICS ***
% Takes Slist: The joint screw axes in the space frame when the manipulator
%              is at the home position, in the format of a matrix with the
%              screw axes as the columns,
%       M: The home configuration of the end-effector,
%       T: The desired end-effector configuration Tsd,
%       thetalist0: An initial guess of joint angles that are close to 
%                   satisfying Tsd,
%       eomg: A small positive tolerance on the end-effector orientation 
%             error. The returned joint angles must give an end-effector 
%             orientation error less than eomg,
%       ev: A small positive tolerance on the end-effector linear position 
%           error. The returned joint angles must give an end-effector 
%           position error less than ev.
% Returns thetalist: Joint angles that achieve T within the specified 
%                    tolerances,
%         success: A logical value where TRUE means that the function found
%                  a solution and FALSE means that it ran through the set 
%                  number of maximum iterations without finding a solution
%                  within the tolerances eomg and ev.
% Uses an iterative Newton-Raphson root-finding method.
% The maximum number of iterations before the algorithm is terminated has 
% been hardcoded in as a variable called maxiterations. It is set to 20 at 
% the start of the function, but can be changed if needed.  
% Example Inputs:
% 
% clear; clc;
% Slist = [[0; 0;  1;  4; 0;    0], ...
%        [0; 0;  0;  0; 1;    0], ...
%        [0; 0; -1; -6; 0; -0.1]];
% M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
% T = [[0, 1, 0, -5]; [1, 0, 0, 4]; [0, 0, -1, 1.6858]; [0, 0, 0, 1]];
% thetalist0 = [1.5; 2.5; 3];
% eomg = 0.01;
% ev = 0.001;
% [thetalist, success] = IKinSpace(Slist, M, T, thetalist0, eomg, ev)
% 
% Output:
% thetalist =
%    1.5707
%    2.9997
%    3.1415
% success =
%     1

thetalist = thetalist0;
i = 0;
maxiterations = 20;
Tsb = FKinSpace(M, Slist, thetalist);
Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6(TransInv(Tsb) * T));
err = norm(Vs(1: 3)) > eomg || norm(Vs(4: 6)) > ev;
while err && i < maxiterations
    thetalist = thetalist + pinv(JacobianSpace(Slist, thetalist)) * Vs;
    i = i + 1;
    Tsb = FKinSpace(M, Slist, thetalist);
    Vs = Adjoint(Tsb) * se3ToVec(MatrixLog6(TransInv(Tsb) * T));
    err = norm(Vs(1: 3)) > eomg || norm(Vs(4: 6)) > ev;
end
success = ~ err;
end
function [S, theta] = AxisAng6(expc6)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a 6-vector of exponential coordinates for rigid-body motion
% S*theta. 
% Returns S: the corresponding normalized screw axis,
%         theta: the distance traveled along/about S.
% Example Input:
% 
% clear; clc;
% expc6 = [1; 0; 0; 1; 2; 3];
% [S, theta] = AxisAng6(expc6)
%  
% Output:
% S =
%     1
%     0
%     0 
%     1
%     2
%     3
% theta =
%     1

theta = norm(expc6(1: 3));
if NearZero(theta)
    theta = norm(expc6(4: 6));
end
S = expc6 / theta;      
end
function judge = TestIfSO3(mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes mat: A 3x3 matrix.
% Check if mat is close to or on the manifold SO(3).
% Example Inputs:
% 
% clear; clc;
% mat = [1.0, 0.0,   0.0;
%        0.0, 0.1, -0.95;
%        0.0, 1.0,   0.1];
% judge = TestIfSO3(mat)
% 
% Output:
% dudge =
%     0

judge = norm(DistanceToSO3(mat)) < 1e-3;
end
function d = DistanceToSO3(mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes mat: A 3x3 matrix.
% Returns the Frobenius norm to describe the distance of mat from the SO(3) 
% manifold.
% Computes the distance from R to the SO(3) manifold using the following 
% method:
% If det(mat) <= 0, return a large number. 
% If det(mat) > 0, return norm(mat' * mat - I).
% Example Inputs:
% 
% clear; clc;
% mat = [1.0, 0.0,   0.0;
%        0.0, 0.1, -0.95;
%        0.0, 1.0,   0.1];
% d = DistanceToSO3(mat)
% 
% Output:
% d =
%     0.0884

if det(mat) > 0
	d = norm(mat' * mat - eye(3), 'fro');
else
    d = 1e+9;
end
end
function taulist = ComputedTorque(thetalist, dthetalist, eint, g, ...
                                  Mlist, Glist, Slist, thetalistd, ...
                                  dthetalistd, ddthetalistd, Kp, Ki, Kd)
% *** CHAPTER 11: ROBOT CONTROL ***
% Takes thetalist: n-vector of joint variables,
%       dthetalist: n-vector of joint rates,
%       eint: n-vector of the time-integral of joint errors,
%       g: Gravity vector g,
%       Mlist: List of link frames {i} relative to {i-1} at the home 
%              position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
%       thetalistd: n-vector of reference joint variables,
%       dthetalistd: n-vector of reference joint velocities,
%       ddthetalistd: n-vector of reference joint accelerations,
%       Kp: The feedback proportional gain (identical for each joint),
%       Ki: The feedback integral gain (identical for each joint),
%       Kd: The feedback derivative gain (identical for each joint).
% Returns taulist: The vector of joint forces/torques computed by the 
%                  feedback linearizing controller at the current instant.
% Example Input:
% 
% clc; clear;
% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3];
% eint = [0.2; 0.2; 0.2];
% g = [0; 0; -9.8];
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% thetalistd = [1; 1; 1];
% dthetalistd = [2; 1.2; 2];
% ddthetalistd = [0.1; 0.1; 0.1];
% Kp = 1.3;
% Ki = 1.2;
% Kd = 1.1;
% taulist ...
% = ComputedTorque(thetalist, dthetalist, eint, g, Mlist, Glist, Slist, ...
%                thetalistd, dthetalistd, ddthetalistd, Kp, Ki, Kd)
% 
% Output:
% taulist =
%  133.0053
%  -29.9422
%   -3.0328

e = thetalistd - thetalist;
taulist ...
= MassMatrix(thetalist, Mlist, Glist, Slist) ...
  * (Kp * e + Ki * (eint + e) + Kd * (dthetalistd - dthetalist)) ...
  + InverseDynamics(thetalist, dthetalist, ddthetalistd, g, ...
                    zeros(6, 1), Mlist, Glist, Slist);
end
function taumat ...
         = InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, ...
                                     g, Ftipmat, Mlist, Glist, Slist)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetamat: An N x n matrix of robot joint variables,
%       dthetamat: An N x n matrix of robot joint velocities,
%       ddthetamat: An N x n matrix of robot joint accelerations,
%       g: Gravity vector g,
%       Ftipmat: An N x 6 matrix of spatial forces applied by the 
%                end-effector (If there are no tip forces, the user should
%                input a zero and a zero matrix will be used),
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns taumat: The N x n matrix of joint forces/torques for the 
%                 specified trajectory, where each of the N rows is the 
%                 vector of joint forces/torques at each time step.
% This function uses InverseDynamics to calculate the joint forces/torques 
% required to move the serial chain along the given trajectory.
% Example Inputs (3 Link Robot)

% clc; clear;
% %Create a trajectory to follow using functions from Chapter 9
% thetastart = [0; 0; 0];
% thetaend = [pi / 2; pi / 2; pi / 2];
% Tf = 3;
% N= 1000;
% method = 5 ;
% traj = JointTrajectory(thetastart, thetaend, Tf, N, method);
% thetamat = traj;
% dthetamat = zeros(1000, 3);
% ddthetamat = zeros(1000, 3);
% dt = Tf / (N - 1);
% for i = 1: N - 1
%   dthetamat(i + 1, :) = (thetamat(i + 1, :) - thetamat(i, :)) / dt;
%   ddthetamat(i + 1, :) = (dthetamat(i + 1, :) - dthetamat(i, :)) / dt;
% end
% %Initialise robot descripstion (Example with 3 links)
% g = [0; 0; -9.8];
% Ftipmat = ones(N, 6); 
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% taumat = InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, ...
%                                  g, Ftipmat, Mlist, Glist, Slist);
% %Output using matplotlib to plot the joint forces/torques
% time=0: dt: Tf;
% plot(time, taumat(:, 1), 'b')
% hold on
% plot(time, taumat(:, 2), 'g')
% plot(time, taumat(:, 3), 'r')
% title('Plot for Torque Trajectories')
% xlabel('Time')
% ylabel('Torque')
% legend('Tau1', 'Tau2', 'Tau3')
% 

thetamat = thetamat';
dthetamat = dthetamat';
ddthetamat = ddthetamat';
Ftipmat = Ftipmat';
taumat = thetamat;
for i = 1: size(thetamat, 2)
   taumat(:, i) ...
   = InverseDynamics(thetamat(:, i), dthetamat(:, i), ddthetamat(:, i), ...
                     g, Ftipmat(:, i), Mlist, Glist, Slist);
end
taumat = taumat';
end
function judge = NearZero(near)
% *** BASIC HELPER FUNCTIONS ***
% Takes a scalar.
% Checks if the scalar is small enough to be neglected.
% Example Input:
%  
% clear; clc;
% near = -1e-7;
% judge = NearZero(near)
% 
% Output:
% judge =
%     1

judge = norm(near) < 1e-6;
end

function omg = so3ToVec(so3mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a 3x3 skew-symmetric matrix (an element of so(3)).
% Returns the corresponding 3-vector (angular velocity).
% Example Input: 
% 
% clear; clc;
% so3mat = [[0, -3, 2]; [3, 0, -1]; [-2, 1, 0]];
% omg = so3ToVec(so3mat)  
% 
% Output:
% omg =
%     1
%     2
%     3

omg = [so3mat(3, 2); so3mat(1, 3); so3mat(2, 1)];
end
function [R, p] = TransToRp(T)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes the transformation matrix T in SE(3) 
% Returns R: the corresponding rotation matrix
%         p: the corresponding position vector .
% Example Input:
% 
% clear; clc;
% T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
% [R, p] = TransToRp(T)
% 
% Output:
% R =
%     1     0     0
%     0     0    -1
%     0     1     0
% p =
%     0
%     0
%     3

R = T(1: 3, 1: 3);
p = T(1: 3, 4);
end
function [taumat, thetamat] ...
         = SimulateControl(thetalist, dthetalist, g, Ftipmat, Mlist, ...
                           Glist, Slist, thetamatd, dthetamatd, ...
                           ddthetamatd, gtilde, Mtildelist, Gtildelist, ...
                           Kp, Ki, Kd, dt, intRes)
% *** CHAPTER 11: ROBOT CONTROL ***                       
% Takes thetalist: n-vector of initial joint variables,
%       dthetalist: n-vector of initial joint velocities,
%       g: Actual gravity vector g,
%       Ftipmat: An N x 6 matrix of spatial forces applied by the
%                end-effector (If there are no tip forces, the user should 
%                input a zero and a zero matrix will be used),
%       Mlist: Actual list of link frames i relative to i? at the home 
%              position,
%       Glist: Actual spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
%       thetamatd: An Nxn matrix of desired joint variables from the 
%                  reference trajectory,
%       dthetamatd: An Nxn matrix of desired joint velocities,
%       ddthetamatd: An Nxn matrix of desired joint accelerations,
%       gtilde: The gravity vector based on the model of the actual robot
%               (actual values given above),
%       Mtildelist: The link frame locations based on the model of the 
%                   actual robot (actual values given above),
%       Gtildelist: The link spatial inertias based on the model of the 
%                   actual robot (actual values given above),
%       Kp: The feedback proportional gain (identical for each joint),
%       Ki: The feedback integral gain (identical for each joint),
%       Kd: The feedback derivative gain (identical for each joint),
%       dt: The timestep between points on the reference trajectory.
%       intRes: Integration resolution is the number of times integration 
%               (Euler) takes places between each time step. Must be an 
%               integer value greater than or equal to 1.
% Returns taumat: An Nxn matrix of the controller commanded joint 
%                 forces/torques, where each row of n forces/torques 
%                 corresponds to a single time instant,
%         thetamat: An Nxn matrix of actual joint angles.
% The end of this function plots all the actual and desired joint angles.
% Example Usage
% 
% clc; clear;
% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3];
% %Initialize robot description (Example with 3 links)
% g = [0; 0; -9.8];
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% dt = 0.01;
% %Create a trajectory to follow
% thetaend =[pi / 2; pi; 1.5 * pi];
% Tf = 1;
% N = Tf / dt;
% method = 5;
% thetamatd = JointTrajectory(thetalist, thetaend, Tf, N, method);
% dthetamatd = zeros(N, 3);
% ddthetamatd = zeros(N, 3);
% dt = Tf / (N - 1);
% for i = 1: N - 1
%   dthetamatd(i + 1, :) = (thetamatd(i + 1, :) - thetamatd(i, :)) / dt;
%   ddthetamatd(i + 1, :) = (dthetamatd(i + 1, :) ...
%                           - dthetamatd(i, :)) / dt;
% end
% %Possibly wrong robot description (Example with 3 links)
% gtilde = [0.8; 0.2; -8.8];
% Mhat01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.1]; [0, 0, 0, 1]];
% Mhat12 = [[0, 0, 1, 0.3]; [0, 1, 0, 0.2]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% Mhat23 = [[1, 0, 0, 0]; [0, 1, 0, -0.2]; [0, 0, 1, 0.4]; [0, 0, 0, 1]];
% Mhat34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.2]; [0, 0, 0, 1]];
% Ghat1 = diag([0.1, 0.1, 0.1, 4, 4, 4]);
% Ghat2 = diag([0.3, 0.3, 0.1, 9, 9, 9]);
% Ghat3 = diag([0.1, 0.1, 0.1, 3, 3, 3]);
% Gtildelist = cat(3, Ghat1, Ghat2, Ghat3);
% Mtildelist = cat(4, Mhat01, Mhat12, Mhat23, Mhat34); 
% Ftipmat = ones(N, 6);
% Kp = 20;
% Ki = 10;
% Kd = 18;
% intRes = 8;
% [taumat, thetamat] ...
% = SimulateControl(thetalist, dthetalist, g, Ftipmat, Mlist, Glist, ...
%                 Slist, thetamatd, dthetamatd, ddthetamatd, gtilde, ...
%                 Mtildelist, Gtildelist, Kp, Ki, Kd, dt, intRes);
% 

Ftipmat = Ftipmat';
thetamatd = thetamatd';
dthetamatd = dthetamatd';
ddthetamatd = ddthetamatd';
n = size(thetamatd, 2);
taumat = zeros(size(thetamatd));
thetamat = zeros(size(thetamatd));
thetacurrent = thetalist;
dthetacurrent = dthetalist;
eint = zeros(size(thetamatd, 1), 1);
for i=1: n
    taulist ...
    = ComputedTorque(thetacurrent, dthetacurrent, eint, gtilde, ...
                     Mtildelist, Gtildelist, Slist, thetamatd(:, i), ...
                     dthetamatd(:, i), ddthetamatd(:, i), Kp, Ki, Kd);
    for j=1: intRes
        ddthetalist ...
        = ForwardDynamics(thetacurrent, dthetacurrent, taulist, g, ...
                          Ftipmat(:, i), Mlist, Glist, Slist);    
        [thetacurrent, dthetacurrent] ...
        = EulerStep(thetacurrent, dthetacurrent, ddthetalist, ...
                    (dt / intRes));
    end
    taumat(:, i) = taulist;
    thetamat(:, i) = thetacurrent;    
    eint = eint + (dt*(thetamatd(:, i) - thetacurrent));
end
%Output using matplotlib
links = size(thetamat, 1);
leg = cell(1, 2 * links);
time=0: dt: dt * n - dt;
timed=0: dt: dt * n - dt;
figure
hold on
for i=1: links
    col = rand(1, 3);
    plot(time, (thetamat(i, :)'), '-', 'Color', col)
    plot(timed, (thetamatd(i, :)'), '.', 'Color', col)
    leg{2 * i - 1} = (strcat('ActualTheta', num2str(i)));
    leg{2 * i} = (strcat('DesiredTheta', num2str(i)));
end
title('Plot of Actual and Desired Joint Angles')
xlabel('Time')
ylabel('Joint Angles')
legend(leg, 'Location', 'NorthWest')
taumat = taumat';
thetamat = thetamat';
end
function so3mat = MatrixLog3(R)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes R (rotation matrix).
% Returns the corresponding so(3) representation of exponential 
% coordinates.
% Example Input:
% 
% clear; clc;
% R = [[0, 0, 1]; [1, 0, 0]; [0, 1, 0]];
% so3mat = MatrixLog3(R)
% 
% Output:
% angvmat =
%         0   -1.2092    1.2092
%    1.2092         0   -1.2092
%   -1.2092    1.2092         0

acosinput = (trace(R) - 1) / 2;
if acosinput >= 1
    so3mat = zeros(3);
elseif acosinput <= -1
    if ~NearZero(1 + R(3, 3))
        omg = (1 / sqrt(2 * (1 + R(3, 3)))) ...
              * [R(1, 3); R(2, 3); 1 + R(3, 3)];
    elseif ~NearZero(1 + R(2, 2))
        omg = (1 / sqrt(2 * (1 + R(2, 2)))) ...
              * [R(1, 2); 1 + R(2, 2); R(3, 2)];
    else
        omg = (1 / sqrt(2 * (1 + R(1, 1)))) ...
              * [1 + R(1, 1); R(2, 1); R(3, 1)];
    end
    so3mat = VecToso3(pi * omg);
else
	theta = acos(acosinput);
    so3mat = theta * (1 / (2 * sin(theta))) * (R - R');
end
end
function T = MatrixExp6(se3mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a se(3) representation of exponential coordinates.
% Returns a T matrix in SE(3) that is achieved by traveling along/about the 
% screw axis S for a distance theta from an initial configuration T = I.
% Example Input:
% 
% clear; clc;
% se3mat = [ 0,      0,       0,      0;
%          0,      0, -1.5708, 2.3562;
%          0, 1.5708,       0, 2.3562;
%          0,      0,       0,      0]
% T = MatrixExp6(se3mat)
% 
% Output:
% T =
%    1.0000         0         0         0
%         0    0.0000   -1.0000   -0.0000
%         0    1.0000    0.0000    3.0000
%         0         0         0    1.0000 

omgtheta = so3ToVec(se3mat(1: 3, 1: 3));
if NearZero(norm(omgtheta))
    T = [eye(3), se3mat(1: 3, 4); 0, 0, 0, 1];
else
    [omghat, theta] = AxisAng3(omgtheta);
    omgmat = se3mat(1: 3, 1: 3) / theta; 
    T = [MatrixExp3(se3mat(1: 3, 1: 3)), ...
         (eye(3) * theta + (1 - cos(theta)) * omgmat ...
          + (theta - sin(theta)) * omgmat * omgmat) ...
            * se3mat(1: 3, 4) / theta;
         0, 0, 0, 1];
end
end
function norm_v = Normalize(V)
% *** BASIC HELPER FUNCTIONS ***
% Takes in a vector.
% Scales it to a unit vector.
% Example Input:
% 
% clear; clc;
% V = [1; 2; 3];
% norm_v = Normalize(V)
% 
% Output:
% norm_v =
%    0.2673
%    0.5345
%    0.8018

norm_v = V / norm(V);
end


function invR = RotInv(R)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a 3x3 rotation matrix.
% Returns the inverse (transpose).
% Example Input:
% 
% clear; clc;
% R = [0, 0, 1; 1, 0, 0; 0, 1, 0];
% invR = RotInv(R)
% 
% Output:
% invR =
%     0     1     0
%     0     0     1
%     1     0     0

invR = R';
end
function JTFtip = EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: A list of joint variables,
%       Ftip: Spatial force applied by the end-effector expressed in frame 
%             {n+1},
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format 
%              of a matrix with screw axes as the columns,
% Returns JTFtip: The joint forces and torques required only to create the 
%                 end-effector force Ftip.
% This function calls InverseDynamics with g = 0, dthetalist = 0, and 
% ddthetalist = 0.
% Example Input (3 Link Robot):
% 
% clear; clc;
% thetalist = [0.1; 0.1; 0.1];
% Ftip = [1; 1; 1; 1; 1; 1];
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% JTFtip = EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist)
% 
% Output:
% JTFtip =
%    1.4095
%    1.8577
%    1.3924

n = size(thetalist, 1);
JTFtip = InverseDynamics(thetalist, zeros(n, 1), zeros(n, 1), ...
                         [0; 0; 0], Ftip, Mlist, Glist, Slist);
end
function judge = TestIfSE3(mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes mat: A 4x4 matrix.
% Check if mat is close to or on the manifold SE(3).
% Example Inputs:
% 
% clear; clc;
% mat = [1.0, 0.0,   0.0,  1.2;
%        0.0, 0.1, -0.95,  1.5;
%        0.0, 1.0,   0.1, -0.9;
%        0.0, 0.0,   0.1, 0.98];
% judge = TestIfSE3(mat)
% 
% Output:
% judge =
%     0

judge = norm(DistanceToSE3(mat)) < 1e-3;
end
function invT = TransInv(T)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a transformation matrix T.
% Returns its inverse. 
% Uses the structure of transformation matrices to avoid taking a matrix
% inverse, for efficiency.
% Example Input:
% 
% clear; clc;
% T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
% invT = TransInv(T)
% 
% Ouput:
% invT =
%     1     0     0     0
%     0     0     1    -3
%     0    -1     0     0
%     0     0     0     1

[R, p] = TransToRp(T);
invT = [R', -R' * p; 0, 0, 0, 1];
end
function s = CubicTimeScaling(Tf, t)
% *** CHAPTER 9: TRAJECTORY GENERATION ***
% Takes Tf: Total time of the motion in seconds from rest to rest,
%       t: The current time t satisfying 0 < t < Tf.
% Returns s: The path parameter s(t) corresponding to a third-order 
%            polynomial motion that begins and ends at zero velocity.
% Example Input: 
% 
% clear; clc;
% Tf = 2;
% t = 0.6;
% s = CubicTimeScaling(Tf,t)
% 
% Output:
% s =
%    0.2160

s = 3 * (t / Tf) ^ 2 - 2 * (t / Tf) ^ 3;
end
function grav = GravityForces(thetalist, g, Mlist, Glist, Slist)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: A list of joint variables,
%       g: 3-vector for gravitational acceleration,
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns grav: The joint forces/torques required to overcome gravity at 
%               thetalist
% This function calls InverseDynamics with Ftip = 0, dthetalist = 0, and 
% ddthetalist = 0.
% Example Input (3 Link Robot):
% 
% clear; clc;
% thetalist = [0.1; 0.1; 0.1];
% g = [0; 0; -9.8];
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% grav = GravityForces(thetalist, g, Mlist, Glist, Slist)
% 
% Output:
% grav =
%   28.4033
%  -37.6409
%   -5.4416

n = size(thetalist, 1);
grav = InverseDynamics(thetalist, zeros(n, 1), zeros(n, 1) ,g, ...
                       [0; 0; 0; 0; 0; 0], Mlist, Glist, Slist);
end
function [thetalistNext, dthetalistNext] ...
         = EulerStep(thetalist, dthetalist, ddthetalist, dt)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: n-vector of joint variables,
%       dthetalist: n-vector of joint rates,
%       ddthetalist: n-vector of joint accelerations,
%       dt: The timestep delta t.
% Returns thetalistNext: Vector of joint variables after dt from first 
%                        order Euler integration,
%         dthetalistNext: Vector of joint rates after dt from first order 
%                         Euler integration.
% Example Inputs (3 Link Robot):
% 
% clear; clc;
% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3];
% ddthetalist = [2; 1.5; 1];
% dt = 0.1;
% [thetalistNext, dthetalistNext] = EulerStep(thetalist, dthetalist, ...
%                                           ddthetalist, dt)
% 
% Output:
% thetalistNext =
%    0.1100
%    0.1200
%    0.1300
% dthetalistNext =
%    0.3000
%    0.3500
%    0.4000

thetalistNext = thetalist + dt * dthetalist;
dthetalistNext = dthetalist + dt * ddthetalist;
end
function d = DistanceToSE3(mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes mat: A 4x4 matrix.
% Returns the Frobenius norm to describe the distance of mat from the SE(3) 
% manifold.
% Compute the determinant of matR, the top 3x3 submatrix of mat. If 
% det(matR) <= 0, return a large number. If det(matR) > 0, replace the top 
% 3x3 submatrix of mat with matR' * matR, and set the first three entries 
% of the fourth column of mat to zero. Then return norm(mat - I).
% Example Inputs:
% 
% clear; clc;
% mat = [1.0, 0.0,   0.0,  1.2;
%        0.0, 0.1, -0.95,  1.5;
%        0.0, 1.0,   0.1, -0.9;
%        0.0, 0.0,   0.1, 0.98];
% d = DistanceToSE3(mat)
% 
% Output:
% d =
%     0.1349

[R, p] = TransToRp(mat);
if det(R) > 0
	d = norm([R' * R, [0; 0; 0]; mat(4, :)] - eye(4), 'fro');
else
    d = 1e+9;
end
end
function V = se3ToVec(se3mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes se3mat a 4x4 se(3) matrix
% Returns the corresponding 6-vector (representing spatial velocity).
% Example Input:
% 
% clear; clc;
% se3mat = [[0, -3, 2, 4]; [3, 0, -1, 5]; [-2, 1, 0, 6]; [0, 0, 0, 0]];
% V = se3ToVec(se3mat)
% 
% Output:
% V =
%     1
%     2
%     3
%     4
%     5
%     6

V = [se3mat(3, 2); se3mat(1, 3); se3mat(2, 1); se3mat(1: 3, 4)];
end
function [thetalist, success] = IKinBody(Blist, M, T, thetalist0, eomg, ev)
% *** CHAPTER 6: INVERSE KINEMATICS ***
% Takes Blist: The joint screw axes in the end-effector frame when the
%              manipulator is at the home position, in the format of a 
%              matrix with the screw axes as the columns,
%       M: The home configuration of the end-effector,
%       T: The desired end-effector configuration Tsd,
%       thetalist0: An initial guess of joint angles that are close to 
%                   satisfying Tsd,
%       eomg: A small positive tolerance on the end-effector orientation
%             error. The returned joint angles must give an end-effector 
%             orientation error less than eomg,
%       ev: A small positive tolerance on the end-effector linear position 
%           error. The returned joint angles must give an end-effector
%           position error less than ev.
% Returns thetalist: Joint angles that achieve T within the specified 
%                    tolerances,
%         success: A logical value where TRUE means that the function found
%                  a solution and FALSE means that it ran through the set 
%                  number of maximum iterations without finding a solution
%                  within the tolerances eomg and ev.
% Uses an iterative Newton-Raphson root-finding method.
% The maximum number of iterations before the algorithm is terminated has 
% been hardcoded in as a variable called maxiterations. It is set to 20 at 
% the start of the function, but can be changed if needed.  
% Example Inputs:
% 
% clear; clc;
% Blist = [[0; 0; -1; 2; 0; 0], [0; 0; 0; 0; 1; 0], [0; 0; 1; 0; 0; 0.1]];
% M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
% T = [[0, 1, 0, -5]; [1, 0, 0, 4]; [0, 0, -1, 1.6858]; [0, 0, 0, 1]];
% thetalist0 = [1.5; 2.5; 3];
% eomg = 0.01;
% ev = 0.001;
% [thetalist, success] = IKinBody(Blist, M, T, thetalist0, eomg, ev)
% 
% Output:
% thetalist =
%    1.5707
%    2.9997
%    3.1415
% success =
%     1

thetalist = thetalist0;
i = 0;
maxiterations = 20;
Vb = se3ToVec(MatrixLog6(TransInv(FKinBody(M, Blist, thetalist)) * T));
err = norm(Vb(1: 3)) > eomg || norm(Vb(4: 6)) > ev;
while err && i < maxiterations
    thetalist = thetalist + pinv(JacobianBody(Blist, thetalist)) * Vb;
    i = i + 1;
    Vb = se3ToVec(MatrixLog6(TransInv(FKinBody(M, Blist, thetalist)) * T));
    err = norm(Vb(1: 3)) > eomg || norm(Vb(4: 6)) > ev;
end
success = ~ err;
end
function s = QuinticTimeScaling(Tf, t)
% *** CHAPTER 9: TRAJECTORY GENERATION ***
% Takes Tf: Total time of the motion in seconds from rest to rest,
%       t: The current time t satisfying 0 < t < Tf.
% Returns s: The path parameter s(t) corresponding to a fifth-order
%            polynomial motion that begins and ends at zero velocity and 
%            zero acceleration.
% Example Input: 
% 
% clear; clc;
% Tf = 2;
% t = 0.6;
% s = QuinticTimeScaling(Tf,t)
% 
% Output:
% s =
%    0.1631

s = 10 * (t / Tf) ^ 3 - 15 * (t / Tf) ^ 4 + 6 * (t / Tf) ^ 5;
end
function ddthetalist = ForwardDynamics(thetalist, dthetalist, taulist, ...
                                       g, Ftip, Mlist, Glist, Slist)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: A list of joint variables,
%       dthetalist: A list of joint rates,
%       taulist: An n-vector of joint forces/torques,
%       g: Gravity vector g,
%       Ftip: Spatial force applied by the end-effector expressed in frame 
%             {n+1},
%       Mlist: List of link frames i relative to i-1 at the home position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns,
% Returns ddthetalist: The resulting joint accelerations.
% This function computes ddthetalist by solving:
% Mlist(thetalist) * ddthetalist = taulist - c(thetalist,dthetalist) ...
%                                  - g(thetalist) - Jtr(thetalist) * Ftip
% Example Input (3 Link Robot):
% 
% clear; clc;
% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3]; 
% taulist = [0.5; 0.6; 0.7];
% g = [0; 0; -9.8];
% Ftip = [1; 1; 1; 1; 1; 1];
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
% G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
% G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
% ddthetalist = ForwardDynamics(thetalist, dthetalist, taulist, g, ...
%                             Ftip, Mlist, Glist, Slist)
% 
% Output:
% ddthetalist =
%   -0.9739
%   25.5847
%  -32.9150

ddthetalist = MassMatrix(thetalist, Mlist, Glist, Slist) ...
              \ (taulist - VelQuadraticForces(thetalist, dthetalist, ...
                                              Mlist, Glist, Slist) ...
                 - GravityForces(thetalist, g, Mlist, Glist, Slist) ...
                 - EndEffectorForces(thetalist, Ftip, Mlist, Glist, ...
                                     Slist));
end
function traj = ScrewTrajectory(Xstart, Xend, Tf, N, method)
% *** CHAPTER 9: TRAJECTORY GENERATION ***
% Takes Xstart: The initial end-effector configuration,
%       Xend: The final end-effector configuration,
%       Tf: Total time of the motion in seconds from rest to rest,
%       N: The number of points N > 1 (Start and stop) in the discrete
%          representation of the trajectory,
%       method: The time-scaling method, where 3 indicates cubic
%               (third-order polynomial) time scaling and 5 indicates 
%               quintic (fifth-order polynomial) time scaling.
% Returns traj: The discretized trajectory as a list of N matrices in SE(3)
%               separated in time by Tf/(N-1). The first in the list is 
%               Xstart and the Nth is Xend .
% This function calculates a trajectory corresponding to the screw motion 
% about a space screw axis.
% Example Input:
% 
% clear; clc;
% Xstart = [[1 ,0, 0, 1]; [0, 1, 0, 0]; [0, 0, 1, 1]; [0, 0, 0, 1]];
% Xend = [[0, 0, 1, 0.1]; [1, 0, 0, 0]; [0, 1, 0, 4.1]; [0, 0, 0, 1]];
% Tf = 5;
% N = 4;
% method = 3;
% traj = ScrewTrajectory(Xstart, Xend, Tf, N, method)
% 
% Output:
% traj =
%    1.0000         0         0    1.0000
%         0    1.0000         0         0
%         0         0    1.0000    1.0000
%         0         0         0    1.0000
%
%    0.9041   -0.2504    0.3463    0.4410
%    0.3463    0.9041   -0.2504    0.5287
%   -0.2504    0.3463    0.9041    1.6007
%         0         0         0    1.0000
%
%    0.3463   -0.2504    0.9041   -0.1171
%    0.9041    0.3463   -0.2504    0.4727
%   -0.2504    0.9041    0.3463    3.2740
%         0         0         0    1.0000
%
%   -0.0000    0.0000    1.0000    0.1000
%    1.0000   -0.0000    0.0000   -0.0000
%    0.0000    1.0000   -0.0000    4.1000
%         0         0         0    1.0000

timegap = Tf / (N - 1);
traj = cell(1, N);
for i = 1: N
    if method == 3
        s = CubicTimeScaling(Tf, timegap * (i - 1));
    else
        s = QuinticTimeScaling(Tf, timegap * (i - 1));
    end
    traj{i} = Xstart * MatrixExp6(MatrixLog6(TransInv(Xstart) * Xend) * s);
end
end


