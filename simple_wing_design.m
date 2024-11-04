
%%
% ================================================================ %
%
% MIT License
%
% Copyright (c) [year:2022] [author:Metehan Yayla]
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
% ================================================================ %
%%

clc
clear all
close all

format compact

set(0,'defaultAxesFontSize',24)

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','Tex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(0,'DefaultTextInterpreter', 'latex')

opengl software


%% constant parameters

grav = 9.80665;     % m/s^2
rho_SL = 1.225;     % kg/m^3


%% Airfoil Data

% SG6042
% Cl_func = @(alpha) (0.0000097678*alpha.^4 - 0.000274243*alpha.^3 - 0.0021*alpha.^2 + 0.119*alpha + 0.73088);
% Cl_aoa_func = @(alpha) (4*0.0000097678*alpha.^3 - 3*0.000274243*alpha.^2 - 2*0.0021*alpha + 0.119)*180/pi;

% SD7043
Cl_func = @(alpha_deg) (-0.0002741870*alpha_deg.^3 + 0.0015541706*alpha_deg.^2 + 0.1095103276*alpha_deg + 0.4222541374);
Cl_aoa_func = @(alpha_deg) (-3*0.0002741870*alpha_deg.^2 + 2*0.0015541706*alpha_deg + 0.1095103276)*180/pi;


%% Wing Geometry

% performance parameters
V_cruise = 208;                 % cruise velocity, [KCAS]
V_cruise = V_cruise/1.94;       % cruise velocity, [m/s]

% aircraft/wing parameters
Wo = 51000;                     % aircraft mass, [kg]
Sref = 116;                     % reference wing area, [m^2]
AR = 11;                        % Aspect Ratio
taper = 0.5;                    % taper ratio
incidence = 3.0;                % wing incidence angle [deg]
twist = -1.0;                   % wing twist angle [deg]

% wingspan [m]
wing_span = sqrt(Sref*AR);      % entire wing span, [m]
half_wing_span = wing_span/2.0;

% chord length
c_root = 2*Sref/((1+taper)*wing_span);
c_tip = c_root*taper;

% Wing area
half_wing_area = Sref/2;

% mean aerodynamic chord
mac = 2/3*c_root*(1+taper+taper^2)/(1+taper);

% Trailing Edge and Leading Edge locations along Fuselage axis
x_rootChord_LE = c_root;
x_rootChord_TE = 0.0;

x_tipChord_LE = x_rootChord_LE - (c_root-c_tip)/4;
x_tipChord_TE = x_rootChord_TE + 3/4*(c_root-c_tip);

% coordinates in body frame (will be used to plot wing geometry)
coords_inFb = [0.0 0.0;
    c_root 0.0;
    x_tipChord_TE half_wing_span;
    x_tipChord_LE half_wing_span];

% coordinates in matlab axis system
coords_inFm = fliplr(coords_inFb);


%% Design Lift coefficient

% required CL to maintain level flight at the given altitude and speed
CL_design = 2*Wo*grav/(rho_SL*Sref*V_cruise^2);

fprintf('Design Lift Coefficient (CL_design): %.4f \n', CL_design)


%% Lifting-line result

% Lifting-line number of segments
N = 1000;

% lifting line segment angles
theta = pi/2:-pi/(2*(N)):pi/(2*N);

% spanwise control point locations
z = half_wing_span*cos(theta);

% mean geometric chord for each segment
c_mean = c_tip.*(1+(half_wing_span-z).*(1-taper)./(taper.*half_wing_span));

% angle of attack for each segment
alpha_root = incidence;
alpha_tip = alpha_root + twist;
alpha = zeros(1,N);
for k=1:N
    slope = (alpha_tip-alpha_root)/half_wing_span;
    alpha(1,k) = alpha_root + slope*z(k);
end


% lift curve slope for each segment
Cl_aoa = Cl_aoa_func(alpha);
Cl = Cl_func(alpha);
mu = c_mean .* Cl_aoa / (4 * wing_span);
LHS = c_mean .* Cl / (4 * wing_span);


% Solving N equations to find coefficients A(i):
B = zeros(N,N);
for k=1:N
    for j=1:N
        B(k,j) = sin((2*j-1) * theta(k)) * (1 + (mu(k) * (2*j-1)) / sin(theta(k)));
    end
end

A = B\transpose(LHS);

sum1 = zeros(1,N);
sum2 = zeros(1,N);
for k = 1:N
    sum1(k) = 0;
    sum2(k) = 0;
    for j = 1 : N
        sum1(k) = sum1(k) + (2*j-1) * A(j)*sin((2*j-1)*theta(k));
        sum2(k) = sum2(k) + A(j)*sin((2*j-1)*theta(k));
    end
end
CL = 4*wing_span*sum2 ./ c_mean;
CL_dist = [CL, 0];
y_s = [z, wing_span/2];
CL_wing = pi * AR * A(1);


L_dist = zeros(1,N+1);
for i=1:N
    % compute lift per unit span
    dS = c_mean(i);
    L_dist(i) = 1/2*rho_SL*V_cruise^2*dS*CL_dist(i);
end

fprintf('Wing Lift Coefficient (CL_wing): %.4f \n', CL_wing)


%%  elliptical lift distribution

L0_elliptic = 4*Wo*grav/(pi*wing_span);
L_dist_elliptical = nan(size(y_s));
for i=2:length(y_s)
    
    a = 2*Wo*grav/pi/half_wing_span;
    delta_y = y_s(i) - y_s(i-1);
    L_dist_elliptical(i) = L0_elliptic*2/wing_span*sqrt(half_wing_span^2-y_s(i)^2);

end
L_dist_elliptical(1) = L_dist_elliptical(2);

% sum section lift over entire wing
dy = diff(y_s); dy = [dy(1) dy];
L_elliptical = 2*sum(L_dist_elliptical.*dy);

CL_elliptic = (2*L_elliptical)/(rho_SL*V_cruise^2*Sref);
fprintf('Elliptic Lift Coefficient (CL_elliptic): %.4f \n', CL_elliptic)


%% results summary

fprintf('\n ------------- Summary ----------------- \n')
fprintf(' Root chord [m]: %.2f \n', c_root)
fprintf(' Tip chord [m]: %.2f \n', c_tip)
fprintf(' Wing Span [m]: %.2f \n', half_wing_span*2)
fprintf('\n')


%% plot results

% plot wing geometry and lift distribution
fh1 = figure(1); clf
fh1.WindowState = 'maximized';

subplot(2,1,1)
hold on
daspect([1 1 1])%, set(gca,'xtick',[]), set(gca,'ytick',[])
plot([0 half_wing_span ],[c_root*3/4 c_root*3/4], '-- k', 'LineWidth', 2)
plot(coords_inFm([1 2], 1), coords_inFm([1 2], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([1 3], 1), coords_inFm([1 3], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([2 4], 1), coords_inFm([2 4], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([3 4], 1), coords_inFm([3 4], 2), '- k', 'LineWidth', 2)
ylabel('Chord [m]')

subplot(2,1,2)
hold on, grid on
plot(y_s, L_dist/1000, 'LineWidth', 2)
plot(y_s, L_dist_elliptical/1000, 'LineWidth', 2)
legend('actual','elliptical')

xlabel('Half wingspan [m]')
ylabel('Lift Dist. [kN/m]')

