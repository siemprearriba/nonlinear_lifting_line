
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
clear
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


%% clean working directory

clear_working_directory;
if exist('auto_generated_aircraft.txt','file') == 2
    delete('auto_generated_aircraft.txt')
end


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
V_cruise = 26;                  % cruise velocity, [KCAS]
V_cruise = V_cruise/1.94;       % cruise velocity, [m/s

% aircraft/wing parameters
Wo = 12;                        % aircraft mass, [kg]
Sref = 1.7564;                  % reference wing area, [m^2]
AR = 13;                        % Aspect Ratio
wing_span = sqrt(Sref*AR);      % entire wing span, [m]
half_wing_span = wing_span/2.0;
inboard_half_span_ratio = 0.45;  % ratio: inboard_half_span/half_span
inboard_taper = 1.0;
outboard_taper = 0.6;

incidence = 4.0;                % wing incidence angle [deg]
inboard_twist = 0.0;            % inboard wing twist angle [deg]
outboard_twist = -6.0;          % outboard wing twist angle [deg]


% inboard wing portion
inboard_half_span = inboard_half_span_ratio*half_wing_span;

% outboard wing portion
outboard_half_span = half_wing_span - inboard_half_span;

% chord length
c_root = Sref/((1+inboard_taper)*inboard_half_span+(1+outboard_taper)*inboard_taper*outboard_half_span);
c_mid = c_root*inboard_taper;
c_tip = c_mid*outboard_taper;

% Wing area
inboard_half_wing_area = (c_root+c_mid)*inboard_half_span/2;
outboard_half_wing_area = (c_mid+c_tip)*outboard_half_span/2;

% element mean aerodynamic chord
inboard_mac = 2/3*c_root*(1+inboard_taper+inboard_taper^2)/(1+inboard_taper);
outboard_mac = 2/3*c_mid*(1+outboard_taper+outboard_taper^2)/(1+outboard_taper);
c_mean = (inboard_mac*inboard_half_wing_area + outboard_mac*outboard_half_wing_area)/(Sref/2);

% TE and LE locations
x_rootChord_LE = c_root;
x_rootChord_TE = 0.0;

x_midChord_LE = 3/4*c_root + c_mid/4;
x_midChord_TE = 3/4*(c_root - c_mid);

x_tipChord_LE = x_midChord_LE - (c_mid-c_tip)/4;
x_tipChord_TE = x_midChord_TE + 3/4*(c_mid-c_tip);

% wing area calculation check
S_calc = (c_root+c_mid)*inboard_half_span + (c_mid+c_tip)*outboard_half_span;
if abs(S_calc-Sref) > 1e-6
    fprintf('Wrong wing area configuration. Check chord length computations. \n')
    return
end

% coordinates in body frame
coords_inFb = [0.0 0.0;
    c_root 0.0;
    x_midChord_TE inboard_half_span;
    x_midChord_LE inboard_half_span;
    x_tipChord_TE inboard_half_span+outboard_half_span;
    x_tipChord_LE inboard_half_span+outboard_half_span];

% coordinates in matlab axis system
coords_inFm = fliplr(coords_inFb);


%% Design Lift coefficient

CL_design = 2*Wo*grav/(rho_SL*Sref*V_cruise^2);

fprintf('Design Lift Coefficient (CL_design): %.4f \n', CL_design)


%%  elliptical lift distribution

% elliptical lift distribution
iter = 1;
CL0 = 0.5; dCL = 0.01;
err = 1.0; err_norm = norm(err); tol = 5e-3;
while err_norm > tol
    
    CL0 = CL0 + sign(err)*dCL;
    
    y = linspace(0,half_wing_span, 300);
    CL_dist_elliptical = CL0*sqrt(1-(y/half_wing_span).^2);
    
    dy = y(2)-y(1);
    CL_elliptical = sum(CL_dist_elliptical*dy)/half_wing_span;
    
    err = (CL_design-CL_elliptical);
    err_norm = norm(err);
    iter = iter + 1;
    if iter > 500
        break
    end
end

% plot wing geometry and lift distribution
fh = figure(1); clf
fh.WindowState = 'maximized';

subplot(3,1,1)
hold on
daspect([1 1 1])%, set(gca,'xtick',[]), set(gca,'ytick',[])
plot([0 half_wing_span ],[c_root*3/4 c_root*3/4], '-- k', 'LineWidth', 2)
plot(coords_inFm([1 2], 1), coords_inFm([1 2], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([1 3], 1), coords_inFm([1 3], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([2 4], 1), coords_inFm([2 4], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([3 4], 1), coords_inFm([3 4], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([4 6], 1), coords_inFm([4 6], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([3 5], 1), coords_inFm([3 5], 2), '- k', 'LineWidth', 2)
plot(coords_inFm([5 6], 1), coords_inFm([5 6], 2), '- k', 'LineWidth', 2)

subplot(3,1,3)
hold on, grid on
plot(y, CL_dist_elliptical, '--r','LineWidth',2)
xlabel('Half wingspan [m]')

fprintf('Elliptic Lift Coefficient (CL_elliptic): %.4f \n', CL_elliptical)


%% Lifting-line result

% Lifting-line number of segments
N = 50;

% lifting line segment angles
theta = pi/2:-pi/(2*(N)):pi/(2*N);

% spanwise control point locations
z = half_wing_span*cos(theta);

figure(1)
subplot(3,1,1)
for k=1:N
    plot([z(k) z(k)], [0 c_root*1.1], '-.r')
end


% inboard and outboard segment separation
is_inboard_segment = (z<=inboard_half_span);
is_outboard_segment = not(is_inboard_segment);

% unified wing geometric variables
taper_combined = is_inboard_segment*inboard_taper + is_outboard_segment*outboard_taper;
c_root_combined = is_inboard_segment*c_root + is_outboard_segment*c_mid;
c_tip_combined = is_inboard_segment*c_mid + is_outboard_segment*c_tip;
half_span_combined = is_inboard_segment*inboard_half_span + is_outboard_segment*outboard_half_span;
y_k = is_inboard_segment.*z + is_outboard_segment.*(z-inboard_half_span);


% mean geometric chord for each segment
c_mean = c_tip_combined.*(1+(half_span_combined-y_k).*(1-taper_combined)./(taper_combined.*half_span_combined));

% angle of attack for each segment
alpha_root = incidence;
alpha_mid = alpha_root + inboard_twist;
alpha_tip = alpha_mid + outboard_twist;
alpha = zeros(1,N);
for k=1:N
    if z(k) <= inboard_half_span
        slope = (alpha_mid-alpha_root)/inboard_half_span;
        alpha(1,k) = alpha_root + slope*z(k);
        last_inboard_idx = k;
    else
        slope = (alpha_tip-alpha(1,last_inboard_idx))/outboard_half_span;
        alpha(1,k) = alpha(1,last_inboard_idx) + slope*(z(k)-inboard_half_span);
    end
end

% plot segment angle of attack
figure(1)
subplot(3,1,2)
hold on, grid on
plot(z, alpha, 'linewidth', 2, 'color', 'k')
ylabel('Section AoA [deg]')

% lift curve slope for each segmen
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


subplot(3,1,3)
hold on
plot(y_s, CL_dist, '-', 'linewidth', 2, 'color', 'k')
ylabel('Section $C_L$')
legend('Ideal Elliptic Distribution','Actual Lift Distribution')
CL_wing = pi * AR * A(1);

fprintf('Wing Lift Coefficient (CL_wing): %.4f \n', CL_wing)


%% results summary

fprintf('\n ------------- Summary ----------------- \n')
fprintf(' Root chord [m]: %.3f \n', c_root)
fprintf(' Mid chord [m]: %.3f \n', c_mid)
fprintf(' Tip chord [m]: %.3f \n', c_tip)
fprintf(' Wing Span [m]: %.2f \n', half_wing_span*2)
fprintf('\n')


subplot(3,1,1)
title(['Number of Elements: ', num2str(N), ', Root Chord: ', num2str(round(c_root,2)), '[m], Mid Chord: ', ...
    num2str(round(c_mid,2)), '[m], Tip Chord: ', num2str(round(c_tip,2)), '[m], AR:', num2str(AR)])

