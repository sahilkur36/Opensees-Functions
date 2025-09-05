%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
clear all
close all
clc

orientation = [];
fc = 6000*1.5;



%% 4story
nStories = 4;
h = 216; % inches  % 264  336  648
b = 16; % inches 
n_bars = 3;
a_bars = 1.56;
n_layers = 10;
l_be = 51;
cc = 3; % clear cover

%% 8story
% nStories = 8;
% b = 20;
% h = 264; % 264  336  648
% n_bars = 4;
% a_bars = 1.56;
% n_layers = 11;
% l_be = 56;
% cc = 3; % clear cover

%% 12story
% nStories = 12;
% h = 288; % inches  % 264  336  648
% b = 22; % inches 
% n_bars = 4;
% a_bars = 1.56;
% n_layers = 12;
% l_be = 61;
% cc = 3; % clear cover

% d = 261;
% As = 4*1.56*ones(1,22);
% As_d = [3+(0:5:50)  264-flip(3+(0:5:50))];
% As = [4*1.56*ones(1,11), 2*0.31*ones(1,12), 4*1.56*ones(1,11)];
% As_d = [3+(0:5:50), 54+(12:12:144), 264-flip(3+(0:5:50))];

% As = [4*1.56*ones(1,12), 4*1.56*ones(1,12)];
% As_d = [3+(0:5:55), 336-flip(3+(0:5:55))];
% As = [4*1.56*ones(1,12), 2*0.31*ones(1,18), 4*1.56*ones(1,12)];
% As_d = [3+(0:5:55), 66+(0:12:204), 336-flip(3+(0:5:55))];

% As = [4*1.56*ones(1,16), 4*1.56*ones(1,16)];
% As_d = [3+(0:5:75),  h-flip(3+(0:5:75))];
% As = [4*1.56*ones(1,16), 2*0.31*ones(1,40), 4*1.56*ones(1,16)];
% As_d = [3+(0:5:75), 78+(12:12:480), h-flip(3+(0:5:75))];

% % 
% b = 22;
% h = 648; % 288  384  648
% As = 4*1.56*ones(1,34);
% % As_d = [3+(0:5:55)  h-flip(3+(0:5:55))];
% As_d = [3+(0:5:80),  h-flip(3+(0:5:80))];

% b = 22;
% h = 456; % 456  
% As = 3*1.56*ones(1,20);
% % As_d = [3+(0:5:55)  h-flip(3+(0:5:55))];
% As_d = [3+(0:7:63),  h-flip(3+(0:7:63))];

s = round(4*(l_be - cc)/n_layers)/4; % round to the nearest 0.25"
A_be = 3+(1:n_layers)*s;
As = n_bars*a_bars*ones(1,2*n_layers);
As_d = [A_be,  h-flip(A_be)];

P = 1000*(59.4*(nStories - 1)+57.71+22.5*(nStories - 1)+9);

fy = 60000*1.25;
Es = 27000000;
% P = 534000*20/8;
slab_depth = 0;
b_eff = 0;


%% Assumptions
% 1. Whitney stress block
% 2. for Phi assumes columns are tied and not spiral

%% Import Packages
import aci_318.*

%% Inital Setup
% As = str2double(strsplit(strrep(strrep(As{1},'[',''),']',''),','));
% As_d = str2double(strsplit(strrep(strrep(As_d{1},'[',''),']',''),','));
% if strcmp(orientation,'neg')
%     As = fliplr(As);
%     As_d = h - fliplr(As_d);
% elseif strcmp(orientation,'oop') % out of plane bending
%     As = [sum(As)/2, sum(As)/2]; % Assume half the steel on each end
%     As_d = [min(3,h/5), max(h-3,h-h/5)]; % Assume location is min of three inches or 1/5 width from the sides
% end
[ beta_1 ] = fn_beta_1( fc );

%% Begin Method
% Find Location of Neutral Axis
y_prev = -h/2;
step = 0.1;
tolerance = 50; % lbs %(0.85*fc*b*0.85*d/2+abs(P))/5000;
balance_found = 0;
count = 0;
while balance_found == 0
    count = count + 1;
    y(count) = y_prev + step;
    y_prev = y(count);
    if y(count) >= h/2
        warning('Nuetral Axis of Concrete Section Not Found, Assume NA is as if there is no compression steel')
        c = As(end)*fy/(beta_1*0.85*fc*b);
        e_s = 0.003*(c-As_d(end))/c; % Positive = compression strain, negative = tension strain
        fy_eff = fy - 0.85*fc; % Effective compresison steel stress which accounts for displaced concrete.
        fs = max(min(e_s*Es,fy_eff),-fy); % Check both positive and negative yield stress
        break
    end
    c = h/2-y(count);
    [ balance_eq(count), fs, e_s ] = fn_calulate_bending_balance( c, P, As, As_d, b, b_eff, slab_depth, fy, Es, fc, beta_1 );
    if abs(balance_eq(count)) < tolerance
        balance_found = 1;
    elseif count > 1 && (sign(balance_eq(count)) ~= sign(balance_eq(count-1)))
        step = 0.00002;
        switch_index = count;
        while balance_found == 0
            count = count + 1;
            y(count) = y_prev - step;
            y_prev = y(count);
            if count > 5000 + switch_index
                warning('Nuetral Axis of Concrete Section Not Found, Assume NA is as if there is no compression steel')
                c = As(end)*fy/(beta_1*0.85*fc*b);
                e_s = 0.003*(c-As_d)/c; % Positive = compression strain, negative = tension strain
                fy_eff = fy - 0.85*fc; % Effective compresison steel stress which accounts for displaced concrete.
                fs = max(min(e_s*Es,fy_eff),-fy); % Check both positive and negative yield stress
                break
            end
            c = h/2-y(count);
            [ balance_eq(count), fs, e_s ] = fn_calulate_bending_balance( c, P, As, As_d, b, b_eff, slab_depth, fy, Es, fc, beta_1 );
            if abs(balance_eq(count)) < tolerance
                balance_found = 1;
            end
        end
    end
end

% Nominal Moment Capacity
a = beta_1*c;
if slab_depth > 0 % T-beam Section
    if a > slab_depth
        Mn = sum(As.*fs.*(h/2 - As_d)) + (0.85*fc*b_eff*slab_depth)*(h-slab_depth)/2 + (0.85*fc*b*(a-slab_depth))*((h - 3*a - slab_depth)/2);
    else
        Mn = sum(As.*fs.*(h/2 - As_d)) + (0.85*fc*b_eff*a)*(h-a)/2;
    end
else % Rectangular Section
    Mn = sum(As.*fs.*(h/2 - As_d)) + (0.85*fc*b*a)*(h-a)/2;
end

%% Phi and Design Capacity
if min(e_s) < -0.005 % check largest tension strain (tension is negative)
    phi = 0.9;
else
    phi = max(min((abs(min(e_s))-0.002)*(0.9-0.65)/(0.005-0.002)+0.65,0.9),0.65); % Interpolate for phi
end
Mu = phi*Mn;

%% Validation Plots 
if 1 == 1
    plot(y,balance_eq)
    grid on
    xlabel('Distance from Center')
    ylabel('Balance Equation')
end


c_approx = h*(0.1+1.2*P/(fc*h*b));

Mn/1000