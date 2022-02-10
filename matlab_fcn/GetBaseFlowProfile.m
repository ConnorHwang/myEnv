%% Baseflow Coefficient from Borronin(2012)
function [fcn0, fcn1, fcn2] = GetBaseFlowProfile(a1,r1,d1,b1,m,d2,BaseFlowType)
%--------------------------------------------------------------------------
% Returns the baseflow profile in syms type equation followin
% Borronin(2012) paper
% Input : ...
% Output: Baseflow functions
%--------------------------------------------------------------------------
global radius

if(BaseFlowType == 3)
    % Note that Delta_g should be determined according to the viscosity
    % ratio of the gas and the liquid phase.
%     a1 = 0.6; r1 = 0; d1 = 0.7; b1 = 0.4; d2 = 0.25; m = 0.018;
    syms r2 radius
    fcn0 = 0; % dummy output;
    fcn1 = a1*exp(-(radius-r1)^2/d1^2)+b1;

    a2 = (a1*exp(-(1-r1)^2/d1^2)+b1)/exp(-(1-r2)^2/d2^2);
    % To get r2
    eq_r2 = a1/d1^2*(1-r1)*exp(-(1-r1)^2/d1^2) == m*a2/d2^2*(1-r2)*exp(-(1-r2)^2/d2^2);
    ans_r2 = double(solve(eq_r2,r2));

    a2_eq2 = double(subs(a2,ans_r2));
    fcn2 = a2_eq2*exp(-(radius-ans_r2)^2/d2^2);

elseif(BaseFlowType == 4)
%     a1 = 0.4; r1 = 0.5; d1 = 0.4; b1 = 0.6; d2 = 0.25; m = 0.018;
    syms r2 radius
    fcn0 = a1*exp(-(radius-r1)^2/d1^2)+b1;   
    fcn1 = 1*heaviside(-(radius-r1))+heaviside(radius-r1)*(a1*exp(-(radius-r1)^2/d1^2)+b1);
    
    a2 = (a1*exp(-(1-r1)^2/d1^2)+b1)/exp(-(1-r2)^2/d2^2);
    eq_r2 = a1/d1^2*(1-r1)*exp(-(1-r1)^2/d1^2) == m*a2/d2^2*(1-r2)*exp(-(1-r2)^2/d2^2);
    ans_r2 = double(solve(eq_r2,r2));

    a2_eq2 = double(subs(a2,ans_r2));
    fcn2 = a2_eq2*exp(-(radius-ans_r2)^2/d2^2);

end