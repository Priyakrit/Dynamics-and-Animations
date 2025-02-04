clc; clear;

global mu
deg = pi/180;
mu = 398600;
%...Input data (angles in degrees):
e = 0.006;
a=6766.4;
h = sqrt(a*(1-e^2)*mu);
RA = 0;
incl = 51.6185;
w = 15.2076;
TA = 0;
%...
coe = [h, e, RA*deg, incl*deg, w*deg, TA*deg];
%...Algorithm 4.2 (requires angular elements be in radians):
[r0, v0] = sv_from_coe(coe);
period =  2*pi()/sqrt(mu)*a^(3/2);

%% Animating trajectory
rad = 6378;
[x,y,z]=sphere;
surf(rad*x,rad*y,rad*z); 
curve = animatedline('LineWidth',2 ,'Color','red');
set(gca,'XLim',[-8000 8000],'YLim',[-8000,8000],'ZLim', [-8000,8000])
view(67,22);
hold on



time = 1:30:5*period;
for i = 1:length(time)
    [R,V] = rv_from_r0v0(r0,v0,time(i));
    addpoints(curve, R(1),R(2),R(3))
    head = scatter3(R(1),R(2),R(3),'filled','red');
    drawnow
    pause(0.05);
    delete(head)

end