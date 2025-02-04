clc;clear;
rc = [4479.044956928165 -694.69214882738 5052.611874981442];
vc = [-0.9406229319779803 7.37792378051726 1.8541283129314754];

rt = [-1259.8546215159733 -5659.279979405802 -3473.4651254907144];
vt = [4.850181803531396 -3.859539589086098 4.536664669044564];

Ht = cross(rt,vt);
Hc = cross(rc,vc);
global mu
LOI = cross(Ht,Hc);
loi = LOI/norm(LOI);
LOI = 8500*loi;
LOI_opp = -8500*loi;
LOIS = [LOI;LOI_opp];
plot3(LOIS(:,1),LOIS(:,2),LOIS(:,3),'LineWidth',2,'Color','black');
hold on
rad = 6378;
[x,y,z]=sphere;
surf(rad*x,rad*y,rad*z);
map = earth_colormap();
colormap (map)
curve_t = animatedline('LineWidth',2 ,'Color','red');
curve_c = animatedline('LineWidth',2 ,'Color','blue');
set(gca,'XLim',[-8000 8000],'YLim',[-8000,8000],'ZLim', [-8000,8000])
view(145,28);
mu = 398600;
COE_C = coe_from_sv(rc,vc,mu);
a_c = COE_C(7);
Pc = 2*pi()/sqrt(mu)*a_c^(3/2);
time = 0:20:Pc+30;

RT = []; RC = [];
for i = 1:length(time)
    [R_C,V_C] = rv_from_r0v0(rc,vc,time(i));
    [R_T,V_T] = rv_from_r0v0(rt,vt,time(i));
    addpoints(curve_t, R_T(1),R_T(2),R_T(3))
    head_t = scatter3(R_T(1),R_T(2),R_T(3),'filled','red'); %% Target position plot
    RT = [RT; time(i), R_T(1), R_T(2), R_T(3), norm(R_T)];
    RC = [RC; time(i), R_C(1), R_C(2), R_C(3), norm(R_C)];
    addpoints(curve_c, R_C(1),R_C(2),R_C(3))
    head_c = scatter3(R_C(1),R_C(2),R_C(3),'filled','blue'); %% Chaser position plot (ISS)
    drawnow
    
    delete(head_t)
    delete(head_c)
end

COE_C = coe_from_sv(R_C,V_C,mu);
a_c = COE_C(7); T_c = 2*pi()/sqrt(mu)*a_c^(3/2);e_c = COE_C(2);
nu_c = COE_C(6);
E_c = 2*atan(sqrt((1-e_c)/(1+e_c))*tan(nu_c/2));
Me_c = E_c-e_c*sin(E_c);
time_c = 0.5*Me_c/pi()*T_c;
COE_T = coe_from_sv(R_T,V_T,mu);
a_t = COE_T(7); T_t = 2*pi()/sqrt(mu)*a_t^(3/2);e_t = COE_T(2);
nu_t = COE_T(6);
E_t = 2*atan2(sqrt(1-e_t)*sin(nu_t/2),sqrt(1+e_t)*cos(nu_t/2));
Me_t = E_t-e_t*sin(E_t);
time_t = 0.5*Me_t/pi()*T_t;

%% second revolution untill chaser is at its perigee
time2 = 1:20:T_c-time_c+20;
rc = R_C; vc = V_C;rt=R_T;vt=V_T;
for i = 1:length(time2)
    [R_C,V_C] = rv_from_r0v0(rc,vc,time2(i));
    [R_T,V_T] = rv_from_r0v0(rt,vt,time2(i));
    addpoints(curve_t, R_T(1),R_T(2),R_T(3))
    head_t = scatter3(R_T(1),R_T(2),R_T(3),'filled','red'); %% Target position plot
    addpoints(curve_c, R_C(1),R_C(2),R_C(3))
    head_c = scatter3(R_C(1),R_C(2),R_C(3),'filled','blue'); %% Chaser position plot (ISS)
    drawnow
    
    delete(head_t)
    delete(head_c)
end
rc = R_C; vc = V_C;rt=R_T;vt=V_T;
%% determining transfer orbit using lambert
COE_C = coe_from_sv(R_C,V_C,mu);
a_c = COE_C(7); T_c = 2*pi()/sqrt(mu)*a_c^(3/2);e_c = COE_C(2);
nu_c = COE_C(6);
E_c = 2*atan2(sqrt(1-e_c)*sin(nu_c/2),sqrt(1+e_c)*cos(nu_c/2));
Me_c = E_c-e_t*sin(E_c);
time_c = 0.5*Me_c/pi()*T_c;
COE_T = coe_from_sv(R_T,V_T,mu);
a_t = COE_T(7); T_t = 2*pi()/sqrt(mu)*a_t^(3/2);e_t = COE_T(2);
nu_t = COE_T(6);
E_t = 2*atan2(sqrt(1-e_t)*sin(nu_t/2),sqrt(1+e_t)*cos(nu_t/2));
Me_t = E_t-e_t*sin(E_t);
time_t = 0.5*Me_t/pi()*T_t;
delta_t = T_c/2;

delta = [];
for j = 1:size(RT,1)
    [R_C180,V_C180] = rv_from_r0v0(rc,vc,delta_t);
    delta = [delta norm(RT(j,2:4)-R_C180)];
    
end
pos = find(delta==min(delta));

% actually using lambert now
R1 = rc; R2 = RT(pos,2:4);
[V1,V2] = lambert(R1,R2,delta_t,'retro');
%% animating orbit change
rc = R1; vc = V1;rt=R_T;vt=V_T;
time_lambert = 1:20:delta_t;
for i = 1:length(time_lambert)
    [R_C,V_C] = rv_from_r0v0(rc,vc,time_lambert(i));
    [R_T,V_T] = rv_from_r0v0(rt,vt,time_lambert(i));
    addpoints(curve_t, R_T(1),R_T(2),R_T(3))
    head_t = scatter3(R_T(1),R_T(2),R_T(3),'filled','red'); %% Target position plot
    addpoints(curve_c, R_C(1),R_C(2),R_C(3))
    head_c = scatter3(R_C(1),R_C(2),R_C(3),'filled','blue'); %% Chaser position plot (ISS)
    drawnow
    
    delete(head_t)
    delete(head_c)
end
head_t = scatter3(R_T(1),R_T(2),R_T(3),'filled','red'); %% Target position plot
addpoints(curve_c, R_C(1),R_C(2),R_C(3))
head_c = scatter3(R_C(1),R_C(2),R_C(3),'filled','blue'); %% Chaser position plot (ISS)
drawnow