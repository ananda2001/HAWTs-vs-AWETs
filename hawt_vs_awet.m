clc
close all
clear all



%%%%%%%%%%%%%%%% wind profile(default) :
rho = 1.225 ; % air density in kg/m³
K = [1 3] ; % air form factor
for i = 1:2 
v_mean = 15 ; % maximum wind-speed in m/s
k = K(i) ;
a = v_mean./gamma(1./k+1) ; % weibull speed parameter
v_w = 0:1:50; % wind speed range
fv(i,:) = (k/a).*(v_w./a).^(k-1).*exp(-(v_w./a).^k);  % probability distribution of wind
end

% turbine parameters (physical) :
R = 150/2 ; % turbine blade radius
H = 5.7267 ; % turbine inertia constant
A = pi*R^2 ; % area of wind swept by blades
v_c = 5 ; % wind cut_in speed (33% of turbine speed) in m/s
v_r = 25 ; % rated wind speed in m/s
v_t = (0:5:v_r)' ; % turbine blade tip speed in m/s
LD_k = 0.5*rho*mean((v_t-v_w).^2) ; % lift-drag parameter
L = 100 ; % maximum lift force needed in N
c_l = L./LD_k ; % lift co-efficient
D = 100 ; % maximum Drag force needed in N
c_d = D./LD_k ; % drag co-efficient
v_ratio = v_t/max(v_w) ;
    
%%%%% PMSG wind turbine characteristic :
pm_pmsg = 30 ; % power rating in Kw
c_p = c_l.sqrt(1+v_ratio.^2).(v_ratio - (c_d./c_l).*v_ratio.^2) ; % power co-efficient   
pl = 0.5*rho*c_p*A ;
pow_max = 1e-6 * (pl*k*a^3/(k+2))*gamma(1/(k+2)+1); % maximum power output in Kw
pow_pmsg = pow_max.*(cumsum(0.5*v_w.*fv(2,:)))/1.5 ;
pow_pmsg(7,:) = pow_pmsg(2,:);
pow_pmsg(2,:) = [];
%pow_pmsg(6,:) = pow_pmsg(2,:);
subplot(221);
plot(v_w,pow_pmsg);
grid on; title('PMSG characteristic')
xlabel('turbine speed in m/s');
ylabel('PMSG output power in kW'); 
leg1 = legend(["6 m/s" "8 m/s" "10 m/s" "12 m/s" "15 m/s"]);
title(leg1, 'cut-off wind speed')

%%%%%%%%%%%%%%% DFIG turbine characteristic :
pm_dfig = 30 ; % power rating in Kw
c_p = c_l.sqrt(1+v_ratio.^2).(v_ratio-(c_d./c_l).*v_ratio.^2) ; % power co-efficient   
pl = rho.c_p.*A.(c_l./c_d).^2 ;
pow_max = 1e-6 * (pl*k*a^3/(k+2))*gamma(1/(k+2)+1) ; % maximum power output in Kw
pow_dfig = pow_max.*(cumsum(0.5*v_w.*fv(1,:)))/1.8 ; % output power distribution
pow_dfig(7,:) = pow_dfig(2,:);
pow_dfig(2,:) = [];
subplot(222); plot(v_w,pow_dfig);
grid on; title('DFIG characteristic')
xlabel('turbine speed in m/s');
ylabel('DFIG output power in kW'); 
leg2 = legend(["6 m/s" "8 m/s" "10 m/s" "12 m/s" "15 m/s"]);
title(leg2, 'cut-off wind speed')



%%%%%%%%%%%%% Levelized Cost Of Energy(LCOE) analysis :
icc_hawt = 10 ; % initial capital costs for HWAT 
crf_hawt = 0.37 ; % capital recovery factor for HWAT
omc_hawt = 1000 ; %operation & maintenance cost for HWAT
incst_hawt = 1000 ; %initial foundation cost
omc_pmsg = 2000 ; %operation & maintenance cost for HWAT
omc_dfig = 1500 ; %operation & maintenance cost for HWAT
Pow_pmsg = max(pow_pmsg') ;
Pow_pmsg(1)=[] ;
Pow_dfig = max(pow_dfig') ;
Pow_dfig(1)=[] ;
lcoe_hawt_pmsg = ((icc_hawt*crf_hawt)+omc_hawt+omc_pmsg)./(Pow_pmsg) ; 
lcoe_hawt_dfig = ((icc_hawt*crf_hawt)+omc_hawt+omc_dfig)./(Pow_dfig) ; 

icc_awe = 1 ; % initial capital costs for AWE
crf_awe = 0.79 ; % capital recovery factor for AWE
omc_awe = 300 ; %operation & maintenance cost for AWE
incst_awe = 800 ; %initial foundation cost
lcoe_awe_pmsg = ((icc_awe*crf_awe)+omc_awe+omc_pmsg)./(Pow_pmsg) ; 
lcoe_awe_dfig = ((icc_awe*crf_awe)+omc_awe+omc_dfig)./(Pow_dfig) ; 

subplot(2,2,3), plot(incst_hawt-lcoe_hawt_pmsg), hold on
plot(incst_hawt-lcoe_hawt_dfig), grid on
title('LCOE comparison for HAWT'), legend(["PMSG" "DFIG"])
xlabel('year'), ylabel('energy cost per kWh')


subplot(224)
plot(incst_awe-lcoe_awe_pmsg), hold on, plot(incst_awe-lcoe_awe_dfig), grid on
title('LCOE comparison for AWE'), legend('PMSG','DFIG')