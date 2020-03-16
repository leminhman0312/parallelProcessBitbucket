clear all;
clc;

nprocs = [1,2,4,8,16];

speed_ori_fixed_n = [1,1.98126185797461,3.80454003342509,7.20914512068038,8.87125262655941];


speed_mod_fixed_n = [1,1.93619775382304,3.71815393891202,7.10188706058281,8.70724688522684];


f1 = figure(1)
plot(nprocs,speed_ori_fixed_n);
title("PROCESSORS VS. SPEEDUP ORIGINAL CODE, FIXED N");
xlabel("NUMBER OF PROCESSORS");
ylabel("SPEEDUP")
set(gca,"FontSize", 14)
print(f1,"fix_ori.png")


f2 = figure(2)
plot(nprocs,speed_mod_fixed_n);
title("PROCESSORS VS. SPEEDUP MODIFIED CODE, FIXED N");
xlabel("NUMBER OF PROCESSORS");
ylabel("SPEEDUP")
set(gca,"FontSize", 14)
print(f2,"fix_mod.png")


f3 = figure(3)
plot(nprocs,speed_ori_fixed_n, '-xr','markersize', 15);
hold on
plot(nprocs,speed_mod_fixed_n,'-sb','markersize',15);
title("PROCESSORS VS. SPEEDUP FIXED N");
xlabel("NUMBER OF PROCESSORS");
ylabel("SPEEDUP")
h = legend({"ORIGINAL","MODIFIED"});      
set (h, "fontsize", 16,"Location","southeast");
set(gca,"FontSize", 14)
print(f3,"fix_combine.png")




scale_ori_varied_n = [1,0.970144162459363,0.989720847556712,1.00140612481951,0.675045763270162];
scale_mod_varied_n = [1,0.969462666761537,0.977430603735728,0.98727591040333,0.665274972526809];


f4 = figure(4)
plot(nprocs,scale_ori_varied_n);
title("PROCESSORS VS. SCALED EFFICIENCIES ORIGINAL CODE, VARIED N");
xlabel("NUMBER OF PROCESSORS");
ylabel("SCALE EFFIENCIES")
set(gca,"FontSize", 14)
print(f4,"varied_ori.png")

f5 = figure(5)
plot(nprocs,scale_mod_varied_n);
title("PROCESSORS VS. SCALED EFFICIENCIES MODIFIED CODE, VARIED N");
xlabel("NUMBER OF PROCESSORS");
ylabel("SCALE EFFIENCIES")
set(gca,"FontSize", 14)
print(f5,"varied_mod.png")


f6 = figure(6)
plot(nprocs,scale_ori_varied_n, '-xr','markersize', 15);
hold on
plot(nprocs,scale_mod_varied_n,'-sb','markersize',15);
title("PROCESSORS VS. SCALED EFFICIENCIES VARIED N");
xlabel("NUMBER OF PROCESSORS");
ylabel("SCALE EFFIENCIES")
h = legend({"ORIGINAL","MODIFIED"});      
set (h, "fontsize", 16);
set(gca,"FontSize", 14)
print(f6,"varied_combine.png")






