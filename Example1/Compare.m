clc,clear
load Data
load referenceData
%%
figure; % slice 0 
hold on;
plot(R_slice0,U_slice0,'LineWidth',1);
plot(R_slice0_10,U_slice0_10,'LineWidth',1);
plot(R_slice0_20,U_slice0_20,'LineWidth',1);
plot(R_slice0_50,U_slice0_50,'LineWidth',1);
legend('reference','10 \times 10 mesh','20 \times 20 mesh','50 \times 50 mesh');
hold off;



figure; % slice 2
hold on;
plot(R_slice2,U_slice2,'LineWidth',1);
plot(R_slice2_10,U_slice2_10,'LineWidth',1);
plot(R_slice2_20,U_slice2_20,'LineWidth',1);
plot(R_slice2_50,U_slice2_50,'LineWidth',1);
legend('reference','10 \times 10 mesh','20 \times 20 mesh','50 \times 50 mesh');
hold off;


figure; % slice 3
hold on;
plot(R_slice3,U_slice3,'LineWidth',1);
plot(R_slice3_10,U_slice3_10,'LineWidth',1);
plot(R_slice3_20,U_slice3_20,'LineWidth',1);
plot(R_slice3_50,U_slice3_50,'LineWidth',1);
legend('reference','10 \times 10 mesh','20 \times 20 mesh','50 \times 50 mesh');
hold off;


%%

figure; % slice 0 
hold on;
plot(R_slice0,U_slice0,'LineWidth',1);
plot(R_slice0_11,U_slice0_11,'LineWidth',1);
plot(R_slice0_21,U_slice0_21,'LineWidth',1);
plot(R_slice0_51,U_slice0_51,'LineWidth',1);
legend('reference','11 \times 11 mesh','21 \times 21 mesh','51 \times 51 mesh');
hold off;



figure; % slice 2
hold on;
plot(R_slice2,U_slice2,'LineWidth',1);
plot(R_slice2_11,U_slice2_11,'LineWidth',1);
plot(R_slice2_21,U_slice2_21,'LineWidth',1);
plot(R_slice2_51,U_slice2_51,'LineWidth',1);
legend('reference','11 \times 11 mesh','21 \times 21 mesh','51 \times 51 mesh');
hold off;


figure; % slice 3
hold on;
plot(R_slice3,U_slice3,'LineWidth',1);
plot(R_slice3_11,U_slice3_11,'LineWidth',1);
plot(R_slice3_21,U_slice3_21,'LineWidth',1);
plot(R_slice3_51,U_slice3_51,'LineWidth',1);
legend('reference','11 \times 11 mesh','21 \times 21 mesh','51 \times 51 mesh');
hold off;

