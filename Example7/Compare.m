%% compare 
clc,clear
load Data2
load referenceData2
figure; % slice 0 
hold on;
plot(R_slice0,U_slice0,'LineWidth',1);
plot(R_slice0_30,U_slice0_30,'LineWidth',1);
plot(R_slice0_60,U_slice0_60,'LineWidth',1);
legend('reference','30 \times 30 mesh','60 \times 60 mesh'); 
hold off;

figure; % slice 1
hold on;
plot(R_slice1,U_slice1,'LineWidth',1);
plot(R_slice1_30,U_slice1_30,'LineWidth',1);
plot(R_slice1_60,U_slice1_60,'LineWidth',1);
legend('reference','30 \times 30 mesh','60 \times 60 mesh');
hold off;


figure; % slice 2
hold on;
plot(R_slice2,U_slice2,'LineWidth',1);
plot(R_slice2_30,U_slice2_30,'LineWidth',1);
plot(R_slice2_60,U_slice2_60,'LineWidth',1);
legend('reference','30 \times 30 mesh','60 \times 60 mesh');ylim([0.4,0.8]);
hold off;


figure; % slice 3
hold on;
plot(R_slice3,U_slice3,'LineWidth',1);
plot(R_slice3_30,U_slice3_30,'LineWidth',1);
plot(R_slice3_60,U_slice3_60,'LineWidth',1);
legend('reference','30 \times 30 mesh','60 \times 60 mesh');ylim([0.4,0.8]);
hold off;

figure; % slice 4
hold on;
plot(R_slice4,U_slice4,'LineWidth',1);
plot(R_slice4_30,U_slice4_30,'LineWidth',1);
plot(R_slice4_60,U_slice4_60,'LineWidth',1);
legend('reference','30 \times 30 mesh','60 \times 60 mesh');ylim([0.4,0.8]);
hold off;


%% compare 
clc,clear
load Data1
load referenceData1
figure; % slice 0 
hold on;
plot(R_slice0,U_slice0,'LineWidth',1);
plot(R_slice0_10,U_slice0_10,'LineWidth',1);
plot(R_slice0_30,U_slice0_30,'LineWidth',1);
plot(R_slice0_80,U_slice0_80,'LineWidth',1);
legend('reference','10 \times 10 mesh','20 \times 20 mesh','80 \times 80 mesh');
hold off;

figure; % slice 1
hold on;
plot(R_slice1,U_slice1,'LineWidth',1);
plot(R_slice1_10,U_slice1_10,'LineWidth',1);
plot(R_slice1_30,U_slice1_30,'LineWidth',1);
plot(R_slice1_80,U_slice1_80,'LineWidth',1);
legend('reference','10 \times 10 mesh','20 \times 20 mesh','80 \times 80 mesh');
hold off;


figure; % slice 2
hold on;
plot(R_slice2,U_slice2,'LineWidth',1);
plot(R_slice2_10,U_slice2_10,'LineWidth',1);
plot(R_slice2_30,U_slice2_30,'LineWidth',1);
plot(R_slice2_80,U_slice2_80,'LineWidth',1);
legend('reference','10 \times 10 mesh','20 \times 20 mesh','80 \times 80 mesh');
hold off;


figure; % slice 3
hold on;
plot(R_slice3,U_slice3,'LineWidth',1);
plot(R_slice3_10,U_slice3_10,'LineWidth',1);
plot(R_slice3_30,U_slice3_30,'LineWidth',1);
plot(R_slice3_80,U_slice3_80,'LineWidth',1);
legend('reference','10 \times 10 mesh','20 \times 20 mesh','80 \times 80 mesh');
hold off;

figure; % slice 4
hold on;
plot(R_slice4,U_slice4,'LineWidth',1);
plot(R_slice4_10,U_slice4_10,'LineWidth',1);
plot(R_slice4_30,U_slice4_30,'LineWidth',1);
plot(R_slice4_80,U_slice4_80,'LineWidth',1);
legend('reference','10 \times 10 mesh','20 \times 20 mesh','80 \times 80 mesh');
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

figure; % slice 1
hold on;
plot(R_slice1,U_slice1,'LineWidth',1);
plot(R_slice1_11,U_slice1_11,'LineWidth',1);
plot(R_slice1_21,U_slice1_21,'LineWidth',1);
plot(R_slice1_51,U_slice1_51,'LineWidth',1);
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

figure; % slice 4
hold on;
plot(R_slice4,U_slice4,'LineWidth',1);
plot(R_slice4_11,U_slice4_11,'LineWidth',1);
plot(R_slice4_21,U_slice4_21,'LineWidth',1);
plot(R_slice4_51,U_slice4_51,'LineWidth',1);
legend('reference','11 \times 11 mesh','21 \times 21 mesh','51 \times 51 mesh');
hold off;

