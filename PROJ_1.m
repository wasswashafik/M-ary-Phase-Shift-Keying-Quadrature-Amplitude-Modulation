clc; close all; 
%################---Set Values for the following----######################%
Cycle_Count = 25;
Var_Start = 0.01; 
Var_Res = 0.1;
Var_End = 1;   
A_off = 1.00; P_off = 0*(pi/180); %P_off in degrees (Amp-1.01,1.05;Phase 1 or 5)

%**********************Modulation Schemes in MPSK*************************%
%-----@@@@@@@@@@@@@----------BPSK---@@@@@@@@@@@@@@@@----------%
K = 1; M = 2;% 10log(M=2)
%Pre-allocate memory
Bit_Err_Rate = zeros(1,1000);Sym_Err_Rate = zeros(1,1000);
Norm_Bit_Energy = zeros(1,1000);Sym_Err_Rate_BPSK_Theory = zeros(1,1000);
Bit_Err_Rate_BPSK_Theory = zeros(1,1000);
j=0; %temp variable

for noise_var = Var_Start : Var_Res : Var_End
    j=j+1;
    [Bit_Err_Rate(j),Sym_Err_Rate(j)]= RF_Transceiver_Fn_BPSK(noise_var,Cycle_Count,A_off,P_off);
    Norm_Bit_Energy(j) = (1/(2*K*noise_var));%(Es/(2*K*noise_var));
    %Theoretical Value for comparison
    Bit_Err_Rate_BPSK_Theory(j) = qfunc(sqrt(1/noise_var));
    Sym_Err_Rate_BPSK_Theory(j) = 2*qfunc(sqrt(1/noise_var) * sin(pi/M));
end
Bit_Err_Rate_BPSK = Bit_Err_Rate;
Sym_Err_Rate_BPSK = Sym_Err_Rate;
Norm_Bit_Energy_BPSK_dB = 10.*log10(Norm_Bit_Energy);

%------@@@@@@@@@@@@@---------QPSK----@@@@@@@@@@@@@@@@---------%
M = 4; K = 2; % 10log(M=4)
%Pre-allocate memory
Bit_Err_Rate = zeros(1,1000);Sym_Err_Rate = zeros(1,1000);
Norm_Bit_Energy = zeros(1,1000);Sym_Err_Rate_QPSK_Theory = zeros(1,1000);
j=0; %temp variable

for noise_var = Var_Start : Var_Res : Var_End
    j=j+1;
    [Bit_Err_Rate(j),Sym_Err_Rate(j)]= RF_Transceiver_Fn_QPSK(noise_var,Cycle_Count,A_off,P_off);
    Norm_Bit_Energy(j) = (1/(2*K*noise_var));%(Es/(2*K*noise_var));
    %Theoretical Value for comparison
    Sym_Err_Rate_QPSK_Theory(j) = 2*qfunc(sqrt(1/noise_var) * sin(pi/M));
end
Bit_Err_Rate_QPSK = Bit_Err_Rate;
Sym_Err_Rate_QPSK = Sym_Err_Rate;
Norm_Bit_Energy_QPSK_dB = 10.*log10(Norm_Bit_Energy);              

%------@@@@@@@@@@@---------8PSK------@@@@@@@@@@@@@@@-------%
K = 3; M = 8;% 10log(M=8)
%Pre-allocate memory
Bit_Err_Rate = zeros(1,1000);Sym_Err_Rate = zeros(1,1000);
Norm_Bit_Energy = zeros(1,1000);Sym_Err_Rate_8PSK_Theory = zeros(1,1000);
j=0; %temp variable

for noise_var = Var_Start : Var_Res : Var_End
    j=j+1;
    [Bit_Err_Rate(j),Sym_Err_Rate(j)]= RF_Transceiver_Fn_8PSK(noise_var,Cycle_Count,A_off,P_off);
    Norm_Bit_Energy(j) = (1/(2*K*noise_var));%(Es/(2*K*noise_var));
    %Theoretical Value for comparison
    Sym_Err_Rate_8PSK_Theory(j) = 2*qfunc(sqrt(1/noise_var) * sin(pi/M));
end
Bit_Err_Rate_8PSK = Bit_Err_Rate;
Sym_Err_Rate_8PSK = Sym_Err_Rate;
Norm_Bit_Energy_8PSK_dB = 10.*log10(Norm_Bit_Energy);

%-------@@@@@@@@@@@@--------16PSK-------@@@@@@@@@@@@@------%
K = 4; M= 16; % 10log(M=16)
%Pre-allocate memory
Bit_Err_Rate = zeros(1,1000);Sym_Err_Rate = zeros(1,1000);
Norm_Bit_Energy = zeros(1,1000);Sym_Err_Rate_16PSK_Theory = zeros(1,1000);
j=0; %temp variable

for noise_var = Var_Start : Var_Res : Var_End
    j=j+1;
    [Bit_Err_Rate(j),Sym_Err_Rate(j)]= RF_Transceiver_Fn_16PSK(noise_var,Cycle_Count,A_off,P_off);
    Norm_Bit_Energy(j) = (1/(2*K*noise_var));%(Es/(2*K*noise_var));
    %Theoretical Value for comparison
    Sym_Err_Rate_16PSK_Theory(j) = 2*qfunc(sqrt(1/noise_var) * sin(pi/M));
end
Bit_Err_Rate_16PSK = Bit_Err_Rate;
Sym_Err_Rate_16PSK = Sym_Err_Rate;
Norm_Bit_Energy_16PSK_dB = 10.*log10(Norm_Bit_Energy);


%#######################---------Plots------##############################%
figure('Name','Bit Error Rate  vs  Bit Energy','NumberTitle','off');
hold on; grid on;
set(gca, 'YScale', 'log')
plot(Norm_Bit_Energy_BPSK_dB,Bit_Err_Rate_BPSK_Theory,'r--o','LineWidth',1);
plot(Norm_Bit_Energy_BPSK_dB,Bit_Err_Rate_BPSK,'-r','LineWidth',2);
plot(Norm_Bit_Energy_QPSK_dB,Bit_Err_Rate_QPSK,'-b','LineWidth',2);
plot(Norm_Bit_Energy_8PSK_dB,Bit_Err_Rate_8PSK,'-g','LineWidth',2);
plot(Norm_Bit_Energy_16PSK_dB,Bit_Err_Rate_16PSK,'-k','LineWidth',2);
xlabel('Bit Energy(Eb/No) in dB)'); ylabel('Bit Error Rate');
legend('Theor BPSK','BPSK','QPSK','8PSK','16PSK');
xlim([0 15]); ylim([10^-5 1]);
hold off;

figure('Name','Symbol Error Rate  vs  Bit Energy','NumberTitle','off');
hold on; grid on;
set(gca, 'YScale', 'log')
plot(Norm_Bit_Energy_BPSK_dB,Sym_Err_Rate_BPSK,'-r','LineWidth',2);
plot(Norm_Bit_Energy_QPSK_dB,Sym_Err_Rate_QPSK,'-b','LineWidth',2);
plot(Norm_Bit_Energy_8PSK_dB,Sym_Err_Rate_8PSK,'-g','LineWidth',2);
plot(Norm_Bit_Energy_16PSK_dB,Sym_Err_Rate_16PSK,'-k','LineWidth',2);
%Theoretical Values
plot(Norm_Bit_Energy_BPSK_dB,Sym_Err_Rate_BPSK_Theory,'c--o','LineWidth',1);
plot(Norm_Bit_Energy_QPSK_dB,Sym_Err_Rate_QPSK_Theory,'b--o','LineWidth',1);
plot(Norm_Bit_Energy_8PSK_dB,Sym_Err_Rate_8PSK_Theory,'g--o','LineWidth',1);
plot(Norm_Bit_Energy_16PSK_dB,Sym_Err_Rate_16PSK_Theory,'k--o','LineWidth',1);
legend('BPSK','QPSK','8PSK','16PSK','Theor-BPSK','Theor-QPSK','Theor-8PSK','Theor-16PSK');
xlabel('Bit Energy(Eb/No)in dB'); ylabel('Symbol Error Rate');
xlim([0 15]); ylim([10^-5 1]);
hold off;

%#######################----------Functions--------######################%
%--------------------********--BPSK----**********-------------------------%
function [Bit_Err_Rate,Sym_Err_Rate]= RF_Transceiver_Fn_BPSK(noise_var,Cycle_Count,A_off,P_off)
%---------Initialisation----------------% 
%Vector Constellation for BPSK; Normalised -> sqrt(Es=1)
%Reference vectors with Amp & Phase offsets
Ref_vect1_wOff = [A_off*cos(0+P_off),A_off*sin(0+P_off)]; 
Ref_vect2_wOff = [A_off*cos(pi+P_off),A_off*sin(pi+P_off)]; 
%Reference vectors without Amp & Phase offsets
A_off = 1; P_off = 0;
Ref_vect1 = [A_off*cos(0+P_off),A_off*sin(0+P_off)]; 
Ref_vect2 = [A_off*cos(pi+P_off),A_off*sin(pi+P_off)]; 

%Bit Definition for 8PSK - Gray Codes
Ref_bits1 = [0,0];Ref_bits2 = [1,0];

%Counters/Variables
BitError_Counter = 0; SymError_Counter = 0;
Bit_Counter= 0; Sym_Counter = 0; i=1; 

%Preallocate
Scatter_x = zeros(1,Cycle_Count);
Scatter_y = zeros(1,Cycle_Count);

%---Start of Transceiver chain for a given sigma value----%
%------------TXR-------------------%
for Tx_Count = 1:Cycle_Count
x = rand;
if(x>0 && x <= 0.5)
    Tx_Vect = Ref_vect1_wOff; Tx_Bits = Ref_bits1;
elseif(x>0.5 && x<= 1)
    Tx_Vect = Ref_vect2_wOff; Tx_Bits = Ref_bits2;
end

I = Tx_Vect(1); Q = Tx_Vect(2);

Bit_Counter = Bit_Counter + 1;
Sym_Counter = Sym_Counter + 1;

%------------Channel------------------%
r1 = I + normrnd(0,sqrt(noise_var)); %AWGN
r2 = Q + normrnd(0,sqrt(noise_var)); %AWGN
%------------RXR----------------------%
Rx_vect = [r1,r2];
dot1 = dot(Rx_vect,Ref_vect1);
dot2 = dot(Rx_vect,Ref_vect2);
if (dot1 == max([dot1,dot2]))
    Rx_Bits = Ref_bits1;
elseif (dot2 == max([dot1,dot2]))
    Rx_Bits = Ref_bits2;
end

%Vars for scatter plot
Scatter_x(i) = Rx_vect(1); Scatter_y(i) = Rx_vect(2); i=i+1;

%Check for signal Integrity
   
   if(Tx_Bits(1) ~= Rx_Bits(1))
       BitError_Counter = BitError_Counter +1;
       SymError_Counter = SymError_Counter +1; 
   end
   
end %---End of Transceiver chain for a given sigma value----%

    %---Error Calculations-------%
    Bit_Err_Rate = (BitError_Counter/Bit_Counter);
    Sym_Err_Rate = (SymError_Counter/Sym_Counter);

    %---Scatter Plot for High/Mid/Low Variance Cases----%
    %{
    figure('Name','Scatter Plot - High Variance','NumberTitle','off');
    subplot(2,2,1);
    scatter(Scatter_x,Scatter_y,'.r');
    xlim([-1.5 1.5]);ylim([-1.5 1.5]);
    hold on;grid on; 
    ezplot('x^2 + y^2 = 1')
    title('BPSK');
    hold off;
    %}
    %---------------End of Scatter Plot-=---------------%
end

%--------------------********--QPSK----**********-------------------------%
function [Bit_Err_Rate,Sym_Err_Rate]= RF_Transceiver_Fn_QPSK(noise_var,Cycle_Count,A_off,P_off)
%---------Initialisation----------------% 
%Vector Constellation for QPSK; Normalised -> sqrt(Es=1
%Reference vectors with Amp & Phase offsets
Ref_vect1_wOff = [A_off*cos(0+P_off),A_off*sin(0+P_off)]; 
Ref_vect2_wOff = [A_off*cos((pi/2)+P_off),A_off*sin((pi/2)+P_off)]; 
Ref_vect3_wOff = [A_off*cos(pi+P_off),A_off*sin(pi+P_off)]; 
Ref_vect4_wOff = [A_off*cos((3*pi/2)+P_off),A_off*sin((3*pi/2)+P_off)];

%Reference vectors without Amp & Phase offsets
A_off =1; P_off = 0;
Ref_vect1 = [A_off*cos(0+P_off),A_off*sin(0+P_off)]; 
Ref_vect2 = [A_off*cos((pi/2)+P_off),A_off*sin((pi/2)+P_off)]; 
Ref_vect3 = [A_off*cos(pi+P_off),A_off*sin(pi+P_off)]; 
Ref_vect4 = [A_off*cos((3*pi/2)+P_off),A_off*sin((3*pi/2)+P_off)];

%Bit Definition for 8PSK - Gray Codes
Ref_bits1 = [0,0];Ref_bits2 = [0,1]; Ref_bits3 = [1,1];Ref_bits4 = [1,0];

%--------- Counters/Variables
BitError_Counter = 0; SymError_Counter = 0;
Bit_Counter= 0; Sym_Counter = 0; i=1;

%Preallocate
Scatter_x = zeros(1,Cycle_Count);
Scatter_y = zeros(1,Cycle_Count);

%---Start of Transceiver chain for a given sigma value----%
%------------TXR-------------------%
for Tx_Count = 1:Cycle_Count
x = rand;
if(x>0 && x <= 0.25)
    Tx_Vect = Ref_vect1_wOff;     Tx_Bits = Ref_bits1;
elseif(x>0.25 && x<= 0.5)
    Tx_Vect = Ref_vect2_wOff;     Tx_Bits = Ref_bits2;
elseif (x> 0.5 && x<= 0.75)
    Tx_Vect = Ref_vect3_wOff;     Tx_Bits = Ref_bits3;
else 
    Tx_Vect = Ref_vect4_wOff;     Tx_Bits = Ref_bits4;
end

I = Tx_Vect(1); Q = Tx_Vect(2);

Bit_Counter = Bit_Counter + 2;
Sym_Counter = Sym_Counter + 1;

%------------Channel------------------%
r1 = I + normrnd(0,sqrt(noise_var)); %AWGN
r2 = Q + normrnd(0,sqrt(noise_var)); %AWGN
%------------RXR----------------------%
Rx_vect = [r1,r2];
dot1 = dot(Rx_vect,Ref_vect1);dot2 = dot(Rx_vect,Ref_vect2);
dot3 = dot(Rx_vect,Ref_vect3);dot4 = dot(Rx_vect,Ref_vect4);
if (dot1 == max([dot1,dot2,dot3,dot4]))
    Rx_Bits = Ref_bits1;
elseif (dot2 == max([dot1,dot2,dot3,dot4]))
    Rx_Bits = Ref_bits2;
elseif (dot3 == max([dot1,dot2,dot3,dot4]))
    Rx_Bits = Ref_bits3;
elseif (dot4 == max([dot1,dot2,dot3,dot4]))
    Rx_Bits = Ref_bits4;
end

%Vars for scatter plot
Scatter_x(i) = Rx_vect(1); Scatter_y(i) = Rx_vect(2); i=i+1;

%Check for signal Integrity 
   Sym_Err_Flag = 0;
   if(Tx_Bits(1) ~= Rx_Bits(1))
       BitError_Counter = BitError_Counter +1; Sym_Err_Flag = 1; end
   if(Tx_Bits(2) ~= Rx_Bits(2))
       BitError_Counter = BitError_Counter +1; Sym_Err_Flag = 1; end
   if(Sym_Err_Flag == 1)
       SymError_Counter = SymError_Counter + 1; end
   
end %---End of Transceiver chain for a given variance value----%

    %---Error Calculations-------%
    Bit_Err_Rate = (BitError_Counter/Bit_Counter);
    Sym_Err_Rate = (SymError_Counter/Sym_Counter);

    %---Scatter Plot for High/Mid/Low Variance Cases----%
    %{
    %figure('Name','Scatter Plot','NumberTitle','off');
    subplot(2,2,2);
    scatter(Scatter_x,Scatter_y,'.r');
    xlim([-1.5 1.5]);ylim([-1.5 1.5]);
    hold on;grid on; 
    ezplot('x^2 + y^2 = 1')
    title('QPSK');
    hold off;
    %}   
    %---------------End of Scatter Plot-=---------------%
end

%--------------------********--8PSK----**********-------------------------%
function [Bit_Err_Rate,Sym_Err_Rate]= RF_Transceiver_Fn_8PSK(noise_var,Cycle_Count,A_off,P_off)
%---------Initialisation----------------% 
%Vector Constellation for 8PSK; %Normalised -> sqrt(Es=1)
%Reference vectors with Amp & Phase offsets
Ref_vect1_wOff = [A_off*cos(0+P_off),A_off*sin(0+P_off)];
Ref_vect2_wOff = [A_off*cos((1*pi/4)+P_off),A_off*sin((1*pi/4)+P_off)];
Ref_vect3_wOff = [A_off*cos((2*pi/4)+P_off),A_off*sin((2*pi/4)+P_off)];
Ref_vect4_wOff = [A_off*cos((3*pi/4)+P_off),A_off*sin((3*pi/4)+P_off)];
Ref_vect5_wOff = [A_off*cos((4*pi/4)+P_off),A_off*sin((4*pi/4)+P_off)];
Ref_vect6_wOff = [A_off*cos((5*pi/4)+P_off),A_off*sin((5*pi/4)+P_off)]; 
Ref_vect7_wOff = [A_off*cos((6*pi/4)+P_off),A_off*sin((6*pi/4)+P_off)];
Ref_vect8_wOff = [A_off*cos((7*pi/4)+P_off),A_off*sin((7*pi/4)+P_off)]; 
%Reference vectors without Amp & Phase offsets
A_off = 1; P_off = 0;
Ref_vect1 = [A_off*cos(0+P_off),A_off*sin(0+P_off)];
Ref_vect2 = [A_off*cos((1*pi/4)+P_off),A_off*sin((1*pi/4)+P_off)];
Ref_vect3 = [A_off*cos((2*pi/4)+P_off),A_off*sin((2*pi/4)+P_off)];
Ref_vect4 = [A_off*cos((3*pi/4)+P_off),A_off*sin((3*pi/4)+P_off)];
Ref_vect5 = [A_off*cos((4*pi/4)+P_off),A_off*sin((4*pi/4)+P_off)];
Ref_vect6 = [A_off*cos((5*pi/4)+P_off),A_off*sin((5*pi/4)+P_off)]; 
Ref_vect7 = [A_off*cos((6*pi/4)+P_off),A_off*sin((6*pi/4)+P_off)];
Ref_vect8 = [A_off*cos((7*pi/4)+P_off),A_off*sin((7*pi/4)+P_off)]; 

%Bit Definition for 8PSK - Gray Codes
Ref_bits1 = [0,0,0]; Ref_bits2 = [0,0,1]; Ref_bits3 = [0,1,1];Ref_bits4 = [0,1,0]; 
Ref_bits5 = [1,1,0]; Ref_bits6 = [1,1,1];Ref_bits7 = [1,0,1]; Ref_bits8 = [1,0,0];

%Counters/Variables
BitError_Counter = 0; SymError_Counter = 0;
Bit_Counter= 0; Sym_Counter = 0; i=1; 

%Preallocate
Scatter_x = zeros(1,Cycle_Count);
Scatter_y = zeros(1,Cycle_Count);

%---Start of Transceiver chain for a given Variance value----%
%------------TXR-------------------%
for Tx_Count = 1:Cycle_Count
x = rand;
if(x>0 && x <= (1/8))
    Tx_Vect = Ref_vect1_wOff;     Tx_Bits = Ref_bits1;
elseif(x>(1/8) && x<= (2/8))
    Tx_Vect = Ref_vect2_wOff;     Tx_Bits = Ref_bits2;
elseif (x> (2/8) && x<= (3/8))
    Tx_Vect = Ref_vect3_wOff;     Tx_Bits = Ref_bits3;
elseif (x> (3/8) && x<= (4/8))
    Tx_Vect = Ref_vect4_wOff;     Tx_Bits = Ref_bits4;
elseif (x> (4/8) && x<= (5/8))
    Tx_Vect = Ref_vect5_wOff;     Tx_Bits = Ref_bits5;
elseif (x> (5/8) && x<= (6/8))
    Tx_Vect = Ref_vect6_wOff;     Tx_Bits = Ref_bits6;
elseif (x> (6/8) && x<= (7/8))
    Tx_Vect = Ref_vect7_wOff;     Tx_Bits = Ref_bits7;
elseif (x> (7/8) && x<= (1))
    Tx_Vect = Ref_vect8_wOff;     Tx_Bits = Ref_bits8;
end

I = Tx_Vect(1); Q = Tx_Vect(2);  

Bit_Counter = Bit_Counter + 3;
Sym_Counter = Sym_Counter + 1;

%------------Channel------------------%
r1 = I + normrnd(0,sqrt(noise_var)); %AWGN
r2 = Q + normrnd(0,sqrt(noise_var)); %AWGN
%------------RXR----------------------%
Rx_vect = [r1,r2];
dot1 = dot(Rx_vect,Ref_vect1);dot2 = dot(Rx_vect,Ref_vect2);dot3 = dot(Rx_vect,Ref_vect3);
dot4 = dot(Rx_vect,Ref_vect4);dot5 = dot(Rx_vect,Ref_vect5);dot6 = dot(Rx_vect,Ref_vect6);
dot7 = dot(Rx_vect,Ref_vect7);dot8 = dot(Rx_vect,Ref_vect8);
if (dot1 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8])) 
    Rx_Bits = Ref_bits1;
elseif (dot2 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8])) 
    Rx_Bits = Ref_bits2;
elseif (dot3 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8])) 
    Rx_Bits = Ref_bits3;
elseif (dot4 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8])) 
    Rx_Bits = Ref_bits4;
elseif (dot5 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8])) 
    Rx_Bits = Ref_bits5;
elseif (dot6 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8]))
    Rx_Bits = Ref_bits6;
elseif (dot7 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8]))
    Rx_Bits = Ref_bits7;
elseif (dot8 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8]))
    Rx_Bits = Ref_bits8;
end
%Vars for scatter plot
Scatter_x(i) = Rx_vect(1); Scatter_y(i) = Rx_vect(2); i=i+1;

%Check for signal Integrity
   Sym_Err_Flag = 0;      
   if(Rx_Bits(1) ~= Tx_Bits(1))
       BitError_Counter = BitError_Counter +1;Sym_Err_Flag = 1;end
   if(Rx_Bits(2) ~= Tx_Bits(2))
       BitError_Counter = BitError_Counter +1;Sym_Err_Flag = 1;end
   if(Rx_Bits(3) ~= Tx_Bits(3))
       BitError_Counter = BitError_Counter +1;Sym_Err_Flag = 1;end
   if(Sym_Err_Flag == 1)
       SymError_Counter = SymError_Counter + 1; end
    
end %---End of Transceiver chain for a given sigma value----%

    %---Error Calculations-------%
    Bit_Err_Rate = (BitError_Counter/Bit_Counter);
    Sym_Err_Rate = (SymError_Counter/Sym_Counter);

    %---Scatter Plot for High/Mid/Low Variance Cases----%
    %{
    %figure('Name','Scatter Plot','NumberTitle','off');
    subplot(2,2,3);
    scatter(Scatter_x,Scatter_y,'.r');
    xlim([-1.5 1.5]);ylim([-1.5 1.5]);
    hold on;grid on; 
    ezplot('x^2 + y^2 = 1')
    title('8PSK');
    hold off;
    %}  
    %---------------End of Scatter Plot-=---------------%
end

%--------------------********--16PSK----**********-------------------------%
function [Bit_Err_Rate,Sym_Err_Rate]= RF_Transceiver_Fn_16PSK(noise_var,Cycle_Count,A_off,P_off)
%---------Initialisation----------------% 
%Vector Constellation for 16PSK (Es=1,normalised)
%Reference vectors without Amp & Phase offsets
Ref_vect1_wOff = [A_off*cos((0*pi/8)+P_off),A_off*sin((0*pi/8)+P_off)];
Ref_vect2_wOff = [A_off*cos((1*pi/8)+P_off),A_off*sin((1*pi/8)+P_off)];
Ref_vect3_wOff = [A_off*cos((2*pi/8)+P_off),A_off*sin((2*pi/8)+P_off)];
Ref_vect4_wOff = [A_off*cos((3*pi/8)+P_off),A_off*sin((3*pi/8)+P_off)];
Ref_vect5_wOff = [A_off*cos((4*pi/8)+P_off),A_off*sin((4*pi/8)+P_off)];
Ref_vect6_wOff = [A_off*cos((5*pi/8)+P_off),A_off*sin((5*pi/8)+P_off)];
Ref_vect7_wOff = [A_off*cos((6*pi/8)+P_off),A_off*sin((6*pi/8)+P_off)]; 
Ref_vect8_wOff = [A_off*cos((7*pi/8)+P_off),A_off*sin((7*pi/8)+P_off)];
Ref_vect9_wOff = [A_off*cos((8*pi/8)+P_off),A_off*sin((8*pi/8)+P_off)];
Ref_vect10_wOff = [A_off*cos((9*pi/8)+P_off),A_off*sin((9*pi/8)+P_off)];
Ref_vect11_wOff = [A_off*cos((10*pi/8)+P_off),A_off*sin((10*pi/8)+P_off)];
Ref_vect12_wOff = [A_off*cos((11*pi/8)+P_off),A_off*sin((11*pi/8)+P_off)];
Ref_vect13_wOff = [A_off*cos((12*pi/8)+P_off),A_off*sin((12*pi/8)+P_off)];
Ref_vect14_wOff = [A_off*cos((13*pi/8)+P_off),A_off*sin((13*pi/8)+P_off)];
Ref_vect15_wOff = [A_off*cos((14*pi/8)+P_off),A_off*sin((14*pi/8)+P_off)];
Ref_vect16_wOff = [A_off*cos((15*pi/8)+P_off),A_off*sin((15*pi/8)+P_off)];
%Reference vectors without Amp & Phase offsets
A_off = 1; P_off = 0;
Ref_vect1 = [A_off*cos((0*pi/8)+P_off),A_off*sin((0*pi/8)+P_off)];
Ref_vect2 = [A_off*cos((1*pi/8)+P_off),A_off*sin((1*pi/8)+P_off)];
Ref_vect3 = [A_off*cos((2*pi/8)+P_off),A_off*sin((2*pi/8)+P_off)];
Ref_vect4 = [A_off*cos((3*pi/8)+P_off),A_off*sin((3*pi/8)+P_off)];
Ref_vect5 = [A_off*cos((4*pi/8)+P_off),A_off*sin((4*pi/8)+P_off)];
Ref_vect6 = [A_off*cos((5*pi/8)+P_off),A_off*sin((5*pi/8)+P_off)];
Ref_vect7 = [A_off*cos((6*pi/8)+P_off),A_off*sin((6*pi/8)+P_off)]; 
Ref_vect8 = [A_off*cos((7*pi/8)+P_off),A_off*sin((7*pi/8)+P_off)];
Ref_vect9 = [A_off*cos((8*pi/8)+P_off),A_off*sin((8*pi/8)+P_off)];
Ref_vect10 = [A_off*cos((9*pi/8)+P_off),A_off*sin((9*pi/8)+P_off)];
Ref_vect11 = [A_off*cos((10*pi/8)+P_off),A_off*sin((10*pi/8)+P_off)];
Ref_vect12 = [A_off*cos((11*pi/8)+P_off),A_off*sin((11*pi/8)+P_off)];
Ref_vect13 = [A_off*cos((12*pi/8)+P_off),A_off*sin((12*pi/8)+P_off)];
Ref_vect14 = [A_off*cos((13*pi/8)+P_off),A_off*sin((13*pi/8)+P_off)];
Ref_vect15 = [A_off*cos((14*pi/8)+P_off),A_off*sin((14*pi/8)+P_off)];
Ref_vect16 = [A_off*cos((15*pi/8)+P_off),A_off*sin((15*pi/8)+P_off)];

%Bit Definition for 16PSK - Gray Codes
Ref_bits1 = [0,0,0,0]; Ref_bits2 = [0,0,0,1]; Ref_bits3 = [0,0,1,1];
Ref_bits4 = [0,0,1,0]; Ref_bits5 = [0,1,1,0]; Ref_bits6 = [0,1,1,1];
Ref_bits7 = [0,1,0,1]; Ref_bits8 = [0,1,0,0]; Ref_bits9 = [1,1,0,0];
Ref_bits10 = [1,1,0,1]; Ref_bits11 = [1,1,1,1]; Ref_bits12 = [1,1,1,0];
Ref_bits13 = [1,0,1,0]; Ref_bits14 = [1,0,1,1]; Ref_bits15 = [1,0,0,1];
Ref_bits16 = [1,0,0,0];

%Counters/Variables
BitError_Counter = 0; SymError_Counter = 0;
Bit_Counter= 0; Sym_Counter = 0; i=1; 

%Preallocate
Scatter_x = zeros(1,Cycle_Count);
Scatter_y = zeros(1,Cycle_Count);

%---Start of Transceiver chain for a given variance value----%
%------------TXR-------------------%
for Tx_Count = 1:Cycle_Count
x = rand;
if(x>0 && x <= (1/16))
    Tx_Vect = Ref_vect1_wOff;     Tx_Bits = Ref_bits1; 
elseif(x>(1/16) && x<= (2/16))
    Tx_Vect = Ref_vect2_wOff;    Tx_Bits = Ref_bits2;
elseif (x> (2/16) && x<= (3/16))
    Tx_Vect = Ref_vect3_wOff;    Tx_Bits = Ref_bits3;
elseif (x> (3/16) && x<= (4/16))
    Tx_Vect = Ref_vect4_wOff;    Tx_Bits = Ref_bits4;
elseif (x> (4/16) && x<= (5/16))
    Tx_Vect = Ref_vect5_wOff;     Tx_Bits = Ref_bits5;
elseif (x> (5/16) && x<= (6/16))
    Tx_Vect = Ref_vect6_wOff;    Tx_Bits = Ref_bits6;
elseif (x> (6/16) && x<= (7/16))
    Tx_Vect = Ref_vect7_wOff;    Tx_Bits = Ref_bits7;
elseif (x> (7/16) && x<= (8/16))
    Tx_Vect = Ref_vect8_wOff;    Tx_Bits = Ref_bits8;
elseif (x> (8/16) && x<= (9/16))
    Tx_Vect = Ref_vect9_wOff;    Tx_Bits = Ref_bits9;
elseif (x> (9/16) && x<= (10/16))
    Tx_Vect = Ref_vect10_wOff;    Tx_Bits = Ref_bits10;
elseif (x> (10/16) && x<= (11/16))
    Tx_Vect = Ref_vect11_wOff;    Tx_Bits = Ref_bits11;
elseif (x> (11/16) && x<= (12/16))
    Tx_Vect = Ref_vect12_wOff;    Tx_Bits = Ref_bits12;
elseif (x> (12/16) && x<= (13/16))
    Tx_Vect = Ref_vect13_wOff;    Tx_Bits = Ref_bits13;
elseif (x> (13/16) && x<= (14/16))
    Tx_Vect = Ref_vect14_wOff;    Tx_Bits = Ref_bits14;
elseif (x> (14/16) && x<= (15/16))
    Tx_Vect = Ref_vect15_wOff;    Tx_Bits = Ref_bits15;
elseif (x> (15/16) && x<= (16/16))
    Tx_Vect = Ref_vect16_wOff;    Tx_Bits = Ref_bits16;
end

I = Tx_Vect(1); Q = Tx_Vect(2);

Bit_Counter = Bit_Counter + 4;
Sym_Counter = Sym_Counter + 1;

%------------Channel------------------%
r1 = I + normrnd(0,sqrt(noise_var)); %AWGN
r2 = Q + normrnd(0,sqrt(noise_var)); %AWGN
%------------RXR----------------------%
Rx_vect = [r1,r2];
dot1 = dot(Rx_vect,Ref_vect1);dot2 = dot(Rx_vect,Ref_vect2);dot3 = dot(Rx_vect,Ref_vect3);
dot4 = dot(Rx_vect,Ref_vect4);dot5 = dot(Rx_vect,Ref_vect5);dot6 = dot(Rx_vect,Ref_vect6);
dot7 = dot(Rx_vect,Ref_vect7);dot8 = dot(Rx_vect,Ref_vect8);dot9 = dot(Rx_vect,Ref_vect9);
dot10 = dot(Rx_vect,Ref_vect10);dot11 = dot(Rx_vect,Ref_vect11);dot12 = dot(Rx_vect,Ref_vect12);
dot13 = dot(Rx_vect,Ref_vect13);dot14 = dot(Rx_vect,Ref_vect14);
dot15 = dot(Rx_vect,Ref_vect15);dot16 = dot(Rx_vect,Ref_vect16);
if (dot1 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits1;
elseif (dot2 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits2;
elseif (dot3 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits3;
elseif (dot4 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits4;
elseif (dot5 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits5;
elseif (dot6 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits6;
elseif (dot7 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits7;
elseif (dot8 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits8;
elseif (dot9 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits9;
elseif (dot10 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits10;
elseif (dot11 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits11;
elseif (dot12 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits12;
elseif (dot13 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits13;
elseif (dot14 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits14;
elseif (dot15 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits15;
elseif (dot16 == max([dot1,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11,dot12,dot13,dot14,dot15,dot16]))
    Rx_Bits = Ref_bits16; 
end
%Vars for scatter plot
Scatter_x(i) = Rx_vect(1); Scatter_y(i) = Rx_vect(2);i=i+1;

%Check for signal Integrity
       Sym_Err_Flag = 0;     
       if(Rx_Bits(1) ~= Tx_Bits(1))
           BitError_Counter = BitError_Counter +1;Sym_Err_Flag = 1; end
       if(Rx_Bits(2) ~= Tx_Bits(2))
           BitError_Counter = BitError_Counter +1;Sym_Err_Flag = 1; end
       if(Rx_Bits(3) ~= Tx_Bits(3))
           BitError_Counter = BitError_Counter +1;Sym_Err_Flag = 1; end
       if(Rx_Bits(4) ~= Tx_Bits(4))
           BitError_Counter = BitError_Counter +1;Sym_Err_Flag = 1; end
       if(Sym_Err_Flag == 1)
           SymError_Counter = SymError_Counter + 1; end
    
end %---End of Transceiver chain for a given sigma value----%

    %---Error Calculations-------%
    Bit_Err_Rate = (BitError_Counter/Bit_Counter);
    Sym_Err_Rate = (SymError_Counter/Sym_Counter);

    %---Scatter Plot for High/Mid/Low Variance Cases----%
    %{
    %figure('Name','Scatter Plot','NumberTitle','off');
    subplot(2,2,4);
    scatter(Scatter_x,Scatter_y,'.r');
    xlim([-1.5 1.5]);ylim([-1.5 1.5]);
    hold on;grid on; 
    ezplot('x^2 + y^2 = 1')
    title('16PSK');
    hold off;
    %}
    %---------------End of Scatter Plot-=---------------%
end


