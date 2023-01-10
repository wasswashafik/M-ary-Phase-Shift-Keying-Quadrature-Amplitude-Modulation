clc; close all; clear all;
%################---Set Values for the following----######################%
Cycle_Count = 2500;
Var_Start = 0.01; 
Var_Res = 0.01; 
Var_End = 0.01;   

%**********************Modulation Schemes in M-QAM*************************%
%-------@@@@@@@@@@@@--------16QAM-------@@@@@@@@@@@@@------%
K = 4; M= 16; % 10log(M=16)
%Pre-allocate memory
Bit_Err_Rate = zeros(1,500);
Sym_Err_Rate = zeros(1,500);
Bit_Energy_16QAM = zeros(1,500);
Sym_Err_Rate_16QAM_Theory = zeros(1,500);

j=0; %temp variable
for noise_var = Var_Start : Var_Res : Var_End
    j=j+1
   [Bit_Err_Rate(j),Sym_Err_Rate(j),Eav]= RF_Transceiver_Fn_16QAM(noise_var,Cycle_Count);
    Bit_Energy_16QAM(j) = (Eav/(2*K*noise_var));%(Es/(2*K*noise_var));
    %Theoretical Value for comparison
    Sym_Err_Rate_16QAM_Theory(j) = 1 - (1 - 2*qfunc(sqrt(3*Eav/(30*noise_var))))^2;
end
Bit_Err_Rate_16QAM = Bit_Err_Rate;
Sym_Err_Rate_16QAM = Sym_Err_Rate;
Bit_Energy_16QAM_dB = 10.*log10(Bit_Energy_16QAM);

%#######################---------Plots------##############################%
figure('Name','BER & SER for QAM','NumberTitle','off');
hold on; grid on;
set(gca, 'YScale', 'log')
semilogy(Bit_Energy_16QAM_dB,Bit_Err_Rate_16QAM,'-r','LineWidth',2);
semilogy(Bit_Energy_16QAM_dB,Sym_Err_Rate_16QAM,'-b','LineWidth',2);
semilogy(Bit_Energy_16QAM_dB,Sym_Err_Rate_16QAM_Theory,'b--o','LineWidth',1);
xlabel('Bit Energy(Eb/No) in dB)'); ylabel('Bit Error Rate & Symbol Error Rate');
legend('BER','SER','Theoretical SER');
xlim([0 25]); ylim([10^-5 1]);
hold off;

%#######################----------Functions--------######################%
%--------------------********--16QAM----**********-------------------------%
function [Bit_Err_Rate,Sym_Err_Rate,Eav]= RF_Transceiver_Fn_16QAM(noise_var,Cycle_Count)
%---------Initialisation----------------% 
%Vector Constellation for 16QAM
%Reference vectors
QAM_Ref_Vect1 = [1,1]; QAM_Ref_Vect2 = [3,1]; QAM_Ref_Vect3 = [1,3];
QAM_Ref_Vect4 = [3,3]; QAM_Ref_Vect5 = [-1,1]; QAM_Ref_Vect6 = [-3,1];
QAM_Ref_Vect7 = [-1,3]; QAM_Ref_Vect8 = [-3,3]; QAM_Ref_Vect9 = [-1,-1];
QAM_Ref_Vect10 = [-3,-1]; QAM_Ref_Vect11 = [-3,-3]; QAM_Ref_Vect12 = [-1,-3];
QAM_Ref_Vect13 = [1,-1]; QAM_Ref_Vect14 = [3,-1]; QAM_Ref_Vect15 = [1,-3];
QAM_Ref_Vect16 = [3,-3];

%Bit Definition for 16PSK - Gray Codes
QAM_Ref_bits1 = [1,1,1,1]; QAM_Ref_bits2 = [1,0,1,1];
QAM_Ref_bits3 = [1,1,1,0]; QAM_Ref_bits4 = [1,0,1,0];
QAM_Ref_bits5 = [0,1,1,1]; QAM_Ref_bits6 = [0,0,1,1];
QAM_Ref_bits7 = [0,1,1,0]; QAM_Ref_bits8 = [0,0,1,0];
QAM_Ref_bits9 = [0,1,0,1]; QAM_Ref_bits10 = [0,0,0,1];
QAM_Ref_bits11 = [0,0,0,0]; QAM_Ref_bits12 = [0,1,0,0];
QAM_Ref_bits13 = [1,1,0,1]; QAM_Ref_bits14 = [1,0,0,1];
QAM_Ref_bits15 = [1,1,0,0]; QAM_Ref_bits16 = [1,0,0,0];

%Preallocate
Scatter_x = zeros(1,Cycle_Count);
Scatter_y = zeros(1,Cycle_Count);
Sym_Energy = zeros(1,Cycle_Count);

%Counters/Variables
BitError_Counter = 0; SymError_Counter = 0;
Bit_Counter= 0; Sym_Counter = 0;  i=1; 

%---Start of Transceiver chain for a given variance value----%
%------------TXR-------------------%
for Tx_Count = 1:Cycle_Count
x = rand;
if(x>0 && x <= (1/16))
    Tx_Vect = QAM_Ref_Vect1;    Tx_Bits = QAM_Ref_bits1; 
elseif(x>(1/16) && x<= (2/16))
    Tx_Vect = QAM_Ref_Vect2;    Tx_Bits = QAM_Ref_bits2;
elseif (x> (2/16) && x<= (3/16))
    Tx_Vect = QAM_Ref_Vect3;    Tx_Bits = QAM_Ref_bits3;
elseif (x> (3/16) && x<= (4/16))
    Tx_Vect = QAM_Ref_Vect4;    Tx_Bits = QAM_Ref_bits4;
elseif (x> (4/16) && x<= (5/16))
    Tx_Vect = QAM_Ref_Vect5;    Tx_Bits = QAM_Ref_bits5;
elseif (x> (5/16) && x<= (6/16))
    Tx_Vect = QAM_Ref_Vect6;    Tx_Bits = QAM_Ref_bits6;
elseif (x> (6/16) && x<= (7/16))
    Tx_Vect = QAM_Ref_Vect7;    Tx_Bits = QAM_Ref_bits7;
elseif (x> (7/16) && x<= (8/16))
    Tx_Vect = QAM_Ref_Vect8;    Tx_Bits = QAM_Ref_bits8;
elseif (x> (8/16) && x<= (9/16))
    Tx_Vect = QAM_Ref_Vect9;    Tx_Bits = QAM_Ref_bits9;
elseif (x> (9/16) && x<= (10/16))
    Tx_Vect = QAM_Ref_Vect10;    Tx_Bits = QAM_Ref_bits10;
elseif (x> (10/16) && x<= (11/16))
    Tx_Vect = QAM_Ref_Vect11;    Tx_Bits = QAM_Ref_bits11;
elseif (x> (11/16) && x<= (12/16))
    Tx_Vect = QAM_Ref_Vect12;    Tx_Bits = QAM_Ref_bits12;
elseif (x> (12/16) && x<= (13/16))
    Tx_Vect = QAM_Ref_Vect13;    Tx_Bits = QAM_Ref_bits13;
elseif (x> (13/16) && x<= (14/16))
    Tx_Vect = QAM_Ref_Vect14;    Tx_Bits = QAM_Ref_bits14;
elseif (x> (14/16) && x<= (15/16))
    Tx_Vect = QAM_Ref_Vect15;    Tx_Bits = QAM_Ref_bits15;
elseif (x> (15/16) && x<= (16/16))
    Tx_Vect = QAM_Ref_Vect16;    Tx_Bits = QAM_Ref_bits16;
end

I = Tx_Vect(1); Q = Tx_Vect(2); 
Sym_Energy(i) = I^2 + Q^2;


Bit_Counter = Bit_Counter + 4;
Sym_Counter = Sym_Counter + 1;

%------------Channel------------------%
r1 = I + normrnd(0,sqrt(noise_var)); %AWGN
r2 = Q + normrnd(0,sqrt(noise_var)); %AWGN
%------------RXR----------------------%
Rx_vect = [r1,r2];

switch true 
    case r1 >=0 && r2>=0 %Quadrant I
        evm1 = (norm( Rx_vect - QAM_Ref_Vect1 ) / norm( QAM_Ref_Vect1 ));
        evm2 = (norm( Rx_vect - QAM_Ref_Vect2 ) / norm( QAM_Ref_Vect2 ));
        evm3 = (norm( Rx_vect - QAM_Ref_Vect3 ) / norm( QAM_Ref_Vect3 ));
        evm4 = (norm( Rx_vect - QAM_Ref_Vect4 ) / norm( QAM_Ref_Vect4 ));
        if (evm1 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits1;
        elseif (evm2 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits2;
        elseif (evm3 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits3;
        elseif (evm4 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits4;
        end
    case r1 <0 && r2>0 %Quadrant II
        evm1 = (norm( Rx_vect - QAM_Ref_Vect5 ) / norm( QAM_Ref_Vect5 ));
        evm2 = (norm( Rx_vect - QAM_Ref_Vect6 ) / norm( QAM_Ref_Vect6 ));
        evm3 = (norm( Rx_vect - QAM_Ref_Vect7 ) / norm( QAM_Ref_Vect7 ));
        evm4 = (norm( Rx_vect - QAM_Ref_Vect8 ) / norm( QAM_Ref_Vect8 ));
        if (evm1 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits5;
        elseif (evm2 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits6;
        elseif (evm3 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits7;
        elseif (evm4 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits8;
        end
    case r1 <0 && r2<0 %Quadrant III
        evm1 = (norm( Rx_vect - QAM_Ref_Vect9 ) / norm( QAM_Ref_Vect9 ));
        evm2 = (norm( Rx_vect - QAM_Ref_Vect10 ) / norm( QAM_Ref_Vect10 ));
        evm3 = (norm( Rx_vect - QAM_Ref_Vect11 ) / norm( QAM_Ref_Vect11 ));
        evm4 = (norm( Rx_vect - QAM_Ref_Vect12 ) / norm( QAM_Ref_Vect12 ));
        if (evm1 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits9;
        elseif (evm2 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits10;
        elseif (evm3 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits11;
        elseif (evm4 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits12;
        end
    case r1 >0 && r2<0 %Quadrant IV
        evm1 = (norm( Rx_vect - QAM_Ref_Vect13 ) / norm( QAM_Ref_Vect13 ));
        evm2 = (norm( Rx_vect - QAM_Ref_Vect14 ) / norm( QAM_Ref_Vect14 ));
        evm3 = (norm( Rx_vect - QAM_Ref_Vect15 ) / norm( QAM_Ref_Vect15 ));
        evm4 = (norm( Rx_vect - QAM_Ref_Vect16 ) / norm( QAM_Ref_Vect16 ));
        if (evm1 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits13;
        elseif (evm2 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits14;
        elseif (evm3 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits15;
        elseif (evm4 == min([evm1,evm2,evm3,evm4]))
            Rx_Bits = QAM_Ref_bits16;
        end
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
   if(Rx_Bits(4) ~= Tx_Bits(4))
       BitError_Counter = BitError_Counter +1;Sym_Err_Flag = 1;end
   if(Sym_Err_Flag == 1)
   SymError_Counter = SymError_Counter + 1; end
    
end %---End of Transceiver chain for a given variance value----%

    %---Error & Energy Calculations-------%
    Bit_Err_Rate = (BitError_Counter/Bit_Counter);
    Sym_Err_Rate = (SymError_Counter/Sym_Counter);
    Eav = mean(Sym_Energy)
    
     %---Scatter Plot for High/Mid/Low Variance Cases----%
    %{
    figure('Name','Scatter Plot - High Variance','NumberTitle','off');
    scatter(Scatter_x,Scatter_y,'.r');
    xlim([-4 4]);ylim([-4 4]);
    hold on;grid on; 
    rectangle('Position',[-3 -3 6 6],'EdgeColor','c')
    title('16 QAM');
    hold off;
    %}    
    %---------------End of Scatter Plot-=---------------%
    
    
end




