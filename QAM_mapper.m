% QAM_mapper.m
function [symbol] = QAM_mapper(bitseq,Nb)
if Nb==2, symbol=QAM4_mapper(bitseq);
elseif Nb==4, symbol=QAM16_mapper(bitseq);
elseif Nb==6, symbol=QAM64_mapper(bitseq);
else disp('When you need M-QAM modulation, you need to provide log2(M) as the 2nd argument');
end

function [symbol] = QAM4_mapper(bitseq)
j =sqrt(-1);
bitseq_length=length(bitseq);
symbol_length=bitseq_length/2;
symbol=zeros(1,symbol_length);
QAM4_const = [1+j -1+j -1-j 1-j ]/sqrt(2.); 
for ii=0: symbol_length-1, symbol(ii+1)=QAM4_const(bitseq(2*ii+1)*2 +bitseq(2*ii+2)+1); end

function [symbol] =QAM16_mapper(bitseq)
j =sqrt(-1);
bitseq_length=length(bitseq);
symbol_length=bitseq_length/4;
symbol=zeros(1,symbol_length);
QAM16_const = [-3-3*j, -3-j, -3+3*j, -3+j, -1-3*j, -1-j, -1+3*j, -1+j, 3-3*j, 3-j, 3+3*j, 3+j, 1-3*j, 1-j, 1+3*j, 1+j]/sqrt(10.);
for ii=0:symbol_length-1 
    symbol(ii+1) = QAM16_const(bitseq(4*ii+1)*8 +bitseq(4*ii+2)*4 +bitseq(4*ii+3)*2 +bitseq(4*ii+4)+1); 
end

function [symbol] = QAM64_mapper(bitseq)
j =sqrt(-1);
bitseq_length=length(bitseq);
symbol_length=bitseq_length/6;
symbol=zeros(1,symbol_length);
QAM64_table = [-7+7*j,-7+5*j, -7+1*j, -7+3*j, -7-7*j,-7-5*j, -7-1*j, -7-3*j,...
              -5+7*j,-5+5*j, -5+1*j, -5+3*j, -5-7*j,-5-5*j, -5-1*j, -5-3*j,...
              -1+7*j,-1+5*j, -1+1*j, -1+3*j, -1-7*j,-1-5*j, -1-1*j, -1-3*j,...
              -3+7*j,-3+5*j, -3+1*j, -3+3*j, -3-7*j,-3-5*j, -3-1*j, -3-3*j,...
               7+7*j, 7+5*j,  7+1*j,  7+3*j,  7-7*j, 7-5*j,  7-1*j,  7-3*j,...
               5+7*j, 5+5*j,  5+1*j,  5+3*j,  5-7*j, 5-5*j,  5-1*j,  5-3*j,...
               1+7*j, 1+5*j,  1+1*j,  1+3*j,  1-7*j, 1-5*j,  1-1*j,  1-3*j,...
               3+7*j, 3+5*j,  3+1*j,  3+3*j,  3-7*j, 3-5*j,  3-1*j,  3-3*j,...
                ]/sqrt(42.);    
for ii=0:symbol_length-1, symbol(ii+1) =QAM64_table(bitseq(6*ii+1)*32 +bitseq(6*ii+2)*16 +bitseq(6*ii+3)*8 +bitseq(6*ii+4)*4+bitseq(6*ii+5)*2+bitseq(6*ii+6)+1); end