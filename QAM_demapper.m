% QAM_mapper.m
function [bitseq] = QAM_demapper(symbol,Nb)
if Nb==2, bitseq=QAM4_demapper(symbol);
elseif Nb==4, bitseq=QAM16_demapper(symbol);
elseif Nb==6, bitseq=QAM64_demapper(symbol);
else disp('When you need M-QAM modulation, you need to provide log2(M) as the 2nd argument');
end

function [bitseq] = QAM4_demapper(symbolseq)
j=sqrt(-1);
QAM4_const = [1+j -1+j -1-j 1-j ]/sqrt(2.); 
symbolseq_length=length(symbolseq);
temp = [];
for ii=1:symbolseq_length, temp=[temp dec2bin(find(QAM4_const==symbolseq(ii))-1,2)]; end
for ii=1:length(temp), bitseq(ii)=bin2dec(temp(ii)); end

function [bitseq] =QAM16_demapper(symbolseq)
j =sqrt(-1);
symbolseq_length = length(symbolseq);
QAM16_const = [-3-3*j, -3-j, -3+3*j, -3+j, -1-3*j, -1-j, -1+3*j, -1+j, 3-3*j, 3-j, 3+3*j, 3+j, 1-3*j, 1-j, 1+3*j, 1+j]/sqrt(10.);
temp=[];
for ii=0:symbolseq_length-1,  temp=[temp dec2bin(find(QAM16_const==symbolseq(ii+1))-1,4)];end
for ii=1:length(temp), bitseq(ii)=bin2dec(temp(ii)); end

function [bitseq] =QAM64_demapper(symbolseq)
j =sqrt(-1);
symbolseq_length = length(symbolseq);
QAM64_const = [-7+7*j,-7+5*j, -7+1*j, -7+3*j, -7-7*j,-7-5*j, -7-1*j, -7-3*j,...
              -5+7*j,-5+5*j, -5+1*j, -5+3*j, -5-7*j,-5-5*j, -5-1*j, -5-3*j,...
              -1+7*j,-1+5*j, -1+1*j, -1+3*j, -1-7*j,-1-5*j, -1-1*j, -1-3*j,...
              -3+7*j,-3+5*j, -3+1*j, -3+3*j, -3-7*j,-3-5*j, -3-1*j, -3-3*j,...
               7+7*j, 7+5*j,  7+1*j,  7+3*j,  7-7*j, 7-5*j,  7-1*j,  7-3*j,...
               5+7*j, 5+5*j,  5+1*j,  5+3*j,  5-7*j, 5-5*j,  5-1*j,  5-3*j,...
               1+7*j, 1+5*j,  1+1*j,  1+3*j,  1-7*j, 1-5*j,  1-1*j,  1-3*j,...
               3+7*j, 3+5*j,  3+1*j,  3+3*j,  3-7*j, 3-5*j,  3-1*j,  3-3*j,...
               ]/sqrt(42.);  
temp = [];
for ii=0:symbolseq_length-1, temp = [temp dec2bin(find(QAM64_const==symbolseq(ii+1))-1,6)]; end 
for ii = 1: length(temp),  bitseq(ii) = bin2dec(temp(ii)); end