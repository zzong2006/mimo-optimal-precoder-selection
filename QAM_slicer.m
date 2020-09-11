% QAM_slicer.m
function [sliced_symbolseq] = QAM_slicer(symbolseq,Nb)
if Nb==2, sliced_symbolseq=QAM4_slicer(symbolseq);
elseif Nb==4, sliced_symbolseq=QAM16_slicer(symbolseq);
elseif Nb==6, sliced_symbolseq=QAM64_slicer(symbolseq);
else disp('When you need M-QAM modulation, you need to provide log2(M) as the 2nd argument');
end

function [sliced_symbolseq] = QAM4_slicer(symbolseq)
j=sqrt(-1);
symbolseq_length=length(symbolseq);
sliced_symbolseq=zeros(1,symbolseq_length);
for ii=1:symbolseq_length
    real_part=real(symbolseq(ii));
    imag_part=imag(symbolseq(ii));
    if (real_part>0)&&(imag_part>0), sliced_symbolseq(ii)= (1+j)/sqrt(2.);
    elseif (real_part<0)&&(imag_part>0), sliced_symbolseq(ii)= (-1+j)/sqrt(2.);
    elseif (real_part<0)&&(imag_part<0), sliced_symbolseq(ii)=(-1-j)/sqrt(2.);
    elseif (real_part>0)&&(imag_part<0), sliced_symbolseq(ii)=(1-j)/sqrt(2.); end
end

function [sliced_symbolseq] =QAM16_slicer(symbolseq)
j=sqrt(-1);
symbolseq_length = length(symbolseq);
for ii=1:symbolseq_length
   real_part=real(symbolseq(ii));
   imag_part=imag(symbolseq(ii));   
   if real_part<-2./sqrt(10.), real_hat=-3./sqrt(10.);
   elseif (-2./sqrt(10.)<=real_part)&&(real_part<0), real_hat=-1./sqrt(10.);
   elseif (0<=real_part)&&(real_part<2./sqrt(10.)), real_hat=1./sqrt(10.);
   elseif 2./sqrt(10.)<=real_part, real_hat=3./sqrt(10.); end  
   if imag_part<-2./sqrt(10.), imag_hat=-3./sqrt(10.);
   elseif (-2./sqrt(10.)<=imag_part)&&(imag_part<0), imag_hat=-1./sqrt(10.);
   elseif (0<=imag_part)&&(imag_part<2./sqrt(10.)), imag_hat=1./sqrt(10.);
   elseif 2./sqrt(10.)<=imag_part, imag_hat=3./sqrt(10.); end   
   sliced_symbolseq(ii)=real_hat+j*imag_hat;  
end

function [sliced_symbolseq]=QAM64_slicer(symbolseq)
j=sqrt(-1);
symbolseq_length = length(symbolseq);
for ii = 1: symbolseq_length
    real_part = real(symbolseq(ii));
    imag_part = imag(symbolseq(ii));
    if(real_part < -6/sqrt(42.)), real_hat = -7/sqrt(42.);
    elseif(real_part < -4/sqrt(42.)), real_hat = -5/sqrt(42.);
    elseif(real_part < -2/sqrt(42.)), real_hat = -3/sqrt(42.);
    elseif(real_part < 0/sqrt(42.)), real_hat = -1/sqrt(42.);
    elseif(real_part < 2/sqrt(42.)), real_hat = 1/sqrt(42.);
    elseif(real_part < 4/sqrt(42.)), real_hat = 3/sqrt(42.);
    elseif(real_part < 6/sqrt(42.)), real_hat = 5/sqrt(42.);
    elseif(real_part >= 6/sqrt(42.)), real_hat = 7/sqrt(42.); end  
    if(imag_part < -6/sqrt(42.)), imag_hat = -7j/sqrt(42.);
    elseif(imag_part < -4/sqrt(42.)), imag_hat = -5j/sqrt(42.);
    elseif(imag_part < -2/sqrt(42.)), imag_hat = -3j/sqrt(42.);
    elseif(imag_part < 0/sqrt(42.)), imag_hat = -j/sqrt(42.);
    elseif(imag_part < 2/sqrt(42.)), imag_hat = j/sqrt(42.);
    elseif(imag_part < 4/sqrt(42.)), imag_hat = 3j/sqrt(42.);
    elseif(imag_part < 6/sqrt(42.)), imag_hat = 5j/sqrt(42.);
    elseif(imag_part >= 6/sqrt(42.)), imag_hat = 7j/sqrt(42.); end    
    sliced_symbolseq(ii) = real_hat + imag_hat;
end