function IPAout=IPA_measure(x,Fs,Dur)
% this is based on Duchowski, A. T., Krejtz, K., Krejtz, I., Biele, C.,
% Niedzielska, A., Kiefer, P., Raubal, M., & Giannopoulos, I. (2018). The Index of Pupillary Activity.
% In Proceedings of the 46th annual ACM Conference
% on Human Factors in Computing Systems (CHI’18).
% ACM. doi: https://doi.org/10.1145/3173574 .3173856

% x:is the pupil dimaeter after remvoing the blinks as instructed by the reference
% Dur: duration of the signal

WvOrd=round(Fs*130/1000/2); % the window of filteration should be about 130 ms (page 4 of the reference)

WvNm=['sym' num2str(WvOrd)];

% wavelet decomposition
[cA1,cD1] = dwt(x,WvNm,'mode','per');
[cA2,cD2] = dwt(cA1,WvNm,'mode','per');

% normalize by 1=2j , j = 2 for 2-level DWT
cD1=cD1/sqrt(2);
cD2=cD2/sqrt(4); cA2=cA2/sqrt(4);

cD2m = modmax(cD2);

% threshold using universal threshold luniv = sˆp(2logn)
% where sˆ is the standard deviation of the noise

Thr=std(cD2m)*sqrt(2.0*log2(length(cD2m)));

% XD = wden(cD2m,Thr,'h','one',2,WvNm); % h: hard thresholding, 'one': no rescaling, 2 level of wavelet decomposition
cD2mF=cD2m; % to make sure that the wmode is also set to per when reconstructing, I replaced the following lines
cD2mF(abs(cD2m)<Thr)=0;

% cD1F=zeros(size(cD1)); cA2F=zeros(size(cA2));
% 
% XDF1 = idwt(cA2F,cD2mF,WvNm,'mode','per');
% XDF=idwt(XDF1,cD1F,WvNm,'mode','per');

IPAout=sum(abs(cD2mF)>0)/Dur;

function  out=modmax(cD2)

m=abs(cD2);

MPrv=[m(1) m(1:end-1)];
MNxt=[m(2:end) m(end)];

Dff1=m-MPrv;
Dff2=m-MNxt;


Ind=((Dff1>=0 & Dff2>=0) & (Dff1>0 | Dff2>0));
out=zeros(size(cD2));

out(Ind)=abs(cD2(Ind));