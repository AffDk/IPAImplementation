function WrpIPA

clc; clear; close all;
SbFldr='C:\Users\afsamani\OneDrive - Aalborg Universitet\WrokResults\Students\Fall2020\S6';
DtFil1=fullfile(SbFldr,'OutputMat1.mat');
DtFil12=fullfile(SbFldr,'OutputMat2.mat'); % 'OutVrbFx','VrbNmsFx','DtPupil','GzDt','BlnkDt'

DtStrct=load(DtFil12);

DtPupil=DtStrct.DtPupil;
BlnkDt=DtStrct.BlnkDt;
OutVrbFx=DtStrct.OutVrbFx;

CnfPup=DtPupil.confidence;

TmPupil=DtPupil.pupil_timestamp;
EyId=DtPupil.eye_id;

TmPupil0=TmPupil(EyId==0);
TmPupil1=TmPupil(EyId==1);

CnfFix=str2double(OutVrbFx(:,3));
FixDur=(str2double(OutVrbFx(:,3))/1000);
FixStr=str2double(OutVrbFx(:,2)); FixStp=FixStr+FixDur;
FixStr(CnfFix<0.6)=nan;

TmBlnk=BlnkDt.start_timestamp;
DurBlnk=BlnkDt.duration;
BlnkStr=TmBlnk-0.2;  BlnkStp=TmBlnk+DurBlnk+0.2;

[PplFxStrInd0,MnVlStr]=FixGzInd(TmPupil0,FixStr);
[PplFxStpInd0,MnVlStp]=FixGzInd(TmPupil0,FixStp);

PupilDiam=DtPupil.diameter_3d;
PupilDiam(CnfPup<0.6)=nan;
PupilDiam0=PupilDiam(EyId==0);
PupilDiam1=PupilDiam(EyId==1);

FxLn=length(PplFxStrInd0);
OutTtl=nan(1,FxLn); IPACnt=0;
for Cnt=1:FxLn
    X=PupilDiam0(PplFxStrInd0(Cnt):PplFxStpInd0(Cnt));
    TmP=TmPupil0(PplFxStrInd0(Cnt):PplFxStpInd0(Cnt));
    
    Ind=(TmP(1)<=BlnkStr & TmP(end)>=BlnkStp);
    
    if ~any(Ind)
        Dur=TmP(end)-TmP(1);
        if Dur>0.6
            IPACnt=IPACnt+1;
            Fs=120;
            NwTmVct=linspace(TmP(1),TmP(end),round(Dur*Fs));
            IndNn=isnan(X); X(IndNn)=[]; TmP(IndNn)=[];
            Xnw=interp1(TmP,X,NwTmVct,'linear','extrap');
            IndNn=isnan(Xnw);
            Xnw(IndNn)=[];
            IPAout=IPA_measure(Xnw,Fs,Dur);
            OutTtl(IPACnt)=IPAout;
        end
    else
        alaki=0;
    end
end
OutTtl(IPACnt+1:end)=[];
alaki=0;


function [GzFxInd,PgMnVl]=FixGzInd(ZrTm,FxTms)

% FxTms=str2double(OutVrbFx); % fixation onset times
if ~isscalar(FxTms)
    LnFx=length(FxTms);
    LnGz=length(ZrTm);
    
    if LnGz>fix(LnFx/2)
        PgLn=fix(LnFx/2);
        PgNm=fix(LnGz/PgLn)+1;
    else
        PgLn=LnGz;
        PgNm=1;
    end
    PgMn=nan(PgNm,LnFx);
    PgMnId=nan(PgNm,LnFx);
    for PgCnt=1:PgNm
        StrIndx=min((PgCnt-1)*PgLn+1,LnGz+1);
        StpInd=min(PgCnt*PgLn,LnGz);
        
        DtPg=ZrTm(StrIndx:StpInd);
        
        if isempty(DtPg)
            continue;
        end
        ZrMx=repmat(DtPg(:),[1 LnFx]);
        FxMtrx=repmat((FxTms(:))',[length(DtPg) 1]);
        
        Dffmtx=abs(FxMtrx-ZrMx);
        [MnVl,MnInd]=min(Dffmtx,[],1);
        
        try
            PgMn(PgCnt,:)=MnVl;
        catch
            alaki=0;
        end
        PgMnId(PgCnt,:)=MnInd;
        
        
    end
    [PgMnVl,PgInd]=min(PgMn,[],1);
    
    GzFxInd=(PgInd-1)*PgLn+PgMnId((0:PgNm:PgNm*(LnFx-1))+PgInd);
else
    [PgMnVl,GzFxInd]=min(abs(ZrTm-FxTms));
    
end
