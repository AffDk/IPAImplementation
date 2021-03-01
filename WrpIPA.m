function WrpIPA

% within the BetweenBreaks
clc; clear; close all;
SbFldr='C:\Users\afsamani\OneDrive - Aalborg Universitet\WrokResults\Students\Fall2020\S6';
DtFil1=fullfile(SbFldr,'OutputMat1.mat');
DtFil12=fullfile(SbFldr,'OutputMat2.mat'); % 'OutVrbFx','VrbNmsFx','DtPupil','GzDt','BlnkDt'

DtStrctGm=load(DtFil1);
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

% [PplFxStrInd0,MnVlStr]=FixGzInd(TmPupil0,FixStr);
% [PplFxStpInd0,MnVlStp]=FixGzInd(TmPupil0,FixStp);

PupilDiam=DtPupil.diameter_3d;
PupilDiam(CnfPup<0.6)=nan;
PupilDiam0=PupilDiam(EyId==0);
PupilDiam1=PupilDiam(EyId==1);

Evnts=DtStrctGm.OutVrbGmEvnt.Event; PplTmClEvnt=str2double(string(DtStrctGm.OutVrbGmEvnt.PupilTime));
StrMolLoc=(Evnts=='Mole Spawned');


SpwnTm=PplTmClEvnt(StrMolLoc);
IdKssTm=FndGpSpwn(SpwnTm);

VrTm=[min(TmPupil0); IdKssTm; max(TmPupil0)+0.001];
SgNm=length(VrTm);

% FxLn=length(PplFxStrInd0);
OutTtl=nan(1,SgNm); IPACnt=0;
for Cnt=1:SgNm-1
    Inx=(TmPupil0>=VrTm(Cnt) & TmPupil0<VrTm(Cnt+1));
    X=PupilDiam0(Inx);
    TmP=TmPupil0(Inx);
    
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


function TmOut=FndGpSpwn(X)

Xn=[X(1); X];
Dx=diff(Xn);
% Id=find(Dx>9);

TmOut=X(Dx>9);


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
