%V1toMTmodel.m
%
%Can also be run in functional form
%function threshold = V1toMTmodel(stimIn, distortionFactor, sigmaForNormalization)

close all;

clear all;

%The different stimuli presented to the model:
possibleStims=[45 90];

%Parameters of the model (affecting distribution of receptive fields):
distortionFactor=[1.2 2];

%Parameters affecting normalization:
sigmaForNormalization=[2 2];

adaptationGain=0.1; %Parameter affecting MT adaptation characteristics:

%For timing:
tic

%Repetitions of stimulus presentation:
numReps=100;

sigmaStart=pi/8; %Bandwidth of the tuning
sizeBankMT=32; %Cells in the MT layer
sizeBankV1=sizeBankMT*12; %Cells in the input layer


aStart=100; %Scaling factor, describing the stimulus induced gain on firing rate
bStart=10; %Spontaneous firing in the cells (in Hz)

%Same for MT layer:
aMT=100;
bMT=10;

%Make the V1 cell bank:

RFdensity=distortionFactor(1)*(1-cos([8*pi/sizeBankV1:8*pi/sizeBankV1:8*pi]))+0.012;
RFdensity=(RFdensity./sum(RFdensity)).*(2*pi);
makeTprefV1=[0];
aConn=aStart;
bConn=bStart;
sigmaV1=((1-cos(8*pi/sizeBankV1:8*pi/sizeBankV1:8*pi)))*distortionFactor(2)+sigmaStart;

for cell=1:sizeBankV1

    %   disp(['making cell # ' num2str(cell)])

    if cell>1
        TprefInc=RFdensity(cell);
        makeTprefV1=[makeTprefV1 makeTprefV1(cell-1)+TprefInc];
    end

    aV1(cell)=aStart;
    bV1(cell)=bStart;
    tuningV1(cell,:)=aV1(cell).*vonMises([2*pi/sizeBankV1:2*pi/sizeBankV1:2*pi],makeTprefV1(cell),sigmaV1(cell))+bV1(cell);
    %tuningV1(cell,:)=tuningV1(cell,:)./sum(tuningV1(cell,:));
end


maxTuneV1=max(max(tuningV1));

for i=1:sizeBankV1
    tuningV1(:,i)=tuningV1(:,i)./sum(tuningV1(:,i));
end

%Normalize:
tuningV1=tuningV1.*(maxTuneV1/max(max(tuningV1)));


%Set up the connectivity:

sigmaConn=ones(1,sizeBankMT).*sigmaStart;

helperFunc=linspace(0,2*pi,sizeBankMT);
whichCellV1toMT=linspace(0,2*pi,sizeBankMT);

%Make the MT cells:
for cell=1:sizeBankMT

    makeTprefMT(cell)=makeTprefV1(ceil(cell*sizeBankV1/sizeBankMT));
    whichCellV1toMT(cell)=helperFunc(cell);

    connV1toMT(cell,:)=vonMises([2*pi/sizeBankV1:2*pi/sizeBankV1:2*pi],whichCellV1toMT(cell),sigmaConn(cell));
    connV1toMT(cell,:)=connV1toMT(cell,:)-min(connV1toMT(cell,:));
    connV1toMT(cell,:)=connV1toMT(cell,:)./max(connV1toMT(cell,:)).*aMT + bMT;

    connMatrixMT(cell,:)=-1*vonMises([2*pi/sizeBankMT:2*pi/sizeBankMT:2*pi],makeTprefMT(cell),sigmaConn(cell));
    connMatrixMT(cell,:)=connMatrixMT(cell,:);%./sum(connMatrixMT(cell,:));

end

maxTuneMT=max(max(connV1toMT));

for i=1:sizeBankV1
    connV1toMT(:,i)=connV1toMT(:,i)./sum(connV1toMT(:,i));
end

%Normalize:
connV1toMT=connV1toMT.*(maxTuneMT/max(max(connV1toMT)));

tuningMT=connV1toMT*tuningV1;
tuningMT=tuningMT.*(maxTuneMT/max(max(tuningMT)));

MTfwhm=[];

%Determine the full width at half maximum for the MT cells:
for k=1:sizeBankMT
    MTfwhm=[MTfwhm find_fwhm(tuningMT(k,:))];
end


%Loop over the possible stimuli:
for kkk=1:2
    stimIn=possibleStims(kkk);

    %For each stimulus, present it at different variances:
    stimSTD=[0 2.8125 5.625 11.25 16.875 22.5 45 90].*((2*pi/360));
    stimDirections=[stimIn].*((2*pi/360));


    vectorStrength=zeros(numReps,length(stimSTD));
    toPlot=zeros(length(stimSTD),sizeBankMT);


    for rep=1:numReps
        %     if ~mod(rep,20)
        %         disp (['Repetition: ' num2str(rep)])
        %     end
        for stimIndex=1:length(stimDirections)
            for stdIndex=1:length(stimSTD)
                for cell=1:sizeBankV1
                    stim=[randn*stimSTD(stdIndex)+stimDirections(stimIndex)];
                    stimIndexInTuning=mod(ceil(circularize(stim)/(2*pi)*sizeBankV1),sizeBankV1+1);
                    rememberStim(stimIndex,stdIndex,cell)=stimIndexInTuning;
                    ratesV1(stdIndex,cell)=tuningV1(cell,stimIndexInTuning);
                end

                activityV1F1(stdIndex,:)=poissrnd(ratesV1(stdIndex,:)); %Firing for the first stimulus
                activityV1F2(stdIndex,:)=poissrnd(ratesV1(stdIndex,:)); %Firing for the second stimulus


                for cell=1:sizeBankV1
                    squaredV1F1(stdIndex,cell)=(activityV1F1(stdIndex,cell))^2;
                    squaredV1F2(stdIndex,cell)=(activityV1F2(stdIndex,cell))^2;
                end

                for cell=1:sizeBankV1
                    outaV1F1(stdIndex,cell)=(squaredV1F1(stdIndex,cell))./sum(squaredV1F1(stdIndex,:));
                    outaV1F2(stdIndex,cell)=(squaredV1F2(stdIndex,cell))./sum(squaredV1F2(stdIndex,:));
                end


                intoMTF1(stdIndex,:)=outaV1F1(stdIndex,:)./(outaV1F1(stdIndex,:)+sigmaForNormalization(1)) ; %nakaRushton(outaV1(stdIndex,:),nakaRushtonA1,nakaRushtonN1,nakaRushtonC1,nakaRushtonB1);
                intoMTF2(stdIndex,:)=outaV1F2(stdIndex,:)./(outaV1F2(stdIndex,:)+sigmaForNormalization(1)) ; %nakaRushton(outaV1(stdIndex,:),nakaRushtonA1,nakaRushtonN1,nakaRushtonC1,nakaRushtonB1);

                ratesMTF1(stdIndex,:)=intoMTF1(stdIndex,:)*connV1toMT';
                ratesMTF2(stdIndex,:)=intoMTF2(stdIndex,:)*connV1toMT';

                activityMTF1(stdIndex,:)=poissrnd(ratesMTF1(stdIndex,:));
                activityMTF2(stdIndex,:)=poissrnd(ratesMTF2(stdIndex,:));


                %Normalization in MT:
                %             for cell=1:sizeBankMT
                %                 squaredMTF1(stdIndex,cell)=(activityMTF1(stdIndex,cell));%^2;
                %                 squaredMTF2(stdIndex,cell)=(activityMTF2(stdIndex,cell));%^2;
                %
                %             end
                %
                %             for cell=1:sizeBankMT
                %                 outaMTF1(stdIndex,cell)=(squaredMTF1(stdIndex,cell));%./sum(squaredMTF1(stdIndex,:));
                %                 outaMTF2(stdIndex,cell)=(squaredMTF2(stdIndex,cell));%./sum(squaredMTF2(stdIndex,:));
                %
                %             end
                %
                %             MTbackF1(stdIndex,:)=outaMTF1(stdIndex,:);%./(outaMTF1(stdIndex,:)+sigmaForNormalization(2)); %nakaRushton(outaMT(stdIndex,:),nakaRushtonA2,nakaRushtonN2,nakaRushtonC2,nakaRushtonB2);
                %             MTbackF2(stdIndex,:)=outaMTF2(stdIndex,:);%./(outaMTF2(stdIndex,:)+sigmaForNormalization(2)); %nakaRushton(outaMT(stdIndex,:),nakaRushtonA2,nakaRushtonN2,nakaRushtonC2,nakaRushtonB2);

                global fit_makeTpref;
                fit_makeTpref=makeTprefMT;
                global fit_sigmaConn;
                fit_sigmaConn=MTfwhm;
                global fit_activityMT;


                fit_activityMT=activityMTF1(stdIndex,:);

                %Fit the direction coded by the population of MT cells:

                thetaF1(stdIndex,rep) = fminsearch(@logLikeTheta, stimDirections(stimIndex));

                fit_activityMT=activityMTF2(stdIndex,:);

                thetaF2(stdIndex,rep) = fminsearch(@logLikeTheta, stimDirections(stimIndex));



            end

        end


    end



    %Determine the threshold for each level of variance (as the median of
    %the difference distribution):
    for stdIndex=1:length(stimSTD)

        threshold(stdIndex)=median(abs(thetaF1(stdIndex,:)-thetaF2(stdIndex,:)))*(360/(2*pi));

    end

    %Plot the results:
    x=[0 2.8125 5.625 11.25 16.875 22.5 45 90];

    if kkk==1
        plot(x,threshold)
        hold on

    else
        plot(x,threshold,'r')
    end

    %Determine how much time it took:
    timeItTook=toc/60


end
