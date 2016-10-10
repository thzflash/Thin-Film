

%TransferMatrix3 calculates transmission, reflection and absorption of a multilayer of planar homogenous films
%inputs:
%angle of incidence: thetai
%wavelength of incident light: lambda
%thicknesses of the layers: h
%refractive index of the layers: n (may be absorbing and dispersive, but this may require a subfunction)
%polarization: s or p (need one calculation for each for unpolarized light)

cd D:\Dropbox\AfterBurner\Matlab\2016_Thin_Film

clear all
close all
clc; 

thetai=45; %angle of incidence (degrees)

h=[NaN,50*1000,NaN]; %film thicknesses in nm, equal in length to n, start and end with NaN
% h=[NaN,150*1000,150*1000,NaN]; %film thicknesses in nm, equal in length to n, start and end with NaN

pol=1; %polarization, 1 for p and 0 for s
% n=[1.52,sqrt(Au(lambda)),1.4,1.33]; %refractive index data, NaN for frequency dependence

% n = [1,1.30,1]; %%% Pellicle from paper ... Png_2015
% n = [1,3.50,1];   %%% Silicon
n = [1,2.38,1];    %%% Diamond
% n = [1,1.54,1];   %%% HDPE High-Density-PolyEthilene

LAMBDA = (1 : 0.1 : 300) * 1000;

FR = zeros(1,length(LAMBDA));
FR = zeros(1,length(LAMBDA));
FA = zeros(1,length(LAMBDA));

for Ind = 1 : 1 : length(LAMBDA); %vacuum wavelength (nm)

 lambda = LAMBDA(Ind);   
    
[FR(Ind),FT(Ind),FA(Ind)]=Fresnel(lambda,thetai,h,n,pol);

% pause(0.01);
% disp([num2str(lambda*100) '% done...'])
%%% plot results:
xlabel('Wavelength [um]')
% ylabel('Fresnel coefficient')
% title('Fresnel coefficients for transmission (blue), reflection (red) and absorption (green)')
end
figure(1)
hold on
plot(LAMBDA/1000,FT,'b','LineWidth',2);
plot(LAMBDA/1000,FA,'g','LineWidth',2);
plot(LAMBDA/1000,FR,'-r','LineWidth',2);

axis tight;
title('Fresnel coefficients for transmission (blue), reflection (red) and absorption (green)')
grid on;
