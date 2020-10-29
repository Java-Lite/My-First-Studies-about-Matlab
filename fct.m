function [ sisi_miring,luas,keliling ] = segitigasikusiku( sisi_x,sisi_y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sisi_miring=sqrt(sisi_x^2+sisi_y^2);
luas=1/2*sisi_x*sisi_y;
keliling=sisi_x+sisi_y+sisi_miring;

