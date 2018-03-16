
clear;
BRWMDI_LOOCV();
overallauc=positiontooverallauc();
save overallauc overallauc
a=mean(overallauc)
b=std(overallauc)