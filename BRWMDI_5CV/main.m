clear;
for cv=1:100
    BRWMDI5cv();
    overallauc(cv)=positiontooverallauc();
end
save overallauc overallauc
a=mean(overallauc)
b=std(overallauc)

