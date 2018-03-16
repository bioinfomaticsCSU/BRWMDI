License
=========

Copyright (C) 2014 Jianxin Wang(jxwang@mail.csu.edu.cn),Cheng Yan(yancheng01@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn),Cheng Yan(yancheng01@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


Predicting microbe-disease interactions
=================
BRWMDI is one novel computational method, which utilizes bi-random walk algorithm in diease and microbe networks to identify potential microbe-disease interactions.

1.DATA.
1) Disease_Number.xlsx store disease ids and names;

2) Microbe_Number.xlsx store microbe ids and names;

3) Disease-Microbe_Number.xlsx store the known microbe-disease interactions;


2.BRWMDI_5CV.
1) knowndiseasemicrobeinteraction.txt store the known microbe-disease interactions;

2) SymtomBased.mat store the symtom-based similarity of diseases;

3) main.m : main function;

4) BRWMDI5cv.m: function of 5-fold cross validation of BRWMDI; 

5) BRWMDI.m: function of BRWMDI;

6) positiontooverallauc.m: function of calculating the AUC value; 

7) SNF.m: function of similarity network fusion;


3.BRWMDI_LOOCV.
1) knowndiseasemicrobeinteraction.txt store the known microbe-disease interactions;

2) SymtomBased.mat store the symtom-based similarity of diseases;

3) main.m : main function;

4) BRWMDI_LOOCV.m: function of LOOCV of BRWMDI;

5) BRWMDI.m: function of BRWMDI;

6) positiontooverallauc.m: function of calculating the AUC value; 

7) SNF.m: function of similarity network fusion;

All files of Data and Code should be stored in the same folder to run BRWMDI.


