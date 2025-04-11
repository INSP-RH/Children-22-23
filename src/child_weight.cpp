//
//  child_weight.cpp
//
//  This is a function that calculates
//  weight change for children using the dynamic
//  weight model by Kevin D. Hall et al. and
//  Runge Kutta method to solve the ODE system.
//
//  Input:
//  age             .-  Years since individual first arrived to Earth
//  sex             .-  Either 1 = "female" or 0 = "male"
//  FFM             .-  Fat Free Mass (kg) of the individual
//  FM              .-  Fat Mass (kg) of the individual
//  input_EIntake   .-  Energy intake (kcal) of individual per day
//  days            .-  Days to model (integer)
//  dt              .-  Time step used to solve the ODE system numerically
//  K               .-  Richardson parameter
//  Q               .-  Richardson parameter
//  A               .-  Richardson parameter
//  B               .-  Richardson parameter
//  nu              .-  Richardson parameter
//  C               .-  Richardson parameter
//  Note:
//  Weight = FFM + FM. No extracellular fluid or glycogen is considered
//  Please see child_weight.hpp for additional information
//
//  Authors:
//  Dalia Camacho-García-Formentí
//  Rodrigo Zepeda-Tello
//
// References:
//
//  Deurenberg, Paul, Jan A Weststrate, and Jaap C Seidell. 1991. “Body Mass Index as a Measure of Body Fatness:
//      Age-and Sex-Specific Prediction Formulas.” British Journal of Nutrition 65 (2). Cambridge University Press: 105–14.
//
//  Ellis, Kenneth J, Roman J Shypailo, Steven A Abrams, and William W Wong. 2000. “The Reference Child and Adolescent Models of
//      Body Composition: A Contemporary Comparison.” Annals of the New York Academy of Sciences 904 (1). Wiley Online Library: 374–82.
//
//  Fomon, Samuel J, Ferdinand Haschke, Ekhard E Ziegler, and Steven E Nelson. 1982.
//      “Body Composition of Reference Children from Birth to Age 10 Years.” The American Journal of
//      Clinical Nutrition 35 (5). Am Soc Nutrition: 1169–75.
//
//  Hall, Kevin D, Nancy F Butte, Boyd A Swinburn, and Carson C Chow. 2013. “Dynamics of Childhood Growth
//      and Obesity: Development and Validation of a Quantitative Mathematical Model.” The Lancet Diabetes & Endocrinology 1 (2). Elsevier: 97–105.
//
//  Haschke, F. 1989. “Body Composition During Adolescence.” Body Composition Measurements in Infants and Children.
//      Ross Laboratories Columbus, OH, 76–83.
//
//  Katan, Martijn B, Janne C De Ruyter, Lothar DJ Kuijper, Carson C Chow, Kevin D Hall, and Margreet R Olthof. 2016.
//      “Impact of Masked Replacement of Sugar-Sweetened with Sugar-Free Beverages on Body Weight Increases with Initial Bmi:
//      Secondary Analysis of Data from an 18 Month Double–Blind Trial in Children.” PloS One 11 (7). Public Library of Science: e0159771.
//
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute,
// sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
// is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies
// or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
// BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//----------------------------------------------------------------------------------------


#include "child_weight.h"

//Default (classic) constructor for energy matrix
Child::Child(NumericVector input_age, NumericVector input_sex, NumericVector input_bmiCat, NumericVector input_FFM, NumericVector input_FM, NumericMatrix input_EIntake,
             double input_dt, bool checkValues, double input_referenceValues){
    age   = input_age;
    sex   = input_sex;
    bmiCat = input_bmiCat;
    FM    = input_FM;
    FFM   = input_FFM;
    dt    = input_dt;
    EIntake = input_EIntake;
    check = checkValues;
    generalized_logistic = false;
    referenceValues = input_referenceValues;
    build();
}

//Constructor which uses Richard's curve with the parameters of https://en.wikipedia.org/wiki/Generalised_logistic_function
Child::Child(NumericVector input_age, NumericVector input_sex,  NumericVector input_bmiCat, NumericVector input_FFM, NumericVector input_FM, double input_K,
             double input_Q, double input_A, double input_B, double input_nu, double input_C, 
             double input_dt, bool checkValues, double input_referenceValues){
    age   = input_age;
    sex   = input_sex;
    bmiCat = input_bmiCat;
    FM    = input_FM;
    FFM   = input_FFM;
    dt    = input_dt;
    K_logistic = input_K;
    A_logistic = input_A;
    Q_logistic = input_Q;
    B_logistic = input_B;
    nu_logistic = input_nu;
    C_logistic = input_C;
    check = checkValues;
    referenceValues = input_referenceValues;
    generalized_logistic = true;
    ;
    build();
}

Child::~Child(void){
    
}

void Child::build(){
    getParameters();
}

//General function for expressing growth and eb terms
NumericVector Child::general_ode(NumericVector t, NumericVector input_A, NumericVector input_B,
                                 NumericVector input_D, NumericVector input_tA,
                                 NumericVector input_tB, NumericVector input_tD,
                                 NumericVector input_tauA, NumericVector input_tauB,
                                 NumericVector input_tauD){
    
    return input_A*exp(-(t-input_tA)/input_tauA ) +
            input_B*exp(-0.5*pow((t-input_tB)/input_tauB,2)) +
            input_D*exp(-0.5*pow((t-input_tD)/input_tauD,2));
}

NumericVector Child::Growth_dynamic(NumericVector t){
    return general_ode(t, A, B, D, tA, tB, tD, tauA, tauB, tauD);
}

NumericVector Child::Growth_impact(NumericVector t){
    return general_ode(t, A1, B1, D1, tA1, tB1, tD1, tauA1, tauB1, tauD1);
}

NumericVector Child::EB_impact(NumericVector t){
    return general_ode(t, A_EB, B_EB, D_EB, tA_EB, tB_EB, tD_EB, tauA_EB, tauB_EB, tauD_EB);
}

NumericVector Child::cRhoFFM(NumericVector input_FFM){
    return 4.3*input_FFM + 837.0;
}

NumericVector Child::cP(NumericVector FFM, NumericVector FM){
    NumericVector rhoFFM = cRhoFFM(FFM);
    NumericVector C      = 10.4 * rhoFFM / rhoFM;
    return C/(C + FM);
}

NumericVector Child::Delta(NumericVector t){
    return deltamin + (deltamax - deltamin)*(1.0 / (1.0 + pow((t / P),h)));
}

NumericVector Child::FFMReference(NumericVector t){ 
  /*  return ffm_beta0 + ffm_beta1*t; */
NumericVector under = ifelse(bmiCat == 1, 1.0, 0.0);
NumericVector normales = ifelse(bmiCat == 2, 1.0, 0.0);
NumericVector over = ifelse(bmiCat == 3, 1.0, 0.0);
NumericVector obese = ifelse(bmiCat == 4, 1.0, 0.0);

NumericMatrix ffm_ref(17,nind);
  if(referenceValues == 0){
  // -------------------------- Mean values
ffm_ref(0,_)   = 10.134*(1-sex)+9.477*sex;       // 2 years old
ffm_ref(1,_)   = 12.099*(1 - sex) + 11.494*sex;    // 3 years old
ffm_ref(2,_)   = 14.0*(1 - sex) + 13.2*sex;        // 4 years old
ffm_ref(3,_)   = 15.72*(1 - sex) + 14.86*sex;      // 5 years old
ffm_ref(4,_)   = under*(12.7942*(1-sex) + 13.7957*sex) + normales*(17.0238*(1-sex) + 15.2337*sex) + over*(19.3070*(1-sex) + 17.7866*sex) + obese*(22.2248*(1-sex) + 21.2170*sex);   // 6 years old
ffm_ref(5,_)   = under*(17.8106*(1-sex) + 18.4835*sex) + normales*(19.0775*(1-sex) + 17.5198*sex) + over*(20.3344*(1-sex) + 18.9406*sex) + obese*(23.1765*(1-sex) + 22.2733*sex);   // 7 years old
ffm_ref(6,_)   = under*(20.3597*(1-sex) + 18.5363*sex) + normales*(20.4774*(1-sex) + 19.6317*sex) + over*(22.1128*(1-sex) + 21.6080*sex) + obese*(25.8151*(1-sex) + 25.1641*sex);   // 8 years old
ffm_ref(7,_)   = under*(19.3668*(1-sex) + 17.0314*sex) + normales*(22.3768*(1-sex) + 21.3680*sex) + over*(26.7714*(1-sex) + 26.1791*sex) + obese*(31.3143*(1-sex) + 30.1484*sex);   // 9 years old
ffm_ref(8,_)   = under*(22.7116*(1-sex) + 23.3608*sex) + normales*(26.3269*(1-sex) + 26.8882*sex) + over*(31.0871*(1-sex) + 32.8504*sex) + obese*(36.9982*(1-sex) + 35.4395*sex);   // 10 years old
ffm_ref(9,_)   = under*(24.9001*(1-sex) + 20.8802*sex) + normales*(28.0114*(1-sex) + 28.9720*sex) + over*(33.4168*(1-sex) + 34.8732*sex) + obese*(39.9324*(1-sex) + 40.2945*sex);   // 11 years old
ffm_ref(10,_)   = under*(27.9345*(1-sex) + 27.7570*sex) + normales*(32.0987*(1-sex) + 33.6664*sex) + over*(38.4380*(1-sex) + 39.1053*sex) + obese*(44.0899*(1-sex) + 44.5438*sex);   // 12 years old
ffm_ref(11,_)   = under*(24.2164*(1-sex) + 26.9376*sex) + normales*(36.2130*(1-sex) + 36.3192*sex) + over*(42.8030*(1-sex) + 41.4690*sex) + obese*(48.5364*(1-sex) + 46.8082*sex);   // 13 years old
ffm_ref(12,_)   = under*(34.5353*(1-sex) + 29.2222*sex) + normales*(40.8515*(1-sex) + 37.3162*sex) + over*(48.0478*(1-sex) + 45.1033*sex) + obese*(56.8601*(1-sex) + 47.9190*sex);   // 14 years old
ffm_ref(13,_)   = under*(36.2855*(1-sex) + 34.1242*sex) + normales*(43.3825*(1-sex) + 40.3495*sex) + over*(51.0714*(1-sex) + 45.8652*sex) + obese*(58.6970*(1-sex) + 50.1799*sex);   // 15 years old
ffm_ref(14,_)   = under*(38.0709*(1-sex) + 36.8249*sex) + normales*(45.8051*(1-sex) + 40.0525*sex) + over*(56.1791*(1-sex) + 45.0247*sex) + obese*(60.6064*(1-sex) + 52.0168*sex);   // 16 years old
ffm_ref(15,_)   = under*(40.4601*(1-sex) + 36.5201*sex) + normales*(48.6039*(1-sex) + 42.0138*sex) + over*(57.0505*(1-sex) + 48.5732*sex) + obese*(60.3816*(1-sex) + 53.1534*sex);   // 17 years old
ffm_ref(16,_)   = under*(42.9297*(1-sex) + 31.0623*sex) + normales*(49.4997*(1-sex) + 41.3371*sex) + over*(58.3226*(1-sex) + 49.6194*sex) + obese*(57.2369*(1-sex) + 50.6037*sex);   // 18 years old

  }


  if(referenceValues == 1){
    // -------------------------- Median values
ffm_ref(0,_)   = 10.134*(1-sex)+9.477*sex;       // 2 years old
ffm_ref(1,_)   = 12.099*(1 - sex) + 11.494*sex;    // 3 years old
ffm_ref(2,_)   = 14.0*(1 - sex) + 13.2*sex;        // 4 years old
ffm_ref(3,_)   = 15.72*(1 - sex) + 14.86*sex;      // 5 years old
ffm_ref(4,_)   = under*(14.4641*(1-sex) + 13.8627*sex) + normales*(17.1430*(1-sex) + 15.1282*sex) + over*(19.2280*(1-sex) + 17.6859*sex) + obese*(21.9501*(1-sex) + 20.4992*sex);   // 6 years old
ffm_ref(5,_)   = under*(16.3729*(1-sex) + 16.6347*sex) + normales*(18.2285*(1-sex) + 17.2507*sex) + over*(21.7099*(1-sex) + 20.0341*sex) + obese*(24.9713*(1-sex) + 23.4162*sex);   // 7 years old
ffm_ref(6,_)   = under*(18.0019*(1-sex) + 17.2583*sex) + normales*(19.9148*(1-sex) + 19.4286*sex) + over*(24.6404*(1-sex) + 22.1758*sex) + obese*(27.4774*(1-sex) + 26.8346*sex);   // 8 years old
ffm_ref(7,_)   = under*(19.2548*(1-sex) + 17.5150*sex) + normales*(21.9058*(1-sex) + 21.2721*sex) + over*(26.5243*(1-sex) + 25.6952*sex) + obese*(30.8636*(1-sex) + 29.2900*sex);   // 9 years old
ffm_ref(8,_)   = under*(22.7116*(1-sex) + 23.3608*sex) + normales*(26.3269*(1-sex) + 26.8882*sex) + over*(31.0871*(1-sex) + 32.8504*sex) + obese*(36.9982*(1-sex) + 35.4395*sex);   // 10 years old
ffm_ref(9,_)   = under*(24.9001*(1-sex) + 20.8802*sex) + normales*(28.0114*(1-sex) + 28.9720*sex) + over*(33.4168*(1-sex) + 34.8732*sex) + obese*(39.9324*(1-sex) + 40.2945*sex);   // 11 years old
ffm_ref(10,_)   = under*(27.9345*(1-sex) + 27.7570*sex) + normales*(32.0987*(1-sex) + 33.6664*sex) + over*(38.4380*(1-sex) + 39.1053*sex) + obese*(44.0899*(1-sex) + 44.5438*sex);   // 12 years old
ffm_ref(11,_)   = under*(24.2164*(1-sex) + 26.9376*sex) + normales*(36.2130*(1-sex) + 36.3192*sex) + over*(42.8030*(1-sex) + 41.4690*sex) + obese*(48.5364*(1-sex) + 46.8082*sex);   // 13 years old
ffm_ref(12,_)   = under*(34.5353*(1-sex) + 29.2222*sex) + normales*(40.8515*(1-sex) + 37.3162*sex) + over*(48.0478*(1-sex) + 45.1033*sex) + obese*(56.8601*(1-sex) + 47.9190*sex);   // 14 years old
ffm_ref(13,_)   = under*(36.2855*(1-sex) + 34.1242*sex) + normales*(43.3825*(1-sex) + 40.3495*sex) + over*(51.0714*(1-sex) + 45.8652*sex) + obese*(58.6970*(1-sex) + 50.1799*sex);   // 15 years old
ffm_ref(14,_)   = under*(38.0709*(1-sex) + 36.8249*sex) + normales*(45.8051*(1-sex) + 40.0525*sex) + over*(56.1791*(1-sex) + 45.0247*sex) + obese*(60.6064*(1-sex) + 52.0168*sex);   // 16 years old
ffm_ref(15,_)   = under*(40.4601*(1-sex) + 36.5201*sex) + normales*(48.6039*(1-sex) + 42.0138*sex) + over*(57.0505*(1-sex) + 48.5732*sex) + obese*(60.3816*(1-sex) + 53.1534*sex);   // 17 years old
ffm_ref(16,_)   = under*(42.9297*(1-sex) + 31.0623*sex) + normales*(49.4997*(1-sex) + 41.3371*sex) + over*(58.3226*(1-sex) + 49.6194*sex) + obese*(57.2369*(1-sex) + 50.6037*sex);   // 18 years old
  }

NumericVector ffm_ref_t(nind);
int jmin;
int jmax;
double diff;
for(int i=0;i<nind;i++){
  if(t(i)>=18.0){
    ffm_ref_t(i)=ffm_ref(16,i);
  }else{
    jmin=floor(t(i));
    jmin=std::max(jmin,2);
    jmin=jmin-2;
    jmax= std::min(jmin+1,17);
    diff= t(i)-floor(t(i));
    ffm_ref_t(i)=ffm_ref(jmin,i)+diff*(ffm_ref(jmax,i)-ffm_ref(jmin,i));
  } 
}
return ffm_ref_t;
}

NumericVector Child::FMReference(NumericVector t){
   /* return fm_beta0 + fm_beta1*t;*/
NumericVector under = ifelse(bmiCat == 1, 1.0, 0.0);
NumericVector normales = ifelse(bmiCat == 2, 1.0, 0.0);
NumericVector over = ifelse(bmiCat == 3, 1.0, 0.0);
NumericVector obese = ifelse(bmiCat == 4, 1.0, 0.0);

NumericMatrix fm_ref(17,nind);
 if(referenceValues == 0){
  // ---------------------------------------- Mean values

fm_ref(0,_)   = 2.456*(1-sex)+ 2.433*sex;       // 2 years old
fm_ref(1,_)   = 2.576*(1 - sex) + 2.606*sex;    // 3 years old
fm_ref(2,_)   = 2.7*(1 - sex) + 2.8*sex;        // 4 years old
fm_ref(3,_)   = 3.66*(1 - sex) + 4.47*sex;      // 5 years old
fm_ref(4,_)   = under*(1.7764*(1-sex) + 2.5951*sex) + normales*(3.4540*(1-sex) + 3.8303*sex) + over*(4.8055*(1-sex) + 5.7014*sex) + obese*(7.9672*(1-sex) + 9.3883*sex);   // 6 years old
fm_ref(5,_)   = under*(2.3398*(1-sex) + 2.8164*sex) + normales*(3.5859*(1-sex) + 4.2782*sex) + over*(5.4625*(1-sex) + 6.5960*sex) + obese*( 8.4350*(1-sex) + 10.4148*sex);   // 7 years old
fm_ref(6,_)   = under*(3.2767*(1-sex) + 3.0828*sex) + normales*(4.1138*(1-sex) + 5.2226*sex) + over*(5.5455*(1-sex) + 7.3667*sex) + obese*( 9.3266*(1-sex) + 12.0550*sex);   // 8 years old
fm_ref(7,_)   = under*(2.3902*(1-sex) + 2.6538*sex) + normales*(4.1705*(1-sex) + 5.0218*sex) + over*(6.6958*(1-sex) + 8.6945*sex) + obese*(11.5896*(1-sex) + 14.1436*sex);   // 9 years old
fm_ref(8,_)   = under*(2.8822*(1-sex) + 3.1199*sex) + normales*(5.0182*(1-sex) + 5.5597*sex) + over*(8.6036*(1-sex) + 9.3996*sex) + obese*(16.5022*(1-sex) + 15.0802*sex);   // 10 years old
fm_ref(9,_)   = under*(3.1730*(1-sex) + 2.5592*sex) + normales*(5.5337*(1-sex) + 6.1030*sex) + over*( 8.9291*(1-sex) + 10.9054*sex) + obese*(17.9139*(1-sex) + 20.1630*sex);   // 11 years old
fm_ref(10,_)   = under*(3.4905*(1-sex) + 3.7443*sex) + normales*(6.4588*(1-sex) + 7.2476*sex) + over*(11.5237*(1-sex) + 12.6802*sex) + obese*(20.4754*(1-sex) + 23.3280*sex);   // 12 years old
fm_ref(11,_)   = under*(2.6404*(1-sex) + 3.2124*sex) + normales*(7.0512*(1-sex) + 8.5389*sex) + over*(12.9732*(1-sex) + 14.3787*sex) + obese*(24.2060*(1-sex) + 24.3422*sex);   // 13 years old
fm_ref(12,_)   = under*(3.9475*(1-sex) + 3.9076*sex) + normales*(8.2349*(1-sex) + 8.7253*sex) + over*(15.7762*(1-sex) + 16.4665*sex) + obese*(30.8151*(1-sex) + 28.9186*sex);   // 14 years old
fm_ref(13,_)   = under*(4.3234*(1-sex) + 3.8050*sex) + normales*(8.6141*(1-sex) + 9.9317*sex) + over*(17.2515*(1-sex) + 17.9368*sex) + obese*(31.7228*(1-sex) + 29.0934*sex);   // 15 years old
fm_ref(14,_)   = under*(4.4563*(1-sex) + 4.9751*sex) + normales*(9.4644*(1-sex) + 9.8079*sex) + over*(20.0891*(1-sex) + 16.2521*sex) + obese*(32.2330*(1-sex) + 31.7539*sex);   // 16 years old
fm_ref(15,_)   = under*(4.2495*(1-sex) + 4.2212*sex) + normales*(10.1763*(1-sex) + 10.1097*sex) + over*(19.9104*(1-sex) + 18.5442*sex) + obese*(30.9616*(1-sex) + 33.4933*sex);   // 17 years old
fm_ref(16,_)   = under*(4.7453*(1-sex) + 3.1877*sex) + normales*(11.1200*(1-sex) +  9.3930*sex) + over*(20.9956*(1-sex) + 19.9462*sex) + obese*(31.1828*(1-sex) + 30.6175*sex);   // 18 years old

 }

 if(referenceValues == 1){
  // ---------------------------------------- Median values

fm_ref(0,_)   = 2.456*(1-sex)+ 2.433*sex;       // 2 years old
fm_ref(1,_)   = 2.576*(1 - sex) + 2.606*sex;    // 3 years old
fm_ref(2,_)   = 2.7*(1 - sex) + 2.8*sex;        // 4 years old
fm_ref(3,_)   = 3.66*(1 - sex) + 4.47*sex;      // 5 years old
fm_ref(4,_)   = under*(2.0359*(1-sex) + 2.5660*sex) + normales*(3.4642*(1-sex) + 3.7042*sex) + over*(4.6220*(1-sex) + 5.6735*sex) + obese*(7.1058*(1-sex) + 8.7339*sex);   // 6 years old
fm_ref(5,_)   = under*(2.3771*(1-sex) + 2.9560*sex) + normales*(3.6030*(1-sex) + 4.1865*sex) + over*(5.5651*(1-sex) + 6.4374*sex) + obese*(8.0501*(1-sex) + 9.3100*sex);   // 7 years old
fm_ref(6,_)   = under*(2.1231*(1-sex) + 3.0917*sex) + normales*(3.6729*(1-sex) + 4.8531*sex) + over*(5.8971*(1-sex) + 7.0172*sex) + obese*( 8.9372*(1-sex) + 11.5469*sex);   // 8 years old
fm_ref(7,_)   = under*(2.4068*(1-sex) + 2.9027*sex) + normales*(4.0597*(1-sex) + 4.8707*sex) + over*(6.5720*(1-sex) + 8.7112*sex) + obese*(10.8084*(1-sex) + 12.7559*sex);   // 9 years old
fm_ref(8,_)   = under*(2.8822*(1-sex) + 3.1199*sex) + normales*(5.0182*(1-sex) + 5.5597*sex) + over*(8.6036*(1-sex) + 9.3996*sex) + obese*(16.5022*(1-sex) + 15.0802*sex);   // 10 years old
fm_ref(9,_)   = under*(3.1730*(1-sex) + 2.5592*sex) + normales*(5.5337*(1-sex) + 6.1030*sex) + over*( 8.9291*(1-sex) + 10.9054*sex) + obese*(17.9139*(1-sex) + 20.1630*sex);   // 11 years old
fm_ref(10,_)   = under*(3.4905*(1-sex) + 3.7443*sex) + normales*(6.4588*(1-sex) + 7.2476*sex) + over*(11.5237*(1-sex) + 12.6802*sex) + obese*(20.4754*(1-sex) + 23.3280*sex);   // 12 years old
fm_ref(11,_)   = under*(2.6404*(1-sex) + 3.2124*sex) + normales*(7.0512*(1-sex) + 8.5389*sex) + over*(12.9732*(1-sex) + 14.3787*sex) + obese*(24.2060*(1-sex) + 24.3422*sex);   // 13 years old
fm_ref(12,_)   = under*(3.9475*(1-sex) + 3.9076*sex) + normales*(8.2349*(1-sex) + 8.7253*sex) + over*(15.7762*(1-sex) + 16.4665*sex) + obese*(30.8151*(1-sex) + 28.9186*sex);   // 14 years old
fm_ref(13,_)   = under*(4.3234*(1-sex) + 3.8050*sex) + normales*(8.6141*(1-sex) + 9.9317*sex) + over*(17.2515*(1-sex) + 17.9368*sex) + obese*(31.7228*(1-sex) + 29.0934*sex);   // 15 years old
fm_ref(14,_)   = under*(4.4563*(1-sex) + 4.9751*sex) + normales*(9.4644*(1-sex) + 9.8079*sex) + over*(20.0891*(1-sex) + 16.2521*sex) + obese*(32.2330*(1-sex) + 31.7539*sex);   // 16 years old
fm_ref(15,_)   = under*(4.2495*(1-sex) + 4.2212*sex) + normales*(10.1763*(1-sex) + 10.1097*sex) + over*(19.9104*(1-sex) + 18.5442*sex) + obese*(30.9616*(1-sex) + 33.4933*sex);   // 17 years old
fm_ref(16,_)   = under*(4.7453*(1-sex) + 3.1877*sex) + normales*(11.1200*(1-sex) +  9.3930*sex) + over*(20.9956*(1-sex) + 19.9462*sex) + obese*(31.1828*(1-sex) + 30.6175*sex);   // 18 years old
 }



  
NumericVector fm_ref_t(nind);
int jmin;
int jmax;
double diff;
for(int i=0;i<nind;i++){
  if(t(i)>=18.0){
    fm_ref_t(i)=fm_ref(16,i);
  }else{
    jmin=floor(t(i));
    jmin=std::max(jmin,2);
    jmin=jmin-2;
    jmax= std::min(jmin+1,17);
    diff= t(i)-floor(t(i));
    fm_ref_t(i)=fm_ref(jmin,i)+diff*(fm_ref(jmax,i)-fm_ref(jmin,i));
  } 
}
return fm_ref_t;
}

NumericVector Child::IntakeReference(NumericVector t){
    NumericVector EB      = EB_impact(t);
    NumericVector FFMref  = FFMReference(t);
    NumericVector FMref   = FMReference(t);
    NumericVector delta   = Delta(t);
    NumericVector growth  = Growth_dynamic(t);
    NumericVector p       = cP(FFMref, FMref);
    NumericVector rhoFFM  = cRhoFFM(FFMref);
    return EB + K + (22.4 + delta)*FFMref + (4.5 + delta)*FMref +
                230.0/rhoFFM*(p*EB + growth) + 180.0/rhoFM*((1-p)*EB-growth);
}

NumericVector Child::Expenditure(NumericVector t, NumericVector FFM, NumericVector FM){
    NumericVector delta     = Delta(t);
    NumericVector Iref      = IntakeReference(t);
    NumericVector Intakeval = Intake(t);
    NumericVector DeltaI    = Intakeval - Iref;
    NumericVector p         = cP(FFM, FM);
    NumericVector rhoFFM    = cRhoFFM(FFM);
    NumericVector growth    = Growth_dynamic(t);
    NumericVector Expend    = K + (22.4 + delta)*FFM + (4.5 + delta)*FM +
                                0.24*DeltaI + (230.0/rhoFFM *p + 180.0/rhoFM*(1.0-p))*Intakeval +
                                growth*(230.0/rhoFFM -180.0/rhoFM);
    return Expend/(1.0+230.0/rhoFFM *p + 180.0/rhoFM*(1.0-p));
}

//Rungue Kutta 4 method for Adult
List Child::rk4 (double days){
    
    //Initial time
    NumericMatrix k1, k2, k3, k4;
    
    //Estimate number of elements to loop into
    int nsims = floor(days/dt);
    
    //Create array of states
    NumericMatrix ModelFFM(nind, nsims + 1); //in rcpp
    NumericMatrix ModelFM(nind, nsims + 1); //in rcpp
    NumericMatrix ModelBW(nind, nsims + 1); //in rcpp
    NumericMatrix AGE(nind, nsims + 1); //in rcpp
    NumericVector TIME(nsims + 1); //in rcpp
    
    //Create initial states
    ModelFFM(_,0) = FFM;
    ModelFM(_,0)  = FM;
    ModelBW(_,0)  = FFM + FM;
    TIME(0)  = 0.0;
    AGE(_,0)  = age;
    
    //Loop through all other states
    bool correctVals = true;
    for (int i = 1; i <= nsims; i++){

        
        //Rungue kutta 4 (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
        k1 = dMass(AGE(_,i-1), ModelFFM(_,i-1), ModelFM(_,i-1));
        k2 = dMass(AGE(_,i-1) + 0.5 * dt/365.0, ModelFFM(_,i-1) + 0.5 * k1(0,_), ModelFM(_,i-1) + 0.5 * k1(1,_));
        k3 = dMass(AGE(_,i-1) + 0.5 * dt/365.0, ModelFFM(_,i-1) + 0.5 * k2(0,_), ModelFM(_,i-1) + 0.5 * k2(1,_));
        k4 = dMass(AGE(_,i-1) + dt/365.0, ModelFFM(_,i-1) + k3(0,_), ModelFM(_,i-1) +  k3(1,_));
        
        //Update of function values
        //Note: The dt is factored from the k1, k2, k3, k4 defined on the Wikipedia page and that is why
        //      it appears here.
        ModelFFM(_,i) = ModelFFM(_,i-1) + dt*(k1(0,_) + 2.0*k2(0,_) + 2.0*k3(0,_) + k4(0,_))/6.0;        //ffm
        ModelFM(_,i)  = ModelFM(_,i-1)  + dt*(k1(1,_) + 2.0*k2(1,_) + 2.0*k3(1,_) + k4(1,_))/6.0;        //fm
        
        //Update weight
        ModelBW(_,i) = ModelFFM(_,i) + ModelFM(_,i);
        
        //Update TIME(i-1)
        TIME(i) = TIME(i-1) + dt; // Currently time counts the time (days) passed since start of model
        
        //Update AGE variable
        AGE(_,i) = AGE(_,i-1) + dt/365.0; //Age is variable in years
    }
    
    return List::create(Named("Time") = TIME,
                        Named("Age") = AGE,
                        Named("Fat_Free_Mass") = ModelFFM,
                        Named("Fat_Mass") = ModelFM,
                        Named("Body_Weight") = ModelBW,
                        Named("Correct_Values")=correctVals,
                        Named("Model_Type")="Children");


}

NumericMatrix  Child::dMass (NumericVector t, NumericVector FFM, NumericVector FM){
    
    NumericMatrix Mass(2, nind); //in rcpp;
    NumericVector rhoFFM    = cRhoFFM(FFM);
    NumericVector p         = cP(FFM, FM);
    NumericVector growth    = Growth_dynamic(t);
    NumericVector expend    = Expenditure(t, FFM, FM);
    Mass(0,_)               = (1.0*p*(Intake(t) - expend) + growth)/rhoFFM;    // dFFM
    Mass(1,_)               = ((1.0 - p)*(Intake(t) - expend) - growth)/rhoFM; //dFM
    return Mass;
    
}

void Child::getParameters(void){
    
    //General constants
    rhoFM    = 9.4*1000.0;
    deltamin = 10.0;
    P        = 12.0;
    h        = 10.0;
    
    //Number of individuals
    nind     = age.size();
    
    //Sex specific constants
    ffm_beta0 = 2.9*(1 - sex)  + 3.8*sex;
    ffm_beta1 = 2.9*(1 - sex)  + 2.3*sex;
    fm_beta0  = 1.2*(1 - sex)  + 0.56*sex;
    fm_beta1  = 0.41*(1 - sex) + 0.74*sex;
    K         = 800*(1 - sex)  + 700*sex;
    deltamax  = 19*(1 - sex)   + 17*sex;
    A         = 3.2*(1 - sex)  + 2.3*sex;
    B         = 9.6*(1 - sex)  + 8.4*sex;
    D         = 10.1*(1 - sex) + 1.1*sex;
    tA        = 4.7*(1 - sex)  + 4.5*sex;       //years
    tB        = 12.5*(1 - sex) + 11.7*sex;      //years
    tD        = 15.0*(1-sex)   + 16.2*sex;      //years
    tauA      = 2.5*(1 - sex)  + 1.0*sex;       //years
    tauB      = 1.0*(1 - sex)  + 0.9*sex;       //years
    tauD      = 1.5*(1 - sex)  + 0.7*sex;       //years
    A_EB      = 7.2*(1 - sex)  + 16.5*sex;
    B_EB      = 30*(1 - sex)   + 47.0*sex;
    D_EB      = 21*(1 - sex)   + 41.0*sex;
    tA_EB     = 5.6*(1 - sex)  + 4.8*sex;
    tB_EB     = 9.8*(1 - sex)  + 9.1*sex;
    tD_EB     = 15.0*(1 - sex) + 13.5*sex;
    tauA_EB   = 15*(1 - sex)   + 7.0*sex;
    tauB_EB   = 1.5*(1 -sex)   + 1.0*sex;
    tauD_EB   = 2.0*(1 - sex)  + 1.5*sex;
    A1        = 3.2*(1 - sex)  + 2.3*sex;
    B1        = 9.6*(1 - sex)  + 8.4*sex;
    D1        = 10.0*(1 - sex) + 1.1*sex;
    tA1       = 4.7*(1 - sex)  + 4.5*sex;
    tB1       = 12.5*(1 - sex) + 11.7*sex;
    tD1       = 15.0*(1 - sex) + 16.0*sex;
    tauA1     = 1.0*(1 - sex)  + 1.0*sex;
    tauB1     = 0.94*(1 - sex) + 0.94*sex;
    tauD1     = 0.69*(1 - sex) + 0.69*sex;
}


//Intake in calories
NumericVector Child::Intake(NumericVector t){
    if (generalized_logistic) {
        return A_logistic + (K_logistic - A_logistic)/pow(C_logistic + Q_logistic*exp(-B_logistic*t), 1/nu_logistic); //t in years
    } else {
        int timeval = floor(365.0*(t(0) - age(0))/dt); //Example: Age: 6 and t: 7.1 => timeval = 401 which corresponds to the 401 entry of matrix
        return EIntake(timeval,_);
    }
    
}
