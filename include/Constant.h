#ifndef Constant_h
#define Constant_h

#define simulation_type 0     // 0: beta simulation(including major bkg simulation), 1: alpha simulation, 2: cosmic ray muon simulation

namespace TPCsystem{
	//Important: Constants changed
  //const int Nfec = 4;
  //const int Nchip = 4;
  //const int Nch =64;

  const int Nfec = 2;
  const int Nchip = 4;
  const int Nch = 68;


  const int Nsp =512; //sampling points
  const int Tchip = Nfec*Nchip;
  const int Tch = Tchip*Nch;
  const int Nconnector =8;
  const int Ncchn=128;

  //variables for simulation digitization
  const int nch = 120;          //total channel numbers for x and y each (for detector)
  const double chnwidth = 1.33;        //diagonal of one pad, unit is mm
  const double E_ion = 26.0;          //ionization energy of Ar, in eV
  const double timestep = 40.0;       //the timestep of one ADC sampling point
  
  //for Ar+3.5% iC4H10 @ 120V/cm
  const double v_drift = 3.72;         //drift velocity in cm/us

  //diffusion coefficients, the relation between diffused distance sigma and drift distance z is :
  // sigma_(t or l) = (dt or dl)* sqrt(z)
  const double dl = 0.037;            //longitudinal diffusion coefficient in Ar+3.5% iC4H10 at ~120V/cm (unit is cm^0.5)
  const double dt = 0.055;            //transverse diffusion coefficient in Ar+3.5% iC4H10 at ~120V/cm (unit is cm^0.5)

  // //for Ar+2.5% iC4H10 @ 100V/cm
  // const double v_drift = 3.46;         //drift velocity in cm/us
  // const double dl = 0.0435;
  // const double dt = 0.0656;

  const bool on_plane_spread = true;    // whether to enable spread sim on the resist layer
  const double sigma_spread = 0.5;      //sigma of the signal spreading due to the resist layer

  //response function, bin width is 40ns
  const double response_func[] = {4.118, 7.618, -2.782, -1.982, 8.318, 7.818, 7.018,
      25.418, 23.418, -5.782, -3.882, 17.918, 12.718, 7.218, 11.718, 17.918, 22.418,
      18.918, 7.718, 17.118, 24.918, 18.918, 17.318, 32.518, 31.818, 10.918, 5.218,
      14.818, 10.818, 7.918, 15.118, 29.818, 128.618, 425.318, 941.018, 1568.718,
      2208.018, 2814.318, 3355.718, 3809.818, 4170.118, 4435.618, 4612.418, 4706.418,
      4718.318, 4672.318, 4582.418, 4445.818, 4267.718, 4070.018, 3840.418, 3587.518,
      3355.418, 3127.018, 2882.218, 2634.518, 2404.418, 2178.018, 1971.918, 1792.618,
      1604.918, 1423.118, 1271.618, 1132.818, 983.518, 850.718, 747.718, 649.018, 560.618,
      495.818, 431.318, 367.018, 315.018, 267.218, 224.018, 187.918, 156.018, 120.718,
      94.618, 86.818, 74.718, 57.418, 49.018, 39.418, 22.518, 16.818, 11.118, -5.482, -11.482,
      -8.082, -9.082, -17.082, -4.382, 1.518, -8.182, -18.782, -29.882, -24.982, -12.482,
      -11.982, -13.582, -10.382, -15.882, -14.882, 2.718, 9.418, -4.382, -11.482, -12.982,
      -23.582, -26.782, -23.382, -23.982, -21.882, -13.082, -18.582, -31.082, -33.782, -19.682,
      -11.782, -18.782, -32.482, -33.282, -25.882, -27.382, -27.482, -29.282, -29.082, -30.782,
      -34.382, -32.482, -28.382, -27.782, -21.982, -18.282, -27.382, -30.782, -21.082, -19.782,
      -28.582, -40.182, -37.282, -22.882, -30.182, -43.882, -42.182, -22.982, -28.882, -50.082, -42.882};
  //maximum value position of the response function
  const int peak_pos = 45;
  //gain of the TPC system
  const double gain = 2300;
  //factor for the response function (for 1pC CSA gain)
  const double factor = 1000;
  //time bin width for the response function
  const double timebinwidth = 2e-9;
  //CSA gain: 10pC/1pC/240fC/120fC
  const double CSAgain = 0.12;

}
#endif
