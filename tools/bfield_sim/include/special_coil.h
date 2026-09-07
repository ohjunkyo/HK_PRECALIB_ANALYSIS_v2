// Calculate B-field from diagonal special coil
void calc_coilD () {
   
      // Define coil parameters
      Int_t n_phi = 360;
      Float_t d_phi = 2*TMath::Pi() / n_phi;
      Float_t theta_earth = b_earth.Theta()-TMath::Pi()/2;
      Float_t theta_coil = TMath::Pi()/2-theta_earth;
      // Float_t curr_coil = 1634.;
      Float_t curr_coil = 1000.;
   
      // Position vectors
      TVector3 v_start(r_tank, 0, -r_tank*tan(theta_coil));
      TVector3 v_prev = v_start;
   
      // Loop for phi segments
      for (Int_t iPhi=1; iPhi<=n_phi; iPhi++) {
   
         // End position of this segment
         Float_t x_this = r_tank * cos(iPhi*d_phi);
         Float_t y_this = r_tank * sin(iPhi*d_phi);
         Float_t z_this = -tan(theta_coil) * x_this;
         TVector3 v_this(x_this, y_this, z_this);
   
         // Middle position of this segment
         TVector3 v_seg = (v_this + v_prev) * 0.5;
               
         // Current vector
         TVector3 v_curr = (v_this - v_prev).Unit() * curr_coil;
   
         // Length of this segment
         Float_t d_seg = (v_this - v_prev).Mag();
   
         // Calculate B-field for this segment
         calc_segment(v_seg, v_curr, d_seg);
   
         // Update previous position
         v_prev = v_this;
   
      }
    
   }
   
   // Calculate B-field from vertical special coil
   void calc_coilV_new () {
   
      // Coil parameters
      Int_t N_COIL_V_NEW = 4;
      Float_t theta_14 = 30. * TMath::Pi() / 180.; // Angle between south and coil 1,4 start point (deg)
      Float_t theta_23 = 60. * TMath::Pi() / 180.; // Angle between south and coil 2,3 start point (deg)
      Float_t x_ver_new[N_COIL_V_NEW] = {+r_tank*cos(theta_14), +r_tank*cos(theta_23), -r_tank*cos(theta_23), -r_tank*cos(theta_14)};
      Float_t z_ver_new[N_COIL_V_NEW] = {21000, 17800, 21000, 17800};
      Float_t curr_ver_new[N_COIL_V_NEW] = {0., 300., 300., 0.};
      // Float_t curr_ver_new[N_COIL_V_NEW] = {150., 150., 150., 150.};
   
      // Loop for vertical coils
      for (Int_t iCoil=0; iCoil<N_COIL_V_NEW; iCoil++) {
   
         printf("\rVertical coil new %d", iCoil);
         fflush(stdout);
         
         // --- Top, bottom ---
   
         // Calculate Y position
         Float_t y_new  = sqrt(pow(r_tank, 2) - pow(x_ver_new[iCoil], 2));
   
         // Start/end position
         TVector3 v_start_top(x_ver_new[iCoil], +y_new, +z_ver_new[iCoil]);
         TVector3 v_end_top  (x_ver_new[iCoil], -y_new, +z_ver_new[iCoil]);
         TVector3 v_start_bot(x_ver_new[iCoil], -y_new, -z_ver_new[iCoil]);
         TVector3 v_end_bot  (x_ver_new[iCoil], +y_new, -z_ver_new[iCoil]);
         
         // Calculate B-field 
         calc_radial(curr_ver_new[iCoil], v_start_top, v_end_top);
         calc_radial(curr_ver_new[iCoil], v_start_bot, v_end_bot);
   
         // --- Side ---
         
         // Coil position
         TVector3 v_start1(x_ver_new[iCoil], +y_new, -z_ver_new[iCoil]);
         TVector3 v_end1  (x_ver_new[iCoil], +y_new, +z_ver_new[iCoil]);
         TVector3 v_start2(x_ver_new[iCoil], -y_new, +z_ver_new[iCoil]);
         TVector3 v_end2  (x_ver_new[iCoil], -y_new, -z_ver_new[iCoil]);
   
         // Calculate B-field
         calc_linear(curr_ver_new[iCoil], v_start1, v_end1);
         calc_linear(curr_ver_new[iCoil], v_start2, v_end2);
         
      } // End of loop for vertical coils
   
      printf("\n");
   
   }
   
   // Calculate B-field from vertical cosine coil
   void calc_coilV_cos () {
      
      // Coil parameters
      Int_t N_THETA_COS = 9;  // Number of theta segments. Total number of coils is N_THETA_COS * 4 + 2
      Float_t d_cos = 0.1;
      // Int_t N_THETA_COS = 4;  // Number of theta segments. Total number of coils is N_THETA_COS * 4 + 2
      // Float_t d_cos = 0.22;
      Float_t curr_cos = 400.;
      // Float_t curr_cos = 800.;
   
      // First coil at +X (south)
      Float_t angle_width = 2.5;
      // Float_t angle_width = 6.0;
      Float_t x_coil0 = r_tank * cos(angle_width * TMath::Pi() / 180.); 
      Float_t y_coil0 = r_tank * sin(angle_width * TMath::Pi() / 180.);
      Float_t z_coil0 = 21000;
      TVector3 v_coil0[4];
      v_coil0[0].SetXYZ(+x_coil0, +y_coil0, +z_coil0);
      v_coil0[1].SetXYZ(+x_coil0, -y_coil0, +z_coil0);
      v_coil0[2].SetXYZ(+x_coil0, -y_coil0, -z_coil0);
      v_coil0[3].SetXYZ(+x_coil0, +y_coil0, -z_coil0);
   
      // Loop for X sides
      for (Int_t iX=0; iX<2; iX++) {
   
         // Define coil position at +X or -X
         TVector3 v_coil0_this[4];
         for (Int_t iCorner=0; iCorner<4; iCorner++) {
            v_coil0_this[iCorner] = v_coil0[iCorner];
            if (iX == 1) {
               // Mirror X
               v_coil0_this[iCorner].SetX(-v_coil0_this[iCorner].X());
            }
         }
   
         // First coil at +X or -X
         calc_radial(curr_cos, v_coil0_this[0], v_coil0_this[1]);
         calc_linear(curr_cos, v_coil0_this[1], v_coil0_this[2]);
         calc_radial(curr_cos, v_coil0_this[2], v_coil0_this[3]);
         calc_linear(curr_cos, v_coil0_this[3], v_coil0_this[0]);
   
         // Loop for Y sides
         for (Int_t iY=0; iY<2; iY++) {
   
            // Rotate direction
            Float_t d_cos_this = +d_cos;
            if (iY == 1) {
               d_cos_this = -d_cos;
            }
   
            // Loop for theta
            for (Int_t iTheta=0; iTheta<N_THETA_COS; iTheta++) {
               
               // Define rotation angle
               Float_t cos_this = d_cos_this * (iTheta+1);
               Float_t theta_this = asin(cos_this); // It is not cos but sin in this coordinate
   
               // Define coil position
               TVector3 v_coil_this[4];
               for (Int_t iCorner=0; iCorner<4; iCorner++) {
                  v_coil_this[iCorner] = v_coil0_this[iCorner];
                  v_coil_this[iCorner].RotateZ(theta_this);
               }
   
               // Calculate B-field
               calc_radial(curr_cos, v_coil_this[0], v_coil_this[1]);
               calc_linear(curr_cos, v_coil_this[1], v_coil_this[2]);
               calc_radial(curr_cos, v_coil_this[2], v_coil_this[3]);
               calc_linear(curr_cos, v_coil_this[3], v_coil_this[0]);
   
               printf("\rVertical coil cos %d-%d-%d", iX, iY, iTheta);
               fflush(stdout);
   
            }
            
         }
         
      }
   
      // printf("\n");
   
   }
   
   // Calculate B-field from vertical cosine coil with current variation
   void calc_coilV_cos2 () {
      
      // Coil parameters
      Int_t N_COIL = 14;
      Float_t curr_cos = 600.;
   
      // First coil at +X (south)
      Float_t angle_width = TMath::Pi()*2 / N_COIL;
      // Float_t x_coil0 = (r_tank-500) * cos(angle_width/2.);
      // Float_t y_coil0 = (r_tank-500) * sin(angle_width/2.);
      Float_t x_coil0 = r_tank * cos(angle_width/2.);
      Float_t y_coil0 = r_tank * sin(angle_width/2.);
      Float_t z_coil0 = 21000;
      TVector3 v_coil0[4];
      v_coil0[0].SetXYZ(+x_coil0, +y_coil0, +z_coil0);
      v_coil0[1].SetXYZ(+x_coil0, -y_coil0, +z_coil0);
      v_coil0[2].SetXYZ(+x_coil0, -y_coil0, -z_coil0);
      v_coil0[3].SetXYZ(+x_coil0, +y_coil0, -z_coil0);
   
      // Loop for coils
      for (Int_t iCoil=0; iCoil<N_COIL; iCoil++) {
         
         printf("\rVertical coil cos2 %d", iCoil);
         fflush(stdout);
   
         // Angle of this coil
         Float_t angle_this = iCoil * angle_width;
   
         // Coil position
         TVector3 v_coil_this[4];
         for (Int_t iCorner=0; iCorner<4; iCorner++) {
            v_coil_this[iCorner] = v_coil0[iCorner];
            v_coil_this[iCorner].RotateZ(angle_this);
         }
   
         // Define current
         Float_t curr_this = curr_cos * cos(angle_this);
   
         // Calculate B-field
         calc_radial(curr_this, v_coil_this[0], v_coil_this[1]);
         calc_linear(curr_this, v_coil_this[1], v_coil_this[2]);
         calc_radial(curr_this, v_coil_this[2], v_coil_this[3]);
         calc_linear(curr_this, v_coil_this[3], v_coil_this[0]);
   
      }
   
      printf("\n");
   
   }
   
   // Calculate B-field from vertical fan-shape coil
   void calc_coilV_fan () {
      
      // Coil parameters
      Int_t N_COIL = 12; // Must be even number
      Float_t angle_fan = 60. * TMath::Pi()/180.; // 60 must be the best number
      Float_t curr_fan = 100.;
      Float_t angle_step = angle_fan*2 / (N_COIL-1);
   
      // Loop for coils
      for (Int_t iCoil=0; iCoil<N_COIL; iCoil++) {
         
         printf("\rVertical coil fan %d", iCoil);
         fflush(stdout);
   
         // Angle of this coil
         Float_t angle1, angle2;
         if (iCoil < N_COIL/2) {
            angle1 =                TMath::Pi()/2. - angle_fan + angle_step * iCoil;
            angle2 = angle1 + angle_step * N_COIL/2 + TMath::Pi();
         } else {
            angle1 = TMath::Pi() - (TMath::Pi()/2. - angle_fan + angle_step * (iCoil-N_COIL/2));
            angle2 = angle1 - angle_step * N_COIL/2 + TMath::Pi();
         }
   
         // Coil position
         Float_t x1 = r_tank * cos(angle1);
         Float_t y1 = r_tank * sin(angle1);
         Float_t x2 = r_tank * cos(angle2);
         Float_t y2 = r_tank * sin(angle2);
         Float_t z1 = +21000;
         TVector3 v_coil_this[4];
         v_coil_this[0].SetXYZ(x1, y1, -z1);
         v_coil_this[1].SetXYZ(x1, y1, +z1);
         v_coil_this[2].SetXYZ(x2, y2, +z1);
         v_coil_this[3].SetXYZ(x2, y2, -z1);
   
         // Calculate B-field
         calc_linear(curr_fan, v_coil_this[0], v_coil_this[1]);
         calc_radial(curr_fan, v_coil_this[1], v_coil_this[2]);
         calc_linear(curr_fan, v_coil_this[2], v_coil_this[3]);
         calc_radial(curr_fan, v_coil_this[3], v_coil_this[0]);
   
         // // Additional top coil
         // // Float_t z_add1 = +15200;
         // // Float_t z_add2 = +21000;
         // Float_t z_add1 = +16600;
         // Float_t z_add2 = +19600;
         // Float_t curr_add = 100.;
         // v_coil_this[0].SetXYZ(x1, y1, +z_add1);
         // v_coil_this[1].SetXYZ(x1, y1, +z_add2);
         // v_coil_this[2].SetXYZ(x2, y2, +z_add2);
         // v_coil_this[3].SetXYZ(x2, y2, +z_add1);
         // calc_linear(curr_add, v_coil_this[0], v_coil_this[1]);
         // calc_radial(curr_add, v_coil_this[1], v_coil_this[2]);
         // calc_linear(curr_add, v_coil_this[2], v_coil_this[3]);
         // calc_radial(curr_add, v_coil_this[3], v_coil_this[0]);
   
         // // Additional bottom coil
         // v_coil_this[0].SetXYZ(x1, y1, -z_add2);
         // v_coil_this[1].SetXYZ(x1, y1, -z_add1);
         // v_coil_this[2].SetXYZ(x2, y2, -z_add1);
         // v_coil_this[3].SetXYZ(x2, y2, -z_add2);
         // calc_linear(curr_add, v_coil_this[0], v_coil_this[1]);
         // calc_radial(curr_add, v_coil_this[1], v_coil_this[2]);
         // calc_linear(curr_add, v_coil_this[2], v_coil_this[3]);
         // calc_radial(curr_add, v_coil_this[3], v_coil_this[0]);
   
      }
   
      // // Cancel out top and bottom part above 60 deg.
      // // Add linear coil instead.
      // Float_t x1 = r_tank * cos((90-60.) * TMath::Pi()/180.);
      // Float_t y1 = r_tank * sin((90-60.) * TMath::Pi()/180.);
      // Float_t z1 = +21000;
      // TVector3 v_coil_this[2];
      // v_coil_this[0].SetXYZ(+x1, +y1, +z1);
      // v_coil_this[1].SetXYZ(+x1, -y1, +z1);
      // calc_radial(-curr_fan*N_COIL/2, v_coil_this[0], v_coil_this[1]);
      // calc_linear(+curr_fan*N_COIL/2, v_coil_this[0], v_coil_this[1]);
      // v_coil_this[0].SetXYZ(+x1, -y1, -z1);
      // v_coil_this[1].SetXYZ(+x1, +y1, -z1);
      // calc_radial(-curr_fan*N_COIL/2, v_coil_this[0], v_coil_this[1]);
      // calc_linear(+curr_fan*N_COIL/2, v_coil_this[0], v_coil_this[1]);
      // v_coil_this[0].SetXYZ(-x1, +y1, +z1);
      // v_coil_this[1].SetXYZ(-x1, -y1, +z1);
      // calc_radial(-curr_fan*N_COIL/2, v_coil_this[0], v_coil_this[1]);
      // calc_linear(+curr_fan*N_COIL/2, v_coil_this[0], v_coil_this[1]);
      // v_coil_this[0].SetXYZ(-x1, -y1, -z1);
      // v_coil_this[1].SetXYZ(-x1, +y1, -z1);
      // calc_radial(-curr_fan*N_COIL/2, v_coil_this[0], v_coil_this[1]);
      // calc_linear(+curr_fan*N_COIL/2, v_coil_this[0], v_coil_this[1]);
         
      printf("\n");
   
   }
   
// Calculate B-field from HH coils
void calc_coilHH () {

   for (Int_t iCoil=0; iCoil<N_COIL_HH; iCoil++) {

#ifdef HH_COIL5_EQUIV
      // Select coil#5 equivalent
      if (iCoil!=2 && iCoil!=3) continue;
#endif

      // Calc z position
      // Float_t z_this = z_hor[1] - (z_hor[1] - z_hor[12]) / 5 * iCoil;
      Float_t z_this = z_hh[iCoil];

      printf("Horizontal HH coil %d (HH%02d) z=%.1f\n", iCoil, iCoil+1, z_this);

      // Coil position
      // TVector3 v_start(19500, 0, z_this);
      TVector3 v_start(19000, 0, z_this);
      TVector3 v_end = v_start;

      // Current
      // Float_t curr_this = 31.35*4*2.0;
      Float_t curr_this = 31.35*4*1.8;
      // Float_t curr_this = 31.35*4*1.6;
      if (iCoil==0 || iCoil==5) {
         curr_this = 28.20*2*4*2.0;
      }

      // Calculate B-field
      calc_radial(curr_this, v_start, v_end);

   }
   
}

// Calculate B-field from VV coil 1
void calc_coilVV01 () {
   
   printf("Vertical VV coil 0 (VV01)\n");

   // VV1 port positions in SK coordinate
   Float_t r_vv1 = 19500 - 300;
   TVector3 xy_vv1[2];
   xy_vv1[0].SetXYZ(r_vv1*cos(18*TMath::DegToRad()),   r_vv1*sin(18*TMath::DegToRad()),   0);
   xy_vv1[1].SetXYZ(r_vv1*cos(18*4*TMath::DegToRad()), r_vv1*sin(18*4*TMath::DegToRad()), 0);

   // Convert to coil coordinate
   for (Int_t iPort=0; iPort<2; iPort++) {
      xy_vv1[iPort].RotateZ(angle_south);
   }

   // Calculate B-field for radial part
   Float_t curr_this = 31.10*4;
   TVector3 v_start1(xy_vv1[1].X(), xy_vv1[1].Y(), +20700);
   TVector3 v_end1  (xy_vv1[0].X(), xy_vv1[0].Y(), +20700);
   TVector3 v_start2(xy_vv1[0].X(), xy_vv1[0].Y(), -18700);
   TVector3 v_end2  (xy_vv1[1].X(), xy_vv1[1].Y(), -18700);
   calc_radial(curr_this, v_start1, v_end1);
   calc_radial(curr_this, v_start2, v_end2);

   // Calculate B-field for linear part
   TVector3 v_start3(xy_vv1[0].X(), xy_vv1[0].Y(), +20700);
   TVector3 v_end3  (xy_vv1[0].X(), xy_vv1[0].Y(), -18700);
   TVector3 v_start4(xy_vv1[1].X(), xy_vv1[1].Y(), -18700);
   TVector3 v_end4  (xy_vv1[1].X(), xy_vv1[1].Y(), +20700);
   calc_linear(curr_this, v_start3, v_end3);
   calc_linear(curr_this, v_start4, v_end4);

}