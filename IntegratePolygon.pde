
//██╗███╗   ██╗████████╗███████╗ ██████╗ ██████╗  █████╗ ████████╗███████╗
//██║████╗  ██║╚══██╔══╝██╔════╝██╔════╝ ██╔══██╗██╔══██╗╚══██╔══╝██╔════╝
//██║██╔██╗ ██║   ██║   █████╗  ██║  ███╗██████╔╝███████║   ██║   █████╗  
//██║██║╚██╗██║   ██║   ██╔══╝  ██║   ██║██╔══██╗██╔══██║   ██║   ██╔══╝  
//██║██║ ╚████║   ██║   ███████╗╚██████╔╝██║  ██║██║  ██║   ██║   ███████╗
//╚═╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝   ╚══════╝
                       
static int A;   /* alpha */
static int B;   /* beta */
static int C;   /* gamma */

/* projection integrals */
static float P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

/* face integrals */
static float Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

/* volume integrals */
static float T0, 
            T1[] = {0,0,0},
            T2[] = {0,0,0}, 
            TP[] = {0,0,0};
      
//██████████
PVector integratePoly(Polygon p){
  float r[] = {0,0,0};  /* center of mass */
  float J[][] = {{0,0,0},{0,0,0},{0,0,0}};                 /* inertia tensor */

  compVolumeIntegrals(p);

  print("\nT1 = ", T0, "(Volume) \n\n");

  print("Tx =   ", T1[0],"\n");
  print("Ty =   ", T1[1],"\n");
  print("Tz =   ", T1[2],"\n\n");
  
  print("Txx =  ", T2[0],"\n");
  print("Tyy =  ", T2[1],"\n");
  print("Tzz =  ", T2[2],"\n\n");

  print("Txy =  ", TP[0],"\n");
  print("Tyz =  ", TP[1],"\n");
  print("Tzx =  ", TP[2],"\n\n");

  //p.density /= T0;
  //p.mass = p.density * T0;
  p.density = p.mass/T0;
  println("mass = ", p.mass);
  println("density = ", p.density,'\n');
  
  /* compute center of mass */
  r[0] = T1[0] / T0;
  r[1] = T1[1] / T0;
  r[2] = T1[2] / T0;

  /* compute inertia tensor */
  J[0][0] = p.density * (T2[1] + T2[2]);
  J[1][1] = p.density * (T2[2] + T2[0]);
  J[2][2] = p.density * (T2[0] + T2[1]);
  J[0][1] = J[1][0] = - p.density * TP[0];
  J[1][2] = J[2][1] = - p.density * TP[1];
  J[2][0] = J[0][2] = - p.density * TP[2];

  /* translate inertia tensor to center of mass */
  J[0][0] -= p.mass * (r[1]*r[1] + r[2]*r[2]);
  J[1][1] -= p.mass * (r[2]*r[2] + r[0]*r[0]);
  J[2][2] -= p.mass * (r[0]*r[0] + r[1]*r[1]);
  J[0][1] = J[1][0] += p.mass * r[0] * r[1]; 
  J[1][2] = J[2][1] += p.mass * r[1] * r[2]; 
  J[2][0] = J[0][2] += p.mass * r[2] * r[0]; 

  print("center of mass: ( ", r[0]," ", r[1]," ", r[2], " ) \n\n");

  p.I = new PMatrix3D(J[0][0], J[0][1], J[0][2], 0,
                      J[1][0], J[1][1], J[Y][2], 0,
                      J[2][0], J[2][1], J[2][2], 0,
                            0,       0,       0, 1);
  
  p.moveShape(new PVector(-r[X], -r[Y], -r[Z]) );
  return new PVector(r[X], r[Y], r[Z]);
}

//██████████
PVector getNormalToFace(PShape f) {
  //get the normal to the face
  PVector normal = new PVector(0,0,0);
    
  //Newell's method for getting the normal
  int numberOfVertices = f.getVertexCount();
  for (int j = 0; j < numberOfVertices; ++j){
    normal.x += (f.getVertex(j).y - f.getVertex((j+1)%numberOfVertices).y) 
              * (f.getVertex(j).z + f.getVertex((j+1)%numberOfVertices).z);
    normal.y += (f.getVertex(j).z - f.getVertex((j+1)%numberOfVertices).z) 
              * (f.getVertex(j).x + f.getVertex((j+1)%numberOfVertices).x);
    normal.z += (f.getVertex(j).x - f.getVertex((j+1)%numberOfVertices).x) 
              * (f.getVertex(j).y + f.getVertex((j+1)%numberOfVertices).y);
  }  
  normal.normalize();
  return normal;
}
     
//██████████
void compVolumeIntegrals(Polygon p0) {
  PShape f;
  float nx, ny, nz;
  int i;
  
  PVector normal;
  
  T0 = 0; T1[0] = 0; T1[1] = 0; T1[2] = 0;
          T2[0] = 0; T2[1] = 0; T2[2] = 0;
          TP[0] = 0; TP[1] = 0; TP[2] = 0;
  
  //do each of the faces
  for (i = 0; i < p0.s.getChildCount() ; i++) {

    f = p.s.getChild(i);
    
    //get the normal to the plane of the face
    normal = getNormalToFace(f);
    normal.mult(1);
    
    nx = abs(normal.x);
    ny = abs(normal.y); 
    nz = abs(normal.z);
    
    float n[] = {normal.x, normal.y, normal.z};
    
    float w = - normal.dot( f.getVertex(0) ); 
    
  /*if (w<0) print('-');
    else print('+');
    if (i%80 == 0) {println();}
  */
    
    /*choose a-b-c as a right-handed permutation of x-y-z that maximizes |n.c|*/
    if      (nx > ny && nx > nz) C = 0; //if x is biggest use c = x = 0
    else if (ny > nz)            C = 1; //if y is biggest use c = y = 1
    else                         C = 2; //if z is biggest use c = z = 2
    A = (C + 1) % 3;
    B = (A + 1) % 3;

    compFaceIntegrals(f, n, w);

    T0 += normal.x * ((A == 0) ? Fa : ((B == 0) ? Fb : Fc));

    T1[A] += n[A] * Faa;
    T1[B] += n[B] * Fbb;
    T1[C] += n[C] * Fcc;
    T2[A] += n[A] * Faaa;
    T2[B] += n[B] * Fbbb;
    T2[C] += n[C] * Fccc;
    TP[A] += n[A] * Faab;
    TP[B] += n[B] * Fbbc;
    TP[C] += n[C] * Fcca;
  }

  T1[0] /= 2; T1[1] /= 2; T1[2] /= 2;
  T2[0] /= 3; T2[1] /= 3; T2[2] /= 3;
  TP[0] /= 2; TP[1] /= 2; TP[2] /= 2;
}

//██████████
void compFaceIntegrals(PShape f, float n[], float w) {
  float k1, k2, k3, k4;
  
  compProjectionIntegrals(f);
  
  k1 = 1 / n[C]; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;
  
  Fa = k1 * Pa;
  Fb = k1 * Pb;
  Fc = -k2 * (n[A]*Pa + n[B]*Pb + w*P1);

  Faa = k1 * Paa;
  Fbb = k1 * Pbb;
  Fcc = k3 * (SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb
   + w*(2*(n[A]*Pa + n[B]*Pb) + w*P1));

  Faaa = k1 * Paaa;
  Fbbb = k1 * Pbbb;
  Fccc = -k4 * (CUBE(n[A])*Paaa + 3*SQR(n[A])*n[B]*Paab 
     + 3*n[A]*SQR(n[B])*Pabb + CUBE(n[B])*Pbbb
     + 3*w*(SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb)
     + w*w*(3*(n[A]*Pa + n[B]*Pb) + w*P1));

  Faab = k1 * Paab;
  Fbbc = -k2 * (n[A]*Pabb + n[B]*Pbbb + w*Pbb);
  Fcca = k3 * (SQR(n[A])*Paaa + 2*n[A]*n[B]*Paab + SQR(n[B])*Pabb
   + w*(2*(n[A]*Paa + n[B]*Pab) + w*Pa));
} 
                                                                        
/*██████████ compute various integrations over projection of face */
void compProjectionIntegrals(PShape f) {
  float a0, a1, da;
  float b0, b1, db;
  float a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
  float a1_2, a1_3, b1_2, b1_3;
  float C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
  float Cab, Kab, Caab, Kaab, Cabb, Kabb;
  int i;

  P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;
  
  float vertexI[] = {0,0,0};
  float vertexII[]= {0,0,0};
  
  for (i = 0; i < f.getVertexCount(); i++) {
    vertexI[0] = f.getVertex(i).x;
    vertexI[1] = f.getVertex(i).y;
    vertexI[2] = f.getVertex(i).z;
    
    vertexII[0] = f.getVertex( (i+1) % f.getVertexCount() ).x;
    vertexII[1] = f.getVertex( (i+1) % f.getVertexCount() ).y;
    vertexII[2] = f.getVertex( (i+1) % f.getVertexCount() ).z;
    
    a0 = vertexI[A];
    b0 = vertexI[B];
    a1 = vertexII[A];
    b1 = vertexII[B];
    
    da = a1 - a0;
    db = b1 - b0;
    a0_2 = a0 * a0;   a0_3 = a0_2 * a0;   a0_4 = a0_3 * a0;
    b0_2 = b0 * b0;   b0_3 = b0_2 * b0;   b0_4 = b0_3 * b0;
    a1_2 = a1 * a1;   a1_3 = a1_2 * a1; 
    b1_2 = b1 * b1;   b1_3 = b1_2 * b1;

    C1 = a1 + a0;
    Ca = a1*C1 + a0_2;         Caa = a1*Ca + a0_3;   Caaa = a1*Caa + a0_4;
    Cb = b1*(b1 + b0) + b0_2;  Cbb = b1*Cb + b0_3;   Cbbb = b1*Cbb + b0_4;
    
    Cab = 3*a1_2 + 2*a1*a0 + a0_2; 
    Kab = a1_2 + 2*a1*a0 + 3*a0_2;
    
    Caab = a0*Cab + 4*a1_3;    
    Kaab = a1*Kab + 4*a0_3;
    Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
    Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

    P1 += db*C1;
    Pa += db*Ca;
    Paa += db*Caa;
    Paaa += db*Caaa;
    
    Pb += da*Cb;
    Pbb += da*Cbb;
    Pbbb += da*Cbbb;
    
    Pab += db*(b1*Cab + b0*Kab);
    Paab += db*(b1*Caab + b0*Kaab);
    Pabb += da*(a1*Cabb + a0*Kabb);
  }

  P1   /=   2.0;
  Pa   /=   6.0;
  Paa  /=  12.0;
  Paaa /=  20.0;
  Pb   /= - 6.0;
  Pbb  /= -12.0;
  Pbbb /= -20.0;
  Pab  /=  24.0;
  Paab /=  60.0;
  Pabb /= -60.0;
}
