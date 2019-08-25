float  SQR(float x)  {return x*x   ;}
float CUBE(float x) {return x*x*x ;}

//unit conversion
float m = 100;

void showMatrix(PMatrix3D M,String s, float x, float y){
  pushMatrix();
  translate( x, y);
  text(s + "\n" + M.m00 + "\n" + M.m10 + "\n" + M.m20,  20, 0 );
  text(    "\n" + M.m01 + "\n" + M.m11 + "\n" + M.m21, 220, 0 );  
  text(    "\n" + M.m02 + "\n" + M.m12 + "\n" + M.m22, 420, 0 );
  popMatrix();
}

void printVector(float dist, PVector prop, String s){
  translate(dist,0);
  text(s + "\n" + prop.x + "\n" + prop.y + "\n" + prop.z + "\n" + prop.mag(), 10, 30 );
  translate(-dist,0);
}

// find the min and max corner of a body
void checkBoundingBox(PVector min, PVector max, PVector suspect) {
  if (min.x > suspect.x) min.x = suspect.x;
  if (min.y > suspect.y) min.y = suspect.y;
  if (min.z > suspect.z) min.z = suspect.z;
  
  if (max.x < suspect.x) max.x = suspect.x;
  if (max.y < suspect.y) max.y = suspect.y;
  if (max.z < suspect.z) max.z = suspect.z;
}

// add one matrix to another
void matrixPlusEquals( PMatrix3D R, PMatrix3D RDot)  {
    R.m00 += RDot.m00;  R.m01 += RDot.m01;  R.m02 += RDot.m02; R.m03 += RDot.m03;
    R.m10 += RDot.m10;  R.m11 += RDot.m11;  R.m12 += RDot.m12; R.m13 += RDot.m13;
    R.m20 += RDot.m20;  R.m21 += RDot.m21;  R.m22 += RDot.m22; R.m23 += RDot.m23;
    R.m30 += RDot.m03;  R.m31 += RDot.m21;  R.m32 += RDot.m32; R.m23 += RDot.m33;
}

//create a copy of a matrix
PMatrix3D matrixClone(PMatrix3D M) {
  return new PMatrix3D(M.m00, M.m01, M.m02, M.m03,
                       M.m10, M.m11, M.m12, M.m13,
                       M.m20, M.m21, M.m22, M.m23,
                       M.m30, M.m31, M.m32, M.m33);
}

// vector a = a + b*c 
void vectorPlusEqualsScale(PVector vec, PVector vecdot, float dt) {
  vec.x += vecdot.x*dt; 
  vec.y += vecdot.y*dt; 
  vec.z += vecdot.z*dt;
}

//draw a line and spheres around a point showing a vector
void drawVector(PVector a, PVector b){  
  pushMatrix();
  translate( a.x, a.y, a.z );
  fill(255,0,0);
  sphere(3);
  translate( b.x - a.x, b.y - a.y, b.z - a.z);
  fill(0,255,0);
  sphere(3);
  popMatrix();
  
  beginShape();
  stroke(0,0,255);
  vertex(a.x,a.y,a.z);
  vertex(b.x,b.y,b.z);
  vertex(a.x,a.y+0.01,a.z);
  endShape();
}

//orthogonalise the top 3x3 submatrix with grahm-shmidt
void orthog(PMatrix3D M) {
  PVector v1 = new PVector(M.m00, M.m10, M.m20);
  PVector v2 = new PVector(M.m01, M.m11, M.m21);
  PVector v3 = new PVector(M.m02, M.m12, M.m22);
  float k1, k2;
  
  k1 = - v1.dot(v2) / v1.mag();
  v2.add( v1.copy().mult(k1) );
  
  k1 = - v1.dot(v3) / v1.mag();
  k2 = - v2.dot(v3) / v2.mag();
  v3.add( v1.copy().mult(k1) );
  v3.add( v2.copy().mult(k2) );
  
  v1.normalize();
  v2.normalize();
  v3.normalize();
  
  M.m00 = v1.x;  M.m01 = v2.x;  M.m02 = v3.x; 
  M.m10 = v1.y;  M.m11 = v2.y;  M.m12 = v3.y; 
  M.m20 = v1.z;  M.m21 = v2.z;  M.m22 = v3.z; 
}
