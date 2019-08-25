//██████╗  ██████╗ ██╗  ██╗   ██╗ ██████╗  ██████╗ ███╗   ██╗
//██╔══██╗██╔═══██╗██║  ╚██╗ ██╔╝██╔════╝ ██╔═══██╗████╗  ██║
//██████╔╝██║   ██║██║   ╚████╔╝ ██║  ███╗██║   ██║██╔██╗ ██║
//██╔═══╝ ██║   ██║██║    ╚██╔╝  ██║   ██║██║   ██║██║╚██╗██║
//██║     ╚██████╔╝███████╗██║   ╚██████╔╝╚██████╔╝██║ ╚████║
//╚═╝      ╚═════╝ ╚══════╝╚═╝    ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝

class Polygon {
  boolean exists;
  float len;
  float wid;
  float dep;
  PShape s,sh;
  PShape sRotated, sBounding, sBoundingRotated;
  PVector minCorner, maxCorner;
  
  float density;
  float mass;
  PVector position;
  PVector velocity;
  PVector acceleration;
  PVector force;
  
  PMatrix3D I, Iinv, ICuboid;
  PMatrix3D R, Rinv, RDot;
  PVector angularVelocity; 
  PMatrix3D angularVelocityCrossMatrix;
  PVector angularAcceleration = new PVector();
  PVector torque;
  
  PVector angularMomentum;
  
  float springConstant;
  
  ArrayList<PVector> vertices = new ArrayList<PVector>();
  ArrayList<PVector> verticesRotated = new ArrayList<PVector>();
              
  //██████████ Setup
  Polygon() { 
    position = new PVector(0,0,0);
    velocity = new PVector(0,0,0);
    acceleration = new PVector(0,0,0);
    angularVelocity = new PVector(0, 0, 0);  //construct base angular momentum vector
    torque = new PVector(0,0,0);
    force = new PVector(0,0,0);
    
    //construct base rotation matrix
    R = new PMatrix3D(1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1);
  
    mass = 1; 
    springConstant = 1;
  } 
  
  //██████████ remove any cuboids saved
  void clearCuboids() {
    vertices = new ArrayList<PVector>();
    verticesRotated = new ArrayList<PVector>();;
  }
  
  //██████████ Reset speed/position
  void resetPosition() {
    position        = new PVector(0,0,0);
    velocity        = new PVector(0,0,0);
    angularVelocity = new PVector(0,0,0);
    R               = new PMatrix3D(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
  }
  
  //██████████ Initialise Cuboid
  void makeCuboid() {
    PVector max = maxCorner;
    PVector min = minCorner;
    
    //make an arraylist of the vertices of the Polygon 
    vertices.add(new PVector(max.x, max.y, max.z));
    vertices.add(new PVector(min.x, max.y, max.z));
    vertices.add(new PVector(min.x, min.y, max.z));
    vertices.add(new PVector(max.x, min.y, max.z));
    
    vertices.add(new PVector(max.x, max.y, min.z));
    vertices.add(new PVector(min.x, max.y, min.z));
    vertices.add(new PVector(min.x, min.y, min.z));
    vertices.add(new PVector(max.x, min.y, min.z));
      
    len = max.x-min.x;
    wid = max.y-min.y;
    dep = max.z-min.z;
      
    //transcribe these to the rotated points arraylist
    for (int i = 0; i < vertices.size() ; ++i)
    {verticesRotated.add(new PVector(0,0,0) );}
     
    //construct base PShape
    sBounding = constructCuboidPShape(vertices);
    
    // get tensor of inertia and it's inverse
    I = new PMatrix3D( mass*(wid*wid + dep*dep)/12, 0, 0, 0,
                       0, mass*(len*len + dep*dep)/12, 0, 0,
                       0, 0, mass*(len*len + wid*wid)/12, 0,
                       0, 0, 0,                          1);
                      
    ICuboid = matrixClone(I);
    Iinv = matrixClone(I);
    Iinv.invert();
  }
  
  //██████████ Initialise Custom Obj Shape
  void makeObj(String dir){
    s = loadShape(dir);
    sRotated = loadShape(dir);
    exists = true;
  }
  
  //██████████ 
  void setupObj(){
    exists = false;
    
    scaleShape(-30); //scales and finds min max corners
  
    clearCuboids();
    makeCuboid();
    
    println("Cuboid Volume = ",  len * wid * dep, "");
    println("density roughly is ", mass/(len * wid * dep),"\n");
    
    //vertexNumber = s.getChildCount()/2;
    centerOfMass = integratePoly(this);
    
    print("Center OF Mass = ", centerOfMass.x, " ", centerOfMass.y, " ", centerOfMass.z, " \n");
    
    //set the Tensor of Inertia inverse used to integrate angular velocity
    Iinv = matrixClone(I);
    Iinv.invert();
    
    exists = true;
  }  
  
  //██████████ Scale the size of the Shape
  void scaleShape(float k) {
    PVector temp;
    int j;
    
    minCorner = new PVector( 255, 255, 255);
    maxCorner = new PVector(-255,-255,-255);
    
    //scale each of the points of the PShape
    for (int i = 0; i < s.getChildCount(); ++i ) {
      for (j = 0; j<3 ; ++j) {
        temp = s.getChild(i).getVertex(j);
        temp.mult(k);
        p.s.getChild(i).setVertex(j, temp.copy() );
      
        checkBoundingBox(minCorner, maxCorner, temp); 
      }
    }
    
    //also move the bounding box
    for (int i = 0; i < vertices.size(); ++i){
      vertices.get(i).mult(k); 
     }
    
    //also alter the inertia tensor
    if (exists) {
      I.scale(k*k);
      ICuboid.scale(k*k);
      Iinv.scale(1.0/(k*k));
    }
  }
  
  //██████████ Scale the size of the Shape
  void scaleMass(float k) {
  if (exists) {
      mass *= k;
      I.scale(k);
      ICuboid.scale(k);
      Iinv.scale(1.0/k);
    }
  }
  
  //██████████ Translate the shape to newPos
  void moveShape(PVector newPos) {
    PVector temp;
    int j;
  
    //move the PShape vertices
    for (int i = 0; i < s.getChildCount(); ++i ) {
      for (j = 0; j<3 ; ++j) {
        temp = s.getChild(i).getVertex(j).copy();
        temp.add(newPos);
        p.s.getChild(i).setVertex(j, temp);
      
        checkBoundingBox(minCorner, maxCorner, temp); 
      }
    }
    
    //also move the bounding box
    for (int i = 0; i < vertices.size(); ++i){
      vertices.get(i).add(newPos); 
     }
  }
  
  //██████████ Construct the shape from the vectors
  PShape constructCuboidPShape(ArrayList<PVector> points) {
    //Join the cuboid vertices together in a PShape
    PShape sh = createShape();
    sh.beginShape();
    
    sh.stroke(255);
    sh.noFill();
    
    int i;
    for (i = 0; i < 4; ++i) {
      int k;
      k =         i; sh.vertex( points.get(k).x,  points.get(k).y,  points.get(k).z );
      k =       4+i; sh.vertex( points.get(k).x,  points.get(k).y,  points.get(k).z );
      k = 4+(i+1)%4; sh.vertex( points.get(k).x,  points.get(k).y,  points.get(k).z );
      k = (i+1)%4  ; sh.vertex( points.get(k).x,  points.get(k).y,  points.get(k).z );
      k =         i; sh.vertex( points.get(k).x,  points.get(k).y,  points.get(k).z );
    }

    sh.endShape();
    return sh;
  }
    
  //██████████
  void rotateR(float dt) {
    //get angular velocity cross matrix
    angularVelocityCrossMatrix 
      = new PMatrix3D(  0                 , -angularVelocity.z, +angularVelocity.y, 0, 
                        +angularVelocity.z, 0                 , -angularVelocity.x, 0, 
                        -angularVelocity.y, +angularVelocity.x, 0                 , 0, 
                        0                 , 0                 , 0                 , 0);
    
    //midpoint integration 
    PMatrix3D midmatrix = matrixClone(angularVelocityCrossMatrix);    
    
    midmatrix.apply(R);    
    midmatrix.scale(dt/2);

    matrixPlusEquals(midmatrix, R);
    
    //Find RDot at the midpoint and scale it to be dt
    RDot = matrixClone(angularVelocityCrossMatrix);
    RDot.apply(R);
    RDot.scale(dt);
    
    //use the midpoint slope to find the final point
    matrixPlusEquals(R, RDot);
    
    //orthogonalise to clean out errors
    orthog(R);
  }
  
  //██████████ Use R Matrix to rotate bodyspace vectors to worldspace
  void rotatePoints() {
    //make some temporary vectors to rotate
    PVector vecOld = new PVector(0,0,0);
    PVector vecNew = new PVector(0,0,0);
    
    //rotate the points and make the rotated vertices as such
    int j;
    for (int i = 0; i < s.getChildCount(); ++i) {
      for (j = 0; j < s.getChild(i).getVertexCount(); ++j) {
        vecOld.set( s.getChild(i).getVertex(j) );        //vecOld = unrotated vertex
        R.mult( vecOld, vecNew );                        //multyiply R*vecOld -> vecNew
        vecNew.add(position);                            //add the CoM Position
        sRotated.getChild(i).setVertex(j, vecNew);       //rotated vertex = vecNew
      }
    }
    
    //also rotate the bounding box points
    for (int i = 0; i < vertices.size(); ++i){
        vecOld.set( vertices.get(i).copy() );         //vecOld = unrotated vertex
        R.mult( vecOld, vecNew );                     //multyiply R*vecOld -> vecNew
        vecNew.add(position);                         //add the CoM Position
        verticesRotated.set(i, vecNew.copy() );       //rotated vertex = vecNew
    } 
  }
  
  
  //██████████ Apply a force on a point of the polygon
  void applyForce(PVector r0, PVector f, float deltaT) {
    //find the dispacement vector between CoM and the point of force
    PVector r = r0.copy();
            r.sub(position);
            
    //convert force -> torque -> angular acceleration
    //and force -> acceleration
    force = f.copy();
    torque = r.cross(f);
    Iinv.mult(torque, angularAcceleration);
    acceleration = new PVector(f.x/mass, f.y/mass, f.z/mass); 
    
    //integrate linear velocity w/ linear acceleration
    vectorPlusEqualsScale(velocity, acceleration, deltaT);
    
    //change torque to body space    
    Rinv = matrixClone(R);
    Rinv.invert();
    Rinv.mult(torque, torque);
    
    //integrate velocity w/ angular acceleration
    angularAcceleration.mult( deltaT );
    angularVelocity.add( angularAcceleration );
  }
  
  //██████████ ""Fix"" any floating point errors
  void dampen(float k) {
    PVector delta;
    
    delta = velocity.copy();
    delta.normalize();
    delta.mult( k*k*velocity.magSq() );
    velocity.sub( delta );
    
    delta = angularVelocity.copy();
    delta.normalize();
    delta.mult( k*angularVelocity.magSq() );
    angularVelocity.sub(delta);
  }
  
  //██████████ 
  void recenter() {
    PVector newPosition = maxCorner.copy();
    newPosition.sub(minCorner);
    newPosition.div(2);
    newPosition.add(minCorner);
    newPosition.mult(-1);
    moveShape(newPosition);
    scaleShape(1);
  }
  
  PVector getSpringForce(){
    PVector spring;
    spring = mouse.copy().sub(forcePoint);
    
    float springLength = spring.mag();
    float restLength = 0.5*m;
    spring.normalize();
    
    if (springLength>restLength)
      spring.mult(springLength - restLength);
    else
      spring.mult(0);
    
    return spring;
  }
  
  //██████████ Carry out what happens each frame
  void doFrame(float DeltaT) { 
   if (exists == true) {
     
    //carry out forces & rotation
    PVector spring = getSpringForce();
    applyForce(position, g, mass*DeltaT);
    applyForce(forcePoint, spring, springConstant*DeltaT);
    
    dampen(0.0005);

    //integrate position
    vectorPlusEqualsScale(position, velocity, DeltaT);
    
    rotateR(DeltaT);
    rotatePoints();
   }
  } 
  
  //██████████ Show the PShape
  void display() {
      drawVector(forcePoint, mouse);
      shape( sRotated );
  }
}
