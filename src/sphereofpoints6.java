import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import processing.core.PApplet;
import processing.core.PVector;
import processing.opengl.PGraphicsOpenGL;

import javax.imageio.ImageIO;
import javax.media.opengl.GL2;
import javax.media.opengl.GLProfile;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;

public class sphereofpoints6 extends PApplet {














float[] projMatrix;
float[] mvMatrix;
String skyboxName = "besiege";   
float a = 0.0f;

public final int[] nodes_vbo = new int[1];
public int vnodelength;

int round =0;
GLSLshader shader;
PGraphicsOpenGL pgl;
GL2 gl;
Texture particleImg, ringImg;
Geodesic rt;
Point3d[] rts;
boolean growcomplete = false;
int radius = 30;
int nodecount = 30252;
int skybox;

boolean showlines = false;

public void setup() {
  
  size(1000, 800, OPENGL);
  hint(DISABLE_DEPTH_TEST);
  float fov = PI/3.0f;
  float cameraZ = (height/2.0f) / tan(fov/2.0f);
  perspective(fov, (PApplet.parseFloat(width)/PApplet.parseFloat(height)), cameraZ/10.0f, cameraZ*500);
  float distance = nodecount/11;
  camera(distance, distance, distance / tan(PI*30.0f / 180.0f), width/2.0f, height/2.0f, 0, 0, 1, 0);
  projMatrix = new float[16];
  mvMatrix = new float[16];
  loadSkybox(skyboxName, ".png");
  rt = new Geodesic(nodecount);
  rts = rt.getPointList();

  pgl = (PGraphicsOpenGL) g;
  gl = pgl.beginPGL().gl.getGL2();
  gl.setSwapInterval(1); 
  initShaders();
  initNodes();
  try {
    particleImg = AWTTextureIO.newTexture(new File(dataPath("particle.png")), true);
    ringImg = AWTTextureIO.newTexture(new File(dataPath("select.png")), true); 
  }
  catch (IOException e) {    
    println("Texture file (particle.png) is missing");
    exit();  // or handle it some other way
  }
   pgl.beginPGL();

  skybox = gl.glGenLists(1);
  gl.glNewList(skybox, GL2.GL_COMPILE);
  
  gl.glFrontFace(GL2.GL_CCW);
  gl.glEnable(GL2.GL_CULL_FACE);
  TexturedCube();
  gl.glDisable(GL2.GL_CULL_FACE);
  gl.glEndList();

  pgl.endPGL();
//  gl.glHint (GL.GL_POINT_SMOOTH_HINT, GL.GL_NICEST);
 
  // writestartlog("Stage 3 - Textures - Point smooth end");
  gl.glPointParameterf(GL2.GL_POINT_SIZE_MAX, 1000.0f);
  gl.glPointParameterf(GL2.GL_POINT_SIZE_MIN, 0.0f );
  gl.glPointParameterf( GL2.GL_POINT_FADE_THRESHOLD_SIZE, 60.0f );
  gl.glTexEnvi(GL2.GL_POINT_SPRITE, GL2.GL_COORD_REPLACE, GL2.GL_TRUE);
  gl.glEnable( GL2.GL_BLEND );
  // frameRate(2000);
   gl = pgl.beginPGL().gl.getGL2();
 // cam.setDistance(3000);
 // cam.setMinimumDistance(75);
 // cam.rotateX(.05);
}
//3,4,12
float rtx = 0;
public void draw() {

  background(0);
  if (frameCount % 30 == 0) {
    println(frameRate);
  }

  translate(width/2, height/2, 50);
  rtx += .003f;
  rotateY(rtx);
  
  loadMatrix();  //EXCLUDE if running newer than 2.0b7
  pgl.beginPGL();
  gl.glCallList( skybox );
  pgl.endPGL();
  

  gl = pgl.beginPGL().gl.getGL2(); 
  gl.glBlendFunc(GL2.GL_SRC_ALPHA, GL2.GL_ONE);
  pgl.beginPGL();
  gl.glDisable(GL2.GL_POINT_SMOOTH);
  gl.glEnable(GL2.GL_POINT_SPRITE);
  gl.glEnableClientState(GL2.GL_VERTEX_ARRAY);
  gl.glEnable(GL2.GL_VERTEX_PROGRAM_POINT_SIZE_ARB);

  gl.glColor3f(1, 1, 1);
  gl.glPointSize(35); 
  //  gl.glBegin(GL.GL_POINTS);
  shader.startShader();

  particleImg.bind(gl);
  particleImg.enable(gl);

  gl.glUniform1f(shader.getUniformLocation("gui"), 0.0f);
  gl.glUniform1f(shader.getUniformLocation("fogval"), 3000.0f);
  gl.glUniform1f(shader.getUniformLocation("pointSize"), 200.0f);
  
  gl.glBindBuffer(GL2.GL_ARRAY_BUFFER, nodes_vbo[0]);
  gl.glVertexPointer(3, GL2.GL_FLOAT, 0, 0);
  gl.glDrawArrays(GL2.GL_POINTS, 0, nodecount);
  gl.glBindBuffer( GL2.GL_ARRAY_BUFFER, 0);
  

  particleImg.disable(gl);
  
  ringImg.bind(gl);
  ringImg.enable(gl);
     gl.glUniform1f(shader.getUniformLocation("pointSize"), 200.0f);
 gl.glBindBuffer(GL2.GL_ARRAY_BUFFER, nodes_vbo[0]);
  gl.glVertexPointer(3, GL2.GL_FLOAT, 0, 0);
  gl.glDrawArrays(GL2.GL_POINTS, 0, nodecount);
  gl.glBindBuffer( GL2.GL_ARRAY_BUFFER, 0);
  
  ringImg.disable(gl);
  
  shader.endShader();

  gl.glDisable(GL2.GL_POINT_SPRITE);
  gl.glDisable(GL2.GL_VERTEX_PROGRAM_POINT_SIZE);
  pgl.endPGL();
 gl.glColor3f(1, 1, 1);

  gl.glBlendFunc(GL2.GL_SRC_ALPHA, GL2.GL_ONE_MINUS_SRC_ALPHA);
  
  //smooth();
 if (showlines) {
  strokeWeight(2);
  stroke(255, 255, 255, 20);
  for (int i=0;i< rts.length; i++) {
    Vector3d ps = new Vector3d(rts[i]);
    line(0.0f, 0.0f, 0.0f, (float)ps.x, (float)ps.y, (float)ps.z);
  }

 }
  if (!growcomplete) {
    Vector3d ps1 = new Vector3d(rts[0]);
    Vector3d ps2 = new Vector3d(rts[rt.getNextNode()]);
    float pdist = dist((float)ps1.x, (float)ps1.y, (float)ps1.z, (float)ps2.x, (float)ps2.y, (float)ps2.z);
    if (pdist < 45.0f) {
      radius = radius + (int) sqrt(rt.getNumberOfPoints())/2;
      //cam.setDistance(radius * 2 + 100, 10);
      rts = rt.getPointList();
      initNodes();
    } 
    else {
      initNodes();
      growcomplete = true;
    }
  }
   
}

public void keyPressed() {
  showlines = !showlines;
}




public void loadMatrix() {
  gl.glMatrixMode(GL2.GL_PROJECTION);
  projMatrix[0] = pgl.projection.m00;
  projMatrix[1] = pgl.projection.m10;
  projMatrix[2] = pgl.projection.m20;
  projMatrix[3] = pgl.projection.m30;
 
  projMatrix[4] = pgl.projection.m01;
  projMatrix[5] = pgl.projection.m11;
  projMatrix[6] = pgl.projection.m21;
  projMatrix[7] = pgl.projection.m31;
 
  projMatrix[8] = pgl.projection.m02;
  projMatrix[9] = pgl.projection.m12;
  projMatrix[10] = pgl.projection.m22;
  projMatrix[11] = pgl.projection.m32;
 
  projMatrix[12] = pgl.projection.m03;
  projMatrix[13] = pgl.projection.m13;
  projMatrix[14] = pgl.projection.m23;
  projMatrix[15] = pgl.projection.m33;
 
  gl.glLoadMatrixf(projMatrix, 0);
 
  gl.glMatrixMode(GL2.GL_MODELVIEW);
  mvMatrix[0] = pgl.modelview.m00;
  mvMatrix[1] = pgl.modelview.m10;
  mvMatrix[2] = pgl.modelview.m20;
  mvMatrix[3] = pgl.modelview.m30;
 
  mvMatrix[4] = pgl.modelview.m01;
  mvMatrix[5] = pgl.modelview.m11;
  mvMatrix[6] = pgl.modelview.m21;
  mvMatrix[7] = pgl.modelview.m31;
 
  mvMatrix[8] = pgl.modelview.m02;
  mvMatrix[9] = pgl.modelview.m12;
  mvMatrix[10] = pgl.modelview.m22;
  mvMatrix[11] = pgl.modelview.m32;
 
  mvMatrix[12] = pgl.modelview.m03;
  mvMatrix[13] = pgl.modelview.m13;
  mvMatrix[14] = pgl.modelview.m23;
  mvMatrix[15] = pgl.modelview.m33;
  gl.glLoadMatrixf(mvMatrix, 0);
}
public class Geodesic {

  private int N; //Tesselation frequency

  private int vertNum; // Number of vertices
  private int next;
  private boolean closecheck;
  /** Creates a new instance of Geodesic */

  public Geodesic(int tesselationFrequency)

  {
    //  N = tesselationFrequency;
    int vertcount = 0;
    next = 0;
    closecheck = false;
    while (vertNum < tesselationFrequency) {
      N = round(sqrt((tesselationFrequency-2)/10)) + vertcount;

      vertNum = 10*N*N+2;
      vertcount = vertcount + 1;
    }
    // println(N);
    //println(vertNum);
  }

  public int getNextNode() {
    return next;
  }

  public Point3d[] getPointList()

  {

    Point3d[] pointList = new Point3d[vertNum];



    for (int i = 0; i < vertNum; i++)

    {

      pointList[i] = createGeo(i);
      if (!closecheck) {
        if (next == 0) {
          next = i;
        } 
        else {
          Vector3d ps1 = new Vector3d(pointList[0]);
          Vector3d ps2 = new Vector3d(pointList[next]);
          Vector3d ps3 = new Vector3d(pointList[i]);
          float pdist1 = dist((float)ps1.x, (float)ps1.y, (float)ps1.z, (float)ps2.x, (float)ps2.y, (float)ps2.z);
          float pdist2 = dist((float)ps1.x, (float)ps1.y, (float)ps1.z, (float)ps3.x, (float)ps3.y, (float)ps3.z);
          if (pdist2 < pdist1) {
            next = i;
          }
        }
      }
    }
    closecheck = true;
    return pointList;
  }



  public int getNumberOfPoints()

  {

    return vertNum;
  }







  // x mod y

  private double mod(double x, double y)

  {

    return (y == 0) ? x : (x - (double)Math.floor(x/y)*y);
  }



  /* createGeo(m,N) outputs spherical coordinates of a geodesic sphere.
   
   * Where m is a single integer between 0 and 10N^2 + 2
   
   *        and values of m correspond to:
   
   *
   
   *                             0<= m < 12         Icosahedron vertex point
   
   *                           12 <= m < 30N - 18   Icosahedron edge point
   
   *                     30N - 18 <= m < 10N^2 + 2  Interior surface point
   
   */



  private Point3d createGeo(int m)

  {         

    SphericalCoordinates sc;



    if (m == 0)

    {

      sc = getPoint(0, 0, 1); //Top Vertex
    }

    else if (m == 1)

    {

      sc = getPoint(0, 0, 20); //Bottom Vertex;
    }

    else if (m < 12)

    {

      sc = getPoint(0, 0, m+4); //Top of one of the middle vertices
    }

    else if (m < 30*N -18) //Icosahedron edges

    {

      int edge  = (int)mod(m-12, 30); 

      int point = (int)Math.floor((m-12)/30);



      if (edge < 5)

      {

        sc = getPoint(0, point + 1, edge + 1);
      }

      else if (edge < 10)

      {

        sc = getPoint(point + 1, N-(point + 1), edge - 4);
      }

      else if (edge < 15)

      {

        sc = getPoint(0, point+1, edge - 4);
      }

      else if (edge < 20)

      {

        sc = getPoint(point + 1, 0, edge - 9);
      }

      else if (edge < 25)

      {

        sc = getPoint(point + 1, N-(point + 1), edge - 4);
      }

      else

      {

        sc = getPoint(0, point + 1, edge - 9);
      }
    }

    else    //Inner vertices

    {

      int face = (int)mod((m-(30*N)+18), 20);

      int point = (int)Math.floor((m-(30*N)+18)/20)+1;

      int offset = N - 2;

      int z;

      for (z = 1; z <=N; z++)

      {

        if (point <= offset)

        {

          break;
        }

        else

        {

          offset = offset + N - (z + 2);
        }
      }



      int y = offset - point + 1;

      int x = N - y - z;

      sc = getPoint(x, y, face + 1);
    }        

    return sc.getCartesian();
  }





  public SphericalCoordinates getPoint(int x, int y, int face)

  {

    SphericalCoordinates sc = topFace(x, y);



    if (face == 1)

    {
    }

    else if (face <= 5)

    {

      rotate(sc, face-1);
    }

    else if (face <= 10)

    {

      shiftDown(sc, face - 5);
    }

    else if (face <= 15)

    {

      shiftDown(sc, face - 10);

      flip(sc);
    }

    else

    {

      rotate(sc, face - 16);

      flip(sc);
    }         

    return sc;
  }          



  private SphericalCoordinates topFace(int x, int y)

  {

    int z = N - x - y;

    double x1 = x*Math.sin(2*Math.PI/5);

    double y1 = y + x*Math.cos(2*Math.PI/5);

    double z1 = .5f*N + (N - x - y)/((1 + Math.sqrt(5))/2);

    double phi = Math.atan2(x1, y1);

    double theta = Math.atan2(Math.sqrt(Math.pow(x1, 2) + Math.pow(y1, 2)), z1);

    SphericalCoordinates sc = new SphericalCoordinates(theta, phi);

    return sc;
  }



  private void rotate(SphericalCoordinates sc, int number)

  {

    sc.phi = sc.phi + number*2*Math.PI/5;
  }



  private void shiftDown(SphericalCoordinates sc, int number)

  {

    rotate(sc, number - 1);

    double phi0 = Math.PI/5 + 2*(number - 1)*Math.PI/5;

    Point3d p3d = new Point3d();

    double r1 = Math.sin(sc.theta)*Math.cos(sc.phi - phi0);

    double r3 = Math.cos(sc.theta);

    double sqrt5 = Math.sqrt(5);



    p3d.x = r1/sqrt5 + 2*r3/sqrt5;

    p3d.y = Math.sin(sc.theta)*Math.sin(sc.phi - phi0);

    p3d.z = r3/sqrt5 - 2*r1/sqrt5;



    sc.phi = phi0 + Math.PI/5 + Math.atan2(p3d.y, p3d.x);

    sc.theta = Math.atan2(Math.sqrt(p3d.x*p3d.x + p3d.y*p3d.y), p3d.z);
  } 



  private void flip(SphericalCoordinates sc)

  {

    sc.phi = sc.phi + Math.PI/5;

    sc.theta = Math.PI - sc.theta;
  }
}

public class SphericalCoordinates {
  double r;
  double theta;
  double phi;


  public SphericalCoordinates(int x, int y)
  {
    r = radius;
    theta = (double) x;
    phi = (double) y;
  }

  public SphericalCoordinates(double x, double y)
  {
    r = radius;
    theta = x;
    phi = y;
  }




  public Point3d getCartesian() {
    final double xyz[] = new double[3];
    xyz[0] = r * Math.sin(theta) * Math.cos(phi);
    xyz[1] = r * Math.sin(theta) * Math.sin(phi);
    xyz[2] = r * Math.cos(theta);
    Point3d ps = new Point3d(xyz);
    return ps ;
  }
}

public int fastDist( int dx, int dy )
{
  int min, max, approx;

  if ( dx < 0 ) dx = -dx;
  if ( dy < 0 ) dy = -dy;

  if ( dx < dy )
  {
    min = dx;
    max = dy;
  } 
  else {
    min = dy;
    max = dx;
  }

  approx = ( max * 1007 ) + ( min * 441 );
  if ( max < ( min << 4 ))
    approx -= ( max * 40 );

  // add 512 for proper rounding
  return (( approx + 512 ) >> 10 );
} 
public void initShaders(){
  PGraphicsOpenGL pgl = (PGraphicsOpenGL) g; 
  GL2 gl = pgl.beginPGL().gl.getGL2(); 
  shader = new GLSLshader(gl);
  shader.loadVertexShader("pglslvs.vert");
  shader.loadFragmentShader("pglslfs.frag");
  shader.useShaders();
  pgl.endPGL();
  }
class GLSLshader
{
  GL2 gl;
  int programObject;
  int vertexShader;
  int fragmentShader;

  GLSLshader(GL2 gl0)
  {
    gl = gl0;
    programObject = gl.glCreateProgramObjectARB();
    vertexShader = -1;
    fragmentShader = -1;
  }

  public void loadVertexShader(String file)
  {
    String shaderSource = join(loadStrings(file), "\n");
    vertexShader = gl.glCreateShaderObjectARB(GL2.GL_VERTEX_SHADER);
    gl.glShaderSourceARB(vertexShader, 1, new String[]{
      shaderSource                    }
    , (int[]) null, 0);
    gl.glCompileShaderARB(vertexShader);
    checkLogInfo(gl, vertexShader);
    gl.glAttachObjectARB(programObject, vertexShader);
  }

  public void loadFragmentShader(String file)
  {
    String shaderSource = join(loadStrings(file), "\n");
    fragmentShader = gl.glCreateShaderObjectARB(GL2.GL_FRAGMENT_SHADER);
    gl.glShaderSourceARB(fragmentShader, 1, new String[]{
      shaderSource                    }
    ,(int[]) null, 0);
    gl.glCompileShaderARB(fragmentShader);
    checkLogInfo(gl, fragmentShader);
    gl.glAttachObjectARB(programObject, fragmentShader);
  }

  public int getAttribLocation(String name)
  {
    return(gl.glGetAttribLocation(programObject, name));
  }

  public int getUniformLocation(String name)
  {
    return(gl.glGetUniformLocation(programObject, name));
  }

  public void useShaders()
  {
    gl.glLinkProgramARB(programObject);
    gl.glValidateProgramARB(programObject);
    checkLogInfo(gl, programObject);
  }

  public void startShader()
  {
    gl = pgl.beginPGL().gl.getGL2();
    gl.glUseProgramObjectARB(programObject);
  }

  public void endShader()
  {
    gl.glUseProgramObjectARB(0);
  }

  public void checkLogInfo(GL2 gl, int obj)
  {
    IntBuffer iVal = ByteBuffer.allocateDirect(1).order(ByteOrder.nativeOrder()).asIntBuffer();
    gl.glGetObjectParameterivARB(obj, GL2.GL_OBJECT_INFO_LOG_LENGTH_ARB, iVal);
    int length = 0;
    try {
     length = iVal.get();
    } catch (Exception e) {}
    if (length <= 1) return;
    ByteBuffer infoLog =  ByteBuffer.allocateDirect(length).order(ByteOrder.nativeOrder());
    iVal.flip();
    gl.glGetInfoLogARB(obj, length, iVal, infoLog);
    byte[] infoBytes = new byte[length];
    infoLog.get(infoBytes);
    println("GLSL Validation: \n" + new String(infoBytes));
  }
}





float p = 40000;   // half skybox size
float m = -p;
// create cube edges
PVector P000 = new PVector (m,m,m);
PVector P010 = new PVector (m,p,m);
PVector P110 = new PVector (p,p,m);
PVector P100 = new PVector (p,m,m);
PVector P001 = new PVector (m,m,p);
PVector P011 = new PVector (m,p,p);
PVector P111 = new PVector (p,p,p);
PVector P101 = new PVector (p,m,p);
Texture tex1,tex2,tex3,tex4,tex5,tex6;   // texture images
// load six skybox images as cube texture
public void loadSkybox(String skyboxName, String fExt) 
{ 
  try {
    tex1 = AWTTextureIO.newTexture(GLProfile.getDefault(), ImageIO.read(new File(dataPath(skyboxName + "_front" + fExt))), true);
    tex2 = AWTTextureIO.newTexture(GLProfile.getDefault(), ImageIO.read(new File(dataPath(skyboxName + "_back" + fExt))), true);
    tex3 = AWTTextureIO.newTexture(GLProfile.getDefault(), ImageIO.read(new File(dataPath(skyboxName + "_left" + fExt))), true);
    tex4 = AWTTextureIO.newTexture(GLProfile.getDefault(), ImageIO.read(new File(dataPath(skyboxName + "_right" + fExt))), true);
    tex5 = AWTTextureIO.newTexture(GLProfile.getDefault(), ImageIO.read(new File(dataPath(skyboxName + "_bottom" + fExt))), true);
    tex6 = AWTTextureIO.newTexture(GLProfile.getDefault(), ImageIO.read(new File(dataPath(skyboxName + "_top" + fExt))), true);
  }
  catch (IOException e) {    
    println( e);
  }

  //textureMode(NORMALIZED);
}
// Assign six texture to the six cube faces
public void TexturedCube() 
{
  TexturedCubeSide (P100, P000, P010, P110, tex1);   // -Z "front" face
  TexturedCubeSide (P001, P101, P111, P011, tex2);   // +Z "back" face
  TexturedCubeSide (P000, P001, P011, P010, tex3);   // -X "left" face
  TexturedCubeSide (P101, P100, P110, P111, tex4);   // +X "right" face
  TexturedCubeSide (P110, P010, P011, P111, tex5);   // +Y "base" face
  TexturedCubeSide (P101, P001, P000, P100, tex6);   // -Y "top" face
}
// create a cube side given by 4 edge vertices and a texture
public void TexturedCubeSide(PVector P1, PVector P2, PVector P3, PVector P4, Texture tex)
{

  tex.enable(gl);
  tex.bind(gl);
  gl.glBegin(GL2.GL_QUADS);
  gl.glTexCoord2f(1.0f, 0.0f);
  gl.glVertex3f(P1.x, P1.y, P1.z);
  gl.glTexCoord2f(0.0f, 0.0f);
  gl.glVertex3f(P2.x, P2.y, P2.z);
  gl.glTexCoord2f(0.0f, 1.0f);
  gl.glVertex3f(P3.x, P3.y, P3.z);
  gl.glTexCoord2f(1.0f, 1.0f);
  gl.glVertex3f(P4.x, P4.y, P4.z);
  gl.glEnd();
  tex.disable(gl);
}

public void initNodes() {

        
        if (nodecount > 0) {
           pgl.beginPGL();
           gl.glGenBuffers( 1, nodes_vbo, 0 );
           gl.glBindBuffer( GL2.GL_ARRAY_BUFFER, nodes_vbo[0] );
           gl.glBufferData( GL2.GL_ARRAY_BUFFER, rts.length * 3 * 4, point3d_to_float_buffer( rts ), GL2.GL_STATIC_DRAW );
           gl.glBindBuffer( GL2.GL_ARRAY_BUFFER, 0);
           pgl.endPGL();
         }
        // iconvbo = null;
    }
    
    
    // buffer converters
public FloatBuffer point3d_to_float_buffer( Point3d[] _vector )
{
  
 FloatBuffer a  = ByteBuffer.allocateDirect(_vector.length * 3 * 4).order(ByteOrder.nativeOrder()).asFloatBuffer();
 
 for( int i = 0; i < _vector.length; i++ )
 {
   Vector3d v = new Vector3d(_vector[i]);
   a.put( (float)v.x );
   a.put( (float)v.y );
   a.put( (float)v.z );
 }
 a.rewind();
 return a;
}


  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "sphereofpoints6" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
