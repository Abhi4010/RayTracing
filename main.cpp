#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include<math.h>
#include<windows.h>
#include<GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define pi 2*acos(0)



class Vector {
public:
	double x, y, z;


	double getVectX() { return x; }
	double getVectY() { return y; }
	double getVectZ() { return z; }

	double magnitude () {
		return sqrt((x*x) + (y*y) + (z*z));
	}

	Vector normalize() {
		double magnitude = sqrt((x*x) + (y*y) + (z*z));
		return Vector (x/magnitude, y/magnitude, z/magnitude);
	}

	Vector negative () {
		return Vector (-x, -y, -z);
	}

	double dotProduct(Vector v) {
		return x*v.getVectX() + y*v.getVectY() + z*v.getVectZ();
	}

	Vector crossProduct(Vector v) {
		return Vector (y*v.getVectZ() - z*v.getVectY(), z*v.getVectX() - x*v.getVectZ(), x*v.getVectY() - y*v.getVectX());
	}

	Vector vectAdd (Vector v) {
		return Vector (x + v.getVectX(), y + v.getVectY(), z + v.getVectZ());
	}
	Vector operator + (Vector v)
	{
		return Vector (x + v.getVectX(), y + v.getVectY(), z + v.getVectZ());

	}

	Vector vectMult (double scalar) {
		return Vector (x*scalar, y*scalar, z*scalar);
	}
	Vector operator * (double scalar)
	{
		return Vector (x*scalar, y*scalar, z*scalar);

	}
	Vector () {
	x = 0;
	y = 0;
	z = 0;
}

Vector (double i, double j, double k) {
	x = i;
	y = j;
	z = k;
}
};



struct Camera{
	Vector campos, camdir, camright, camdown;
	Camera ()
	{
		campos = Vector(0,0,0);
		camdir = Vector(0,0,1);
		camright = Vector(0,0,0);
		camdown = Vector(0,0,0);
	}
	Camera (Vector campos, Vector camdir, Vector camright, Vector camdown): campos(campos),camdir(camdir),
	camright(camright), camdown(camdown)
	{

	}

};
struct Color
{
	/* data */
	double red, green, blue, special;
	Color()
	{
		red = 0.5;
	green = 0.5;
	blue = 0.5;

	}
	Color(double red, double green, double blue, double special) :
	red(red), green(green), blue(blue),special(special) {}
	Color operator + (Color color) {
		return Color (red + color.red, green + color.green, blue + color.blue, special);
	}

	Color operator * (Color color) {
		return Color (red*color.red, green*color.green, blue*color.blue, special);
	}

	Color operator %(Color color) {
		return Color ((red + color.red)/2, (green + color.green)/2, (blue + color.blue)/2, special);
	}

};


struct  Ray
{

	/* data */
	Vector origin, direction;
	Ray()
	{
		origin = Vector(0,0,0);
		direction = Vector(1,0,0);
	}
  Ray(Vector origin, Vector direction) : origin(origin), direction(direction) {}
};

Color clip(Color myColor) {
		double alllight = myColor.red + myColor.green + myColor.blue;
		double excesslight = alllight - 3;
		if (excesslight > 0) {
			myColor.red = myColor.red + excesslight*(myColor.red/alllight);
			myColor.green = myColor.green + excesslight*(myColor.green/alllight);
			myColor.blue = myColor.blue + excesslight*(myColor.blue/alllight);
		}
		if (myColor.red > 1) {myColor.red = 1;}
		if (myColor.green > 1) {myColor.green = 1;}
		if (myColor.blue > 1) {myColor.blue = 1;}
		if (myColor.red < 0) {myColor.red = 0;}
		if (myColor.green < 0) {myColor.green = 0;}
		if (myColor.blue < 0) {myColor.blue = 0;}

		return Color (myColor.red, myColor.green, myColor.blue, myColor.special);
	}

	Color colorScalar (Color myColor, double scalar) {
		return Color (myColor.red*scalar, myColor.green *scalar, myColor.blue *scalar, myColor.special);
	}



class Src {
	public:

	virtual Vector getLightPosition() {return Vector(0, 0, 0);}
	virtual Color getLightColor() {return Color(1, 1, 1, 0);}
	Src() {}

};



class Light : public Src {
	Vector position;
	Color color;

	public:

	Light () {
	position = Vector(0,0,0);
	color = Color(1,1,1, 0);
}

Light (Vector p, Color c) {
	position = p;
	color = c;
}
	virtual Vector getLightPosition () { return position; }
	virtual Color getLightColor () { return color; }

};
class Object {
	public:

	Object ()
	{

	}

	virtual Color getColor () { return Color (0.0, 0.0, 0.0, 0); }

	virtual Vector getNormalAt(Vector pos_intersect) {
		return Vector (0, 0, 0);
	}

	virtual double findIntersection(Ray ray) {
		return 0;
	}

};

class ObjSphere : public Object {
	Vector center;
	double radius;
	Color color;

	public:

	ObjSphere ();

	ObjSphere (Vector, double, Color);

	// method functions
	Vector getSphereCenter () { return center; }
	double getSphereRadius () { return radius; }
	virtual Color getColor () { return color; }

	virtual Vector getNormalAt(Vector point) {
		// normal always points away from the center of a sphere
		Vector normal_Vect = point.vectAdd(center.negative()).normalize();
		return normal_Vect;
	}

	virtual double findIntersection(Ray ray) {
		Vector ray_origin = ray.origin;
		double ray_origin_x = ray_origin.getVectX();
		double ray_origin_y = ray_origin.getVectY();
		double ray_origin_z = ray_origin.getVectZ();

		Vector ray_direction = ray.direction;
		double ray_direction_x = ray_direction.getVectX();
		double ray_direction_y = ray_direction.getVectY();
		double ray_direction_z = ray_direction.getVectZ();

		Vector s_c = center;
		double s_c_x = s_c.getVectX() ;
		double s_c_y = s_c.getVectY();
		double s_c_z = s_c.getVectZ();

		double a = 1; // normalized
		double b = (2*(ray_origin_x - s_c_x)*ray_direction_x) + (2*(ray_origin_y - s_c_y)*ray_direction_y) + (2*(ray_origin_z - s_c_z)*ray_direction_z);
		double c = pow(ray_origin_x - s_c_x, 2) + pow(ray_origin_y - s_c_y, 2) + pow(ray_origin_z - s_c_z, 2) - (radius*radius);

		double discriminant = b*b - 4*c;

		if (discriminant > 0) {

			double root_1 = ((-1*b - sqrt(discriminant))/2) - 0.000001;

			if (root_1 > 0) {
				return root_1;
			}
			else {
				double root_2 = ((sqrt(discriminant) - b)/2) - 0.000001;
				return root_2;
			}
		}
		else {
			return -1;
		}
	}

};

ObjSphere::ObjSphere () {
	center = Vector(0,0,0);
	radius = 1.0;
	color = Color(0.5,0.5,0.5, 0);
}

ObjSphere::ObjSphere (Vector centerValue, double radiusValue, Color colorValue) {
	center = centerValue;
	radius = radiusValue;
	color = colorValue;
}
class ObjPlane : public Object {
	Vector normal;
	double distance;
	Color color;

	public:

	ObjPlane ();

	ObjPlane (Vector, double, Color);

	// method functions
	Vector getPlaneNormal () { return normal; }
	double getPlaneDistance () { return distance; }
	virtual Color getColor () { return color; }

	virtual Vector getNormalAt(Vector point) {
		return normal;
	}

	virtual double findIntersection(Ray ray) {
		Vector ray_direction = ray.direction;

		double a = ray_direction.dotProduct(normal);

		if (a == 0) {
			// ray is parallel to the plane
			return -1;
		}
		else {
			double b = normal.dotProduct(ray.origin + (normal*distance).negative());
			return -1*b/a;
		}
	}

};

ObjPlane::ObjPlane () {
	normal = Vector(1,0,0);
	distance = 0;
	color = Color(0.5,0.5,0.5, 0);
}

ObjPlane::ObjPlane (Vector normalValue, double distanceValue, Color colorValue) {
	normal = normalValue;
	distance = distanceValue;
	color = colorValue;
}



bool openGLDisplay;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;




using namespace std;



struct RGB {
	double r;
	double g;
	double b;
};

	int depth_per_inch = 72;
	int windowWidth = 700;

	int windowHeight = 480;
	int n = windowWidth*windowHeight;
	RGB *pixels = new RGB[n];
	unsigned int imgdta[640][480][3];
	Color colors1,colors2,colors3;
	Color colorcb1,colorcb2;
	double sp1ReflectionIndex, sp2ReflectionIndex, sp3ReflectionIndex;
	Color ColorIntersectObj;


void generateBMP (const char *filename, int w, int h, int depth_per_inch, RGB *data) {
	FILE *f;
	int k = w*h;
	int s = 4*k;
	int filesize = 54 + s;

	double factor = 39.375;
	int m = static_cast<int>(factor);

	int ppm = depth_per_inch*m;

	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0,0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};

	bmpfileheader[ 2] = (unsigned char)(filesize);
	bmpfileheader[ 3] = (unsigned char)(filesize>>8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);

	bmpinfoheader[ 4] = (unsigned char)(w);
	bmpinfoheader[ 5] = (unsigned char)(w>>8);
	bmpinfoheader[ 6] = (unsigned char)(w>>16);
	bmpinfoheader[ 7] = (unsigned char)(w>>24);

	bmpinfoheader[ 8] = (unsigned char)(h);
	bmpinfoheader[ 9] = (unsigned char)(h>>8);
	bmpinfoheader[10] = (unsigned char)(h>>16);
	bmpinfoheader[11] = (unsigned char)(h>>24);

	bmpinfoheader[21] = (unsigned char)(s);
	bmpinfoheader[22] = (unsigned char)(s>>8);
	bmpinfoheader[23] = (unsigned char)(s>>16);
	bmpinfoheader[24] = (unsigned char)(s>>24);

	bmpinfoheader[25] = (unsigned char)(ppm);
	bmpinfoheader[26] = (unsigned char)(ppm>>8);
	bmpinfoheader[27] = (unsigned char)(ppm>>16);
	bmpinfoheader[28] = (unsigned char)(ppm>>24);

	bmpinfoheader[29] = (unsigned char)(ppm);
	bmpinfoheader[30] = (unsigned char)(ppm>>8);
	bmpinfoheader[31] = (unsigned char)(ppm>>16);
	bmpinfoheader[32] = (unsigned char)(ppm>>24);

	f = fopen(filename,"wb");

	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);



	for (int i = 0; i < k; i++) {
		RGB rgb = data[i];

		double red = (data[i].r)*255;
		double green = (data[i].g)*255;
		double blue = (data[i].b)*255;

		unsigned char color[3] = {(int)floor(blue),(int)floor(green),(int)floor(red)};

		fwrite(color,1,3,f);
	}

	fclose(f);

}


int getInteractedObject(int index_of_minimum_value,vector<double> Obj_intersections )
{
    // prevent unnessary calculations
	if (Obj_intersections.size() == 0) {
		// if there are no intersections
		return -1;
	}
	else if (Obj_intersections.size() == 1) {
		if (Obj_intersections.at(0) > 0) {
			// if that intersection is greater than zero then its our index of minimum value
			return 0;
		}
		else {
			// otherwise the only intersection value is negative
			return -1;
		}
	}
	else {
		// otherwise there is more than one intersection
		// first find the maximum value

		double max = 0;
		for (int i = 0; i < Obj_intersections.size(); i++) {
			if (max < Obj_intersections.at(i)) {
				max = Obj_intersections.at(i);
			}
		}

		// then starting from the maximum value find the minimum positive value
		if (max > 0) {
			// we only want positive intersections
			for (int index = 0; index < Obj_intersections.size(); index++) {
				if (Obj_intersections.at(index) > 0 && Obj_intersections.at(index) <= max) {
					max = Obj_intersections.at(index);
					index_of_minimum_value = index;
				}
			}

			return index_of_minimum_value;
		}
		else {
			// all the intersections were negative
			return -1;
		}
	}
}

int getFirstIntersectedObj(vector<double> Obj_intersections) {
	// return the index of the winning intersection
	int index_of_minimum_value;

	int p = getInteractedObject(index_of_minimum_value,Obj_intersections);
	return p;


}



void colorSetOfPlane( Vector pos_intersect, bool ifCheckerboard)
{

    if(ifCheckerboard)
    {
       int square = (int)floor(pos_intersect.getVectX()) +
		 (int)floor(pos_intersect.getVectZ());


		if ((square % 2) == 0) {
			ColorIntersectObj.blue = colorcb1.red;
			ColorIntersectObj.green = colorcb1.green;
			ColorIntersectObj.blue = colorcb1.blue;
		}
		else {
			ColorIntersectObj.red = colorcb2.red;
			ColorIntersectObj.green = colorcb2.green;
			ColorIntersectObj.blue = colorcb2.blue;
		}




}
}


void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			drawgrid=1-drawgrid;
			break;

		default:
			break;
	}
}
int current;

Color renderColorValue(Vector pos_intersect, Vector rayDirIntersecting, vector<Object*>
                 s_b, int indexFirstIntersected, vector<Src*> light_sources,
                 double accuracy, double ambientlight) {

	Color ColorIntersectObj = s_b.at(indexFirstIntersected)->getColor();
	Vector normalIntersectingObj = s_b.at(indexFirstIntersected)->getNormalAt
	(pos_intersect);
	bool val =true;

	if (ColorIntersectObj.special == 2) {

		colorSetOfPlane(pos_intersect,val);
		int square = (int)floor(pos_intersect.getVectX()) +
		 (int)floor(pos_intersect.getVectZ());


		if ((square % 2) == 0) {
			ColorIntersectObj.blue = colorcb1.red;
			ColorIntersectObj.green = colorcb1.green;
			ColorIntersectObj.blue = colorcb1.blue;
		}
		else {
			// white tile
			ColorIntersectObj.red = colorcb2.red;
			ColorIntersectObj.green = colorcb2.green;
			ColorIntersectObj.blue = colorcb2.blue;
		}



	}

	Color colorFinal = colorScalar(ColorIntersectObj,ambientlight);

	if (ColorIntersectObj.special > 0 && ColorIntersectObj.special <= 1) {
		double dot1 = normalIntersectingObj.dotProduct(rayDirIntersecting.negative());
		Vector scalar1 = normalIntersectingObj.vectMult(dot1);
		Vector add1 = scalar1 + (rayDirIntersecting);
		Vector scalar2 = add1.vectMult(2);
		Vector add2 = rayDirIntersecting.negative()+ (scalar2);
		Vector reflection_direction = add2.normalize();

		Ray reflection_ray (pos_intersect, reflection_direction);

		vector<double> reflection_intersections;

		for (int reflection_index = 0; reflection_index < s_b.size(); reflection_index++) {
			reflection_intersections.push_back(s_b.at(reflection_index)->findIntersection(reflection_ray));
		}

		int indexFirstIntersected_with_reflection = getFirstIntersectedObj(reflection_intersections);

		if (indexFirstIntersected_with_reflection != -1) {
			if (reflection_intersections.at(indexFirstIntersected_with_reflection) > accuracy) {

				Vector reflection_intersection_position = pos_intersect.vectAdd(reflection_direction.vectMult(reflection_intersections.at(indexFirstIntersected_with_reflection)));
				Vector reflection_intersection_ray_direction = reflection_direction;

				Color reflection_intersection_color = renderColorValue(reflection_intersection_position, reflection_intersection_ray_direction, s_b, indexFirstIntersected_with_reflection, light_sources, accuracy, ambientlight);

				colorFinal = colorFinal + colorScalar(reflection_intersection_color, ColorIntersectObj.special);
			}
		}
	}

	for (int light_index = 0; light_index < light_sources.size(); light_index++) {
		Vector light_direction = light_sources.at(light_index)->getLightPosition().
		vectAdd(pos_intersect.negative()).normalize();

		float cosine_angle = normalIntersectingObj.dotProduct(light_direction);


		if (cosine_angle > 0) {
			bool shadowed = false;
			Vector distance_to_light = light_sources.at(light_index)->getLightPosition() + (pos_intersect.negative()).normalize();
			float distance_to_light_magnitude = distance_to_light.magnitude();

			Ray shadow_ray (pos_intersect, light_sources.at(light_index)->
                   getLightPosition().vectAdd(pos_intersect.negative()).normalize());

			vector<double> secondary_intersections;

			for (int Obj_index = 0; Obj_index < s_b.size() && shadowed == false;
			 Obj_index++)
            {
				secondary_intersections.push_back(s_b.at(Obj_index)
                                      ->findIntersection(shadow_ray));
			}

			for (int c = 0; c < secondary_intersections.size(); c++) {
				if (secondary_intersections.at(c) > accuracy) {
					if (secondary_intersections.at(c) <= distance_to_light_magnitude) {
						shadowed = true;
					}
					break;
				}

			}

			if (shadowed == false) {
				colorFinal = colorFinal + (ColorIntersectObj * colorScalar(light_sources.at(light_index)->getLightColor(),cosine_angle));

				if (ColorIntersectObj.special > 0 && ColorIntersectObj.special <= 1) {
					double dot1 = normalIntersectingObj.dotProduct(rayDirIntersecting.negative());
					Vector scalar1 = normalIntersectingObj * dot1;
					Vector add1 = scalar1 + rayDirIntersecting;
					Vector scalar2 = add1 *2;
					Vector add2 = rayDirIntersecting.negative()+ scalar2;
					Vector reflection_direction = add2.normalize();

					double specular = reflection_direction.dotProduct(light_direction);
					if (specular > 0) {
						specular = pow(specular, 10);
						colorFinal = colorFinal+(colorScalar(light_sources.at(light_index)->getLightColor(),specular*ColorIntersectObj.special));
					}
				}

			}

		}
	}

	return clip(colorFinal);
}



void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			cameraHeight -= 3.0;
			break;
		case GLUT_KEY_UP:		// up arrow key
			cameraHeight += 3.0;
			break;

		case GLUT_KEY_RIGHT:
			cameraAngle += 0.03;
			break;
		case GLUT_KEY_LEFT:
			cameraAngle -= 0.03;
			break;

		case GLUT_KEY_PAGE_UP:
			break;
		case GLUT_KEY_PAGE_DOWN:
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}
const unsigned int W = 640;
const unsigned int H = 480;
 unsigned int data[W][H][3];

void display(){


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

	glClearColor( 0, 0, 0, 1 );
    glClear( GL_COLOR_BUFFER_BIT );

    glDrawPixels( windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_INT, imgdta );

	gluLookAt(100*cos(cameraAngle), 100*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	glutSwapBuffers();
}

void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	drawgrid=0;
	drawaxes=1;
	cameraHeight=100.0;
	cameraAngle=1.0;
	angle=0;

	glClearColor(0,0,0,0);

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

	gluPerspective(80,	1,	1,	10000.0);
}







void initRendering()
{

    freopen("C:\\Users\\Ashikee Abhi\\Google Drive\\ray_tracer_input.txt","r",stdin);
    double sp1rad, sp2rad, sp3rad;

    cin>>sp1rad;
    cout<<"sphere 1 rad:"<<sp1rad<<endl;
    double t1, t2, t3;
    cin>>t1>>t2>>t3;
    cout<<"sphere 1 pos: "<<t1<<" "<<t2<<" "<<t3<<endl;

    cin>>sp1ReflectionIndex;
    Color colors1(t1,t2,t3,sp1ReflectionIndex);
    cin>>t1>>t2>>t3;
    Vector new_sphere_location_one (t1,t2, t3);

    cin>>sp2rad;
    cin>>t1>>t2>>t3;
    cin>>sp2ReflectionIndex;
    Color colors2(t1,t2,t3,sp2ReflectionIndex);
    cin>>t1>>t2>>t3;
    cout<<"sphere 2 pos: "<<t1<<" "<<t2<<" "<<t3<<endl;

    Vector new_sphere_location_two (t1,t2, t3);
    cin>>sp3rad;
    cin>>t1>>t2>>t3;
    cin>>sp3ReflectionIndex;
    Color colors3(t1,t2,t3,sp3ReflectionIndex);
    cin>>t1>>t2>>t3;
    cout<<"sphere 3 pos: "<<t1<<" "<<t2<<" "<<t3<<endl;

    Vector new_sphere_location_three (t1,t2, t3);

    //Camera position
    cin>>t1>>t2>>t3;
    Vector campos (t1, t2, t3);
    cout<<"cam pos pos: "<<t1<<" "<<t2<<" "<<t3<<endl;


    //Camera Lookat
    cin>>t1>>t2>>t3;
    Vector look_at (t1, t2, t3);
    cout<<"cam look at : "<<t1<<" "<<t2<<" "<<t3<<endl;


    Vector O (0,0,0);
	Vector X (1,0,0);
	Vector Y (0,1,0);
	Vector Z (0,0,1);

	// Get other camera vectors
	Vector diff_btw (campos.getVectX() - look_at.getVectX(), campos.getVectY() - look_at.getVectY(), campos.getVectZ() - look_at.getVectZ());

	Vector camdir = diff_btw.negative().normalize();
	Vector camright = Y.crossProduct(camdir).normalize();
	Vector camdown = camright.crossProduct(camdir);
	Camera scene_cam (campos, camdir, camright, camdown);


    //Light input

    cin>>t1>>t2>>t3;

    Color lightColor(t1,t2,t3,0);

    cin>>t1>>t2>>t3;
    Vector light_position (t1,t2,t3);
    cout<<"light pos : "<<t1<<" "<<t2<<" "<<t3<<endl;


    cin>>t1>>t2>>t3;
    colorcb1.red = t1;
    colorcb1.green  = t2;
    colorcb1.blue  = t3;

    cin>>t1>>t2>>t3;

   colorcb2.red = t1;
    colorcb2.green  = t2;
    colorcb2.blue  = t3;




    int depth;
    cin>>depth;
    cout<<depth<<endl;
	double aathreshold = 0.1;
	double aspectratio = (double)windowWidth/(double)windowHeight;
	double ambientlight = 0.5;
	double accuracy = 0.00000001;

    ObjSphere s_s (new_sphere_location_one, sp1rad, colors1);
	ObjSphere s_s2 (new_sphere_location_two, sp2rad, colors2);
	ObjSphere s_s3( new_sphere_location_three, sp3rad, colors3);
	ObjPlane scene_plane (Y, -1, Color(1, 1, 1, 2));


	//Vector light_position (-7,10,-10);
	Light scene_light (light_position, lightColor);
	vector<Src*> light_sources;
	light_sources.push_back(dynamic_cast<Src*>(&scene_light));

	// scene Objs



	vector<Object*> s_b;
	s_b.push_back(dynamic_cast<Object*>(&s_s));
	s_b.push_back(dynamic_cast<Object*>(&s_s2));
	s_b.push_back(dynamic_cast<Object*>(&s_s3));
	s_b.push_back(dynamic_cast<Object*>(&scene_plane));


	int current, aa_index;
	double xamnt, yamnt;
	double tempRed, tempGreen, tempBlue;


	for (int x = 0; x < windowWidth; x++) {
		for (int y = 0; y < windowHeight; y++) {
			current = y*windowWidth + x;

			double tempRed[depth*depth];
			double tempGreen[depth*depth];
			double tempBlue[depth*depth];

			for (int aax = 0; aax < depth; aax++) {
				for (int aay = 0; aay < depth; aay++) {

					aa_index = aay*depth + aax;

					srand(time(0));

					if (depth == 1) {

						if (windowWidth > windowHeight) {
							xamnt = ((x+0.5)/windowWidth)*aspectratio - (((windowWidth-windowHeight)/(double)windowHeight)/2);
							yamnt = ((windowHeight - y) + 0.5)/windowHeight;
						}
						else if (windowHeight > windowWidth) {
							xamnt = (x + 0.5)/ windowWidth;
							yamnt = (((windowHeight - y) + 0.5)/windowHeight)/aspectratio - (((windowHeight - windowWidth)/(double)windowWidth)/2);
						}
						else {
							xamnt = (x + 0.5)/windowWidth;
							yamnt = ((windowHeight - y) + 0.5)/windowHeight;
						}
					}
					else {
						if (windowWidth > windowHeight) {
							xamnt = ((x + (double)aax/((double)depth - 1))/windowWidth)*aspectratio - (((windowWidth-windowHeight)/(double)windowHeight)/2);
							yamnt = ((windowHeight - y) + (double)aax/((double)depth - 1))/windowHeight;
						}
						else if (windowHeight > windowWidth) {
                			xamnt = (x + (double)aax/((double)depth - 1))/ windowWidth;
							yamnt = (((windowHeight - y) + (double)aax/((double)depth - 1))/windowHeight)/aspectratio - (((windowHeight - windowWidth)/(double)windowWidth)/2);
						}
						else {
							xamnt = (x + (double)aax/((double)depth - 1))/windowWidth;
							yamnt = ((windowHeight - y) + (double)aax/((double)depth - 1))/windowHeight;
						}
					}

                    Vector cam_ray_origin = scene_cam.campos;
					Vector cam_ray_direction = (camdir+ ( (camright*(xamnt - 0.5) ) + (camdown*(yamnt - 0.5))  )).normalize();

					Ray cam_ray (cam_ray_origin, cam_ray_direction);

					vector<double> intersections;

					for (int index = 0; index < s_b.size(); index++) {
						intersections.push_back(s_b.at(index)->findIntersection(cam_ray));
					}

					int indexFirstIntersected = getFirstIntersectedObj(intersections);

					if (indexFirstIntersected == -1) {
						tempRed[aa_index] = 0;
						tempGreen[aa_index] = 0;
						tempBlue[aa_index] = 0;
					}
					else{
						if (intersections.at(indexFirstIntersected) > accuracy) {

							Vector pos_intersect = cam_ray_origin + (cam_ray_direction * (intersections.at(indexFirstIntersected)));
							Vector rayDirIntersecting = cam_ray_direction;

							Color intersection_color = renderColorValue(pos_intersect, rayDirIntersecting, s_b, indexFirstIntersected, light_sources, accuracy, ambientlight);

							tempRed[aa_index] = intersection_color.red;
							tempGreen[aa_index] = intersection_color.green;
							tempBlue[aa_index] = intersection_color.blue;
						}
					}
				}
			}

        	double totalRed = 0;
			double totalGreen = 0;
			double totalBlue = 0;

			for (int iRed = 0; iRed < depth*depth; iRed++) {
				totalRed = totalRed + tempRed[iRed];
			}
			for (int iGreen = 0; iGreen < depth*depth; iGreen++) {
				totalGreen = totalGreen + tempGreen[iGreen];
			}
			for (int iBlue = 0; iBlue < depth*depth; iBlue++) {
				totalBlue = totalBlue + tempBlue[iBlue];
			}

			double avgRed = totalRed/(depth*depth);
			double avgGreen = totalGreen/(depth*depth);
			double avgBlue = totalBlue/(depth*depth);

			pixels[current].r = avgRed;
			pixels[current].g = avgGreen;
			pixels[current].b = avgBlue;
		}
	}


}

void ImageGLDisplay(bool isdisplaying)
{

	if(isdisplaying)
	{
	    int cc=0;
		for(int i=0;i<windowWidth;i++)
	        for(int j=0;j<windowHeight;j++)
	    {
	        imgdta[i][j][0]=pixels[cc].r*255*255*255*255;
	        imgdta[i][j][1]=pixels[cc].g*255*255*255*255;
	        imgdta[i][j][2]=pixels[cc].b*255*255*255*255;
	        cc++;

	    }
	}
}

int main (int argc, char *argv[]) {


	initRendering();
    	generateBMP("output.bmp",windowWidth,windowHeight,depth_per_inch,pixels);

	glutInit(&argc,argv);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(200, 300);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    openGLDisplay = true;
	ImageGLDisplay(openGLDisplay);

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL
	return 0;
}

