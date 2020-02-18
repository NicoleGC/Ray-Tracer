/*            PURPOSE : Simple framework for ray-tracing

		PREREQUISITES : matrix.h
		Author: Prof.beauchemin and Nicole Garcia
 */


#include <iostream>
#include "GL/glut.h"
#include "GL/freeglut.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include "matrix.h"



#define X                1
#define Y                2
#define Z                3

#define INFINITE_PLANE   0
#define PLANE            1
#define SPHERE           2
#define CONE		 3
#define BASE             4

#define EPSILON          0.000001
#define N_OBJECTS        1024
#define MAX_INTENSITY    255.0

#define Ex               8.0
#define Ey               8.0
#define Ez              - 7.0

#define Gx               0.0
#define Gy               0.0
#define Gz               -1

#define UPx              0.0
#define UPy              0.0
#define UPz              1.0

#define Lx              7.0
#define Ly              4.0
#define Lz             -7.0

#define Near             1.0
#define Far              25.0

#define THETA            45.0
#define ASPECT           1.5

#define H               200

//#define M_PI 3.14159265

struct window_t {
	int width, height;
};

struct camera_t {
	dmatrix_t UP;
	dmatrix_t E;
	dmatrix_t G;
	dmatrix_t u, v, n;
};

struct color_t {
	double r, g, b;

	color_t() {}
	color_t(double _r, double _g, double _b) {
		r = _r;
		g = _g;
		b = _b;
	}
	color_t(const color_t& initializer) {
		r = initializer.r;
		g = initializer.g;
		b = initializer.b;
	}

	color_t operator*(double scalar) {
		return color_t(r * scalar, g * scalar, b * scalar);
	}

	color_t operator+(const color_t& col) {
		return color_t(r + col.r, g + col.g, b + col.b);
	}
};

struct object_t {
	int type;
	double(*intersection)(dmatrix_t *, dmatrix_t *);
	dmatrix_t *(*normal)(dmatrix_t *);
	dmatrix_t M, Minv;
	color_t specular_color, diffuse_color, ambient_color;
	double density, reflectivity, specular_coeff, diffuse_coeff, ambient_coeff, f;
};

struct light_t {
	dmatrix_t position;
	color_t color;
	color_t intensity;
};

static object_t object[N_OBJECTS];
int nobjects = 0;
object_t *build_object(int object_type, dmatrix_t *M, color_t ambient_color, color_t diffuse_color, color_t specular_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double f, double reflectivity) ;
void OnDisplay();
void OnKeyboard(unsigned char key, int x, int y);
void Draw();

const int nChars = ((int)(H * ASPECT)) * H * 3;

unsigned char frame[nChars];

color_t foregroundColor;

void initGLUT(int argc, char** argv, window_t& window) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(window.width, window.height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Assignment 4");

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glShadeModel(GL_FLAT);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glutDisplayFunc(OnDisplay);
	glutKeyboardFunc(OnKeyboard);
}

void SetCurrentColorX(double r, double g, double b) {

	foregroundColor.r = r;
	foregroundColor.g = g;
	foregroundColor.b = b;
}

void SetPixelX(window_t& window, int i, int j) {
	if (i >= window.width || j >= window.height)
		return;

	unsigned int index = 3 * (j * window.width + i);
	frame[index] = (int)(255*foregroundColor.r);
	frame[index + 1] = (int)(255*foregroundColor.g);
	frame[index + 2] = (int)(255*foregroundColor.b);
}

void OnDisplay() {
	memset(frame, 255, nChars);
	Draw();

	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels((int)(H * ASPECT), H, GL_RGB, GL_UNSIGNED_BYTE, (GLubyte*)frame);
	glutSwapBuffers();
	glFlush();
}

void QuitX() {
	exit(0);
}

// Allocates and creates a rotation matrix
dmatrix_t *rotate(double Vx, double Vy, double Vz, double angle)

{
	dmatrix_t *I, *J, *V;
	dmatrix_t *temp1,*temp2,*temp3,*temp4, *temp5;

	I = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	J = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	V = (dmatrix_t *)malloc(sizeof(dmatrix_t));

	dmat_alloc(I, 3, 3);
	dmat_alloc(J, 3, 3);
	dmat_alloc(V, 3, 1);

	I = dmat_identity(I);
	J = dmat_init(J, 0.0);

	(*V).m[1][1] = Vx;
	(*V).m[2][1] = Vy;
	(*V).m[3][1] = Vz;

	V = dmat_normalize(V);

	(*J).m[2][3] = -(*V).m[1][1];
	(*J).m[3][2] = (*V).m[1][1];

	(*J).m[1][3] = (*V).m[2][1];
	(*J).m[3][1] = -(*V).m[2][1];

	(*J).m[1][2] = -(*V).m[3][1];
	(*J).m[2][1] = (*V).m[3][1];
	
	temp1=dmat_mult(J, J);
	temp2=dmat_scalar_mult(temp1, 1.0 - cos(angle));
	temp3=dmat_scalar_mult(J, sin(angle));
	temp4=dmat_add(temp3, temp2);
	temp5=dmat_add(I, temp4);
	

	dmatrix_t* ret = to_homogeneous(temp5, 1.0);
	
	delete_dmatrix(I);
	delete_dmatrix(J);
	delete_dmatrix(V);
	delete_dmatrix(temp1);
	delete_dmatrix(temp2);
	delete_dmatrix(temp3);
	delete_dmatrix(temp4);
	delete_dmatrix(temp5);
	return ret;
}

// Allocates and creates a translation matrix
dmatrix_t *translate(double Tx, double Ty, double Tz)

{
	dmatrix_t *T;

	T = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(T, 4, 4);

	T = dmat_identity(T);

	(*T).m[1][4] = Tx;
	(*T).m[2][4] = Ty;
	(*T).m[3][4] = Tz;

	return T;
}

// Allocates and creates a scale matrix
dmatrix_t *scale(double Sx, double Sy, double Sz)

{
	dmatrix_t *S;

	S = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(S, 4, 4);

	S = dmat_identity(S);

	(*S).m[1][1] = Sx;
	(*S).m[2][2] = Sy;
	(*S).m[3][3] = Sz;

	return S;
}

color_t color_init(double r, double g, double b) {

	color_t s;

	s.r = r;
	s.g = g;
	s.b = b;

	return s;
}

color_t color_mult(double a, color_t c) {

	color_t s;

	s.r = a * c.r;
	s.g = a * c.g;
	s.b = a * c.b;

	return s;
}

color_t color_add(color_t c1, color_t c2) {

	color_t s;

	s.r = c1.r + c2.r;
	s.g = c1.g + c2.g;
	s.b = c1.b + c2.b;

	return s;
}

// Create/allocate a light
light_t *build_light(light_t *light, dmatrix_t *position, color_t color, color_t intensity) {

	dmat_alloc(&light->position, 4, 1);

	light->position = *position;
	light->color.r = color.r;
	light->color.g = color.g;
	light->color.b = color.b;
	light->intensity.r = intensity.r;
	light->intensity.g = intensity.g;
	light->intensity.b = intensity.b;
	return light;
}

window_t *build_window(window_t *Window, int height, float aspect) {

	Window->height = height;
	Window->width = (int)(aspect * height);

	return(Window);
}

camera_t *build_camera(camera_t *Camera, window_t *Window) {

	dmat_alloc(&Camera->E, 4, 1);

	Camera->E.m[X][1] = Ex;
	Camera->E.m[Y][1] = Ey;
	Camera->E.m[Z][1] = Ez;
	Camera->E.m[4][1] = 1.0;

	dmat_alloc(&Camera->G, 4, 1);

	Camera->G.m[X][1] = Gx;
	Camera->G.m[Y][1] = Gy;
	Camera->G.m[Z][1] = Gz;
	Camera->G.m[4][1] = 1.0;

	dmat_alloc(&Camera->n, 4, 1);
	Camera->n = *dmat_normalize(dmat_sub(&Camera->E, &Camera->G));
	Camera->n.l = 3;

	dmat_alloc(&Camera->UP, 4, 1);

	Camera->UP.l = 3;

	Camera->UP.m[X][1] = UPx;
	Camera->UP.m[Y][1] = UPy;
	Camera->UP.m[Z][1] = UPz;
	Camera->UP.m[4][1] = 1.0;

	dmat_alloc(&Camera->u, 4, 1);

	Camera->u = *dmat_normalize(dcross_product(&Camera->UP, &Camera->n));
	Camera->v = *dmat_normalize(dcross_product(&Camera->n, &Camera->u));

	return(Camera);
}

dmatrix_t *intersection_coordinates(dmatrix_t *e, dmatrix_t *direction, double t) {

	dmatrix_t *intersection;

	intersection = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(intersection, 4, 1);

	intersection->m[X][1] = e->m[X][1] + direction->m[X][1] * t;
	intersection->m[Y][1] = e->m[Y][1] + direction->m[Y][1] * t;
	intersection->m[Z][1] = e->m[Z][1] + direction->m[Z][1] * t;
	intersection->m[4][1] = 1.0;

	return intersection;
}

double infinite_plane_intersection(dmatrix_t *e, dmatrix_t *d) {

	double t;

	if (fabs(d->m[Z][1]) < EPSILON) {
		t = -1.0;
	}
	else {
		t = -e->m[Z][1] / d->m[Z][1];
		if (t <= 0.0) {
			t = -1.0;
		}
		else {
			t = -1.0*e->m[Z][1] / d->m[Z][1];
		}
	}
	return t;
}

double plane_intersection(dmatrix_t *e, dmatrix_t *d) {

	double t;
	dmatrix_t *intersection;

	if (fabs(d->m[Z][1]) < EPSILON) {
		t = -1.0;
	}
	else {
		t = -e->m[Z][1] / d->m[Z][1];
		if (t <= 0.0) {
			t = -1.0;
		}
		else {
			intersection = intersection_coordinates(e, d, t);
			if ((fabs(intersection->m[X][1]) > 1.0) || (fabs(intersection->m[Y][1]) > 1.0)) {
				t = -1.0;
			}
			delete_dmatrix(intersection);
		}
	}
	return t;
}



double solve_quadratic(double a, double b, double c) {

	double discriminant, t1, t2, min;

	discriminant = b * b - a * c;
	if (discriminant < 0.0) {
		return -1.0;
	}
	else {
		if (discriminant < EPSILON) {
			return -b / a;
		}
		else {
			t1 = -b / a - sqrtl(discriminant) / a;
			t2 = -b / a + sqrtl(discriminant) / a;

			if (t1 < t2) {
				min = t1;
			}
			else {
				min = t2;
			}

			if (min > EPSILON) {
				return min;
			}
			else {
				return -1.0;
			}
		}
	}
}

double sphere_intersection(dmatrix_t *e, dmatrix_t *d) {
		double a = ddot_product(d, d);
		double b = ddot_product(from_homogeneous(e), from_homogeneous(d));
		double c = ddot_product(from_homogeneous(e), from_homogeneous(e)) - 1.0;
	
	return solve_quadratic(a, b, c);
}


/*FUNCTION: cone_intersection
AUTHOR: Rocio Nicole Garcia
DATE: 2019-04-06 
INPUTS: Point e, ray directio d
OUTPUT: cone instersection t */
double cone_intersection(dmatrix_t *e, dmatrix_t *d) {
	double t1,t2;
	dmatrix_t * intersection;

	double a = (d->m[X][1] * d->m[X][1]) + (d->m[Y][1] * d->m[Y][1]) - (d->m[Z][1] * d->m[Z][1]);
	double b = (e->m[X][1] * d->m[X][1]) + (e->m[Y][1] * d->m[Y][1]) + (d->m[Z][1]) - (e->m[Z][1] * d->m[Z][1]);
	double c = (e->m[X][1] * e->m[X][1]) + (e->m[Y][1] * e->m[Y][1])  - (1)+ (2 * e->m[Z][1]) - (e->m[Z][1] * e->m[Z][1]);
		
	t1=solve_quadratic(a,b,c);

	intersection = intersection_coordinates(e, d, t1);
		
		if ((intersection->m[Z][1] < 0.0) || intersection->m[Z][1]>1.0){
			t1=-1;
			delete_dmatrix(intersection);
			return t1;
		}
		
		return  t1;
	
}

/*FUNCTION: base_intersection
AUTHOR: Rocio Nicole Garcia
DATE: 2019-04-06 
INPUTS: Point e, ray directio d
OUTPUT: cone's base instersection t */
double base_intersection(dmatrix_t *e, dmatrix_t *d){

	double t;
	dmatrix_t *intersection;



	if(fabs(d->m[Z][1])<EPSILON){
		t=-1.0;
	}
	else {
		t=-e->m[Z][1]/d->m[Z][1];
		if(t<=0.0){
			t=-1.0;
		}
		else{
			intersection=intersection_coordinates(e,d,t);

			if(((intersection->m[X][1]*intersection->m[X][1])+(intersection->m[Y][1]*intersection->m[Y][1]))>1.0){
					t=-1.0;
			}
		delete_dmatrix(intersection);
		}
	}
		return t;	
}


/*FUNCTION: base_normal
AUTHOR: Rocio Nicole Garcia
DATE: 2019-04-06 
INPUTS:  Intersection point
OUTPUT: cone's base normal vector*/
dmatrix_t *base_normal(dmatrix_t *intersection){
	
	dmatrix_t *normal;
	normal = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(normal, 4, 1);

	normal->m[X][1] = 0.0;
	normal->m[Y][1] = 0.0;
	normal->m[Z][1] = 1.0;
	normal->m[4][1] = 0.0;

	return dmat_normalize(normal);


}

dmatrix_t *sphere_normal(dmatrix_t *intersection) {

	dmatrix_t *normal;

	normal = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(normal, 4, 1);

	normal->m[X][1] = intersection->m[X][1];
	normal->m[Y][1] = intersection->m[Y][1];
	normal->m[Z][1] = intersection->m[Z][1];
	normal->m[4][1] = 0.0;

	return dmat_normalize(normal);
}

dmatrix_t *plane_normal(dmatrix_t *intersection) {

	dmatrix_t *normal;

	normal = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(normal, 4, 1);

	normal->m[X][1] = 0.0;
	normal->m[Y][1] = 0.0;
	normal->m[Z][1] = -1.0;
	normal->m[4][1] = 0.0;

	return dmat_normalize(normal);
}



/*FUNCTION: cone_normal
AUTHOR: Rocio Nicole Garcia
DATE: 2019-04-03 
INPUTS: intersection point
OUTPUT: cone'snormal vector*/
dmatrix_t *cone_normal(dmatrix_t *intersection) {
	dmatrix_t *normal;
	normal = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(normal, 4, 1);

         
	normal->m[X][1] = intersection->m[X][1];
	normal->m[Y][1] = intersection->m[Y][1];
	normal->m[Z][1] = 1-intersection->m[Z][1];
	normal->m[4][1] = 0.0;

	return dmat_normalize(normal);



}

int find_min_hit_time(double t0[N_OBJECTS]) {

	double min_t = std::numeric_limits<double>::max();
	int position = -1;

	for (int i = 0; i < nobjects; i++) {
		if (t0[i] != -1.0) {
			if (t0[i] < min_t) {
				min_t = t0[i];
				position = i;
			}
		}
	}
	
	return position;
}

dmatrix_t *ray_direction(camera_t *Camera, window_t *Window, double height, double width, double i, double j) {

	int k;
	dmatrix_t *d;

	d = (dmatrix_t *)malloc(sizeof(dmatrix_t));
	dmat_alloc(d, 3, 1);

	for (k = 1; k <= 3; k++) {
		d->m[k][1] = -1.0*Near*Camera->n.m[k][1] + width * (2.0*i / Window->width - 1.0)*Camera->u.m[k][1] + height * (2.0*j / Window->height - 1.0)*Camera->v.m[k][1];
	}

	dmatrix_t* ret = to_homogeneous(d, 0.0);
	delete_dmatrix(d);
	return ret;
}

dmatrix_t *vector_to_light_source(dmatrix_t *intersection, dmatrix_t *light_position) {

	return dmat_normalize(dmat_sub(light_position, intersection));
}

dmatrix_t *vector_to_center_of_projection(dmatrix_t *intersection, dmatrix_t *e) {

	return dmat_normalize(dmat_sub(e, intersection));
}

dmatrix_t *vector_to_specular_reflection(dmatrix_t *N, dmatrix_t *S) {

	return dmat_normalize(dmat_add(dmat_scalar_mult(S, -1.0), dmat_scalar_mult(N, 2.0*ddot_product(N, S))));
}

int shadowed(dmatrix_t *e, dmatrix_t *d) {

	int h, k;
	double t0[N_OBJECTS];

	dmatrix_t *Te, *Td;
	for (k = 0; k < nobjects; k++) {
		Te = dmat_mult(&object[k].Minv, e);
		Td = dmat_mult(&object[k].Minv, d);
		t0[k] = (object[k].intersection)(Te, Td);
		delete_dmatrix(Te);
		delete_dmatrix(Td);

	}
	h = find_min_hit_time(t0);
	return h != -1;
}



/*FUNCTION: shade
AUTHOR: Rocio Nicole Garcia
DATE: 2019-04-02
INPUTS: light, object, camera location coords, ray direction d, color, background color, level of reflectivity
OUTPUT: objects' pixels shade*/
color_t shade(light_t *light, object_t *object, dmatrix_t *e, dmatrix_t *d, color_t color, color_t background, int level) {

	int h,shadows;
	color.r = 0.0;
	color.b = 0.0;
	color.g = 0.0;
	double t0[N_OBJECTS];
	
	dmatrix_t *Te, *Td, *Ts, *normal, *intersection_point, *s_vector, *r_vector, *v_vector, *light_transformed, *intersection_point_transformed,*temp,*temp2;
	double Id, Is;
	
	for (int k = 0; k < nobjects; k++) {
		Te = dmat_mult(&object[k].Minv, e);
		Td = dmat_mult(&object[k].Minv, d);
		t0[k] = (object[k].intersection)(Te,Td);
		delete_dmatrix(Te);
		delete_dmatrix(Td);
		
	}
	
	h = find_min_hit_time(t0); //intersection t


	if (h != -1) {//if it intersected the object do this
		Ts = dmat_mult(&object[h].Minv, &light->position);
		Te = dmat_mult(&object[h].Minv, e);
		Td = dmat_mult(&object[h].Minv, d);
		intersection_point = intersection_coordinates(Te, Td, t0[h]);//find intersection coords
	
		
		//calculate normal
		temp=from_homogeneous(intersection_point);
		temp2 = dmat_normalize(temp);
		normal = to_homogeneous(temp2,0);
		
	
		
		//calculate shading vectors
		s_vector = vector_to_light_source(intersection_point, Ts);
		r_vector = vector_to_specular_reflection(normal, s_vector);
		v_vector = vector_to_center_of_projection(intersection_point, Te);

		
		Id = ddot_product(normal, s_vector);//diffuse intensity
		Is = ddot_product(r_vector, v_vector); //specular intensity
		
		
		//calculate existence of shadows
		light_transformed = dmat_mult(&object[h].M, s_vector);
		intersection_point_transformed = dmat_mult(&object[h].M, intersection_point);
		shadows = shadowed(intersection_point_transformed, light_transformed);
		

		//delete all matrices
		delete_dmatrix(Te);
		delete_dmatrix(Td);
		delete_dmatrix(Ts);
		delete_dmatrix(normal);
		delete_dmatrix(intersection_point);
		delete_dmatrix(s_vector);
		delete_dmatrix(r_vector);
		delete_dmatrix(v_vector);
		delete_dmatrix(light_transformed);
		delete_dmatrix(intersection_point_transformed);
		delete_dmatrix(temp);
                delete_dmatrix(temp2);


		if (Id < 0.0) {
			Id = 0.0;
		}
		if (Is < 0.0) {
			Is = 0.0;
		}
		else  if (Is > 0.0) {
			Is = pow(Is, object[h].f);
		}
			

		//check for shadows, if shadowed only shade with ambient lighting
		if (shadows == 1.0) {

			color.r = light->color.r *light->intensity.r*(object[h].ambient_coeff*object[h].ambient_color.r);
			color.g = light->color.g *light->intensity.g*(object[h].ambient_coeff*object[h].ambient_color.g);
			color.b = light->color.b *light->intensity.b*(object[h].ambient_coeff*object[h].ambient_color.b);
			
		}
		
		//else if not shadowed then compute colour normally
		else {
			//red components
			double ambientRed = light->color.r *light->intensity.r*(object[h].ambient_coeff*object[h].ambient_color.r);
			double diffuseRed = Id * object[h].diffuse_coeff*object[h].diffuse_color.r;
			double specularRed = Is * object[h].specular_coeff*object[h].specular_color.r;

			color.r = ambientRed + diffuseRed + specularRed;

			//green components
			double ambientGreen = light->color.g *light->intensity.g*(object[h].ambient_coeff*object[h].ambient_color.g);
			double diffuseGreen = Id * object[h].diffuse_coeff*object[h].diffuse_color.g;
			double specularGreen = Is * object[h].specular_coeff*object[h].specular_color.g;

			color.g = ambientGreen + diffuseGreen + specularGreen;

			//blue components
			double ambientBlue = light->color.b *light->intensity.b*(object[h].ambient_coeff*object[h].ambient_color.b);
			double diffuseBlue = Id * object[h].diffuse_coeff*object[h].diffuse_color.b;
			double specularBlue = Is * object[h].specular_coeff*object[h].specular_color.b;

			color.b = ambientBlue + diffuseBlue + specularBlue;
		}

	}
		
	return color;
}

object_t *build_object(int object_type, dmatrix_t *M, color_t ambient_color, color_t diffuse_color, color_t specular_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double f, double reflectivity) {

	object_t *object;

	object = (object_t*)malloc(sizeof(*object));
	object->type = object_type;

	dmat_alloc(&object->M, 4, 4);
	object->M = *dmat_duplicate(M);

	dmat_alloc(&object->Minv, 4, 4);
	object->Minv = *dmat_inverse(&object->M);

	object->specular_color = color_init(specular_color.r, specular_color.g, specular_color.b);
	object->diffuse_color = color_init(diffuse_color.r, diffuse_color.g, diffuse_color.b);
	object->ambient_color = color_init(ambient_color.r, ambient_color.g, ambient_color.b);

	object->specular_coeff = specular_coeff;
	object->diffuse_coeff = diffuse_coeff;
	object->ambient_coeff = ambient_coeff;

	object->f = f;
	object->reflectivity = reflectivity;

	switch (object_type) {

	case SPHERE:

		object->intersection = &sphere_intersection;
		object->normal = &sphere_normal;
		break;

	case PLANE:
		object->intersection = &plane_intersection;
		object->normal = &plane_normal;
		break;

	case INFINITE_PLANE:

		object->intersection = &infinite_plane_intersection;
		object->normal = &plane_normal;
		break;
	
	case CONE:
		object->intersection = &cone_intersection;
		object->normal = &cone_normal;
		break;
	case BASE:
		object->intersection = &base_intersection;
		object->normal=&base_normal;

	default:
		break;

	}
	//nobjects++;
	return(object);
}

camera_t Camera;
window_t Window;
light_t light;
color_t background;

int main(int argc, char** argv) {
	/* Set the background color */

	background = color_init(0.0, 0.0, 0.0);

	/* Set up light position, intensity, and color */

	dmatrix_t M, light_position;
	dmat_alloc(&light_position, 4, 1);

	light_position.m[X][1] = Lx;
	light_position.m[Y][1] = Ly;
	light_position.m[Z][1] = Lz;
	light_position.m[4][1] = 1.0;

	color_t light_intensity = color_init(1.0, 1.0, 1.0);
	color_t light_color = color_init(1.0, 1.0, 1.0);
	light = *build_light(&light, &light_position, light_color, light_intensity);

	/* Build display window and synthetic camera */

	Window = *build_window(&Window, H, ASPECT);
	Camera = *build_camera(&Camera, &Window);


	/*shared lighting properties*/
	double specular_coeff = 0.2;
	double diffuse_coeff = 0.5;
	double ambient_coeff = 0.3;

	double f = 10.0;
	double reflectivity = 0.0;

	color_t specular_color = color_init(1.0, 1.0, 1.0);//this property will be the same for all shapes	
	
	/*Build 2 cones*/

	//cone 1
	
	
	color_t diffuse_color = color_init(0.7, 0.4, 1.0); //cone 1's color info
	color_t ambient_color = color_init(0.7, 0.4, 1.0);

	
	dmatrix_t Mc = *translate(6.0,2.0,-4.5);//transformation matrix for cone 1
	
	object_t base = *build_object(BASE, &Mc, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity); //build cone 1's cap
	object[nobjects] = base;
	nobjects++;

	object_t cone = *build_object(CONE, &Mc, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity); //build cone 1
	object[nobjects] = cone;
	nobjects++;

	//cone 2


	diffuse_color = color_init(1.0, 1.0, 0.0);  //cone 2's colour info
	ambient_color = color_init(1.0, 1.0, 0.0);

	dmatrix_t Mc0 =*scale(0.7,0.7,0.7); //build transformation matrix for cone 2
	dmatrix_t Mc1 =*rotate(1.0,1.0,2.0,180);
	dmatrix_t Mc2=*translate(6.0,6.5,-4.0);
	dmatrix_t Mc3=*dmat_mult(&Mc2,&Mc1);
	Mc= *dmat_mult(&Mc0,&Mc3);

	object_t base2 = *build_object(BASE, &Mc, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity); //build cone 2's cap
	object[nobjects] = base2;
	nobjects++;

	object_t cone2 = *build_object(CONE, &Mc, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity); //build cone 2
	object[nobjects] = cone2;
	nobjects++;


	/* Build 4 spheres */
	
	//sphere 1

	diffuse_color = color_init(1.0, 0.4, 0.7);//sphere 1's colour info
	ambient_color = color_init(1.0, 0.4, 0.7);
	
	dmatrix_t Msa=*scale(0.5,0.5,0.5); //build transformation matrix for sphere 1
	dmatrix_t Msb = *translate(6.0, 6.0, -5.0);
	M=*dmat_mult(&Msa,&Msb);

	object_t sphere1 = *build_object(SPHERE, &M, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity);//build sphere 1
	object[nobjects] = sphere1;
	nobjects++;

	//sphere 2

	diffuse_color = color_init(0.0, 1.0, 1.0);//sphere 2's colour info
	ambient_color = color_init(0.0, 1.0, 1.0);
	

	Msa= *scale(0.5,1.0,1.0);
	Msb = *translate(0.0, 4.0, -4.0);
	M = *dmat_mult(&Msb, &Msa); //build transformation matrix for sphere 2
	


	object_t sphere2 = *build_object(SPHERE, &M, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity);//build sphere 2
	object[nobjects] = sphere2;
	nobjects++;

	//sphere 3

	diffuse_color = color_init(0.4, 1.0, 0.4);//sphere 3's colour info
	ambient_color = color_init(0.4, 1.0, 0.4);
	
	Msa= *scale(2.0,0.5,1.0);
	Msb = *translate(3.0, 2.0, -4.0);
	M = *dmat_mult(&Msb, &Msa);//build transformation matrix for sphere 3
	
	

	object_t sphere3 = *build_object(SPHERE, &M, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity);//build pshere 3
	object[nobjects] = sphere3;
	nobjects++;

	//sphere 4

	diffuse_color = color_init(1.0, 0.7, 0.0);//sphere 4's colour info
	ambient_color = color_init(1.0, 0.7, 0.0);
	
	Msa=*scale(1.0,1.0,0.5);
	Msb=*translate(2.0,4.5,-6.0);
	M=*dmat_mult(&Msa,&Msb);//build transformation matrix for sphere 4
	
	object_t sphere4 = *build_object(SPHERE, &M, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity);//build sphere 4
	object[nobjects] = sphere4;
	nobjects++;

	
	/*Build Floor*/
	
	specular_coeff = 0.2;//special lighting properties for floor
	diffuse_coeff = 0.4;
	ambient_coeff = 0.4;
	f=30;
	
	diffuse_color = color_init(1.0, 1.0, 1.0);//floor's colour info
	ambient_color = color_init(1.0, 1.0, 1.0);

	dmatrix_t Mp1=*scale(7.0,7.0,1.0);
	dmatrix_t Mp2 = *translate(0.0,0.0,-2.0);
	dmatrix_t Mp=*dmat_mult(&Mp1, &Mp2);//build floor's transformation matrix
	
	object_t plane = *build_object(PLANE, &Mp, ambient_color, diffuse_color, specular_color, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity);//build floor
	object[nobjects] = plane;
	nobjects++;


	initGLUT(argc, argv, Window);
	glutMainLoop();
	return 0;
}

void Draw() {
	double aspect = ASPECT; /* Set near plane dimensions */
	double height = Near * tan(M_PI / 180.0 * THETA / 2.0);
	double width = height * aspect;

	dmatrix_t *direction;
	int i, j;
	color_t pixel;

	for (i = 0; i < Window.width; i++) {
		for (j = 0; j < Window.height; j++) {
			direction = ray_direction(&Camera, &Window, height, width, (double)i, (double)j);
			pixel = shade(&light, object, &Camera.E, direction, pixel, background, 3);
			SetCurrentColorX(pixel.r, pixel.g, pixel.b);
			SetPixelX(Window, i, Window.height - (j + 1));
			delete_dmatrix(direction);
		}
	}
}

void OnKeyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 'q':
		QuitX();
		break;
	}
}
