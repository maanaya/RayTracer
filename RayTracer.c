#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#define PI 3.142857

//Structs

typedef struct {
    double x, y, z;
    double dx, dy, dz;
} Raytype;

typedef struct{
    double r, g, b;
}ColorType;

typedef struct{
    double x, y, z;
}Vec3;

typedef struct{
    double x, y;
}Vec2;

typedef struct{
	double x, y, z, w, r, g, b;
}light;

typedef struct{
	double odr, odg, odb, osr, osg, osb, ka, kd, ks, n, alpha, ior;
}material;

typedef struct{
	int type;//1 is sphere, 2 is triangle
    double x, y, z;
    double rad;
    material mat;
    int texNum;
    Vec3 v1;
    Vec3 v2;
    Vec3 v3;
    Vec3 vn1;
    Vec3 vn2;
    Vec3 vn3;
    Vec2 vt1;
    Vec2 vt2;
    Vec2 vt3;
}ShapeType;

typedef struct{
	int width, height;
	Vec3 *pixels;
}texture;

typedef struct{
	int top;
	double stackArray[11];
}stack;

/*
Math functions
*/
void multiplyScalar(double num, Vec3 a){
	a.x = a.x * num;
	a.y = a.y * num;
	a.z = a.z * num;
}

Vec3 multiplyScalarNew(double num, Vec3 a){
	Vec3 vec = {a.x * num, a.y * num, a.z * num};
	return vec;
}

Vec3 addNew(Vec3 a, Vec3 b){
	Vec3 vec = {a.x + b.x, a.y + b.y, a.z + b.z};
	return vec;
}

Vec3 subtractNew(Vec3 a, Vec3 b){
	Vec3 vec = {a.x - b.x, a.y - b.y, a.z - b.z};
	return vec;
}


Vec3 cross(Vec3 a, Vec3 b){
	Vec3 vec;
	vec.x = (a.y * b.z) - (a.z * b.y);
	vec.y = (a.z * b.x) - (a.x * b.z);
	vec.z = (a.x * b.y) - (a.y * b.x);
	return vec;
}

double dot(Vec3 a, Vec3 b){
	double n;
	double x = a.x * b.x;
	double y = a.y * b.y;
	double z = a.z * b.z;
	n = x + y + z;
	return n;
}

Vec3 normalize(Vec3 a){
	Vec3 vec;
	double length = sqrt((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
	vec.x = a.x/length;
	vec.y = a.y/length;
	vec.z = a.z/length;
	return vec;
}

double max(double a, double b){
	if(a > b){
		return a;
	}
	else{
		return b;
	}
}

int sphereEqual(ShapeType a, ShapeType b){
	if(a.x == b.x && a.y == b.y && a.z == b.z && a.rad == b.rad){
		return 1;
	}
	return 0;
}

//Global vars
ColorType bkgColor = {0,0,0};
double bkgIndex = 0;
int shapeCounter = 0;
int lightCounter = 0;
Vec3 viewDir = {0,0,0};
light lights[100];
Vec3 eye = {0,0,0};
texture textures[10];

/*
Helper function to calculate 
ray intersections
for various geometry
*/
int allocateTexture(char *file, texture *t){
	FILE *fp;
	fp = fopen(file, "r");
	char line[256];
	fgets(line, sizeof(line), fp);
	char word[10];
	int w;
    int h;
    int rgb;
	if(!sscanf(line, "%s %d %d %d", word, &w, &h, &rgb) == 4){
		fclose(fp);
        return 0;	
    }
    if(strcmp(word, "P3")){
    	fclose(fp);
    	return 0;
    }
    //allocate here
    t->width = w;
    t->height = h;
    t->pixels = (Vec3*)malloc((t->width * t->height) * sizeof(Vec3));
    w = 0;
    h = 0;
    double r;
	double g;
	double b;
	int count = 0;
    while(fscanf(fp, "%lf %lf %lf", &r, &g, &b) == 3){
    	Vec3 pixel = (Vec3) {r,g,b};
    	t->pixels[count] = pixel;
    	count++;
    	w++;
    }
	fclose(fp);
	return 1;
}

int freeTexture(texture *t){
    free(t->pixels);
}

void findThetas(ShapeType shape, Vec3 intersection, double *regTheta, double *polTheta){
	double nz = (intersection.z - shape.z)/shape.rad;
	double ny = (intersection.y - shape.y)/shape.rad;
	double nx = (intersection.x - shape.x)/shape.rad;
	
	*polTheta = acos(nz);
	*regTheta = atan2(ny,nx);
}

void findIntersection(ShapeType shape, Raytype ray, double *tpos, double *tneg){
	if(shape.type == 1){
		double b = 2*(ray.dx *(ray.x-shape.x) + ray.dy * (ray.y - shape.y) + ray.dz *(ray.z-shape.z));
		double c = ((ray.x - shape.x) * (ray.x - shape.x)) + ((ray.y - shape.y) * (ray.y - shape.y)) + ((ray.z - shape.z) * (ray.z - shape.z)) - (shape.rad * shape.rad);
		if(!isnan(sqrt(b*b - 4*c))){
			*tneg = (-b - sqrt(b*b - 4*c))/2;
			*tpos = (-b + sqrt(b*b - 4*c))/2;
		}
	}
	else if(shape.type == 2){
		Vec3 e1 = subtractNew(shape.v2, shape.v1);
		Vec3 e2 = subtractNew(shape.v3, shape.v1);
		Vec3 n = cross(e1, e2);
		double a = n.x;
		double b = n.y;
		double c = n.z;
		
		double d = -1*(a*shape.v1.x + b*shape.v1.y + c*shape.v1.z);
		
		if(a*ray.dx + b*ray.dy + c*ray.dz != 0){
			double placeholder = (-1 * (a*ray.x + b*ray.y + c*ray.z + d))/(a*ray.dx + b*ray.dy + c*ray.dz);
			Vec3 rayDir = {ray.dx, ray.dy, ray.dz};
			Vec3 intersection = multiplyScalarNew(placeholder, rayDir);
			Vec3 ep = subtractNew(intersection, shape.v1);
			
			//linear algebra math to solve for baycentric coordinates
			double d11 = dot(e1, e1);
			double d12 = dot(e1, e2);
			double d22 = dot(e2, e2);
			double d1p = dot(e1, ep);
			double d2p = dot(e2, ep);
			
			double det = (d11*d22 - d12*d12);
			
			if(det == 0){
				*tpos = -1;
			}
			else{
				double b = (d22*d1p - d12*d2p)/det;
				double y = (d11*d2p - d12*d1p)/det;
				if(b < 1 && b > 0 && y < 1 && y > 0 && b + y < 1){
					*tpos = placeholder;
				}
			}
		}
		else{
			*tpos = -1;
		}
	}
}

double Shadow_Checking(ShapeType shape, ShapeType shapes[], Raytype shadowRay, double t){
	double closest = 100000.0;
	double shadowNum = 0;
	int closestShape = -1;
	for(int i = 0; i < shapeCounter; i++){
		ShapeType test = shapes[i];
		double tpos = 100000.0;
		double tneg = 100000.0;
		//if may not be necessary
		//may need to check for something other than 10000
		if(!sphereEqual(shape, test)){
			findIntersection(test, shadowRay, &tpos, &tneg);
			if(tneg < closest && tneg > t){
				closest = tneg;
				closestShape = i;
			}
			if(tpos < closest && tpos > t){
				closest = tpos;
				closestShape = i;
			}
		}
	}
	if(closest == 100000.0){
		return shadowNum;
	}
	else{
		shadowNum = shapes[closestShape].mat.alpha;
		double recVal = Shadow_Checking(shape, shapes, shadowRay, closest);
		return shadowNum + recVal;
	}
}

ColorType Trace_Ray(Raytype ray, ShapeType shapes[], int depth, int oori);

/*
Function to calculate
what color to return to ray713
based on material and
other variables
Will be more in depth later
*/
ColorType Shade_Ray(ShapeType shape, ShapeType shapes[], Raytype ray, double t, int depth, int oori){
	//find i for each rgb and then multiply by 255
	double r = 0.0;
	double g = 0.0;
	double b = 0.0;
	double rlspecdiff = 0.0;
	double glspecdiff = 0.0;
	double blspecdiff = 0.0;
	double bayb = 0;
	double bayy = 0;
	Vec3 intersection = {ray.x + t * ray.dx, ray.y + t * ray.dy, ray.z + t * ray.dz};
	Vec3 rayDir = {ray.dx, ray.dy, ray.dz};	
	Vec3 I = normalize(multiplyScalarNew(-1, rayDir));
	Vec3 n = {1,1,1};
	if(shape.type == 1){
		Vec3 center = {shape.x, shape.y, shape.z};
		n = multiplyScalarNew(1/shape.rad, subtractNew(intersection, center));
	}
	else{
		/*
		if smooth  shading, use baycentric coordinates and 
		vertex normals to calculate normal vector
		if regular triangle, just use cross e1,e2
		*/
		Vec3 e1 = subtractNew(shape.v2, shape.v1);
		Vec3 e2 = subtractNew(shape.v3, shape.v1);
		Vec3 bayN = cross(e1, e2);
		double a = bayN.x;
		double b = bayN.y;
		double c = bayN.z;			
		double d = -1*(a*shape.v1.x + b*shape.v1.y + c*shape.v1.z);	
		double det = 0;		
		if(a*ray.dx + b*ray.dy + c*ray.dz != 0){
			double placeholder = (-1 * (a*ray.x + b*ray.y + c*ray.z + d))/(a*ray.dx + b*ray.dy + c*ray.dz);
			Vec3 rayDir = {ray.dx, ray.dy, ray.dz};
			Vec3 intersection = multiplyScalarNew(placeholder, rayDir);
			Vec3 ep = subtractNew(intersection, shape.v1);				
			//linear algebra math to solve for baycentric coordinates
			double d11 = dot(e1, e1);
			double d12 = dot(e1, e2);
			double d22 = dot(e2, e2);
			double d1p = dot(e1, ep);
			double d2p = dot(e2, ep);
			
			det = (d11*d22 - d12*d12);				
			if(det != 0){
				bayb = (d22*d1p - d12*d2p)/det;
				bayy = (d11*d2p - d12*d1p)/det;
			}
		}
		if(shape.vn1.x != 0 && shape.vn1.y != 0 && shape.vn1.z != 0){
			//baycentric for normal					
			if(det == 0){
				n = normalize(cross(e1, e2));
			}
			else{
				if(bayb < 1 && bayb > 0 && bayy < 1 && bayy > 0 && bayb + bayy < 1){
					double a = 1 - (bayb + bayy);
					Vec3 normVn1 = normalize(shape.vn1);
					Vec3 normVn2 = normalize(shape.vn2);
					Vec3 normVn3 = normalize(shape.vn3);
					
					Vec3 aVn1 = multiplyScalarNew(a, normVn1);
					Vec3 bVn2 = multiplyScalarNew(bayb, normVn2);
					Vec3 yVn3 = multiplyScalarNew(bayy, normVn3);
					
					Vec3 aPlusb = addNew(aVn1, bVn2);
					n = normalize(addNew(aPlusb, yVn3));
				}
				else{
					n = normalize(cross(e1, e2));
				}
			}
		}
		else{
			n = normalize(cross(e1, e2));
		}
	}
	Vec3 rayOrigin = {ray.x, ray.y, ray.z};
	double cosThetai = dot(I,n);
	
	//check if inside or not, when not supposed to be
	if(cosThetai < 0 && oori){
		ColorType zeros = {0,0,0};
		return zeros;
	}
	if(cosThetai < 0){
		n = multiplyScalarNew(-1, n);
		cosThetai = dot(I,n);
	}
	
	Vec3 v = normalize(subtractNew(rayOrigin, intersection));
	//v should be -rayDir now?
	material m = shape.mat;
	//calculate texture coords for sphere
	if(!(shape.texNum < 0) && shape.type == 1){
		//calculate for thetas and change m.od
		double regTheta = 0;
		double polTheta = 0;
		findThetas(shape, intersection, &regTheta, &polTheta);
		double u = 0.5 + regTheta/(2*PI);
		double v = polTheta/PI;
		
		int texNum = shape.texNum;
		int w = textures[texNum].width;
		int h = textures[texNum].height;
		
		int wu = (int) w * u;
		int hv = (int) h * v;
		
		m.odr = textures[texNum].pixels[wu + hv * w].x/255;
		m.odg = textures[texNum].pixels[wu + hv * w].y/255;
		m.odb = textures[texNum].pixels[wu + hv * w].z/255;		
	}
	//calculate texture coords for triangle
	if(!(shape.texNum < 0) && shape.type == 2){
		double a = 1 - (bayb + bayy);
		double u = (a * shape.vt1.x + bayb * shape.vt2.x + bayy * shape.vt3.x);
		double v = (a * shape.vt1.y + bayb * shape.vt2.y + bayy * shape.vt3.y);
		
		int texNum = shape.texNum;
		int w = textures[texNum].width;
		int h = textures[texNum].height;
		
		int wu = (int) w * u;
		int hv = (int) h * v;
		m.odr = textures[texNum].pixels[wu + hv * w].x/255;
		m.odg = textures[texNum].pixels[wu + hv * w].y/255;
		m.odb = textures[texNum].pixels[wu + hv * w].z/255;
	}
	for(int i = 0; i < lightCounter; i++){
		light tlight = lights[i];
		if(tlight.w == 1){
			//point
			Vec3 p = {tlight.x, tlight.y, tlight.z};
			Vec3 pminusx = subtractNew(p,intersection);
			Vec3 l = normalize(pminusx);
			Vec3 h = normalize(addNew(l,v));

			double rdiffuse = m.kd*m.odr*max(0,dot(n,l));
			double gdiffuse = m.kd*m.odg*max(0,dot(n,l));
			double bdiffuse = m.kd*m.odb*max(0,dot(n,l));
			
			double rspecular = m.ks*m.osr*pow(max(0,dot(n,h)),m.n);
			double gspecular = m.ks*m.osg*pow(max(0,dot(n,h)),m.n);
			double bspecular = m.ks*m.osb*pow(max(0,dot(n,h)),m.n);
			
			//define ray to check for shadow, loop through each object and see if it intersects
			Raytype shadowray = {intersection.x, intersection.y, intersection.z, l.x, l.y, l.z};
			double closest = pminusx.x/l.x;
			for(int i = 0; i < shapeCounter; i++){
				ShapeType test = shapes[i];
				double tpos = 100000.0;
				double tneg = 100000.0;
				if(!sphereEqual(shape, test)){
					findIntersection(test, shadowray, &tpos, &tneg);
					if(tneg < closest && tneg > 0){
						closest = tneg;
					}
					if(tpos < closest && tpos > 0){
						closest = tpos;
					}
				}
			}
			//check to see if there was object in way
			if(closest == pminusx.x/l.x){
				rlspecdiff += tlight.r * (rdiffuse + rspecular);
				glspecdiff += tlight.g * (gdiffuse + gspecular);
				blspecdiff += tlight.b * (bdiffuse + bspecular);
			}
		}
		if(tlight.w == 0){
			//directional			
			Vec3 dl = {tlight.x, tlight.y, tlight.z};
			Vec3 negdl = multiplyScalarNew(-1, dl);
			Vec3 l = normalize(negdl);
			Vec3 h = normalize(addNew(l,v));

			double rdiffuse = m.kd*m.odr*max(0,dot(n,l));
			double gdiffuse = m.kd*m.odg*max(0,dot(n,l));
			double bdiffuse = m.kd*m.odb*max(0,dot(n,l));
			
			double rspecular = m.ks*m.osr*pow(max(0,dot(n,h)),m.n);
			double gspecular = m.ks*m.osg*pow(max(0,dot(n,h)),m.n);
			double bspecular = m.ks*m.osb*pow(max(0,dot(n,h)),m.n);
			
			//define ray to check for shadow, loop through each object and see if it intersects
			Raytype shadowray = {intersection.x, intersection.y, intersection.z, l.x, l.y, l.z};
			double closest = 100000.0;
			for(int i = 0; i < shapeCounter; i++){
				ShapeType test = shapes[i];
				double tpos = 100000.0;
				double tneg = 100000.0;
				if(!sphereEqual(shape, test)){
					findIntersection(test, shadowray, &tpos, &tneg);
					if(tneg < closest && tneg > 0){
						closest = tneg;
					}
					if(tpos < closest && tpos > 0){
						closest = tpos;
					}
				}
			}
			//check to see if there was object in way
			if(closest == 100000.0){
				rlspecdiff += tlight.r * (rdiffuse + rspecular);
				glspecdiff += tlight.g * (gdiffuse + gspecular);
				blspecdiff += tlight.b * (bdiffuse + bspecular);
			}
		}
	}
	double rambient = m.ka * m.odr;
	double gambient = m.ka * m.odg;
	double bambient = m.ka * m.odb;
	
	r = 255*(rambient + rlspecdiff);
	g = 255*(gambient + glspecdiff);
	b = 255*(bambient + blspecdiff);
    ColorType computed_color = {r, g, b};
    
    
    //specular reflection calculations (f and r values)
	double aSpecReflect = dot(n, I);
	Vec3 A = multiplyScalarNew(aSpecReflect, n);
	Vec3 S = subtractNew(A, I);
	Vec3 R = addNew(A, S);
	
	//Fresnel Relfectance
	double ior = m.ior;
	double f0 = pow(((ior-1)/(ior+1)), 2);
	double fr = f0 + (1-f0)*pow(1-cosThetai, 5);
	
	//recurcive call for reflectance
	Raytype Rray = {intersection.x + R.x/100, intersection.y + R.y/100, intersection.z+ R.z/100, R.x, R.y, R.z};
	int nextDepth = depth + 1;
	ColorType added_color = Trace_Ray(Rray, shapes, nextDepth, 1);
	//add to computed color (do more stuff here)
	computed_color.r = computed_color.r + (fr * added_color.r);
	computed_color.g = computed_color.g + (fr * added_color.g);
	computed_color.b = computed_color.b + (fr * added_color.b);
	
    return(computed_color);
}

/*
Function to check if
ray intersects with any
of the in world objects
returns color of the 
closest object
*/
ColorType Trace_Ray(Raytype ray, ShapeType shapes[], int depth, int oori){
	ColorType return_color;
	if(depth == 10){
		return_color = (ColorType) {0,0,0};
		return return_color;
	}
    return_color = bkgColor;
    double closest = 100000.0;
    for(int i = 0; i < shapeCounter; i++){
		ShapeType test = shapes[i];
		double tpos = 100000.0;
		double tneg = 100000.0;
		findIntersection(test, ray, &tpos, &tneg);
		if(tneg < closest && tneg > 0){
			closest = tneg;
			ColorType closestColor = Shade_Ray(test, shapes, ray, closest, depth, oori);
			return_color = closestColor;
		}
		if(tpos < closest && tpos > 0){
			closest = tpos;
			ColorType closestColor = Shade_Ray(test, shapes, ray, closest, depth, oori);
			return_color = closestColor;
		}
	}
    return (return_color);
}

/*
Main function controls file input
and sets up ray tracing formulas
*/
int main(int argc, char *argv[]){
	if(argc != 2){
		printf("Arg error");
		return 1;
	}
    char *filename = argv[1];
    int exit = 1;
    FILE *fp;

    //important variable from file
    Vec2 imSize;
    double hfov;
    Vec3 upDir;
    
    double d = 5.0;
     /*
    Read from command line for filename
    Open filename
    Returns if fopen returns NULL
    */        	
    fp = fopen(filename, "r");
    if(fp == NULL){
        printf("Not a valid filename\n");
        return 1;
    }
    else{ 
    	int imSizePres = 0;
    	int eyePres = 0;
    	int viewDirPres = 0;
    	int hfovPres = 0;
    	int upDirPres = 0;
    	int bkgColorPres = 0;
    	int matCounter = 0;
    	int vCounter = 0;
    	int vnCounter = 0;
    	int vtCounter = 0;
    	int texCounter = 0;
    	shapeCounter = 0;
    	lightCounter = 0;
    	
    	ShapeType shapes[100];
    	material mats[100];
    	Vec3 vertices[100];
    	Vec3 vertexNorms[100];
    	Vec2 vertexTextures[100];
        /*
        Check for keywords to assign values to important variables
        Make Array that lists keyword and just use lots of if statements?
        */
        char line[256];
        while (fgets(line, sizeof(line), fp) && exit){
            char word[10];
            double first;
            double second;
            double third;
            double fourth;
            double fifth;
            double sixth;
            double seventh;
            double eigth;
            double ninth;
            double tenth;
            double eleventh;
            double twelfth;
            int numRead = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", word, &first, &second, &third, &fourth, &fifth, &sixth, &seventh, &eigth, &ninth, &tenth, &eleventh, &twelfth);
            if(numRead >= 1 && numRead <= 13){
            	//printf("numred: %d\n", numRead);
                if(!strcmp(word, "imsize") && numRead == 3 && first > 0 && second > 0){
                	imSize = (Vec2) {first, second};
                	imSizePres = 1;
                }
                if(!strcmp(word, "eye") && numRead == 4){
                	eye = (Vec3) {first, second, third};
                	eyePres = 1;
                }
                if(!strcmp(word, "viewdir") && numRead == 4){
                	viewDir = (Vec3) {first, second, third};
                	viewDirPres = 1;
                }
                if(!strcmp(word, "hfov") && numRead == 2 && first > 0){
                	hfov = first;
                	hfovPres = 1;
                }
                if(!strcmp(word, "updir") && numRead == 4){
                	upDir = (Vec3) {first, second, third};
                	upDirPres = 1;
                }
                if(!strcmp(word, "bkgcolor") && numRead == 5 && first >= 0 && first <= 1.0 && second >= 0 && second <= 1.0 && third >= 0 && third <= 1.0){
                	bkgColor = (ColorType) {first * 255, second * 255, third * 255};
                	bkgIndex = fourth;
                	bkgColorPres = 1;
                }
                if(!strcmp(word, "light") && numRead == 8){
                	lights[lightCounter] = (light) {first, second, third, fourth, fifth, sixth, seventh};
                	lightCounter++;
                }
                //might need to be numRead = 13
                if(!strcmp(word, "mtlcolor") && numRead == 13 && first >= 0 && first <= 1.0 && second >= 0 && second <= 1.0 && third >= 0 && third <= 1.0 && fourth >= 0 && fourth <= 1.0 && fifth >= 0 && fifth <= 1.0 && sixth >= 0 && sixth <= 1.0 && seventh >= 0 && seventh <= 1.0 && eigth >= 0 && eigth <= 1.0 && ninth >= 0 && ninth <= 1.0){
                	mats[matCounter] = (material) {first, second, third, fourth, fifth, sixth, seventh, eigth, ninth, tenth, eleventh, twelfth};
                	matCounter++;
                }
                if(!strncmp(word, "texture", 8)){
                	char holder[10];
                	char texFile[20];
                	if(sscanf(line, "%s %s", holder, texFile) == 2){
                		texture t = (texture){0,0};
                		allocateTexture(texFile, &t);
                		textures[texCounter] = t;
                		texCounter++;
                	}
                }
                if(!strcmp(word, "sphere") && numRead == 5 && fourth > 0){
                	if(matCounter > 0){
                		material m = mats[matCounter-1];
                    	shapes[shapeCounter] = (ShapeType){1,first, second, third, fourth, m, texCounter - 1};
                    	shapeCounter++;
                	}
                }
                if(!strcmp(word, "v") && numRead == 4){
                	vCounter++;
                	vertices[vCounter] = (Vec3) {first, second, third};
                }
                if(!strcmp(word, "vn") && numRead == 4){
                	vnCounter++;
                	Vec3 holder = (Vec3) {first, second, third};
                	vertexNorms[vnCounter] = normalize(holder);
                }
                if(!strcmp(word, "vt") && numRead == 3){
                	vtCounter++;
                	vertexTextures[vtCounter] = (Vec2) {first, second};
                }
                if(!strcmp(word, "f")){
                	char holder[10];
                	if(sscanf(line, "%s %lf/%lf/%lf %lf/%lf/%lf %lf/%lf/%lf", holder, &first, &second, &third, &fourth, &fifth, &sixth, &seventh, &eigth, &ninth) == 10){
                		material m = mats[matCounter-1];
                    	shapes[shapeCounter] = (ShapeType){2, 0, 0, 0, 0, m, texCounter - 1,vertices[(int) first], vertices[(int) fourth], vertices[(int) seventh], vertexNorms[(int) second], vertexNorms[(int) fifth], vertexNorms[(int) eigth], vertexTextures[(int) third], vertexTextures[(int) sixth], vertexTextures[(int) ninth]};
                    	shapeCounter++;
                	}
                	else if(sscanf(line, "%s %lf//%lf %lf//%lf %lf//%lf", holder, &first, &second, &third, &fourth, &fifth, &sixth) == 7){
                		material m = mats[matCounter-1];
                    	shapes[shapeCounter] = (ShapeType){2, 0, 0, 0, 0, m, texCounter - 1,vertices[(int) first], vertices[(int) third], vertices[(int) fifth], vertexNorms[(int) second], vertexNorms[(int) fourth], vertexNorms[(int) sixth]};
                    	shapeCounter++;
                	}
                	else if(sscanf(line, "%s %lf/%lf %lf/%lf %lf/%lf", holder, &first, &second, &third, &fourth, &fifth, &sixth) == 7){
                		material m = mats[matCounter-1];
                		Vec3 zero = (Vec3) {0,0,0};
                    	shapes[shapeCounter] = (ShapeType){2, 0, 0, 0, 0, m, texCounter - 1,vertices[(int) first], vertices[(int) third], vertices[(int) fifth], zero, zero, zero, vertexTextures[(int) second], vertexTextures[(int) fourth], vertexTextures[(int) sixth]};
                    	shapeCounter++;
                	}
                	else if(sscanf(line, "%s %lf %lf %lf", holder, &first, &second, &third) == 4){
                		material m = mats[matCounter-1];
                    	shapes[shapeCounter] = (ShapeType){2, 0, 0, 0, 0, m, texCounter - 1, vertices[(int) first], vertices[(int) second], vertices[(int) third]};
                    	shapeCounter++;
                	}
                }
            }
        }
        fclose(fp);
        //do a logic check to make sure we got everything
        if(!imSizePres || !eyePres || !viewDirPres || !hfovPres || !upDirPres || !bkgColorPres){
        	//printf("imSize: %d\neyePres: %d\nviewDir: %d\nhfov: %d\nupDir: %d\nbkg: %d\n", imSizePres, eyePres, viewDirPres, hfovPres, upDirPres, bkgColorPres);
        	printf("Incorrect formatting\n");
        }
        else{
        	//pixel array
        	//ColorType pixels[(int)(imSize.x * imSize.y)];
        	ColorType *pixels = (ColorType*)malloc((imSize.x * imSize.y) * sizeof(ColorType));
        
            //calculate u, v, and w values here
            Vec3 w = {-viewDir.x, -viewDir.y, -viewDir.z};
            Vec3 uPrime = cross(viewDir, upDir);
            Vec3 u = normalize(uPrime);;
            Vec3 vPrime = cross(u, viewDir);
            Vec3 v = normalize(vPrime);
            
            //calculating constants to find four corners
            Vec3 n = normalize(viewDir);
            Vec3 dn =  multiplyScalarNew(d, n);
            double aspectRatio = imSize.x/imSize.y;
            double hfovRad = hfov * (PI /180);
            //weird thing happening with tan, not giving whole, so need to cast it as int
            double width = (int)(2*d * tan(hfovRad/2));
            double height = (int)(width/aspectRatio);
            
            Vec3 wTwo = multiplyScalarNew(width/2, u);
            Vec3 hTwo = multiplyScalarNew(height/2, v);
            Vec3 vPlusDN = addNew(eye, dn);
            
            //calculate corners
            Vec3 ur = addNew(vPlusDN, addNew(wTwo, hTwo));
            Vec3 ul = addNew(vPlusDN, addNew(multiplyScalarNew(-1,wTwo), hTwo));
            Vec3 lr = addNew(vPlusDN, addNew(wTwo, multiplyScalarNew(-1,hTwo)));
            Vec3 ll = addNew(vPlusDN, addNew(multiplyScalarNew(-1,wTwo), multiplyScalarNew(-1,hTwo)));

            Vec3 deltah = multiplyScalarNew(1/(imSize.x - 1), subtractNew(ur, ul));
            Vec3 deltav = multiplyScalarNew(1/(imSize.y - 1), subtractNew(ll, ul));

            int pCounter = 0;
            for(int i = 0; i < imSize.y; i++){
            	//htracker = vtracker;
            	for(int j = 0; j < imSize.x; j++){
            		//calculate ray and then put in pixel array
            		Vec3 jxdh = multiplyScalarNew(j, deltah);
            		Vec3 ixdv = multiplyScalarNew(i, deltav);
            		Vec3 pt = addNew(addNew(ul, jxdh), ixdv);
            		Vec3 rayDir = normalize(subtractNew(pt, eye));
            		Raytype ting = {eye.x, eye.y, eye.z,rayDir.x,rayDir.y,rayDir.z};
            		ColorType pColor = Trace_Ray(ting, shapes, 1, 1);
            		pixels[pCounter] = pColor;
            		pCounter++;
            	}
            }
            //write to file
            char *pp = strtok(filename, ".");
            char ppm[100];
            strcpy(ppm,pp);
            strcat(ppm, ".ppm");
            fp = fopen(ppm, "w+");
            fprintf(fp, "P3\n%d %d\n255\n", (int)imSize.x, (int)imSize.y);
            for(int i = 0; i < imSize.x * imSize.y; i++){
            	fprintf(fp, "%d %d %d\n", (int)pixels[i].r, (int)pixels[i].g, (int)pixels[i].b);
            }
            for(int i = 0; i < texCounter; i++){
            	freeTexture(&textures[i]);
            }
            fclose(fp);
            free(pixels);
        }
    } 
	return 0;
}
