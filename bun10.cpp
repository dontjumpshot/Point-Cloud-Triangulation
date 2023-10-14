#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <queue>
#include <list>
#include <vector>
//1363680307
//1363682010
//1363682772
//1363682932
//1363683031
#define EPS .00001

int maxdepth = -1;

float xcoor[200000];
float ycoor[200000];
float zcoor[200000];

float minx  =  100000;
float maxx  = -100000;
float miny  =  100000;
float maxy  = -100000;
float minz  =  100000;
float maxz  = -100000;
	
	
//用于显示文件中的点
FILE *fp;
int n  = 200000;
float zoom = 2.5;
char etext[1000];

typedef struct mypoint{
	float x, y, z;
}mypoint;

#define MAXEDGECOUNT 20000
char edgematrix[40000][40000] = {0};
int *xcolist;
int *ycolist;
int tcount=0;


char s[100];

int theta = 45, phi = 45;
int gluLookAt_On = 1;
int width, height;
int mycount = 0;

const double radianFactor = 2 * 3.1415926535 / 360;  //弧度的因素

class Point2d{
	public:
		int x, y;						// 顶点的索引
		Point2d();
		Point2d(int xa, int ya);
};

Point2d::Point2d(){
	x = 0;
	y = 0;
}

Point2d::Point2d(int xa, int ya){
	x = xa;
	y = ya;
}

std::list<int> myoriglist;
std::list<int> mydestlist;
std::list<int> myoriglistc;
std::list<int> mydestlistc;

int myoriglistint[50000];
int mydestlistint[50000];
int edgecount = 0;
//#define MAXEDGECOUNT 20000

class triangle;  //前向声明

class Edge {
	public:
		Edge *nextedge;			// 在逆时针方向上同一个面上的下一条边的指针
		Edge *dualedge;			// 与该边相连的另一条边的指针
		int   origin;
		int   dest;				// 边上的起点和终点
		int num;
		triangle *leftface;		
		triangle *rightface;	// 指向左右两侧面的指针
		Edge();
};

class triangle {
	public:
		triangle *child[3];  //三角形的子三角形数组
		Edge *startingEdge;  //三角形的起始边的指针
		triangle *Locate(int);
		int a, b, c;  //三角形的顶点索引
		triangle(int a, int b, int c);
		triangle();
		void InsertSite(int p);  //用于在三角形中插入一个新的点
};

triangle *start;
Edge *globalstart;			// 最新添加的边

// 由三条边组成的三联体它们将被插入到一个新的位置
class triplet{
	public:
		Edge *e[3];
};

//创建Edge对象时初始化其成员变量的值
inline Edge::Edge(){
	num 	  = 0;
	origin 	  = 0;
	dest   	  = 0;
	nextedge  = NULL;
	dualedge  = NULL;
	leftface  = NULL;
	rightface = NULL;
}

//用于计算三个点组成的三角形的面积
inline float TriArea(int a, int b, int c){
	return (xcoor[b] - xcoor[a])*(ycoor[c] - ycoor[a]) - (ycoor[b] - ycoor[a])*(xcoor[c] - xcoor[a]);	
}

// 如果点a, b, c是逆时针顺序，则返回TRUE
//如果面积为正，则说明abc是逆时针方向的，返回值为1；如果面积为负，则说明abc是顺时针方向的，返回值为0。
int ccw(int a, int b, int c){
	return (TriArea(a, b, c) > 0);
}

/*计算点a、b、c、d到原点的距离的平方，即xcoor[i]*xcoor[i] + ycoor[i]*ycoor[i]
计算三个小三角形的面积TriArea，并用这些面积进行线性组合
判断线性组合是否大于0，如果大于0，则点d位于外接圆内，返回1，否则返回0。*/
int InCircle(triangle *t, int d){
	int a = t->a;	int b = t->b;	int c = t->c;
	return (xcoor[a]*xcoor[a] + ycoor[a]*ycoor[a]) * TriArea(b, c, d) -
			(xcoor[b]*xcoor[b] + ycoor[b]*ycoor[b]) * TriArea(a, c, d) +
			(xcoor[c]*xcoor[c] + ycoor[c]*ycoor[c]) * TriArea(a, b, d) -
			(xcoor[d]*xcoor[d] + ycoor[d]*ycoor[d]) * TriArea(a, b, c) > 0;	
}

/*InTriangle函数使用ccw函数来判断点p与三角形的三个顶点之间的关系。
如果点p与三个顶点的顺序都是逆时针方向（即ccw函数返回值为真），则表示点p在三角形内部，函数返回1；否则，表示点p不在三角形内部，函数返回0。*/
int InTriangle(int p, triangle *node){
	if(ccw(node->a, node->b, p) && ccw(node->b, node->c, p) && ccw(node->c, node->a, p))
		return 1;
	else
		return 0;
}

//检查点x是否在边e的右边
/*RightOf函数使用ccw函数来判断点x相对于边e所组成的有向线段的位置关系。
如果点x在边e的右侧（即边e的起点是顺时针方向，ccw函数返回值为真），则函数返回1；否则，表示点x不在边e的右侧，函数返回0。*/
int RightOf(int x, Edge* e){
	return ccw(x, e->dest, e->origin);
}

//检查点x是否在边e的左边
/*LeftOf函数使用ccw函数来判断点x相对于边e所组成的有向线段的位置关系。
如果点x在边e的左侧（即边e的起点是逆时针方向，ccw函数返回值为真），则函数返回1；否则，表示点x不在边e的左侧，函数返回0。*/
int LeftOf(int x, Edge* e){
	return ccw(x, e->origin, e->dest);
}

/*mylocate函数使用递归的方式进行三角形定位。首先判断起始三角形start是否为叶子节点（即没有子三角形），如果是，则直接返回该三角形。
如果不是叶子节点，则遍历起始三角形的子三角形。

对于每个子三角形，调用InTriangle函数判断点p是否在该子三角形内部。如果是，则递归调用mylocate函数，将该子三角形作为起始三角形，并传递点索引p作为参数。
递归过程会继续在子三角形中寻找点所在的三角形。如果找到匹配的三角形，则将其赋值给e，并返回。

如果在起始三角形和它的所有子三角形中都没有找到点所在的三角形，则返回e（未初始化的变量，可能需要初始化）。*/
triangle *mylocate(triangle *start, int p){
	if(start->child[0] == NULL && start->child[1] == NULL && start->child[2] == NULL){
		return start;
	}

	triangle *e;
	triangle *current;
	current = start;

	for(int i=0;i<3;++i){
		if(current->child[i] != NULL){
			if(InTriangle(p, current->child[i])){
				e = mylocate(current->child[i], p);
				return e;
			}
		}
	}
	return e;
}

/*Locate函数将起始节点start作为参数传递给mylocate函数，同时传递点索引x作为参数。
mylocate函数会遍历三角形网格中的所有三角形，直到找到包含点x的三角形为止。如果找到了这个三角形，mylocate函数会将其返回，然后Locate函数也将它返回。
如果在三角形网格中没有找到包含点x的三角形，mylocate函数会返回一个未初始化的变量，而Locate函数会将其直接返回。*/
triangle* triangle::Locate(int x){
	return mylocate(start, x);
}

/*函数首先将对象的三个顶点a、b、c分别赋值为0，即未初始化的状态。然后将该三角形的起始边startingEdge赋值为空指针NULL，表示该三角形还没有边。
最后，将子三角形指针数组child[3]中的所有元素都赋值为空指针NULL，表示该三角形还没有子三角形。*/
triangle::triangle(){
	a = 0;	b = 0; c = 0;
	startingEdge = NULL;
	for(int i=0;i<3;++i)
			child[i] = NULL;
}

/*带参构造函数首先通过调用ccw函数检查顶点a1、b1和c1是否满足逆时针（CCW）的顺序。如果不满足，则交换顶点b1和c1的值，确保顶点按照逆时针的顺序排列。
接着，构造函数将参数a1、b1和c1分别赋值给三角形对象的成员变量a、b和c。
然后，创建三个新的边对象e1、e2和e3，并通过设置nextedge、dualedge、leftface和rightface等属性，将它们连接成一个环形链表，并与当前的三角形对象关联起来。
最后，构造函数将边e1设置为三角形对象的起始边startingEdge，并将子三角形指针数组child[3]中的所有元素都赋值为空指针NULL。*/
triangle::triangle(int a1, int b1, int c1){
	if(!ccw(a1, b1, c1)){	// 检查这三点是否在CCW意义上
		int temp;
		temp = b1;
		b1   = c1;
		c1   = temp;		
	}

	a = a1;
	b = b1;
	c = c1;

	Edge *e1, *e2, *e3;
	e1 = new Edge();
	e2 = new Edge();
	e3 = new Edge();

	e1->nextedge = e2;		e2->nextedge = e3;		e3->nextedge = e1;
	e1->dualedge = NULL;	e2->dualedge = NULL;	e3->dualedge = NULL;
	e1->leftface = this;	e2->leftface = this;	e3->leftface = this;	
	e1->rightface= NULL;	e2->rightface= NULL;	e3->rightface= NULL;
	
	e1->origin = a1; e1->dest = b1;
	e2->origin = b1; e2->dest = c1;
	e3->origin = c1; e3->dest = a1;

	startingEdge = e1;

	for(int i=0;i<3;++i)
		child[i] = NULL;
}

/*函数通过循环遍历triplet结构体中的三个边对象，并对它们对应的对偶边进行修正。
具体地，对于第i条边，将其对偶边的起始点origin设置为该边的终点dest，将其终点dest设置为该边的起始点origin，
将其下一条边nextedge设置为另外两条边中排在中间位置的那一条（即(edges->e[i]->num+2)%3）。*/
void correctdual(triplet *edges){
	for(int i=0;i<3;++i){
		edges->e[i]->dualedge->origin   = edges->e[i]->dest;
		edges->e[i]->dualedge->dest     = edges->e[i]->origin;
		edges->e[i]->dualedge->nextedge = edges->e[(edges->e[i]->num+2)%3];
	}
}

/*
该函数主要用于在进行三角剖分算法时，通过递归修正某条边pipj及其相邻的三角形，使其满足Delaunay三角形的条件。

具体来说，函数首先检查pipj的右边界(rightface)是否为空，如果为空则返回1。这表示pipj是边界边，不需要修正。

如果pipj的右边界不为空，就根据边pipj和它的下一条边(temp)以及对偶边(dual)的顶点信息计算四个顶点v1、v2、v3和v4。然后通过调用InCircle函数判断顶点v4是否在以pipj的左边界(leftface)所代表的三角形内。如果在内部，则进行修正。

在修正过程中，首先更新edgematrix数组，将pipj的起始点和终点之间的边设置为无效。

然后创建两个新的边对象e和edual，并设置它们的顶点信息以及双向连接关系。

接着创建两个新的三角形f1和f2，并将它们设为父子关系。

设置f1和f2的顶点信息，以及其起始边信息。

设置边e和edual的左右边界信息，以及其他相关边的左边界信息。

最后，递归调用legalizeedge函数修正新生成的三角形f1和f2中的边。

如果顶点v4不在pipj的左边界内，函数返回1表示无需修正。

这个函数的作用是确保进行三角剖分时，每条边的对偶边满足Delaunay三角形的要求，以保证最终得到的三角网格质量较高。*/
int legalizeedge(int pr, Edge *pipj, triangle* t){
	if(pipj->rightface == NULL){
		return 1;						
	}	
	
	Edge *temp = pipj->nextedge;
	
	int v1 = pipj->dest;
	int v2 = pipj->nextedge->dest;
	int v3 = pipj->origin;
	int v4 = pipj->dualedge->nextedge->dest;		
	
	if(InCircle(pipj->leftface, v4)){
		edgematrix[pipj->origin][pipj->dest] = '0';
		edgematrix[pipj->dest][pipj->origin] = '0';
		
		Edge *e1 = pipj->nextedge;
		Edge *e2 = pipj->nextedge->nextedge;
		Edge *e3 = pipj->dualedge->nextedge;
		Edge *e4 = pipj->dualedge->nextedge->nextedge;
		
		Edge *e 	= new Edge();
		Edge *edual = new Edge();
		
		e->dualedge 	= edual;
		edual->dualedge = e;
		
		e->origin     = pr;
		e->dest       = v4;
		edual->origin = e->dest;
		edual->dest   = e->origin;
		
		edgematrix[e->origin][e->dest] = '1';
		edgematrix[e->dest][e->origin] = '1';
		
		// make two new faces and set the pointers
		triangle *f1 = new triangle();
		triangle *f2 = new triangle();		 	
		
		t->child[0]  			  = f1;
		t->child[1]  			  = f2;	
		pipj->rightface->child[0] = f1;
		pipj->rightface->child[1] = f2;
		
		// set the points of the faces
		f2->a = pr;		f2->b = v3;		f2->c = v4;
		f1->a = pr;		f1->b = v4;		f1->c = v1;				
		
		f1->startingEdge = e;
		f2->startingEdge = edual;
		
		// set the faces
		e->leftface      = f1;
		e->rightface	 = f2;
		edual->leftface  = f2;
		edual->rightface = f1;
		
		e1->leftface    = f1;	if(e1->dualedge != NULL)	e1->dualedge->rightface = f1;
		e2->leftface	= f2;	if(e2->dualedge != NULL)	e2->dualedge->rightface = f2;
		e3->leftface	= f2;	if(e3->dualedge != NULL)	e3->dualedge->rightface = f2;
		e4->leftface	= f1;	if(e4->dualedge != NULL) 	e4->dualedge->rightface = f1;		
		
		// set the nextedges
		e->nextedge 	= e4;
		e4->nextedge    = e1;
		e1->nextedge    = e;
		edual->nextedge = e2;
		e2->nextedge    = e3;
		e3->nextedge    = edual;				
		
		legalizeedge(pr, e3, f2);
		legalizeedge(pr, e4, f1);
		return 0;
	}
	else{		
		return 1;
	}
}

/*这段代码定义了一个名为InsertSite的函数，它接受一个整数参数x。该函数用于在三角剖分中插入一个新的点x。

函数首先调用Locate函数定位到包含点x的三角形t。

然后创建两个triplet对象edges和dualedges，用于存储新生成的边和对偶边。

接下来进入一个循环，循环3次，用于创建3条边和对偶边，并设置它们的起始点、终点、对偶边的指向关系以及序号。

然后获取三角形t的起始边和其后续两条边，并将它们存储在se数组中。

接着循环遍历edges中的每条边，设置它们的目标点和下一条边，同时更新edgematrix数组，表示这些边的起始点和终点之间的边已经存在。

然后通过修改se数组中的每条边的nextedge属性，修正原来三角形t内的边的循环顺序。

接着调用correctdual函数两次来修正对偶边的属性。

接下来创建三个新的triangle对象face，并分别为它们设置startingEdge、a、b和c属性。

然后连接三角形t和新生成的三角形face，设置边的leftface和rightface属性。

最后，更新全局变量globalstart，增加边的计数edgecount，并调用legalizeedge函数将新生成的三角形face中的边进行修正。

这个函数的作用是向三角剖分中插入一个新的点，并根据该点更新相邻的边和三角形，使得剖分仍然满足Delaunay三角形的要求。*/
void triangle::InsertSite(int x){
	triangle* t = Locate(x);
		
	triplet *edges 	   =  new triplet();
	triplet *dualedges =  new triplet();

	for(int i=0;i<3;++i){
		edges->e[i]     		  = new Edge();
		dualedges->e[i] 		  = new Edge();
		edges->e[i]->dualedge     = dualedges->e[i];
		dualedges->e[i]->dualedge = edges->e[i];
		edges->e[i]->origin   	  = x;
		edges->e[i]->num      	  = i;
	}

	Edge *se[3];
	se[0] = t->startingEdge;
	se[1] = se[0]->nextedge;
	se[2] = se[1]->nextedge;
	
	// setting the destination of new edges
	for(int i=0;i<3;++i){
		edges->e[i]->dest 	  = se[i]->origin;	
		edges->e[i]->nextedge = se[i];
		if(edgematrix[edges->e[i]->dest][edges->e[i]->origin] != '1'){
			edgematrix[edges->e[i]->dest][edges->e[i]->origin]  = '1';
			edgematrix[edges->e[i]->origin][edges->e[i]->dest]	 = '1';	
		}
	}	

	// correct the face cycles
	se[0]->nextedge = edges->e[1]->dualedge;
	se[1]->nextedge = edges->e[2]->dualedge;
	se[2]->nextedge = edges->e[0]->dualedge;

	// correct the properties of the dual edges
	correctdual(edges);
	correctdual(edges);
	correctdual(edges);

	triangle *face[3];
	for(int i=0;i<3;++i){
		face[i] = new triangle();
		face[i]->startingEdge = edges->e[i];
		face[i]->a = x;
		face[i]->b = edges->e[i]->dest;
		face[i]->c = edges->e[i]->nextedge->dest;
		t->child[i] = face[i];
		edges->e[i]->leftface = face[i];
		se[i]->leftface       = face[i];
		if(se[i]->dualedge !=  NULL)
			se[i]->dualedge->rightface  = face[i];
	}
	
	// correcting the right faces
	for(int i=0;i<3;++i){
		edges->e[i]->rightface     = face[(i+2)%3];
		dualedges->e[i]->leftface  = edges->e[i]->rightface;
		dualedges->e[i]->rightface = edges->e[i]->leftface;
	}

	globalstart = edges->e[2];	
	
	edgecount = edgecount+3;
	for(int i=0;i<3;++i)
		legalizeedge(x, se[i], face[i]);	
}

/*这段代码定义了一个名为addnodes的函数，它接受两个参数：一个指向开始的三角形节点的指针start和当前深度depth。

该函数首先检查当前深度是否大于之前记录的最大深度maxdepth，并进行更新。

然后检查当前节点的子节点child[0]是否为NULL，如果是，则说明当前三角形不可再分割，将该三角形的三条边添加到edgematrix数组中，并增加边数计数器edgecount。

如果子节点child[0]不为NULL，则进入一个循环遍历三个子节点，依次递归调用这个函数并将depth加1。

最后，函数返回。这个函数的作用是用于在三角剖分的不同深度节点上添加更多的子节点，直到达到给定的深度为止。这样可以实现更加细致的剖分。*/
// making a list of edges
void addnodes(triangle * start, int depth){
	printf("in addnodes depth = %d\n", depth);
	if(maxdepth<depth)
		maxdepth = depth;
	if(start->child[0] == NULL){
		Edge *et1 = start->startingEdge;				
		
		if(edgematrix[et1->origin][et1->dest] != '1'){
			edgematrix[et1->origin][et1->dest] 				= '1';
			edgematrix[et1->dest][et1->origin] 				= '1';
			++edgecount;
		}
		if(edgematrix[et1->origin][et1->dest] != '1'){
			edgematrix[et1->dest][et1->nextedge->dest] 		= '1';
			edgematrix[et1->nextedge->dest][et1->dest] 		= '1';
			++edgecount;
		}
		if(edgematrix[et1->origin][et1->dest] != '1'){
			edgematrix[et1->nextedge->dest][et1->origin] 	= '1';
			edgematrix[et1->origin][et1->nextedge->dest] 	= '1';
			++edgecount;
		}		
		return;
	}
	
	for(int i=0;i<3;++i){
		if(start->child[i] != NULL){
			addnodes(start->child[i], depth+1);
		}
	}
	return;
}

void *font = GLUT_BITMAP_TIMES_ROMAN_24;  //这段代码定义了一个名为font的指针变量，并将其指向GLUT库中的字体数据——大小为24的TIMES_ROMAN字体，用于在OpenGL应用程序中进行文本渲染。

/*这段代码定义了一个名为outputCharacter的函数，它接受四个参数：x、y和z分别表示文本渲染的位置坐标，string是要渲染的字符串。

函数通过调用glRasterPos3f函数将渲染位置设置为(x, y, z)。然后使用strlen函数获取字符串的长度，并将其转换为整数类型。

接下来进入一个循环，从字符串的第一个字符开始遍历到最后一个字符。在每次循环中，使用glutBitmapCharacter函数将字体font和当前字符string[i]渲染到指定位置。

这个函数的作用是在OpenGL应用程序中绘制指定字符串的文本。根据给定的位置和字体，将字符串逐个字符地渲染到屏幕上。*/
void outputCharacter(float x, float y, float z, char *string) {
  int len, i;
  glRasterPos3f(x, y, z);
  len = (int) strlen(string);
  for (i = 0; i < len; i++) {
    glutBitmapCharacter(font, string[i]);
  }
}

/*这段代码定义了一个名为changeSize的函数，它接受两个参数w和h，分别表示窗口的宽度和高度。

函数将传入的宽度w和高度h分别赋值给全局变量width和height。

然后，通过判断h是否为0，避免除以0的错误，如果h为0，则将其设置为1。计算宽高比ratio，即w除以h的结果。

接下来，通过调用glMatrixMode和glLoadIdentity函数将当前矩阵模式设置为GL_PROJECTION，并将当前矩阵重置为单位矩阵。

使用glViewport函数将视口设置为(0, 0, w, h)，即将整个窗口作为视口。

调用gluPerspective函数来设置透视投影矩阵，参数包括视场角为45度，宽高比为ratio，近裁剪面为1，远裁剪面为1000。

接下来，根据当前的相机参数计算眼睛的位置eyeX、eyeY、eyeZ，中心点的位置centerX、centerY、centerZ，以及上方向的向量upX、upY、upZ。

如果启用了gluLookAt_On标志，调用gluLookAt函数设置相机视图矩阵，参数包括眼睛位置、中心点位置和上方向向量。

调用glScalef函数根据zoom缩放因子对模型视图矩阵进行缩放变换。

最后，将矩阵模式设置为GL_MODELVIEW，表示之后的操作将对模型视图矩阵进行变换。*/
void changeSize(int w, int h){
    width = w;
    height = h;
    
    if(h == 0) h = 1;
    float ratio = 1.0 * w / h;
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();    
    
    glViewport(0, 0, w, h);
	
	gluPerspective(45, ratio, 1, 1000);	
	
	float r = 5.0f;
	float eyeX = r * sin(theta * radianFactor) * cos(phi * radianFactor);
	float eyeY = r * sin(theta * radianFactor) * sin(phi * radianFactor);
	float eyeZ = r * cos(radianFactor * theta);

	float centerX = 0, centerY = 0, centerZ = 0;
	float upX = 0, upY = 1.0f, upZ = 0;
	
	if(gluLookAt_On) {
		gluLookAt(eyeX, eyeY, eyeZ, 
				  centerX, centerY, centerZ,
				  upX, upY, upZ); 
	}

	glScalef(zoom, zoom, zoom);
	glMatrixMode(GL_MODELVIEW);		
}

/*这段代码是一个渲染场景的函数。它使用OpenGL绘制了一些图形元素，并将其显示在屏幕上。

首先，调用glClear函数清除颜色缓冲区和深度缓冲区，以准备开始新的渲染。

然后，调用glLoadIdentity函数将当前矩阵重置为单位矩阵，以确保没有之前的变换操作影响到当前的绘制。

接下来，注释掉了一段绘制三个坐标轴的代码（这部分代码被注释掉了）。

之后，调用glColor3f函数设置绘制图形的颜色为白色。

在使用glBegin(GL_POINTS)和glEnd()之间，使用一个循环来绘制一系列点。循环中的代码根据xcoor、ycoor和zcoor数组中的数据，依次设置每个点的坐标，并使用glVertex3f函数将点绘制出来。

接着，调用glColor3f函数再次设置绘制图形的颜色为白色。

在使用glBegin(GL_LINES)和glEnd()之间，使用一个循环来绘制一系列线段。循环中的代码计算每条线段的长度，并通过判断长度是否小于0.25来确定是否绘制该线段。如果满足条件，则使用glVertex3f函数分别将线段的两个端点绘制出来。

最后，调用glutSwapBuffers函数交换前后缓冲区，以显示绘制的图形。

这个函数通常作为绘制回调函数，在主循环中被调用以重复渲染场景。*/
void renderScene(void){
	int temporig, tempdest;
	
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glLoadIdentity(); 	
	
	/*glColor3f(1.0,1.0,1.0);  
    glBegin(GL_LINES);
        glVertex3f(-0.8f, 0.0f,0.0f);
        glVertex3f(0.8f, 0.0f,0.0f);
    glEnd();
    glBegin(GL_LINES);
        glVertex3f(0.0f, -0.8f,0.0f);
        glVertex3f(0.0f, 0.8f,0.0f);
    glEnd();
    glBegin(GL_LINES);
        glVertex3f(0.0f, 0.0f, -0.8f);
        glVertex3f(0.0f, 0.0f, 0.8f);
    glEnd();*/
        
    glColor3f (1.0, 1.0, 1.0);    
	glBegin(GL_POINTS);
		for(int u=0; u<40000;u = u+5){			
			glVertex3f(xcoor[u], ycoor[u], zcoor[u]);
		}
    glEnd();     
    
    glColor3f(1.0, 1.0, 1.0);	
	glBegin(GL_LINES);
    
    for(int i=0;i<tcount;++i){
		int ti = xcolist[i];
		int tj = ycolist[i];
		float dist = 0;
		dist = sqrt((xcoor[ti]-xcoor[tj])*(xcoor[ti]-xcoor[tj])+(ycoor[ti]-ycoor[tj])*(ycoor[ti]-ycoor[tj])+(zcoor[ti]-zcoor[tj])*(zcoor[ti]-zcoor[tj]));
		if(dist < 0.25){
			glVertex3f(xcoor[ti], ycoor[ti], zcoor[ti]);
			glVertex3f(xcoor[tj], ycoor[tj], zcoor[tj]);
		}
	}	
	glEnd();
	   
    glutSwapBuffers();
}

/*这段代码是一个键盘输入的回调函数，用于处理按键事件。根据按下的键盘字符，执行不同的操作。

首先，使用switch语句根据按下的按键字符进行分支判断。

如果按下的是'i'，则将zoom变量增加0.5，用于实现放大视角的效果。

如果按下的是'o'，则将zoom变量减少0.5，用于实现缩小视角的效果。

如果按下的是't'，则将theta变量增加1，如果超过360，则重置为1。该变量用于控制绕Y轴旋转的角度。

如果按下的是'p'，则将phi变量增加1，如果超过360，则重置为1。该变量用于控制绕X轴旋转的角度。

如果按下的是'T'，则将theta变量减少1，如果小于0，则重置为359。

如果按下的是'P'，则将phi变量减少1，如果小于0，则重置为359。

如果按下的是'g'，则将gluLookAt_On变量取反。gluLookAt_On用于控制是否使用gluLookAt函数来设置视角。

最后，调用changeSize函数，将当前窗口的宽度和高度作为参数传递给该函数，以便在视角改变后更新场景的大小。

这个函数通常作为键盘输入回调函数，在主循环中被注册并调用以处理键盘输入。*/
void inputKey(unsigned char c, int x, int y){
    switch (c) {			
			case 'i' : zoom = zoom+ 0.5; break;
			case 'o' : zoom = zoom-0.5; break;
            case 't' : theta++; if(theta > 360) theta = 1; break;
            case 'p' : phi++; if(phi > 360) phi = 1; break;
            case 'T' : theta--; if(theta < 0) theta = 359; break;
            case 'P' : phi--; if(phi < 0) phi = 359; break;
            case 'g' : gluLookAt_On = !gluLookAt_On;; break;
    }
        changeSize(width, height);
}

/*这段代码是一个初始化函数，用于设置OpenGL的一些初始状态。

首先，调用glClearColor函数，将清除颜色设置为黑色。这意味着在每次绘制之前，屏幕会被清除为黑色。

接着，调用glMatrixMode函数，指定当前操作的矩阵模式为GL_PROJECTION。这表示后续的矩阵操作将影响到投影矩阵。

然后，调用glLoadIdentity函数，将当前的投影矩阵设置为单位矩阵。单位矩阵表示没有任何变换操作被应用。

再次调用glMatrixMode函数，将当前操作的矩阵模式切换回GL_MODELVIEW。这表示后续的矩阵操作将影响到模型视图矩阵。

该函数通常在程序初始化时被调用，用于设置OpenGL的初始状态，比如清除颜色、选择操作的矩阵模式等。*/
void init(){
   glClearColor(0.0, 0.0, 0.0, 0.0);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char **argv){
	int mycount = 0;
	
	for(int i=0;i<11;++i){				// 忽略前七行
		gets(etext);
	}	
	
	float temp1, temp2, temp3;	
	
	for(int i=0;i<n;++i){
		gets(etext);
		if(etext[0] == 'n')	
			continue;		
		
		sscanf(etext, "%f %f %f", &temp1, &temp2, &temp3);
		
		xcoor[mycount] = temp1;
		ycoor[mycount] = temp2;
		zcoor[mycount] = temp3;
		
		if(minx>temp1)
			minx = temp1;
		if(miny>temp2)
			miny = temp2;
		if(minz>temp3)
			minz = temp3;	
		if(maxx<temp1)
			maxx = temp1;
		if(maxy<temp2)
			maxy = temp2;
		if(maxz<temp3)
			maxz = temp3;	
				
		++mycount;
	}	
	
	float meanx = (maxx+minx)/2;
	float meany = (maxy+miny)/2;
	float meanz = (maxz+minz)/2;
	
	printf("PRINTING VALUES mycount =  %d\n", mycount);
	for(int i=0;i< mycount;++i){		
		xcoor[i] = (xcoor[i]-meanx)*1.0/(maxx-minx);
		ycoor[i] = (ycoor[i]-meany)*1.0/(maxy-miny);
		zcoor[i] = (zcoor[i]-meanz)*1.0/(maxz-minz);		
	}
	
	xcoor[0] = -5;	ycoor[0] = -5;
	xcoor[1] = 5; 	ycoor[1] = -5;
	xcoor[2] = 0;  	ycoor[2] = 5;
	
	start = new triangle(0, 1, 2);
	
	unsigned int myseed = time(NULL);
	srand(myseed);
	int thresh = 30;
	int pointsin = 0;
	for(int i=3;i<40000; i = i+1){
		int tt = rand()%10000;
		if(ycoor[i]>0.3)
			thresh = 400;
		else
			thresh = 320;
		if(tt< thresh){
			start->InsertSite(i);
			++pointsin;
		}
	}
	printf("TRIANGULATION DONE seed = %d\n", myseed);
	//addnodes(start, 0);
	
	xcolist = (int *)malloc(sizeof(int)*(edgecount+100));
	ycolist = (int *)malloc(sizeof(int)*(edgecount+100));
	
	
	for(int i=0;i<40000;++i){
		edgematrix[0][i] = '0';
		edgematrix[i][0] = '0';
		edgematrix[1][i] = '0';
		edgematrix[i][1] = '0';
		edgematrix[2][i] = '0';
		edgematrix[i][2] = '0';
	}
	for(int i=0;i<40000;++i){
		for(int j=i+1;j<40000;++j){
			if(edgematrix[i][j] == '1'){
				xcolist[tcount] = i;
				ycolist[tcount] = j;
				++tcount;
			}
		}
	}
	printf("NODES ADDED EDGECOUNT = %d maxdepth = %d\n", tcount, maxdepth);
	printf("NO OF POINTS INSETED %d\n", ++pointsin);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(800,600);
    glutCreateWindow("");
    init();
    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);
    glutKeyboardFunc(inputKey);
    glutReshapeFunc(changeSize);

    glutMainLoop();
    return 0;
}
