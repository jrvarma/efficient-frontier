 /*
    efront.cpp for computation of efficient frontier & efficient portfolios
    Copyright (C) 2001, 2013  Prof. Jayanth R. Varma, jrvarma@iimahd.ernet.in,
    Indian Institute of Management, Ahmedabad 380 015, INDIA

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (see file COPYING); if not, write to the 
    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
    Boston, MA  02111-1307  USA
*/

/*
    This software uses the Eigen  template library for linear algebra 
    which is primarily licenced  under the MPL2
*/

/*
    The algorithm is that described in Chapter VIII of the original
    book on the subject:
    Harry M. Markowitz, 1959,
    Portfolio Selection; Efficient Diversification of Investments
    New York, John Wiley.
    The notation also closely follows that of Markowitz 
*/


#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>
#define RANGECHECK 1
#include "Eigen/Dense"

#define small 1e-6
#define  infinity 1e30

using namespace std;
using namespace Eigen;


//function prototypes
void outputs(void);
int initlz(void);
void iterate(int &j0);
void read(int argc, char *argv[]);
void graphs(void);
void graph_frontier(void);
void graph_composition(void);
void axis(const double max, const double min, const bool yaxis);
void fill_region(ostream& os, const MatrixXd& oldpath, const MatrixXd& path, 
		 const char* fill_color, const char* stroke_color, const int width);
void line(ostream& os, const double x1, const double y1, const double x2, const double y2, 
	  const char* color, const int width);
void outtextxy(ostream& os, const double x, const double y, const char* s,
	       const char* color, const int fontsize, const char* attrib=NULL);
int find_peak(const int j);
void maxmin(double& max, double& min, const double array[], const int num);
void doc_header(ostream& os);
void doc_close(ostream& os);
void text_header(ostream& os);
void text_close(ostream& os);
void svg_header(ostream& os, const bool first_graph);
void svg_close(ostream& os);

//global variables

static char helpscr[] =
"Usage Efront [options] (infile) \n"
"Valid Options Are:\n"
"-t     : Covariance/Correlation Matrix only lower TRIANGULAR half is entered\n"
"-p     : Means, covariances, correlations etc. are of Returns in PERCENT\n"
"-r     : Correlation Matrix instead of Covariance Matrix\n"
"-s     : Standard deviations instead of variances (Meaningful only if -r is used)\n"
"-n     : No SVG output. Produces plain text output without graphs\n"
"-h     : Produces http header for use if running on web server\n"
"input file contains\n"
"N             : Number of securities\n"
"MEANS         : Means of each security\n"
"COVAR         : The covariance matrix\n"
"                   OR (if -r)\n"
"CORREL        : The correlation matrix and\n"
"VARs/SIGMAs   : The variances (or if -s standard deviations\n"
"Blanks and blank lines may be used freely. Break lines freely\n"
"if no infile specified, reads from standard input\n";

const int MaxIter = 100;
int n;                     // number of securities
#define nplus1 (n+1)       // n + 1
int  num;                  // number of iterations
VectorXd Mu;            // Mu[n] vector of means
MatrixXd  C;            // C[n][n], covariance matrix
MatrixXd  M;            // M[nplus1][nplus1], C bordered with 1's
MatrixXd  N;            // N[n][n],
VectorXd  R;            // R[n],  (Zero  | 1)' where Zero is a vector of n zeros
VectorXd  S;            // S[n], (Mu | n)'
bool    *In;               // In[n] which securities are in the efficient portfolio
double lambdaE;
VectorXd Y[MaxIter];    // portfolio weights at each efficient portfolio
double  meanret[MaxIter];  // mean returns          ..
double  variance[MaxIter]; // variances             ..
double  std_dev[MaxIter];  // std deviations               .. 
double  lambdas[MaxIter];  // lambdas               .. 
bool HttpHeaderRequired;   // whether this is running on a web server 
bool NoSVG;                // no SVG output, only plain text output

int main(int argc, char *argv[])
{
  try{
    read(argc, argv); // process command line arguments and read input
    // return value is whether 
  }
  catch(const char *s){
    cout << s << endl << "Exiting ..." << endl;
    exit(1);
  }
  int j0 = initlz(); //initialisation before first iteration
  //j0 is the security to enter the efficient portfolio
  //it is passed on to iterate
  do {
    iterate(j0); // iterate until lambdaE is zero 
  }while ( lambdaE > small && num < MaxIter);
  doc_header(cout); // write content header and opening tags for document
  outputs();        // write outputs
  if(!NoSVG)        // graphs can be drawn only with SVG
    graphs();       // draw graphs
  doc_close(cout);  // write closing tags for document
  /* 
     Ideally, we should have an HTML/XML file in which the textual output and the graphs 
     are both embedded as separate elements
     Unfortunately, many SVG readers do not render HTML
     Hence if SVG is enabled, we write out the textual output as comments in a pure SVG file 
     in which only the graphs are rendered
  */
}

//initialisation before first iteration
int initlz()
{
  int j0;
  //  M =      C  1        as in Markowitz
  //           1  0
  M = MatrixXd(nplus1, nplus1);
  {
    M.topLeftCorner(n, n) = C;
    M.rightCols(1) = VectorXd::Constant(nplus1, 1.0);
    M.bottomRows(1) = MatrixXd::Constant(1, nplus1, 1.0);
    M(n,n) = 0.0;         // M[n+1,n+1] = 0
  }
  //  R =      (Zero  | 1.0)' where Zero is a vector of n zeros
  //           as in Markowitz
  R= VectorXd::Zero(nplus1);
  R(n) = 1; 
  //  S        (Mu | 0.0)'        as in Markowitz
  S = VectorXd(nplus1);
  S.head(n) = Mu;
  S(n) = 0.0;
  // find the first security to be IN (the one with highest return)
  // This is Markowitz' step 1
  double maxMu = 0;
  for (int i = 0; i < n; i++)
    if (Mu(i) > maxMu){
      j0 = i;
      maxMu = Mu(i);
    }
  for(int i =0; i< nplus1; i++)
    In[i] = false;
  In[j0] = true;
  N = MatrixXd::Zero(nplus1, nplus1);
  // As in Markowitz, N is the inverse of M~ with zero crosses 
  // replacing unit crosses. M~ is M with the rows and columns
  // corresponding to variables that are out replaced by unit
  // crosses.
  // This is Markowitz' step 2
  N(j0,n) = 1;
  N(n,j0) = 1;
  N(n,n) = -M(j0,j0);
  lambdaE = infinity;
  num = 0;
  return j0;
}

//do one iteration of the Markowitz algorithm
void iterate(int& j0)
{
  num = num + 1;
  VectorXd NR(nplus1), NS(nplus1), MNR(nplus1), MNS(nplus1), Lambda_j(n), scalar(1);
  NR = N*R;  // same as Markowitz' T
  NS = N*S;  // same as Markowitz' U
  MNR = M*NR;// same as Markowitz' V
  MNS = M*NS;// Markowitz defines W to be MNS - S
  int oldj0 = j0;
  for (int i=0; i<n; i++){
    if (In[i] && fabs(NS(i)) > small )
      // at what lambda does security i leave the efficient portfolio
      // see Markowitz' step 5
      Lambda_j(i) = -NR(i)/NS(i);
    else  if ( In[i] == 0 && fabs(MNS(i) - Mu(i)) > small)
      // at what lambda does security i enter the efficient portfolio
      // see Markowitz' step 3
      Lambda_j(i) = MNR(i)/(Mu(i)  - MNS(i));
    else
      // negative lambda is the same as never
      Lambda_j(i) = -1;
    if(Lambda_j(i) >= lambdaE) Lambda_j(i) = -1;
  }
  Lambda_j(oldj0) = -1;
  double maxLambdaj = -1;
  for (int i = 0; i<n; i++)
    if (Lambda_j(i) > maxLambdaj){
      j0 = i;
      maxLambdaj = Lambda_j(i);
    }
  if (Lambda_j(j0) < 0 )
    lambdaE = 0;
  else
    lambdaE = Lambda_j(j0);
  // Y[num] = VectorXd(n);
  Y[num] = (NR + lambdaE*NS).head(n);
  meanret[num] = Y[num].dot(Mu); 
  variance[num] = Y[num].transpose() * C *  Y[num];
  lambdas[num] = lambdaE;
  if(lambdaE > 0){
    In[j0] = ! (In[j0]);
    if (In[j0]){  // This is Markowitz' step 4
      VectorXd NMh(nplus1);
      VectorXd Mh(nplus1);
      // Mh is the j0'th column of M same as Markowitz' Cj0
      Mh = M.col(j0);
      // NMh 
      NMh= N*Mh; // same as Markowitz' B
      double b = Mh.dot(NMh); // same as Markowitz
      double c = M(j0,j0) - b;   // same as Markowitz
      for (int j=0; j<nplus1; j++)
        for (int i=0; i<nplus1; i++)
          if ( i!=j0 && j!=j0)
            N(j,i) = N(j,i) + (NMh(i)*NMh(j))/c;
      for (int j=0; j<nplus1; j++)
        if(j!=j0)
          N(j,j0) = N(j0,j) = -(NMh(j) / c);
      N(j0,j0) = 1/c;
    }else{ // i.e. not In[j0]  This is Markowitz' step 8
      for (int j=0; j<nplus1; j++)
        for (int i=0; i<nplus1; i++)
          if ( j !=j0 && i != j0)
            N(j,i) = N(j,i) - ((N(j,j0) * N(j0,i))/ N(j0,j0));
      for (int j=0; j<nplus1; j++){
	N(j,j0) = 0; N(j0,j) = 0;
      }
    } // if In[j0]
  } // if(lambdaE > 0)
}

// process command line arguments and read input
#define READ(x) ((argc==argk) ? (cin >> x) : infile >> x)
void read(int argc, char *argv[])
{
  bool Percent = false, Triangular = false, Correl = false, Sigma = false;
  HttpHeaderRequired = false;
  NoSVG = false;
  //process command line option -t -p -w
  int argk;
  for(argk = 1; argk < argc && *argv[argk]=='-'; argk++)
    switch (argv[argk][1]){
    case 't' : Triangular = true; break;
    case 'p' : Percent = true; break;
    case 'n' : NoSVG = true; break;
    case 'r' : Correl = true; break;
    case 's' : Sigma = true; break;
    case 'h' : HttpHeaderRequired = true; break;
    default  : throw helpscr; break; //any other option is an error
    }
  ifstream infile;
  if (argc > argk){
    infile.open(argv[argk]);
    if(!infile.is_open())
      throw "Cannot open input file";
  }
  if (!READ(n))               // read no of securities
    throw "Error during read";
  Mu = VectorXd(n);                       // set dimension of these vectors
  C = MatrixXd(n, n);
  In = new bool[n];
  for (int i = 0; i < n;i++)
    if(!READ(Mu(i)))
      throw "Error during read";
  for (int i = 0; i < n; i++){
    for (int j = 0; j < (Triangular ? i+1 : n); j++){
      if(!READ(C(i,j)))
	throw "Error during read";
      if(Triangular)
	C(j,i)=C(i,j);
    }
  }
  if(Correl){
    
    MatrixXd Sigmas = MatrixXd::Zero(n, n);
    for (int j = 0; j < n; j++){
      if(!READ(Sigmas(j, j)))
        throw "Error during read";
      if(!Sigma)
        Sigmas(j, j) = sqrt(Sigmas(j, j));
    }
    C = Sigmas*C*Sigmas;
  }
  if(Percent){
    Mu = 0.01*Mu;
    C = 0.0001*C;
  }
}

// produce outputs
void outputs(void)
{
  text_header(cout);  // write tags for plain text
  cout << setprecision(4);
  for(int i = 1; i <= num; i++){
    //write out the means and variances
    cout << "Iteration: " << i  << ". Mean = " << 100*meanret[i] 
	 << "%. Variance = " << 100*variance[i] << "%. Sigma = " << 100*sqrt(variance[i]) 
	 << "%" << endl;
    //write out the portfolio weights
    for(int j = 0; j < n; j++)
      cout << "Sec" << j+1 << " = " << 100*Y[i](j) << "% ";
    cout << endl ;
  }
  text_close(cout);
}

void graphs(void)
{
  if(num==1)
    return; // clearly there is nothing to graph
  // convert all means and std dev to percentages for graphing
  for (int i=1;i<=num; i++){
    meanret[i] = meanret[i]*100;
    std_dev[i] = sqrt(variance[i])*100;
  }
  svg_header(cout, true);      // SVG header for first graph
  graph_frontier();            // draw the efficient frontier
  svg_close(cout);             // close the SVG tag 
  svg_header(cout, false);     // SVG header for other graph
  graph_composition();         // draw the portfolio composition graph
  svg_close(cout);             // close the SVG tag 
}
//various SVG coordinates
#define maxx 1000
#define maxy 1000
#define ROW1 (maxy-60)
#define ROW2 (maxy-40)
#define ROW3 (maxy-20)
#define BOTTOM (maxy-80)
#define LEFT 100
#define RIGHT (maxx-LEFT)
#define TOP 80
#define TITLEROW (TOP/2)
#define CENTREX maxx/2
#define CENTREY maxy/2
#define COL1 20
#define COL2 40
#define COL3 60
#define OFFSET 20

//various SVG text attributes
#define VERTICAL " writing-mode=\"tb\""
#define RIGHTALIGN " text-anchor=\"end\""
#define CENTREALIGN " text-anchor=\"middle\""

//define a symbol for the quotation mark
#define qt "\""

//thse definitions allow us to use XAXIS and YAXIS instead of the boolean variable
//required by the function axis
#define YAXIS true
#define XAXIS false

class Scale
{
public:
  double xmin, xmax, ymax, ymin;
  //simple linear scaling
  //xmax and xmin are the two endpoints in user coordinates
  //RIGHT and LEFT are the two endpoints in SVG coordinates 
  //similarly ymax and ymin are the two endpoints in user coordinates
  //TOP and BOTTOM are the two endpoints in SVG coordinates 
  //but the direction of the positive Y axis is from TOP to BOTTOM
  inline double X(const double x) {return LEFT + (x-xmin)*(RIGHT-LEFT)/(xmax-xmin);}
  inline double Y(const double y) {return BOTTOM - (y-ymin)*(BOTTOM-TOP)/(ymax-ymin);}
};
Scale scale;

//draw the efficient frontier
void graph_frontier()
{
  //find max and min of standard deviation and draw x axis accordingly
  double max, min;
  maxmin(max, min, std_dev, num);
  double sd_range = max - min;
  axis(max, min, XAXIS);
  //find max and min of mean return and draw y axis accordingly
  maxmin(max, min, meanret, num);
  axis(max, min, YAXIS);
  //chart title
  outtextxy(cout, CENTREX, TITLEROW, "The Efficient Frontier", "black", 48, CENTREALIGN);
  //x axis title
  outtextxy(cout, CENTREX, ROW2, "Standard Deviation %", "black", 24);
  //y axis title
  outtextxy(cout, COL1, CENTREY, "Mean %",  "black", 24, VERTICAL CENTREALIGN);
  for(int i=1; i<num; i++){
    //since the efficient frontier is nonlinear, we need to plot a number of points within
    //each segment. we plan for 50 points in all and divide this up between the various
    //segments in proportion to the change in std deviation in each segment
    double KMAX = ceil(50*(std_dev[i]-std_dev[i+1])/sd_range);
    for(int k=0; k<KMAX; k++){
      //linear interpolation of portfolio weights
      VectorXd x(n);
      x = (1-k/KMAX) * Y[i] + (k/KMAX) * Y[i+1];
      //compute mean and variance
      cout << "mean1...\n";
      double mean1 = x.dot(Mu);
      cout << "...mean2\n";
      double var1 = x.transpose() * C * x;
      //same for the next point
      x = (1-(k+1)/KMAX) * Y[i] + ((k+1)/KMAX) * Y[i+1];
      cout << "mean2...\n";
      double mean2 = x.dot(Mu);
      cout << "...mean2\n";
      double var2 = x.transpose() * C * x;
      //draw a line between these points
      line(cout, 
	   scale.X(100*sqrt(var1)), 
	   scale.Y(100*mean1), 
	   scale.X(100*sqrt(var2)), 
	   scale.Y(100*mean2), 
	   "red", 
	   1);
    } // k loop
  } // i loop
}

//this class manages storage for the Matrices path and oldpath
//on construction the oldpath is a straight line from x1,y1 to x2,y2
//every time the class is incremented, path and oldpath are swapped
//the major reason for using this class is that the actual matrices are
//not swapped, only pointers to them are swapped
class PathArrays
{
private:
  MatrixXd _path;
  MatrixXd _oldpath;
  MatrixXd start;
  MatrixXd *path, *oldpath;
  int j;
public:
  PathArrays(const int n, const double x1, const double y1, const double x2, const double y2) : 
    _path(2,num), _oldpath(2,num), start(2,2) {
    start(0,0) = x1;
    start(1,0) = y1;
    start(0,1) = x2;
    start(1,1) = y2;
    j = 0;
    oldpath = &start;
    path = &_path;
  }
  void operator ++(){
    j++;
    if(j == 1){
      oldpath = path;
      path = &_oldpath;
    }else{
      MatrixXd* temp = oldpath;
      oldpath = path;
      path = temp;
    }
  }
  MatrixXd& Path(){return *path;};
  MatrixXd& Oldpath(){return *oldpath;};
};


//draw the graph for portfolio composition
void graph_composition(void)
{
  static const char* color[] ={"pink", "yellow", "cyan", "lightgreen", "grey", 
			 "red", "blue", "green", "brown"};

  double max, min;
  maxmin(max, min, std_dev, num);
  axis(max, min, XAXIS);
  axis(99, 0, YAXIS);
  //find and store the iterations  at which secuirty has highest weightage
  VectorXi peak(n);
  for(int j=0; j<n; j++){
    peak(j) = find_peak(j);
  }
  //convert portfolio weights into cumulative wieghts for graphing purposes
  for(int i=1; i<=num; i++)
    for(int j=0; j<n; j++){
      Y[i](j) = Y[i](j)*100;
      if(j > 0) Y[i](j) += Y[i](j-1);
    }
  //chart title
  outtextxy(cout, CENTREX, TITLEROW, "Composition of Efficient Portfolios", 
	    "black", 48, CENTREALIGN);
  //x axis title
  outtextxy(cout, CENTREX, ROW2, "Standard Deviation %", "black", 24);
  //y axis title
  outtextxy(cout, COL1, CENTREY, "Composition %",  "black", 24, VERTICAL CENTREALIGN);
  int NCOLORS = sizeof(color)/sizeof(char*);
  //the class PathArrays manages the Matrices path and oldpath 
  //for the first path the oldpath is a portion of the X axis
  PathArrays P(num, scale.X(std_dev[1]), BOTTOM, scale.X(std_dev[num]), BOTTOM);
  for (int j=0; j<n; j++, ++P){
    for (int i=1; i<=num; i++){
      P.Path()(0, i-1) = scale.X(std_dev[i]); 
      P.Path()(1, i-1) = scale.Y(Y[i](j));
    }
    fill_region(cout, P.Oldpath(), P.Path(), color[j % NCOLORS], "black", 1);
    //we now write the security number at the point where the security has max wieghtage
    //and therefore there is maximum space to write it
    int k = peak(j);
    char s[5];
    sprintf(s, "%d", j+1);
    double offset = (k==1) ? OFFSET : 0; 
    outtextxy(cout, 
	      scale.X(std_dev[k])-offset,
	      (j > 0) ? scale.Y(0.5*(Y[k](j)+Y[k](j-1))) : scale.Y(0.5*Y[k](j)),
	      s, "black", 24);
  }
}

//decide on the scale for an axis, draw the axis and the tick-mark labels
//max and min are the maximum and minimum values of the data
//the actual maximum and minimum points on the axis are different
//since these have to be round values
//yaxis determines whether we draw the xaxis or the yaxis
//the chosen scale is stored in the global variable scale
void axis(const double max, const double min, const bool yaxis)
{
  static double grid[]= {2.5, 1, 0.5, 0.25};
  static char s[10];

  if(yaxis)
    line(cout, LEFT, BOTTOM, LEFT, TOP, "black", 5);
  else
    line(cout, LEFT, BOTTOM, RIGHT, BOTTOM, "black", 5);
  //stretch out the given max and min a little
  double M = max + 0.01*fabs(max);
  double m = min - 0.01*fabs(min);
  //we now multiply or divide by a power of ten to get the msd (Most Significant Digit)
  //lying between 1 and 10. msd is a double and not an integer
  double msd = M - m; 
  if(msd == 0)
    msd = 0.01; 
  double Exp = 1;
  while(msd > 10){
    msd /= 10; Exp *= 10;
  }
  while(msd < 1){
    msd *= 10; Exp /= 10;
  }
  //we now decide whether the distance between two ticks 
  //depending on how large msd is 
  //for example if msd is more than 5 the step is grid[0] which is 2.5
  int type = (msd >5) ? 0 : (msd > 3) ? 1 : (msd > 1.5) ? 2 : 3;
  //the tick values at the two ends of the axis
  double maxtick = ceil(M/Exp/grid[type])*grid[type];
  double mintick = floor(m/Exp/grid[type])*grid[type];
  //store the chosen scale in the global variable scale
  if(yaxis){
    scale.ymax = maxtick*Exp;
    scale.ymin = mintick*Exp;
  }else{
    scale.xmax = maxtick*Exp;
    scale.xmin = mintick*Exp;
  }
  //we write out the tick mark labels
  for(double tick = mintick; tick <= maxtick; tick += grid[type]){
    sprintf(s, "%5.1f", tick*Exp);
    if(yaxis)
      outtextxy(cout, LEFT, scale.Y(tick*Exp), s, "black", 24, RIGHTALIGN);
    else
      //if(tick > mintick)
      outtextxy(cout, scale.X(tick*Exp), ROW1, s, "black", 24, CENTREALIGN);
  }
}

//write SVG tags for a line between the points x1,y1 and x2,y2
//the given color and width are used for the line
void line(ostream& os, const double x1, const double y1, const double x2, const double y2, 
	  const char* color, const int width) 
{
  os << "<path d=" << qt << "M " << x1 << " " << y1 << " L " << x2 << " " << y2 << qt 
     << endl << "fill=" << qt << "none" << qt << " stroke=" << qt << color << qt 
     << " stroke-width=" << qt << width << qt << "/>" << endl;
}

//write the SVG tags for constructing the filled region lying between two paths
//the region is bounded by the two paths and the lines joining the two paths
//at either end
//the filled region is constructed using the path element
//stroke-color and width govern the color and width of the boundary line
void fill_region(ostream& os, const MatrixXd& oldpath, const MatrixXd& path, 
		 const char* fill_color, const char* stroke_color, const int width) 
{
  int cols = oldpath.cols();
  os << "<path d=" << qt << "M " << oldpath(0, cols-1) << " " << oldpath(1, cols-1) << " L ";
  for(int i = cols-2; i >= 0; i--)
    os << oldpath(0, i) << " " << oldpath(1, i) << " ";
  for(int i = 0; i < path.cols(); i++)
    os << path(0, i) << " " << path(1, i)+1 << " ";
  os << "z" << qt
     << endl << "fill=" << qt << fill_color << qt << " stroke=" << qt << stroke_color << qt 
     << " stroke-width=" << qt << width << qt << "/>" << endl;
}

//find the iteration number at which the security has the highest weightage
int find_peak(const int j)
{
  double max = -infinity;
  int index;
  for(int i = 1; i <= num; i++){
    if(Y[i](j) > max){
      max = Y[i](j);
      index = i;
    }
  }
  return index;
}

//write SVG tags for text at coordinates x, y
//qt is the quotation mark
//s is the string to be written as an SVG text element
//the color and fontsize parameters are written out as attributes for the text element
//the final (optional) parameter is an SVG text atribute that is written out unchanged
void outtextxy(ostream& os, const double x, const double y, const char* s,
	       const char* color, const int fontsize, const char* attrib /*=NULL*/) 
{
  os << "<text x="<< qt << x << qt << " y=" << qt << y << qt << " " 
     << " color=" << qt << color << qt << " font-size=" << qt << fontsize << qt;
  if(attrib != NULL)
    os << attrib;
  os << ">" << endl << s << "</text>" << endl;
}

//returns maximum and minimum elements of an array
//num is number of elements in the array
void maxmin(double& max, double& min, const double array[], const int num)
{
  max = -infinity;
  min =  infinity;
  for (int i = 1; i<=num; i++){
    if (array[i] > max)
      max = array[i];
    if (array[i] < min)
      min = array[i];
  }
}


/* 
   Ideally, we should have an HTML/XML file in which the textual output and the graphs 
   are both embedded as separate elements
   Unfortunately, many SVG readers do not render HTML
   Hence we write out the textual output as comments in a pure SVG file in which only the 
   graphs are rendered
*/

//write out document header and opening tags
void doc_header(ostream& os)
{
  static char svg_type[] = "Content-type: image/svg+xml\n\n";
  static char text_type[] = "Content-type: text/plain\n\n";
  if (HttpHeaderRequired)
    os << ((NoSVG) ? text_type : svg_type);
  if(!NoSVG)
  os << "<?xml version="<< qt << "1.0" << qt << "?>\n" 
     << "<svg xmlns=" << qt << "http://www.w3.org/2000/svg" << qt << endl
     << "xmlns:xlink=" << qt << "http://www.w3.org/1999/xlink" << qt  << endl
     << "id=" << qt << "canvas" << qt << " width=" << qt << "100%" << qt
	  << " height=" << qt << "100%" << qt << endl
	  << "viewBox=" << qt << "0 0 1000 1000" << qt 
     << ">" << endl;
}

//write out end tags to close document
void doc_close(ostream& os)
{
  if(!NoSVG)
    os << endl << "</svg>" << endl;
}

//write out SVG tags for a graph occupying 50% of the screen
//the first_graph occupies the left 50% (x = 0)
//the other occupies the right 50% (x = CENTREX)
void svg_header(ostream& os, const bool first_graph)
{

  os << "<svg width=" << qt << "50%" << qt << " x=" << qt << (first_graph ? 0 : CENTREX) << qt
     << " viewBox=" << qt << "0 0 " << maxx << " " << maxy << qt 
     << " preserveAspectRatio=" << qt << "none" << qt
     << ">" << endl;
}

//write out closing SVG tags for the graph
void svg_close(ostream& os)
{
  os << endl << "</svg>" << endl;
}

//in SVG file text is written out as a comment, so we open comment
void text_header(ostream& os)
{
  if(!NoSVG)
    os << "<!--" << endl;
}

//in SVG file text is written out as a comment so we close comment
void text_close(ostream& os)
{
  if(!NoSVG)
    os << "-->" << endl;
}

