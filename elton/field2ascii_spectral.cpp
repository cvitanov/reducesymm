#include <iostream>
#include <fstream>
#include <iomanip>
#include "channelflow/flowfield.h"
#include "channelflow/dns.h"
#include "utilities.h"

using namespace std;
using namespace channelflow;

int main(int argc, char* argv[]) {

  string purpose("save FlowField as ascii file of gridpoint values");
  ArgList args(argc, argv, purpose);

  const string uname = args.getstr(2, "<fieldname>",      "flowfield input file");
  const string afile = args.getstr(1, "<asciifile>",      "ascii output file");

  args.check();

  FlowField u(uname);


  fieldstate xzstate  = u.xzstate();
  fieldstate ystate = u.ystate();
  cout << "xz state = " << xzstate << endl;
  cout << "y state = " << ystate << endl;


   int Mx = u.Mx();
   int My = u.My();
   int Mz = u.Mz();
  
   cout << "Mx = " << Mx << endl;
   cout << "My = " << My << endl;
   cout << "Mz = " << Mz << endl;
  
  const int Nx = u.Nx();
  const int Ny = u.Ny();
  const int Nz = u.Nz();
  const int Nd = u.Nd();

  cout << "Nx = " << Nx << endl;
  cout << "Ny = " << Ny << endl;
  cout << "Nz = " << Nz << endl;
  cout << "Nd = " << Nd << endl;
  cout << "Lx = " << u.Lx() << endl;
  cout << "Lz = " << u.Lz() << endl;
  cout << "a = " << u.a() << endl;
  cout << "b = " << u.b() << endl;

  
   int p = 16;  // digits of precision
   int w = p+7; // width of field

  ofstream os(appendSuffix(afile, ".asc").c_str());
  os << setprecision(p) << scientific;
   for  (int i=0; i<Nd  ; ++i)
    for  (int my=0; my<My  ; ++my)
      for  (int mx=0; mx<Mx  ; ++mx)
       	for (int mz=0; mz<Mz ; ++mz)
	  os << setw(w) << u.cmplx(mx,my,mz,i) << '\n';


      ofstream gs(appendSuffix(removeSuffix(afile, ".asc"), "_geom.asc").c_str());
  //  ofstream gs(appendSuffix(removeSuffix(afile, ".asc"), ".geom").c_str());
  gs << setprecision(p);
  gs.setf(ios::left);
  gs << setw(w) << u.Mx() << " %Mx\n";
  gs << setw(w) << u.My() << " %My\n";
  gs << setw(w) << u.Mz() << " %Mz\n";
  gs << setw(w) << u.Nx() << " %Nx\n";
  gs << setw(w) << u.Ny() << " %Ny\n";
  gs << setw(w) << u.Nz() << " %Nz\n";
  gs << setw(w) << u.Nd() << " %Nd\n";
  gs << setw(w) << u.Lx() << " %Lx\n";
  gs << setw(w) << u.Lz() << " %Lz\n";
  gs << setw(w) << u.a()  << " %a\n";
  gs << setw(w) << u.b()  << " %b\n";
  gs << setw(w) << u.Lx()/(2*pi) << " %lx=Lx/(2pi)\n";
  gs << setw(w) << u.Lz()/(2*pi) << " %lz=Lz/(2pi)\n";
  gs << setw(w) << 2*pi/u.Lx()   << " %alpha=2pi/Lx\n";
  gs << setw(w) << 2*pi/u.Lz()   << " %gamma=2pi/Lz\n";

   }

