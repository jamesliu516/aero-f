#ifndef _MATCH_NODE_H_
#define _MATCH_NODE_H_

class BinFileHandler;

//------------------------------------------------------------------------------

class MatchNodeSet {

  int numNodes;

  int (*index)[3];
  double (*gap)[3];

public:

  MatchNodeSet(int);
  ~MatchNodeSet();

  void read(BinFileHandler &, int, int (*)[2]);

  template<class NodeMap>
  void renumberNodes(NodeMap &);

  void exportInfo(int, int (*)[3]);
  void setBufferPosition(int, int);
  void getDisplacement(int, double, double, double, bool *, double (*)[2][3], double (*)[3], 
		       double (*)[3], double (*)[3], double (*)[3], double *);
  double getTemperature(int, double, double, bool*, double*, double*);
  template<int dim>
  void send(double, double (*)[dim], double (*)[dim]);
  double (*getGap(int, int *))[3];

  int size() const { return numNodes; }
  int subMatchNode(int i) { return(index[i][0]); } //HB: returns the subdomain node of the ith matched node 	
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <MatchNode.C>
#endif

#endif
