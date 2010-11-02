#include <MatchNode.h>

#include <BinFileHandler.h>

//------------------------------------------------------------------------------

MatchNodeSet::MatchNodeSet(int value)
{

  numNodes = value;

  // index[][0] -> local number in sub
  // index[][1] -> global match number
  // index[][2] -> position in the MPI buffer

  index = new int[value][3];
  gap = new double[value][3];

}

//------------------------------------------------------------------------------

MatchNodeSet::~MatchNodeSet()
{

  if (index) delete [] index;
  if (gap) delete [] gap;

}

//------------------------------------------------------------------------------

void MatchNodeSet::read(BinFileHandler &file, int numRanges, int (*ranges)[2]) 
{

  // read in number of nodes in cluster (not used)
  int numClusNodes;
  file.read(&numClusNodes, 1);

  // read in the offset for the first node
  BinFileHandler::OffType start;
  file.read(&start, 1);

  int count = 0;

  // read in ranges
  for (int iRange = 0; iRange < numRanges; ++iRange)  {
    // compute number of nodes in range
    int nNodes = ranges[iRange][1] - ranges[iRange][0] + 1;

    // seek to correct position in file
    file.seek(start + ranges[iRange][0] * (2 * sizeof(int) + 3 * sizeof(double)));
    
    for (int i = 0; i < nNodes; i++) {
      // read in cluster node number and global match node number
      file.read(index[count], 2);   
      // read in node position  
      file.read(gap[count], 3);
      count++;
    }

  }

  if (count != numNodes) {
    fprintf(stderr, "*** Error: wrong number of nodes read (%d instead of %d)\n",
	    count, numNodes);
    exit(1);
  }

}

//------------------------------------------------------------------------------

void MatchNodeSet::autoInit(int nNodes) 
{
  numNodes = nNodes;
  index = new int[numNodes][3];
  gap = new double[numNodes][3];

  for(int i=0; i<numNodes; i++) {
    index[i][0] = index[i][1] = i;
    gap[i][0] = gap[i][1] = gap[i][2] = 0.0;
  }
}

//------------------------------------------------------------------------------

MatchNodeSet::MatchNodeSet(const char *name) {
  FILE *fp = fopen(name, "r");

  if (!fp)  {
    fprintf(stderr, "*** ERROR: match file \'%s\' file does not exist\n", name);
    exit(-1);
  }

  char line[MAXLINE];
  fgets(line, MAXLINE, fp);
  fgets(line, MAXLINE, fp);
  sscanf(line, "%*s %d", &numNodes);

  index = new int[numNodes][3];
  gap = new double[numNodes][3];

  // read match points
  int i;
  for (i = 0; i < numNodes; i++) {
    fgets(line, MAXLINE, fp);
    sscanf(line, "%d", index[i]);
    index[i][0]--;
    index[i][1] = index[i][0];
  }

  // read gap vectors
  i = 0;
  while (fgets(line, MAXLINE, fp) != 0) {
    sscanf(line, "%lf %lf %lf", gap[i], gap[i]+1, gap[i]+2);
    i++;
  }

  // check that # of gap vectors match number of matched points
  if (i == 0) {
    fprintf(stderr, " *** WARNING: setting all gap vectors to zero\n");
    for (i = 0; i < numNodes; i++) {
      gap[i][0] = 0.0;
      gap[i][1] = 0.0;
      gap[i][2] = 0.0;
    }
  }
  else if (i != numNodes)  {
    fprintf(stderr, " *** ERROR: incorrect number of gap vectors\n");
    exit(-1);
  }

  fclose(fp);
}

//------------------------------------------------------------------------------

void MatchNodeSet::exportInfo(int iSub, int (*list)[3])
{

  for (int i=0; i<numNodes; ++i) {
    list[i][0] = i;
    list[i][1] = index[i][1];
    list[i][2] = iSub;
  }

}

//------------------------------------------------------------------------------

void MatchNodeSet::setBufferPosition(int i, int pos)
{

  index[i][2] = pos;

}

//------------------------------------------------------------------------------

void MatchNodeSet::getDisplacement(int algNum, double dt, double lscale, double uscale, 
				   bool *flag, double (*disp)[2][3], double (*x0)[3], 
				   double (*x)[3], double (*xdot)[3], double (*dx)[3],
				   double *norms)
{

  norms[0] = 0.0;
  norms[1] = 0.0;
  for (int i=0; i<numNodes; ++i) {
    for (int k=0; k<3; ++k) {
      if (flag==0||flag[ index[i][0] ]) {
	norms[0] += disp[ index[i][2] ][0][k] * disp[ index[i][2] ][0][k];
	norms[1] += disp[ index[i][2] ][1][k] * disp[ index[i][2] ][1][k];
      }

      double dx0 = x[ index[i][0] ][k] - x0[ index[i][0] ][k];

      dx[ index[i][0] ][k] = lscale * disp[ index[i][2] ][0][k] - dx0;
      xdot[ index[i][0] ][k] = uscale * disp[ index[i][2] ][1][k];

      if (algNum == 6)
	dx[ index[i][0] ][k] += 0.5 * dt * xdot[ index[i][0] ][k];
      else if (algNum == 7)
	dx[ index[i][0] ][k] += dt * xdot[ index[i][0] ][k];
    }
  }
}

//------------------------------------------------------------------------------

double MatchNodeSet::getTemperature(int algNum, double dt, double scale, 
				    bool* flag, double* buffer, double* temp)
{

  double norm = 0.0;

  for (int i=0; i<numNodes; ++i) {
    if (flag[ index[i][0] ])
      norm += buffer[ index[i][2] ] * buffer[ index[i][2] ];
    temp[ index[i][0] ] = scale * buffer[ index[i][2] ];
  }

  return norm;

}

//------------------------------------------------------------------------------

double (*MatchNodeSet::getGap(int size, int *list))[3]
{

  if (size != numNodes) {
    fprintf(stderr, "*** Error: wrong number of interface nodes\n");
    exit(1);
  }

  for (int i=0; i<size; ++i) {
    if (list[i] != index[i][0]) {
      fprintf(stderr, "*** Error: inconsistent match data %d vs %d\n", 
	      list[i], index[i][0]);
      exit(1);
    }
  }

  return gap;

}

//------------------------------------------------------------------------------
