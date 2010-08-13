#ifndef _TS_INPUT_H_
#define _TS_INPUT_H_

class IoData;

//------------------------------------------------------------------------------

struct TsInput {

  char *solutions;
  char *positions;
  char *levelsets;
  char *podFile;
	char *snapFile;

// Included
  char *shapederivatives;

  TsInput(IoData &);
  ~TsInput();

};

//------------------------------------------------------------------------------

#endif
