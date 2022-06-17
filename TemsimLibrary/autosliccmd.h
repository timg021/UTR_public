//@@@@@@ start TEG code
#ifndef AUTOSLICCMD_H   // only include this file if its not already
#define AUTOSLICCMD_H   // remember that this has been included

#include <string>
#include <vector>
#include "XA_file.h"

#define TEG_MULTITHREADED 1 //if defined, multithreaded execution is used

// TEG class for keeping a count across multiple threads (e.g. for counting the number of active threads)
class Counter_Obj
{
	static int counter;
	static bool isUpdated;
public:
	Counter_Obj() { counter++; }
	~Counter_Obj() { counter--; }
	int GetCount() { return counter;  }
	void SetUpdated(bool status) { isUpdated = status; }
	bool GetUpdated() { return isUpdated; }
	void SetTerminate() { counter = -100000; }
};

int autosliccmd(std::vector<std::string> params, std::vector<xar::Pair> defocus, xar::Pair astigm, xar::Pair xyshifts, std::vector<std::string> fileout);
int autosliccmd1(std::vector<std::string> params);

#endif //AUTOSLICCMD_H
//@@@@@ end TEG code

