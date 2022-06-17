//Header XA_file.h
//
//	HEADER FILE TITLE:
//
//		Common file facilities for XArray and related classes
//
/*!
	\file		XA_file.h
	\brief		Common file facilities for XArray and related classes
	\par		Description:
		This header contains the basic auxiliary file facilities required by most XArray and related classes.
*/

#if !defined XA_FILE_H
#define XA_FILE_H

#define WINDOWS_OS_TEG 1 // switches on the code using low-level Windows file I/O functions

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_ini.h"
#include <list>
#include <fcntl.h>
#include <sys/stat.h>
#include <filesystem>
#include <cstring>

#if defined WINDOWS_OS_TEG
#include <io.h>
#endif


#undef _WIN32
//#ifdef _WIN32
//#include <windows.h>
//#endif

namespace xar
{

	using std::list;
//---------------------------------------------------------------------------
//	MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
//---------------------------------------------------------------------------
//	CONSTANT DEFINITIONS
//
//---------------------------------------------------------------------------
//	ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	STRUCTURE DEFINITIONS
//
	struct Pair{ double a; double b; }; // pair of double values
//
//---------------------------------------------------------------------------
//	IN-LINE FUNCTION DEFINITIONS
//
//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//Class FilePtr
//
//	FilePtr class for support of the 'resource acquisition is initialization' technique
//
/*!
	\brief		Supports the 'resource acquisition is initialization' technique for files
	\par		Description:
		This class allows one to avoid resource leakage due to a failure to close a file, i.e. when
		an exception is thrown. The file will be automaitcally closed when the FilePtr variable
		goes out of scope. FilePtr is written after B.Stroustrup, The C++ programming language,
		2nd ed., Sect.9.4.
	\par		Example:
\verbatim
short a = 5;
FilePtr pFile("afile.txt", "wt"); // Open the text file 'afile.txt' for writing
fprintf("a=%d", a);
// ...
// The file will be automatically closed when 'pFile' variable goes out of scope
\endverbatim
*/
	class FilePtr
	{
	// Enumerators
	// Structures
	// Constructors
	public: 
		//! Default constructor
		FilePtr() { m_pFile = 0; }
		//! Constructor opening a file  and storing its FILE pointer in the member pointer
		FilePtr(const char* filename, const char* mode) { m_pFile = 0; Open(filename, mode); }
		//! Constructor accepting FILE pointer structure for an opened file
		FilePtr(FILE* p1) { m_pFile = p1; }
		//! Copy constructor
		FilePtr(const FilePtr& other) { m_pFile = other.m_pFile; } // closing file twice is OK
	public:
		//! Destructor, closes the file
		~FilePtr() { Close(); }
	// Operators
	public:
		//! Performs implicit conversion into a FILE pointer
		operator FILE*() { return m_pFile; }
		//! Assignment operator (declared private to prohibit copying)
		FilePtr& operator=(const FilePtr& other) { m_pFile = other.m_pFile; return *this; } // closing file twice is OK
	// Attributes
	public:
	// Operations
	public:
		//! Opens a file and stores its FILE pointer in the member pointer
		void Open(const char* filename, const char* mode) 
		{ 
			Close();
			m_pFile = fopen(filename, mode); 
			if (!m_pFile)
			{
				string aaa("runtime_error in FilePtr::Open(cannot open file ");
				aaa += filename;
				aaa += ")";
				throw std::runtime_error(aaa.c_str());
			}
		}
		//! Closes the file and resets the member pointer
		void Close() { if (m_pFile) fclose(m_pFile); m_pFile = 0; }
		//! Surrenders the owneship of the file pointer
		void Detach() { m_pFile = 0; }
	private:
		// Member variables
		//! Member pointer to the controlled file 
		FILE* m_pFile;
		// Member functions
	};

	#if defined WINDOWS_OS_TEG
	class LLFilePtr
	// Variant of FilePtr using low-level I/O routines	
	{
	public: 
		LLFilePtr() { fd = -1; }
		LLFilePtr(const char* filename, bool b4Writing) { Open(filename, b4Writing); }
		LLFilePtr(int fd1) { fd = fd1; }
		LLFilePtr(const LLFilePtr& other) { fd = other.fd; } // closing file twice is OK
		LLFilePtr& operator=(const LLFilePtr& other) { fd = other.fd; return *this; } // closing file twice is OK
		~LLFilePtr() { Close(); }

		void Open(const char* filename, bool b4Writing) 
		{ 
			if (b4Writing)
				fd = _open(filename, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, _S_IREAD | _S_IWRITE);
			else
				fd = _open(filename, _O_BINARY | _O_RDONLY, _S_IREAD);
			if (fd == -1) 
				throw std::runtime_error("runtime error in LLFilePtr::Open() (cannot open file)");
		}
		void Close() { if (fd != -1) _close(fd); fd = -1; }

		operator int() { return fd; }

	private:
		int fd; // file descriptor
	};
	#endif

	#ifdef _WIN32
		class Win32FilePtr
		// Variant of FilePtr using Win32 routines directly
		{
		public: 
			Win32FilePtr() { p = INVALID_HANDLE_VALUE; }
			Win32FilePtr(const char* filename, bool b4Writing) { Open(filename, b4Writing); }
			Win32FilePtr(HANDLE p1) { p = p1; }
			Win32FilePtr(const Win32FilePtr& other) { p = other.p; } // closing file twice is OK
			Win32FilePtr& operator=(const Win32FilePtr& other) { p = other.p; return *this; } // closing file twice is OK
			~Win32FilePtr() { Close(); }

			void Open(const char* filename, bool b4Writing) 
			{ 
				if (b4Writing)
					p = ::CreateFile(filename, GENERIC_WRITE, 0, NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
				else
					p = ::CreateFile(filename, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
				if (p == INVALID_HANDLE_VALUE) 
					throw std::runtime_error("runtime error in Win32FilePtr::Open() (cannot open file)");
			}
			void Close() { if (p != INVALID_HANDLE_VALUE) CloseHandle(p); p = INVALID_HANDLE_VALUE; }

			operator HANDLE() { return p; }

		private:
			HANDLE p;
		};
	#endif

//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//---------------------------------------------------------------------------
//
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//

	string GetFileExtension(const string& filename, bool Convert2Upper = true);
	string GetFilenameFromPath(const string& strPath, bool Convert2Upper = true);
	string GetPathFromFilename(const string& strFilename, bool Convert2Upper = true);
	bool DoesFileExist(const char* filename);
	string SortFiles(const string& infiles);
	index_t NumberOfFiles(const string& infiles);
	XA_API index_t IOFileLists(const string& infiles, const string& outfile, list<string>& listInputFiles, list<string>& listOutputFiles, int* piNoProblem);
	XA_API index_t IOFileSequence(const string& infiles, const string& outbase, index_t nOutFiles, list<string>& listInputFiles, list<string>& listOutputFiles, bool bIndexFormatInOut, int* piNoProblem);
	XA_API index_t IOFileSet(const string& infiles, const string& outbase, vector<index_t> vIndexes, list<string>& listInputFiles, list<string>& listOutputFiles, bool bIndexFormatInOut, int* piNoProblem);
	XA_API index_t IOFileString2List(const string& infiles, list<string>& listInputFiles, int* piNoProblem);
	index_t FindRefractiveIndices(const string strMaterialLegendFilename, const vector<float> vEnergies, vector<vector<float> >& vOutDelta, vector<vector<float> >& vOutBeta);
	index_t ReadSpectrumFile(const string strSpectrumFilename, vector<float>& vOutEnergies, vector<float>& vOutCounts);
	void ReadDefocusParamsFile(string difile, vector<Pair>& v2angles, vector<vector<Pair> >& vvdefocus, bool bVerboseOutput = true);
	void ReadRelionDefocusParamsFile(string difile, vector<Pair>& v2angles, vector<vector<Pair> >& vvdefocus, vector<Pair>& vastigm, vector<Pair>& v2shifts, bool bVerboseOutput = true);
	void FileNames(index_t nangles, index_t ndefocus, string filenamebase, vector<string>& output);
	void FileNames2(vector<index_t> vndefocus, string filenamebase, vector<string>& output);

} // namespace xar closed


#endif	// XA_FILE_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
