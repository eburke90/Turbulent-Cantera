/**
 * @file ct2ctml.cpp
 * Driver for the system call to the python executable that converts
 * cti files to ctml files (see \ref inputfiles).
 */
// Copyright 2001-2005  California Institute of Technology

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"
#include "../../ext/libexecstream/exec-stream.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <ctime>

using namespace Cantera;
using namespace std;

namespace ctml
{

//! return the full path to the Python interpreter.
/*!
 * Use the environment variable PYTHON_CMD if it is set. If not, return
 * the string 'python'.
 *
 * Note, there are hidden problems here that really direct us to use
 * a full pathname for the location of python. Basically the system
 * call will use the shell /bin/sh, in order to launch python.
 * This default shell may not be the shell that the user is employing.
 * Therefore, the default path to python may be different during
 * a system call than during the default user shell environment.
 * This is quite a headache. The answer is to always set the
 * PYTHON_CMD environmental variable in the user environment to
 * an absolute path to locate the python executable. Then this
 * issue goes away.
 */
static string pypath()
{
    string s = "python";
    const char* py = getenv("PYTHON_CMD");

    if (py) {
        string sp = stripws(string(py));
        if (sp.size() > 0) {
            s = sp;
        }
    }
    return s;
}

// Convert a cti file into a ctml file
/*
 *
 *  @param   file    Pointer to the file
 *  @param   debug   Turn on debug printing
 *
 *  @ingroup inputfiles
 */
void ct2ctml(const char* file, const int debug)
{
#ifdef HAS_NO_PYTHON
    /*
     *  Section to bomb out if python is not
     *  present in the computation environment.
     */
    string ppath = file;
    throw CanteraError("ct2ctml",
                       "python cti to ctml conversion requested for file, " + ppath +
                       ", but not available in this computational environment");
#endif

    string python_output;
    int python_exit_code;
    try {
        exec_stream_t python;
        python.set_wait_timeout(exec_stream_t::s_child, 10000);
        python.start(pypath(), "-i");
        stringstream output_stream;
        python.in() <<
            "if True:\n" << // Use this so that the rest is a single block
            "    import sys\n" <<
            "    sys.stderr = sys.stdout\n" <<
            "    import ctml_writer\n" <<
            "    ctml_writer.convert(r'" << file << "')\n";
        python.close_in();
        std::string line;
        while (std::getline(python.out(), line).good()) {
            output_stream << line << std::endl;;
        }
        python.close();
        python_exit_code = python.exit_code();
        python_output = stripws(output_stream.str());
    } catch (std::exception& err) {
        // Report failure to execute the python
        stringstream message;
        message << "Error executing python while converting input file:\n";
        message << "Python command was: '" << pypath() << "'\n";
        message << err.what() << std::endl;
        throw CanteraError("ct2ctml", message.str());
    }

    if (python_exit_code != 0) {
        // Report a failure in the conversion process
        stringstream message;
        message << "Error converting input file \"" << file << "\" to CTML.\n";
        message << "Python command was: '" << pypath() << "'\n";
        if (python_output.size() > 0) {
            message << "-------------- start of converter log --------------\n";
            message << python_output << std::endl;
            message << "--------------- end of converter log ---------------";
        }
        throw CanteraError("ct2ctml", message.str());
    }

    if (python_output.size() > 0) {
        // Warn if there was any output from the conversion process
        stringstream message;
        message << "Warning: Unexpected output from CTI converter\n";
        message << "-------------- start of converter log --------------\n";
        message << python_output << std::endl;
        message << "--------------- end of converter log ---------------\n";
        writelog(message.str());
    }
}


// Read an ctml file from a file and fill up an XML tree
/*
 *  This is the main routine that reads a ctml file and puts it into
 *  an XML_Node tree
 *
 *  @param node    Root of the tree
 *  @param file    Name of the file
 *  @param debug   Turn on debugging printing
 */
void get_CTML_Tree(Cantera::XML_Node* rootPtr, const std::string file, const int debug)
{
    std::string ff, ext = "";

    // find the input file on the Cantera search path
    std::string inname = findInputFile(file);
    if (debug > 0) {
        writelog("Found file: "+inname+"\n");
    }

    if (inname == "") {
        throw CanteraError("get_CTML_Tree", "file "+file+" not found");
    }

    /*
     * Check whether or not the file is XML. If not, it will be first
     * processed with the preprocessor.
     */
    std::string::size_type idot = inname.rfind('.');
    if (idot != string::npos) {
        ext = inname.substr(idot, inname.size());
    }
    if (ext != ".xml" && ext != ".ctml") {
        try {
            ctml::ct2ctml(inname.c_str(), debug);
        } catch (std::exception& err) {
            writelog("get_CTML_Tree: caught an exception:\n");
            writelog(err.what());
        }
        string ffull = inname.substr(0,idot) + ".xml";
        ff = "./" + getBaseName(ffull) + ".xml";
        if (debug > 0) {
            writelogf("ffull name = %s\n", ffull.c_str());
            writelogf("ff name = %s\n", ff.c_str());
        }
    } else {
        ff = inname;
    }
    if (debug > 0) {
        writelog("Attempting to parse xml file " + ff + "\n");
    }
    ifstream fin(ff.c_str());
    if (!fin) {
        throw
        CanteraError("get_CTML_Tree",
                     "XML file " + ff + " not found");
    }
    rootPtr->build(fin);
    fin.close();
}
}
