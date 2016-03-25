
/****************************************************************************** 
 * 
 *  file:  OrHandler.h
 * 
 *  Copyright (c) 2003, Michael E. Smoot .
 *  Copyright (c) 2004, Michael E. Smoot, Daniel Aarno.
 *  All rights reverved.
 * 
 *  See the file COPYING in the top directory of this distribution for
 *  more information.
 *  
 *  THE SOFTWARE IS PROVIDED _AS IS_, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 *  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
 *  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 *  DEALINGS IN THE SOFTWARE.  
 *  
 *****************************************************************************/ 

#ifndef TCLAP_XORHANDLER_H
#define TCLAP_XORHANDLER_H

#include <tclap/Arg.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>


//ToDo finish. constuctor is set to protected at the moment


class LogicHandler {

};

namespace TCLAP {

class LogicHandler
{
	protected:

		/**
		 * The list of of lists of Arg's to be or'd together.
		 */
		std::vector< std::vector<Arg*> > _argList;


		/**
		 * Constructor.  Does nothing.
		 */
		LogicHandler( ) : _argList(std::vector< std::vector<Arg*> >()) {}

public:
		/**
		 * Add a list of Arg*'s that undergo a logical operation.
		 * \param args - list of Arg*.
		 */
		void add( std::vector<Arg*>& args );

		/**
		 * Checks whether the specified Arg is in one of the arg lists and
		 * if it does match one, returns the size of the xor list that the
		 * Arg matched.  If the Arg matches, then it also sets the rest of
		 * the Arg's in the list. You shouldn't use this.
		 * \param a - The Arg to be checked.
		 */
		virtual int check( const Arg* a ) = 0;

		/**
		 * Returns the short usage of the logical operation.
		 */
		virtual std::string shortUsage() = 0;

		/**
		 * Prints the long usage of the logical operation
		 * \param os - Stream to print to.
		 */
		virtual void printLongUsage(std::ostream& os) = 0;

		/**
		 * Simply checks whether the Arg is contained in one of the arg
		 * lists.
		 * \param a - The Arg to be checked.
		 */
		bool contains( const Arg* a );

		std::vector< std::vector<Arg*> >& getArgList();
};

inline std::vector< std::vector<Arg*> >& LogicHandler::getArgList()
{
	return _argList;
}

inline bool LogicHandler::contains( const Arg* a )
{
	for ( int i = 0; static_cast<unsigned int>(i) < _orList.size(); i++ )
		for ( ArgVectorIterator it = _argList[i].begin();
			  it != _argList[i].end();
			  it++ )
			if ( a == (*it) )
				return true;

	return false;
}

inline void LogicHandler::add( std::vector<Arg*>& ors )
{
	_argList.push_back( ors );
}

/**
 * This class handles lists of Arg's that are to be XOR'd on the command
 * line.  This is used by CmdLine and you shouldn't ever use it.
 */
class OrHandler : public LogicHandler
{

	public:

		/**
		 * Constructor.  Does nothing.
		 */
		OrHandler( ) : _orList(std::vector< std::vector<Arg*> >()) {}

		/**
		 * Checks whether the specified Arg is in one of the xor lists and
		 * if it does match one, returns the size of the xor list that the
		 * Arg matched.  If the Arg matches, then it also sets the rest of
		 * the Arg's in the list. You shouldn't use this.  
		 * \param a - The Arg to be checked.
		 */
		virtual int check( const Arg* a );

		/**
		 * Returns the OR specific short usage.
		 */
		virtual std::string shortUsage();

		/**
		 * Prints the OR specific long usage.
		 * \param os - Stream to print to.
		 */
		virtual void printLongUsage(std::ostream& os);

};


//////////////////////////////////////////////////////////////////////
//BEGIN XOR.cpp
//////////////////////////////////////////////////////////////////////


inline int OrHandler::check( const Arg* a )
{
	// iterate over each XOR list
	for ( int i = 0; static_cast<unsigned int>(i) < _orList.size(); i++ )
	{
		// if the XOR list contains the arg..
		ArgVectorIterator ait = std::find( _orList[i].begin(), 
		                                   _orList[i].end(), a );
		if ( ait != _orList[i].end() )
		{
			// first check to see if a mutually exclusive switch
			// has not already been set
			for ( ArgVectorIterator it = _orList[i].begin(); 
				  it != _orList[i].end(); 
				  it++ )
				if ( a != (*it) && (*it)->isSet() )
					throw(CmdLineParseException(
					      "Mutually exclusive argument already set!",
					      (*it)->toString()));

			// go through and set each arg that is not a
			for ( ArgVectorIterator it = _orList[i].begin(); 
				  it != _orList[i].end(); 
				  it++ )
				if ( a != (*it) )
					(*it)->xorSet();

			// return the number of required args that have now been set
			if ( (*ait)->allowMore() )
				return 0;
			else
				return static_cast<int>(_orList[i].size());
		}
	}

	if ( a->isRequired() )
		return 1;
	else
		return 0;
}




//////////////////////////////////////////////////////////////////////
//END XOR.cpp
//////////////////////////////////////////////////////////////////////

} //namespace TCLAP

#endif 
