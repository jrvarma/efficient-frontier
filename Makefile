##################################################
#    Makefile for efront for Computing Efficient Frontier/Portfolios
#    Copyright (C) 2001  Prof. Jayanth R. Varma, jrvarma@iimahd.ernet.in,
#    Indian Institute of Management, Ahmedabad 380 015, INDIA
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program (see file COPYING); if not, write to the 
#    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
#    Boston, MA  02111-1307  USA
##################################################

##################################################
#To compile with MINGW in a CYGWIN environment invoke this makefile with
#make MINGW=1
#Then the make uses the MINGW include files and the MINGW libraries
#defined below. You must set these paths correctly.
#To compile in any other environment, these flags are not necessary.
#The MINGW headers do not compile under -pedantic so we set this flag
#only under non MINGW 
ifdef MINGW
MINGWINCLUDE = /usr/local/mingw/include/g++-3
MINGLIB1 = /usr/local/mingw/lib
MINGLIB2 = /usr/local/mingw/lib/gcc-lib/mingw32/2.95.3-5
MINGWFLAG = -mno-cygwin -I$(MINGWINCLUDE) -L$(MINGLIB1) -L$(MINGLIB2)
else
MINGWFLAG =
CPPFLAGS = -pedantic
# -Wno-ignored-attributes -Wno-deprecated-declarations
endif
##################################################

CXX = g++

%.o : %.cpp
	$(CXX) -g -c  $(CPPFLAGS) $(CXXFLAGS)  $(MINGWFLAG) $< -o $@

name = efront

objects = ${name}.o 

binary = ${name}

library =  

$(binary) : $(objects)
	$(CXX) -g $(MINGWFLAG)  $(objects) -o $(binary) $(library)

${name}.o : ${name}.cpp 

clean: 
	rm -f $(objects)
