Using the Windows binary
========================

The windows binary is `efront.exe`. This is the only file that you need if you plan to use efront under Windows. All other files in this distribution may be deleted. 

*The Windows binary included here was built using an earlier version of the source code that used Masakatsu Ito's matrix class* `MET`. *The old version was compiled with gcc 2.95 under Cygwin in Windows 98. It was built to link with the MINGW library to produce a binary that runs as a native Windows application without any other dll's. The current source code uses the matrix library Eigen partly because MET does not compile easily with new versions of gcc. There is no change in functionality between the two versions of the code.* 

`efront` takes only one argument which is the input file. If this argument is omitted `efront` reads input from standard input so that you can pipe input into `efront` from some other program.

By default, `efront` produces a graph of the efficient frontier and of the composition of the efficient portfolios along the frontier. These graphs are produced in `SVG` (Scalable Vector Graphics) format, which is the standard laid down by the W3C. All modern web browsers and most image viewers can display an SVG file.

Efront can also produce a plain text output giving the efficient portfolios at each turning point in the efficient frontier. Between these points the efficient portfolios are linear combinations of the portfolios at the two ends. To get plain text output, you must use the -n option.

    Usage Efront [options] (infile)
    Valid Options Are:
    -t     : Covariance/Correlation Matrix only lower TRIANGULAR half is entered
    -p     : Means, covariances, correlations etc. are of Returns in PERCENT
    -r     : Correlation Matrix instead of Covariance Matrix
    -s     : Standard deviations instead of variances (Meaningful only if -r is used)
    -n     : No SVG output. Produces plain text output without graphs
    -h     : Produces http header for use if running on web server
    input file contains
    N             : Number of securities
    MEANS         : Means of each security
    COVAR         : The covariance matrix
                       OR (if -r)
    CORREL        : The correlation matrix and
    VARs/SIGMAs   : The variances (or if -s standard deviations
    Blanks and blank lines may be used freely. Break lines freely

Please see the two examples files in the samples directory.

`samples/markowitz.dat` is the example in the original book on the subject: **Harry M. Markowitz (1959),** *Portfolio Selection; Efficient Diversification of Investments,* New York, John Wiley.

`efront -n samples/markowitz.dat` will give a text output (no SVG) while `efront samples/markowitz.dat` will generate a SVG output.

`samples/example.dat` is another example in which only lower triangular matrix is given and all data are in percent. The `-t` and `-p` switches must be used with this file.

`efront -t -p -n samples/example.dat` will give a text output (no SVG) while `efront -t -p samples/example.dat` will generate SVG output.


Building efront
===============

Efront can be built under Linux (or a Linux-like environment under Windows) using `gcc`. As explained earlier, the Windows binary in this repository is based on an older version of the code and a different matrix library and was built using `gcc` 2.95 under Cygwin in Windows 98. 

I believe that efront does not use any Unix specific or Windows-specific features and should build smoothly in any Operating System. I would imagine that if you have a modern C++ compiler, you should be able to build efront in any environment without difficulty. It compiles with the `-pedantic` switch set in `gcc` implying that it conforms strictly to ANSI C++ and does not use any gcc-specific features. It should therefore compile with any other C++ compiler that conforms to ANSI C++. It does however need a modern C++ compiler as it uses C++ templates and namespaces extensively. It also uses the Standard Template Library (STL).

If you are in Linux or a Unix clone that supports `make`, you can use the `Makefile`: run `make` to build `efront` and run `make clean` to remove all the unwanted object files.  (The Windows binary in this distribution was created using `make MINGW=1` to compile in the `MINGW` libraries while running `gcc` under `cygwin`.  But I see no reason to use this switch except to produce a Windows binary for distribution.)

If you want to build under MS Visual C++ or some other development environment which does not support `make`, just compile `efront.cpp`. Of course, you need a modern compiler that supports templates and STL.

COPYRIGHT
---------

The program `efront` is copyrighted and distributed under GNU GPL. Copyright (C) 2001 Prof. Jayanth R. Varma, jrvarma@iima.ac.in, Indian Institute of Management, Ahmedabad 380 015, INDIA

`efront` uses `Eigen`, a lightweight C++ template library for linear algebra which is subject to the terms of the Mozilla Public License v. 2.0..  Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr> Copyright (C) 2007-2011 Benoit Jacob <jacob.benoit.1@gmail.com>


