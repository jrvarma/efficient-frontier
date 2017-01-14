Building efront
===============

Run `make` to build `efront` and run `make clean` to remove all the unwanted object files. Alternatively, simply compile `efront.cpp`.

Using efront
============

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

COPYRIGHT
---------

The program `efront` is copyrighted and distributed under GNU GPL. Copyright (C) 2001 Prof. Jayanth R. Varma, jrvarma@iima.ac.in, Indian Institute of Management, Ahmedabad 380 015, INDIA

`efront` uses `Eigen`, a lightweight C++ template library for linear algebra which is subject to the terms of the Mozilla Public License v. 2.0..  Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr> Copyright (C) 2007-2011 Benoit Jacob <jacob.benoit.1@gmail.com>


