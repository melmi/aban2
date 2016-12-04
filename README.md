# ABAN II

A hydrodynamic free-surface solver.

## Prerequisites

Aban ii uses [JsonCpp](http://jsoncpp.sourceforge.net/‎) and [LASPack]
(http://www.mgnet.org/mgnet/Codes/laspack/html/laspack.html) libraries. They
should be installed before compiling Aban ii.  You can install `JsonCpp` in
ubuntu using the following command:

    sudo apt-get install libjsoncpp-dev

Installation of `LASPack` is described in itsown documentation. We just mention
it for completeness.

1. Download it from [this link](http://www.netlib.org/linalg/laspack.tgz).

2. Extract it:

    tar xvf laspack.tgz

3. Install it:

    bash install local
