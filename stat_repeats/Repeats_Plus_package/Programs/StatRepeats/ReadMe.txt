Setup Instructions
In this archive all cpp and h files for StatRepeats are provided.

Linux

To compile StatRepeats execute the following command in the directory where you unpacked the archive:
	./StatRepeats.compile.sh for program with interface to relational database or
        ./StatRepeatsNoDB.compile.sh for program without interface to relational database.
These commands produce the StatRepeats.version.release or StatRepeatsNoDB.version.release binaries. In order for the above command to succeed you will need gcc 4.8.4 or later and libodbc (for version with relational database
interface).
In case you are using Debian or a derived distribution the typical way to install ODBC support which is needed to compile StatRepeats is by running the following commands:
	sudo apt-get install odbcinst
        sudo apt-get install unixodbc-dev
For other distributions please check your package manager documentation or build from source.

Windows

StatRepeats executable version (StatRepeats.version.release.exe) is provided in this archive in folder Windows. You will need to install the Visual C++ 2019 runtime in order to run it. It can be downloaded from Microsoft cite (https://www.microsoft.com/en-us/download/details.aspx?id=48145).
To compile from source you will need the free Visual Studio 2019 Community Edition. Open StatRepeats.sln which is included in the archive Windows and build the solution. 