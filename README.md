<img src=https://github.com/motian12ps/FaRSA/blob/master/farsa.jpg?raw=true width=135/> Fast Reduced Space Algorithm
=====
---
[Documentation]()|
[License]()|
[Installation]()|


##Installation Guide:

Download source folder by

	git clone https://github.com/motian12ps/FaRSA

or 

	wget url.FaRSA.zip
	

#### Build on OSX 

On OSX, make sure gcc has been installed. If not, to install gcc compiler on Mac OS X, you need to download and install "Command Line Tools for Xcode", which is available in Apple’s developer page.

Then jump to FaRSA source directory, run the below command in terminal:

	make

If compiling failed due to permission issue, please compile by 

	sudo make

To test whether installation is successful, run the below command:

	farsa -t




## Documentation:

Parameters can be modified in a profile. A demo of such profile has been offered named "FaRSA1.profile".

With customerized profile, run the below command line to call FaRSA for solving aimed problem:

	farsa -p path_to_profile

e.g.

	farsa -p FaRSA.profile









