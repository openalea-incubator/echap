
How to use/contribute to this documentation
###########################################


  * This documentation is generated by a program ("Sphinx"), from text files (*.rst)  shared among Echap members
  * rst files are located on an inria subversion server (https://scm.gforge.inria.fr/svn/openaleapkg/trunk/echap) , and can be shared among participants using a subversion SVN client (eg Tortoise SVN on windows)
  * SVN allows you to get a copy of rst files on your computer ("SVN Update"), that you can edit and then send the modifications back to the server ("SVN commit")

Here is a litle memo of the different steps that you should go through to become  echap contributor. For further assistance, don't hesitate to contact a person that participated to a modelling-sprint !

Installation
============

To contribute one will thus have to install: 
  * Python/Openalea, using the OpenAlea installer
  * Sphinx, eg by typing in a command line : easy_install sphinx (easy_install tool is shiped with openalea)
  * a SVN client (like Tortoise), and initialise it

SVN initialisation
==================

  * To be able to get files from the inria server, you will first need to create an account on the inria server
  * After a few hours, your account will be operational, and you can, using the "SVN checkout" command

Edtion and Generation of the documentation
==========================================

  * Once you have downloaded the files, you can edit them with any editor. Rst files are in echap/doc/user.
  * to see your changes, open a prompt in the echap/doc directory, and type 'make html'.
  * generated pages are in echap/doc/build/html (double-click on 'contents.html' to see the first page)
