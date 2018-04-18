~~~ ARC_decon README ~~~

~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~

TABLE OF CONTENTS

I. Intro;

II. Database Specs;

III. Classes;

IV. Spectral Deconvolution;

V. Spatial Deconvolution;

VI. Notes;

VII. A Basic PSQL Tutorial;

~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~

I. Intro

This program is intended to deconvolve spectral and spatial data read-in from the
  ARC robot system.

~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~

II. Database Specs

The database name being assumed by this program is "det_db", both for input and
  output.

The tables are as follows:

Detector output table: "det_out"

Heatmap table: "heatmap"

Spectral decon results table: "spect_res"

Spatial decon results table: "spatial_res"

~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~

III. Classes

This program uses 6 classes, each with its own header and source file (except the
  main class, of course). CMake is used to manage the build process.

decon_main.cc is the main process. It distributes control to the other five classes
  at various points.

gold_decon.cc and spatial_decon.cc are the two main processing classes. They are
  covered in more detail below.

data_read.cc reads the psql database output from the detector into the program.

histo.cc is a compact class that organizes the energy bin data into a histogram
  format fit for the spectral algorithm.

resp_read.cc reads the experimental response files from the execution directory
  and organizes their data into a format usable by the processing.

~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~

IV. Spectral Deconvolution



~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~

V. Spatial Deconvolution



~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~

VI. Notes

~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~

VII. A Basic PSQL Tutorial

Step 1: PSQL install

If on current linux, we can find postgreSQL on synaptic package manager. Simply
  search under "name" for "postgresql". Everything should work out nicely, no
  dpkg required.

Step 2: Creating a database

Open a shell (I like bash). We're gonna change user to root, then change to user
  "postgres" from there. To accomplish this, we do:

    sudo su

    su - postgres

Now, we can create psql users and databases. For now, we'll stick with user
  "postgres". To create a database, we do:

    createdb -O postgres testdb

  This creates a database called "testdb" with owner "postgres".

Now, before we go any further, we have to take a look at the psql network
  settings. Find the "pg_hba.conf" file that was installed with psql. By
  default, it should be located at /etc/postgresql/<version-number>/main/. Open
  this file, and take a look at the line under the heading "IPv4 local
  connections". Make a note of the IP address, then close the file. This is the
  local address we will use to make a connection to a local psql database.
  Usually, this address will be 127.0.0.1.

Step 3: Getting into the database

We can now make alterations to our database. But before that, we have to set a
  password for our new user "postgres". We can do this in another shell (with
  the default user). Open a shell, and do:

    sudo -u postgres psql postgres

  This will provide access into the psql shell layer. From there, we can do

    \password postgres

  By default, no password is set. Set the password to whatever you want, just
  don't forget it!

Now that we've done that, go back to the shell in which you have changed to
  user "postgres". To get into the database we created above, we can do:

    psql postgres -h <IP-address> -d testdb

  This will prompt for a password, which you should have now. We operating
  within the database and can make edits.

Step 4: Adding data

Tables are used to hold data within psql databases. To create a table while
  operating within a table, we can do:

~~~~~~~~~~~~~~~~~~~~~~~
  create table <table-name> (
    <data-spec1>    <data-type>
    <data-spec2>    <data-type>
    ...             ...
    )
~~~~~~~~~~~~~~~~~~~~~~~

  In this instantiation, the data-spec refers to the name of the data column.
  some examples of this could be "names", "phonenumbers" or "addresses". The
  data-type signifies what form the variables within the data columns will take.
  Examples are "integer", varchar(20) (which indicates a string of 20 chars or
  less) and char(10) (which indicates a numeral string of 10 digits or less).

  ~~~~~~~~~~~~~~~~~~~~~
