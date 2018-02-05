////////////////////////////////////
//                                //
//                                //
//    PSQL - A basic tutorial     //
//                                //
//                                //
////////////////////////////////////

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

  create table <table-name> (
    <data-spec1>    <data-type>
    <data-spec2>    <data-type>
    ...             ...
    )

  In this instantiation, the data-spec refers to the name of the data column.
  some examples of this could be "names", "phonenumbers" or "addresses". The
  data-type signifies what form the variables within the data columns will take.
  Examples are "integer", varchar(20) (which indicates a string of 20 chars or
  less) and char(10) (which indicates a numeral string of 10 digits or less).

  
