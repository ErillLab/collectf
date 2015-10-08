This directory contains all one-time scripts that are run for things like DB
schema changes, pulling some data from DB, RefSeq to UniProt acccession number
migraiton, etc.

To run a script accessing the database through Django models, you need to import
necessary modules/settings, which can be troublesome. The way to run scripts is
to put the script in this directory, define a `run()` function inside the file
which will be run when the following command is executed.

```
python manage.py runscript <script-filename>
```
