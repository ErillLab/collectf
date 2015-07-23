We occasionally need to make some changes in the database.

In rare cases, we may need to add some new tables to the database. If that is
the case, Django's own `syncdb` utility should be enough to do the job.

    ./manage.py syncdb

Most of the time, however, schema changes are simpler, such as adding a new
field to an existing field or modifying an existing field in a table. Since we
are using Django 1.6, we rely on a third-party tool for migrations, called
[South](http://south.readthedocs.org/en/latest/index.html#). (Django 1.7 comes
with a built-in migration tool).

To use South,
- make sure it is installed: `sudo pip install south` should work.
- add south to `INSTALLED_APPS` list,

The procedure to change a model is as follows:

1. South complains if the migration history in the database and `migrations`
   directory (in the app directory) agrees. If you don't care about migration
   history,
   - remove the `migrations` folder
   - drop the migration-history table
     - go to Django database shell: `manage.py dbshell`
     - drop the table: `drop table south_migrationhistory`

2. If this is the first time you're using south or if you dropped the table in
   the previous step, you need to create migration-history table using
   `manage.py syncdb`.

3. Before adding the new field description in `models.py`, run `manage.py
   schemamigration --initial`. Since there already tables existing in the
   database, the first migration should be a "fake" one: `manage.py migrate
   <app-name> --fake`.

4. Finally, make the change in `models.py`, run `manage.py schemamigration
   --auto` followed by `manage.py migrate <app-name>`.

These steps should be sufficient to make a change in an existing table. For more
information, see [South documentation](https://south.readthedocs.org/en/latest/).
   

