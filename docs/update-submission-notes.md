
This file explains the way that submission process (and pretty much everything
else) was updated. It will be used as a guide when deploying all changes to the
server.

1. Update database schema.
   1. Use South migration tool.
   2. Before any changes in the model, migrate DB (with --initial)

   3. Make the changes and run south schemamigration again (with --auto)

   4. Apply changes between two migrations

2. Run the script to complete changes in DB scheme
   1. For curations, make TF-instance many-to-many
   2. For curation-site-instances, populate experimental-techniques field
   3. For curation-site-instances, populate TF-function field

3. Now, you can delete unwanted fields.
