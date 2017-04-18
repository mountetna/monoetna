### This is the command you want to run initialize the migrations.

sequel -m ./migrations postgres://[USER]:[PASS]@[HOST]/metis?search_path=private

Also have a look at file `002_add_default_user.rb`, the user info gets replaced
by our Chef deploy so you want to manual add your own user before you run the 
migration.
